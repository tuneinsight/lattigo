package drlwe

import (
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/rlwe"
	"github.com/ldsec/lattigo/v2/utils"
)

// KeySwitchingProtocol is an interface describing the local steps of a generic RLWE CKS protocol
type KeySwitchingProtocol interface {
	AllocateShare(level int) *CKSShare
	GenShare(skInput, skOutput *rlwe.SecretKey, ct *rlwe.Ciphertext, shareOut *CKSShare)
	AggregateShares(share1, share2, shareOut *CKSShare)
	KeySwitch(combined *CKSShare, ct, ctOut *rlwe.Ciphertext)
}

// CKSProtocol is the structure storing the parameters and and precomputations for the collective key-switching protocol.
type CKSProtocol struct {
	params          rlwe.Parameters
	gaussianSampler *ring.GaussianSampler
	baseconverter   *ring.FastBasisExtender
	tmpQP           rlwe.PolyQP
	tmpDelta        *ring.Poly
}

// CKSShare is a type for the CKS protocol shares.
type CKSShare struct {
	Value *ring.Poly
}

// MarshalBinary encodes a CKS share on a slice of bytes.
func (ckss *CKSShare) MarshalBinary() (data []byte, err error) {
	return ckss.Value.MarshalBinary()
}

// UnmarshalBinary decodes marshaled CKS share on the target CKS share.
func (ckss *CKSShare) UnmarshalBinary(data []byte) (err error) {
	ckss.Value = new(ring.Poly)
	return ckss.Value.UnmarshalBinary(data)
}

// NewCKSProtocol creates a new CKSProtocol that will be used to operate a collective key-switching on a ciphertext encrypted under a collective public-key, whose
// secret-shares are distributed among j parties, re-encrypting the ciphertext under another public-key, whose secret-shares are also known to the
// parties.
func NewCKSProtocol(params rlwe.Parameters, sigmaSmudging float64) *CKSProtocol {
	cks := new(CKSProtocol)
	cks.params = params
	prng, err := utils.NewPRNG()
	if err != nil {
		panic(err)
	}
	cks.gaussianSampler = ring.NewGaussianSampler(prng, params.RingQ(), sigmaSmudging, int(6*sigmaSmudging))
	cks.baseconverter = ring.NewFastBasisExtender(params.RingQ(), params.RingP())
	cks.tmpQP = params.RingQP().NewPoly()
	cks.tmpDelta = params.RingQ().NewPoly()
	return cks
}

// AllocateShare allocates the shares of the CKSProtocol
func (cks *CKSProtocol) AllocateShare(level int) *CKSShare {
	return &CKSShare{cks.params.RingQ().NewPolyLvl(level)}
}

// SampleCRP samples a common random polynomial to be used in the CKS protocol from the provided
// common reference string.
func (cks *CKSProtocol) SampleCRP(level int, crs CRS) CRP {
	return NewCRPAtLvl(cks.params, 1, level, -1, crs)
}

// GenShare computes a party's share in the CKS protocol.
// ct.Value[0] can be nil, computations are only done using ct.Value[1]
// NTT flag for ct.Value[1] is expected to be set correctly
func (cks *CKSProtocol) GenShare(skInput, skOutput *rlwe.SecretKey, ct *rlwe.Ciphertext, shareOut *CKSShare) {

	el := ct.RLWEElement()

	ringQ := cks.params.RingQ()
	ringP := cks.params.RingP()
	ringQP := cks.params.RingQP()

	level := utils.MinInt(len(ringQ.Modulus)-1, el.Value[1].Level())
	levelP := cks.params.PCount() - 1

	ringQ.SubLvl(level, skInput.Value.Q, skOutput.Value.Q, cks.tmpDelta)

	ct1 := el.Value[1]
	if !el.Value[1].IsNTT {
		ringQ.NTTLazyLvl(level, el.Value[1], cks.tmpQP.Q)
		ct1 = cks.tmpQP.Q
	}

	// a * (skIn - skOut) mod Q
	ringQ.MulCoeffsMontgomeryConstantLvl(level, ct1, cks.tmpDelta, shareOut.Value)

	// P * a * (skIn - skOut) mod QP (mod P = 0)
	ringQ.MulScalarBigintLvl(level, shareOut.Value, ringP.ModulusBigint, shareOut.Value)

	if !el.Value[1].IsNTT {
		// InvNTT(P * a * (skIn - skOut)) mod QP (mod P = 0)
		ringQ.InvNTTLazyLvl(level, shareOut.Value, shareOut.Value)

		// Samples e in Q
		cks.gaussianSampler.ReadLvl(level, cks.tmpQP.Q)

		// Extend e to P (assumed to have norm < qi)
		ringQP.ExtendBasisSmallNormAndCenter(cks.tmpQP.Q, levelP, nil, cks.tmpQP.P)

		// InvNTT(P * a * (skIn - skOut) + e) mod QP (mod P = e)
		ringQ.AddNoModLvl(level, shareOut.Value, cks.tmpQP.Q, shareOut.Value)

		// InvNTT(P * a * (skIn - skOut) + e) * (1/P) mod QP (mod P = e)
		cks.baseconverter.ModDownQPtoQ(level, levelP, shareOut.Value, cks.tmpQP.P, shareOut.Value)

	} else {
		// Sample e in Q
		cks.gaussianSampler.ReadLvl(level, cks.tmpQP.Q)

		// Extend e to P (assumed to have norm < qi)
		ringQP.ExtendBasisSmallNormAndCenter(cks.tmpQP.Q, levelP, nil, cks.tmpQP.P)

		// Takes the error to the NTT domain
		ringQ.InvNTTLvl(level, shareOut.Value, shareOut.Value)

		// P * a * (skIn - skOut) + e mod Q (mod P = 0, so P = e)
		ringQ.AddLvl(level, shareOut.Value, cks.tmpQP.Q, shareOut.Value)

		// (P * a * (skIn - skOut) + e) * (1/P) mod QP (mod P = e)
		cks.baseconverter.ModDownQPtoQ(level, levelP, shareOut.Value, cks.tmpQP.P, shareOut.Value)

		ringQ.NTTLvl(level, shareOut.Value, shareOut.Value)
	}

	shareOut.Value.Coeffs = shareOut.Value.Coeffs[:level+1]
}

// AggregateShares is the second part of the unique round of the CKSProtocol protocol. Upon receiving the j-1 elements each party computes :
//
// [ctx[0] + sum((skInput_i - skOutput_i) * ctx[0] + e_i), ctx[1]]
func (cks *CKSProtocol) AggregateShares(share1, share2, shareOut *CKSShare) {
	cks.params.RingQ().AddLvl(share1.Value.Level(), share1.Value, share2.Value, shareOut.Value)
}

// KeySwitch performs the actual keyswitching operation on a ciphertext ct and put the result in ctOut
func (cks *CKSProtocol) KeySwitch(combined *CKSShare, ct, ctOut *rlwe.Ciphertext) {
	el, elOut := ct.RLWEElement(), ctOut.RLWEElement()
	cks.params.RingQ().AddLvl(el.Level(), el.Value[0], combined.Value, elOut.Value[0])
	ring.CopyValuesLvl(el.Level(), el.Value[1], elOut.Value[1])
}
