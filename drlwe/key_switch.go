package drlwe

import (
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/rlwe"
	"github.com/ldsec/lattigo/v2/utils"
)

// KeySwitchingProtocol is an interface describing the local steps of a generic RLWE CKS protocol
type KeySwitchingProtocol interface {
	AllocateShare() *CKSShare
	GenShare(skInput, skOutput *rlwe.SecretKey, ct rlwe.Ciphertext, shareOut *CKSShare)
	AggregateShares(share1, share2, shareOut *CKSShare)
	KeySwitch(combined *CKSShare, ct rlwe.Ciphertext, ctOut rlwe.Ciphertext)
}

// CKSProtocol is the structure storing the parameters and and precomputations for the collective key-switching protocol.
type CKSProtocol struct {
	ringQ           *ring.Ring
	ringP           *ring.Ring
	ringQP          *ring.Ring
	gaussianSampler *ring.GaussianSampler
	baseconverter   *ring.FastBasisExtender

	tmpP     *ring.Poly
	tmpQ     *ring.Poly
	tmpDelta *ring.Poly
	tmpQP    *ring.Poly
	hP       *ring.Poly
}

// CKSShare is a type for the CKS protocol shares.
type CKSShare struct {
	Value *ring.Poly
}

// MarshalBinary encodes a CKS share on a slice of bytes.
func (ckss *CKSShare) MarshalBinary() ([]byte, error) {
	return ckss.Value.MarshalBinary()
}

// UnmarshalBinary decodes marshaled CKS share on the target CKS share.
func (ckss *CKSShare) UnmarshalBinary(data []byte) error {
	if ckss.Value == nil {
		ckss.Value = new(ring.Poly)
	}
	return ckss.Value.UnmarshalBinary(data)
}

// NewCKSProtocol creates a new CKSProtocol that will be used to operate a collective key-switching on a ciphertext encrypted under a collective public-key, whose
// secret-shares are distributed among j parties, re-encrypting the ciphertext under another public-key, whose secret-shares are also known to the
// parties.
func NewCKSProtocol(params rlwe.Parameters, sigmaSmudging float64) *CKSProtocol {
	cks := new(CKSProtocol)
	cks.ringQ = params.RingQ()
	cks.ringP = params.RingP() // TODO this assumes that P larger than 1
	cks.ringQP = params.RingQP()

	prng, err := utils.NewPRNG()
	if err != nil {
		panic(err)
	}
	cks.gaussianSampler = ring.NewGaussianSampler(prng, cks.ringQP, sigmaSmudging, int(6*sigmaSmudging))

	cks.baseconverter = ring.NewFastBasisExtender(cks.ringQ, cks.ringP)

	cks.tmpQ = cks.ringQ.NewPoly()
	cks.tmpP = cks.ringP.NewPoly()
	cks.tmpDelta = cks.ringQ.NewPoly()
	cks.tmpQP = cks.ringQP.NewPoly()
	cks.hP = cks.ringP.NewPoly()

	return cks
}

// AllocateShare allocates the shares of the CKSProtocol
func (cks *CKSProtocol) AllocateShare() *CKSShare {
	return &CKSShare{cks.ringQ.NewPoly()}
}

// GenShare computes a party's share in the CKS protocol.
func (cks *CKSProtocol) GenShare(skInput, skOutput *rlwe.SecretKey, ct rlwe.Ciphertext, shareOut *CKSShare) {

	cks.ringQ.Sub(skInput.Value, skOutput.Value, cks.tmpDelta)

	el := ct.RLWEElement()
	ct1 := el.Value[1]
	if !el.IsNTT {
		cks.ringQ.NTTLazy(el.Value[1], cks.tmpQP)
		ct1 = cks.tmpQP
	}

	cks.ringQ.MulCoeffsMontgomeryConstantLvl(el.Level(), ct1, cks.tmpDelta, shareOut.Value)
	cks.ringQ.MulScalarBigintLvl(el.Level(), shareOut.Value, cks.ringP.ModulusBigint, shareOut.Value)

	if !el.IsNTT {
		cks.ringQ.InvNTTLazy(shareOut.Value, shareOut.Value)

		cks.gaussianSampler.ReadLvl(len(cks.ringQP.Modulus)-1, cks.tmpQP)
		cks.ringQ.AddNoMod(shareOut.Value, cks.tmpQP, shareOut.Value)

		for x, i := 0, uint64(len(cks.ringQ.Modulus)); i < uint64(len(cks.ringQP.Modulus)); x, i = x+1, i+1 {
			tmphP := cks.hP.Coeffs[x]
			tmpNTT := cks.tmpQP.Coeffs[i]
			for j := 0; j < cks.ringQ.N; j++ {
				tmphP[j] += tmpNTT[j]
			}
		}

		cks.baseconverter.ModDownSplitPQ(el.Level(), shareOut.Value, cks.hP, shareOut.Value)
		cks.hP.Zero() // hP is assumed 0 above
	} else {
		cks.gaussianSampler.ReadLvl(el.Level(), cks.tmpQ)

		extendBasisSmallNormAndCenter(cks.ringQ, cks.ringP, cks.tmpQ, cks.tmpP)

		cks.ringQ.NTTLvl(el.Level(), cks.tmpQ, cks.tmpQ)
		cks.ringP.NTT(cks.tmpP, cks.tmpP)

		cks.ringQ.AddLvl(el.Level(), shareOut.Value, cks.tmpQ, shareOut.Value)

		cks.baseconverter.ModDownSplitNTTPQ(el.Level(), shareOut.Value, cks.tmpP, shareOut.Value)
	}
}

// AggregateShares is the second part of the unique round of the CKSProtocol protocol. Upon receiving the j-1 elements each party computes :
//
// [ctx[0] + sum((skInput_i - skOutput_i) * ctx[0] + e_i), ctx[1]]
func (cks *CKSProtocol) AggregateShares(share1, share2, shareOut *CKSShare) {
	cks.ringQ.AddLvl(len(share1.Value.Coeffs)-1, share1.Value, share2.Value, shareOut.Value)
}

// KeySwitch performs the actual keyswitching operation on a ciphertext ct and put the result in ctOut
func (cks *CKSProtocol) KeySwitch(combined *CKSShare, ct rlwe.Ciphertext, ctOut rlwe.Ciphertext) {
	el, elOut := ct.RLWEElement(), ctOut.RLWEElement()
	cks.ringQ.AddLvl(el.Level(), el.Value[0], combined.Value, elOut.Value[0])
	cks.ringQ.CopyLvl(el.Level(), el.Value[1], elOut.Value[1])
}

func extendBasisSmallNormAndCenter(ringQ, ringP *ring.Ring, polQ, polP *ring.Poly) {
	var coeff, Q, QHalf, sign uint64
	Q = ringQ.Modulus[0]
	QHalf = Q >> 1

	for j := 0; j < ringQ.N; j++ {

		coeff = polQ.Coeffs[0][j]

		sign = 1
		if coeff > QHalf {
			coeff = Q - coeff
			sign = 0
		}

		for i, pi := range ringP.Modulus {
			polP.Coeffs[i][j] = (coeff * sign) | (pi-coeff)*(sign^1)
		}
	}
}
