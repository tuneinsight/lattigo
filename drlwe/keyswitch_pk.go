package drlwe

import (
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/rlwe"
	"github.com/ldsec/lattigo/v2/utils"
)

// PublicKeySwitchingProtocol is an interface describing the local steps of a generic RLWE PCKS protocol.
type PublicKeySwitchingProtocol interface {
	AllocateShare(levelQ int) *PCKSShare
	GenShare(skInput *rlwe.SecretKey, pkOutput *rlwe.PublicKey, ct *rlwe.Ciphertext, shareOut *PCKSShare)
	AggregateShares(share1, share2, shareOut *PCKSShare)
	KeySwitch(combined *PCKSShare, ct *rlwe.Ciphertext, ctOut *rlwe.Ciphertext)
}

// PCKSShare represents a party's share in the PCKS protocol.
type PCKSShare struct {
	Value [2]*ring.Poly
}

// PCKSProtocol is the structure storing the parameters for the collective public key-switching.
type PCKSProtocol struct {
	params        rlwe.Parameters
	sigmaSmudging float64

	tmpQP rlwe.PolyQP
	tmpP  [2]*ring.Poly

	baseconverter             *ring.FastBasisExtender
	gaussianSampler           *ring.GaussianSampler
	ternarySamplerMontgomeryQ *ring.TernarySampler
}

// NewPCKSProtocol creates a new PCKSProtocol object and will be used to re-encrypt a ciphertext ctx encrypted under a secret-shared key among j parties under a new
// collective public-key.
func NewPCKSProtocol(params rlwe.Parameters, sigmaSmudging float64) (pcks *PCKSProtocol) {
	pcks = new(PCKSProtocol)
	pcks.params = params
	pcks.sigmaSmudging = sigmaSmudging

	pcks.tmpQP = params.RingQP().NewPoly()
	pcks.tmpP = [2]*ring.Poly{params.RingP().NewPoly(), params.RingP().NewPoly()}

	pcks.baseconverter = ring.NewFastBasisExtender(params.RingQ(), params.RingP())
	prng, err := utils.NewPRNG()
	if err != nil {
		panic(err)
	}
	pcks.gaussianSampler = ring.NewGaussianSampler(prng, params.RingQ(), sigmaSmudging, int(6*sigmaSmudging))
	pcks.ternarySamplerMontgomeryQ = ring.NewTernarySampler(prng, params.RingQ(), 0.5, false)

	return pcks
}

// AllocateShare allocates the shares of the PCKS protocol
func (pcks *PCKSProtocol) AllocateShare(levelQ int) (s *PCKSShare) {
	return &PCKSShare{[2]*ring.Poly{pcks.params.RingQ().NewPolyLvl(levelQ), pcks.params.RingQ().NewPolyLvl(levelQ)}}
}

// GenShare is the first part of the unique round of the PCKSProtocol protocol. Each party computes the following :
//
// [s_i * ctx[0] + (u_i * pk[0] + e_0i)/P, (u_i * pk[1] + e_1i)/P]
//
// and broadcasts the result to the other j-1 parties.
func (pcks *PCKSProtocol) GenShare(sk *rlwe.SecretKey, pk *rlwe.PublicKey, ct *rlwe.Ciphertext, shareOut *PCKSShare) {

	el := ct.RLWEElement()

	ringQ := pcks.params.RingQ()
	ringP := pcks.params.RingP()
	ringQP := pcks.params.RingQP()

	levelQ := el.Level()
	levelP := len(ringP.Modulus) - 1

	// samples MForm(u_i) in Q and P separately
	pcks.ternarySamplerMontgomeryQ.ReadLvl(levelQ, pcks.tmpQP.Q)
	ringQP.ExtendBasisSmallNormAndCenter(pcks.tmpQP.Q, levelP, nil, pcks.tmpQP.P)
	ringQP.MFormLvl(levelQ, levelP, pcks.tmpQP, pcks.tmpQP)
	ringQP.NTTLvl(levelQ, levelP, pcks.tmpQP, pcks.tmpQP)

	shareOutQP0 := rlwe.PolyQP{shareOut.Value[0], pcks.tmpP[0]}
	shareOutQP1 := rlwe.PolyQP{shareOut.Value[1], pcks.tmpP[1]}

	// h_0 = u_i * pk_0
	// h_1 = u_i * pk_1
	ringQP.MulCoeffsMontgomeryLvl(levelQ, levelP, pcks.tmpQP, pk.Value[0], shareOutQP0)
	ringQP.MulCoeffsMontgomeryLvl(levelQ, levelP, pcks.tmpQP, pk.Value[1], shareOutQP1)

	ringQP.InvNTTLvl(levelQ, levelP, shareOutQP0, shareOutQP0)
	ringQP.InvNTTLvl(levelQ, levelP, shareOutQP1, shareOutQP1)

	// h_0 = u_i * pk_0
	pcks.gaussianSampler.ReadLvl(levelQ, pcks.tmpQP.Q)
	ringQP.ExtendBasisSmallNormAndCenter(pcks.tmpQP.Q, levelP, nil, pcks.tmpQP.P)
	ringQP.AddLvl(levelQ, levelP, shareOutQP0, pcks.tmpQP, shareOutQP0)

	// h_1 = u_i * pk_1 + e1
	pcks.gaussianSampler.ReadLvl(levelQ, pcks.tmpQP.Q)
	ringQP.ExtendBasisSmallNormAndCenter(pcks.tmpQP.Q, levelP, nil, pcks.tmpQP.P)
	ringQP.AddLvl(levelQ, levelP, shareOutQP1, pcks.tmpQP, shareOutQP1)

	// h_0 = (u_i * pk_0 + e0)/P
	pcks.baseconverter.ModDownQPtoQ(levelQ, levelP, shareOutQP0.Q, shareOutQP0.P, shareOutQP0.Q)

	// h_1 = (u_i * pk_1 + e1)/P
	pcks.baseconverter.ModDownQPtoQ(levelQ, levelP, shareOutQP1.Q, shareOutQP1.P, shareOutQP1.Q)

	// h_0 = s_i*c_1 + (u_i * pk_0 + e0)/P
	if el.Value[0].IsNTT {
		ringQ.NTTLvl(levelQ, shareOut.Value[0], shareOut.Value[0])
		ringQ.NTTLvl(levelQ, shareOut.Value[1], shareOut.Value[1])
		ringQ.MulCoeffsMontgomeryAndAddLvl(levelQ, el.Value[1], sk.Value.Q, shareOut.Value[0])
	} else {
		// tmp = s_i*c_1
		ringQ.NTTLazyLvl(levelQ, el.Value[1], pcks.tmpQP.Q)
		ringQ.MulCoeffsMontgomeryConstantLvl(levelQ, pcks.tmpQP.Q, sk.Value.Q, pcks.tmpQP.Q)
		ringQ.InvNTTLvl(levelQ, pcks.tmpQP.Q, pcks.tmpQP.Q)

		// h_0 = s_i*c_1 + (u_i * pk_0 + e0)/P
		ringQ.AddLvl(levelQ, shareOut.Value[0], pcks.tmpQP.Q, shareOut.Value[0])
	}
}

// AggregateShares is the second part of the first and unique round of the PCKSProtocol protocol. Each party uppon receiving the j-1 elements from the
// other parties computes :
//
// [ctx[0] + sum(s_i * ctx[0] + u_i * pk[0] + e_0i), sum(u_i * pk[1] + e_1i)]
func (pcks *PCKSProtocol) AggregateShares(share1, share2, shareOut *PCKSShare) {
	levelQ1, levelQ2 := len(share1.Value[0].Coeffs)-1, len(share2.Value[1].Coeffs)-1
	if levelQ1 != levelQ2 {
		panic("cannot aggreate two shares at different levelQs.")
	}
	pcks.params.RingQ().AddLvl(levelQ1, share1.Value[0], share2.Value[0], shareOut.Value[0])
	pcks.params.RingQ().AddLvl(levelQ1, share1.Value[1], share2.Value[1], shareOut.Value[1])

}

// KeySwitch performs the actual keyswitching operation on a ciphertext ct and put the result in ctOut
func (pcks *PCKSProtocol) KeySwitch(combined *PCKSShare, ct, ctOut *rlwe.Ciphertext) {
	el, elOut := ct.RLWEElement(), ctOut.RLWEElement()
	pcks.params.RingQ().AddLvl(el.Level(), el.Value[0], combined.Value[0], elOut.Value[0])
	ring.CopyValuesLvl(el.Level(), combined.Value[1], elOut.Value[1])
}

// MarshalBinary encodes a PCKS share on a slice of bytes.
func (share *PCKSShare) MarshalBinary() (data []byte, err error) {
	data = make([]byte, share.Value[0].GetDataLen(true)+share.Value[1].GetDataLen(true))
	var inc, pt int
	if inc, err = share.Value[0].WriteTo(data[pt:]); err != nil {
		return nil, err
	}
	pt += inc

	if _, err = share.Value[1].WriteTo(data[pt:]); err != nil {
		return nil, err
	}
	return
}

// UnmarshalBinary decodes marshaled PCKS share on the target PCKS share.
func (share *PCKSShare) UnmarshalBinary(data []byte) (err error) {
	var pt, inc int
	share.Value[0] = new(ring.Poly)
	if inc, err = share.Value[0].DecodePolyNew(data[pt:]); err != nil {
		return
	}
	pt += inc

	share.Value[1] = new(ring.Poly)
	if _, err = share.Value[1].DecodePolyNew(data[pt:]); err != nil {
		return
	}
	return
}
