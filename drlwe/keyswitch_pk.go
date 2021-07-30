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
	ringQ         *ring.Ring
	ringP         *ring.Ring
	sigmaSmudging float64

	tmpQ       *ring.Poly
	tmpP       *ring.Poly
	share0tmpQ *ring.Poly
	share1tmpQ *ring.Poly
	share0tmpP *ring.Poly
	share1tmpP *ring.Poly

	baseconverter             *ring.FastBasisExtender
	gaussianSampler           *ring.GaussianSampler
	ternarySamplerMontgomeryQ *ring.TernarySampler
}

// NewPCKSProtocol creates a new PCKSProtocol object and will be used to re-encrypt a ciphertext ctx encrypted under a secret-shared key among j parties under a new
// collective public-key.
func NewPCKSProtocol(params rlwe.Parameters, sigmaSmudging float64) *PCKSProtocol {
	pcks := new(PCKSProtocol)
	pcks.ringQ = params.RingQ()
	pcks.ringP = params.RingP()
	pcks.sigmaSmudging = sigmaSmudging

	pcks.tmpQ = pcks.ringQ.NewPoly()
	pcks.tmpP = pcks.ringP.NewPoly()
	pcks.share0tmpQ = pcks.ringQ.NewPoly()
	pcks.share1tmpQ = pcks.ringQ.NewPoly()
	pcks.share0tmpP = pcks.ringP.NewPoly()
	pcks.share1tmpP = pcks.ringP.NewPoly()

	pcks.baseconverter = ring.NewFastBasisExtender(pcks.ringQ, pcks.ringP)
	prng, err := utils.NewPRNG()
	if err != nil {
		panic(err)
	}
	pcks.gaussianSampler = ring.NewGaussianSampler(prng, pcks.ringQ, sigmaSmudging, int(6*sigmaSmudging))
	pcks.ternarySamplerMontgomeryQ = ring.NewTernarySampler(prng, pcks.ringQ, 0.5, false)

	return pcks
}

// AllocateShare allocates the shares of the PCKS protocol
func (pcks *PCKSProtocol) AllocateShare(levelQ int) (s *PCKSShare) {
	return &PCKSShare{[2]*ring.Poly{pcks.ringQ.NewPolyLvl(levelQ), pcks.ringQ.NewPolyLvl(levelQ)}}
}

// GenShare is the first part of the unique round of the PCKSProtocol protocol. Each party computes the following :
//
// [s_i * ctx[0] + (u_i * pk[0] + e_0i)/P, (u_i * pk[1] + e_1i)/P]
//
// and broadcasts the result to the other j-1 parties.
func (pcks *PCKSProtocol) GenShare(sk *rlwe.SecretKey, pk *rlwe.PublicKey, ct *rlwe.Ciphertext, shareOut *PCKSShare) {

	el := ct.RLWEElement()

	ringQ := pcks.ringQ
	ringP := pcks.ringP

	levelQ := el.Level()
	levelP := len(pcks.ringP.Modulus) - 1

	// samples MForm(u_i) in Q and P separately
	pcks.ternarySamplerMontgomeryQ.ReadLvl(levelQ, pcks.tmpQ)
	extendBasisSmallNormAndCenter(ringQ, ringP, pcks.tmpQ, pcks.tmpP)
	ringQ.MFormLvl(levelQ, pcks.tmpQ, pcks.tmpQ)
	ringP.MFormLvl(levelP, pcks.tmpP, pcks.tmpP)
	ringQ.NTTLvl(levelQ, pcks.tmpQ, pcks.tmpQ)
	ringP.NTTLvl(levelP, pcks.tmpP, pcks.tmpP)

	// h_0 = NTT(u_i * pk_0)
	ringQ.MulCoeffsMontgomeryLvl(levelQ, pcks.tmpQ, pk.Value[0][0], pcks.share0tmpQ)
	ringP.MulCoeffsMontgomeryLvl(levelP, pcks.tmpP, pk.Value[0][1], pcks.share0tmpP)
	ringQ.InvNTTLvl(levelQ, pcks.share0tmpQ, pcks.share0tmpQ)
	ringP.InvNTTLvl(levelP, pcks.share0tmpP, pcks.share0tmpP)

	// h_1 = NTT(u_i * pk_1)
	ringQ.MulCoeffsMontgomeryLvl(levelQ, pcks.tmpQ, pk.Value[1][0], pcks.share1tmpQ)
	ringP.MulCoeffsMontgomeryLvl(levelP, pcks.tmpP, pk.Value[1][1], pcks.share1tmpP)
	ringQ.InvNTTLvl(levelQ, pcks.share1tmpQ, pcks.share1tmpQ)
	ringP.InvNTTLvl(levelP, pcks.share1tmpP, pcks.share1tmpP)

	// h_0 = u_i * pk_0 + e0
	pcks.gaussianSampler.ReadLvl(levelQ, pcks.tmpQ)
	extendBasisSmallNormAndCenter(ringQ, ringP, pcks.tmpQ, pcks.tmpP)
	ringQ.AddLvl(levelQ, pcks.share0tmpQ, pcks.tmpQ, pcks.share0tmpQ)
	ringP.AddLvl(levelP, pcks.share0tmpP, pcks.tmpP, pcks.share0tmpP)

	// h_1 = u_i * pk_1 + e1
	pcks.gaussianSampler.ReadLvl(levelQ, pcks.tmpQ)
	extendBasisSmallNormAndCenter(ringQ, ringP, pcks.tmpQ, pcks.tmpP)
	ringQ.AddLvl(levelQ, pcks.share1tmpQ, pcks.tmpQ, pcks.share1tmpQ)
	ringP.AddLvl(levelP, pcks.share1tmpP, pcks.tmpP, pcks.share1tmpP)

	// h_0 = (u_i * pk_0 + e0)/P
	pcks.baseconverter.ModDownQPtoQ(levelQ, levelP, pcks.share0tmpQ, pcks.share0tmpP, shareOut.Value[0])

	// h_1 = (u_i * pk_1 + e1)/P
	pcks.baseconverter.ModDownQPtoQ(levelQ, levelP, pcks.share1tmpQ, pcks.share1tmpP, shareOut.Value[1])

	// h_0 = s_i*c_1 + (u_i * pk_0 + e0)/P
	if el.Value[0].IsNTT {
		ringQ.NTTLvl(levelQ, shareOut.Value[0], shareOut.Value[0])
		ringQ.NTTLvl(levelQ, shareOut.Value[1], shareOut.Value[1])
		ringQ.MulCoeffsMontgomeryAndAddLvl(levelQ, el.Value[1], sk.Value[0], shareOut.Value[0])
	} else {
		// tmp = s_i*c_1
		ringQ.NTTLazyLvl(levelQ, el.Value[1], pcks.tmpQ)
		ringQ.MulCoeffsMontgomeryConstantLvl(levelQ, pcks.tmpQ, sk.Value[0], pcks.tmpQ)
		ringQ.InvNTTLvl(levelQ, pcks.tmpQ, pcks.tmpQ)

		// h_0 = s_i*c_1 + (u_i * pk_0 + e0)/P
		ringQ.AddLvl(levelQ, shareOut.Value[0], pcks.tmpQ, shareOut.Value[0])
	}
}

// AggregateShares is the second part of the first and unique round of the PCKSProtocol protocol. Each party uppon receiving the j-1 elements from the
// other parties computes :
//
// [ctx[0] + sum(s_i * ctx[0] + u_i * pk[0] + e_0i), sum(u_i * pk[1] + e_1i)]
func (pcks *PCKSProtocol) AggregateShares(share1, share2, shareOut *PCKSShare) {
	levelQ1, levelQ2 := len(share1.Value[0].Coeffs)-1, len(share2.Value[0].Coeffs)-1
	if levelQ1 != levelQ2 {
		panic("cannot aggreate two shares at different levelQs.")
	}
	pcks.ringQ.AddLvl(levelQ1, share1.Value[0], share2.Value[0], shareOut.Value[0])
	pcks.ringQ.AddLvl(levelQ1, share1.Value[1], share2.Value[1], shareOut.Value[1])
}

// KeySwitch performs the actual keyswitching operation on a ciphertext ct and put the result in ctOut
func (pcks *PCKSProtocol) KeySwitch(combined *PCKSShare, ct, ctOut *rlwe.Ciphertext) {
	el, elOut := ct.RLWEElement(), ctOut.RLWEElement()
	pcks.ringQ.AddLvl(el.Level(), el.Value[0], combined.Value[0], elOut.Value[0])
	ring.CopyValuesLvl(el.Level(), combined.Value[1], elOut.Value[1])
}

// MarshalBinary encodes a PCKS share on a slice of bytes.
func (share *PCKSShare) MarshalBinary() ([]byte, error) {
	lenR1 := share.Value[0].GetDataLen(true)
	lenR2 := share.Value[1].GetDataLen(true)

	data := make([]byte, lenR1+lenR2)
	_, err := share.Value[0].WriteTo(data[0:lenR1])
	if err != nil {
		return []byte{}, err
	}

	_, err = share.Value[1].WriteTo(data[lenR1 : lenR1+lenR2])
	if err != nil {
		return []byte{}, err
	}

	return data, nil
}

// UnmarshalBinary decodes marshaled PCKS share on the target PCKS share.
func (share *PCKSShare) UnmarshalBinary(data []byte) error {

	if share.Value[0] == nil {
		share.Value[0] = new(ring.Poly)
	}

	if share.Value[1] == nil {
		share.Value[1] = new(ring.Poly)
	}

	err := share.Value[0].UnmarshalBinary(data[0 : len(data)/2])
	if err != nil {
		return err
	}

	err = share.Value[1].UnmarshalBinary(data[len(data)/2:])
	if err != nil {
		return err
	}

	return nil
}
