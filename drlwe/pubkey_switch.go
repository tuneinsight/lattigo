package drlwe

import (
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/rlwe"
	"github.com/ldsec/lattigo/v2/utils"
)

// PublicKeySwitchingProtocol is an interface describing the local steps of a generic RLWE PCKS protocol.
type PublicKeySwitchingProtocol interface {
	AllocateShare(level uint64) *PCKSShare
	GenShare(skInput *rlwe.SecretKey, pkOutput *rlwe.PublicKey, ct rlwe.Ciphertext, shareOut *PCKSShare)
	AggregateShares(share1, share2, shareOut *PCKSShare)
	KeySwitch(combined *PCKSShare, ct rlwe.Ciphertext, ctOut rlwe.Ciphertext)
}

// PCKSShare represents a party's share in the PCKS protocol.
type PCKSShare struct {
	Value [2]*ring.Poly
}

// PCKSProtocol is the structure storing the parameters for the collective public key-switching.
type PCKSProtocol struct {
	ringQ         *ring.Ring
	ringP         *ring.Ring
	ringQP        *ring.Ring
	sigmaSmudging float64

	tmp       *ring.Poly
	share0tmp *ring.Poly
	share1tmp *ring.Poly

	baseconverter            *ring.FastBasisExtender
	gaussianSampler          *ring.GaussianSampler
	ternarySamplerMontgomery *ring.TernarySampler
}

// NewPCKSProtocol creates a new PCKSProtocol object and will be used to re-encrypt a ciphertext ctx encrypted under a secret-shared key among j parties under a new
// collective public-key.
func NewPCKSProtocol(params rlwe.Parameters, sigmaSmudging float64) *PCKSProtocol {
	pcks := new(PCKSProtocol)
	pcks.ringQ = params.RingQ()
	pcks.ringP = params.RingP()
	pcks.ringQP = params.RingQP()
	pcks.sigmaSmudging = sigmaSmudging

	pcks.tmp = pcks.ringQP.NewPoly()
	pcks.share0tmp = pcks.ringQP.NewPoly()
	pcks.share1tmp = pcks.ringQP.NewPoly()

	pcks.baseconverter = ring.NewFastBasisExtender(pcks.ringQ, pcks.ringP)
	prng, err := utils.NewPRNG()
	if err != nil {
		panic(err)
	}
	pcks.gaussianSampler = ring.NewGaussianSampler(prng, pcks.ringQP, sigmaSmudging, int(6*sigmaSmudging))
	pcks.ternarySamplerMontgomery = ring.NewTernarySampler(prng, pcks.ringQP, 0.5, true)

	return pcks
}

// AllocateShares allocates the shares of the PCKS protocol
func (pcks *PCKSProtocol) AllocateShares(level int) (s *PCKSShare) {
	s = new(PCKSShare)
	s.Value[0] = pcks.ringQ.NewPolyLvl(level)
	s.Value[1] = pcks.ringQ.NewPolyLvl(level)
	return
}

// GenShare is the first part of the unique round of the PCKSProtocol protocol. Each party computes the following :
//
// [s_i * ctx[0] + (u_i * pk[0] + e_0i)/P, (u_i * pk[1] + e_1i)/P]
//
// and broadcasts the result to the other j-1 parties.
func (pcks *PCKSProtocol) GenShare(sk *rlwe.SecretKey, pk *rlwe.PublicKey, ct rlwe.Ciphertext, shareOut *PCKSShare) {

	el := ct.RLWEElement()

	// samples u_i
	pcks.ternarySamplerMontgomery.Read(pcks.tmp)
	pcks.ringQP.NTTLazy(pcks.tmp, pcks.tmp)
	// h_0 = u_i * pk_0
	pcks.ringQP.MulCoeffsMontgomeryConstant(pcks.tmp, pk.Value[0], pcks.share0tmp)
	// h_1 = u_i * pk_1
	pcks.ringQP.MulCoeffsMontgomeryConstant(pcks.tmp, pk.Value[1], pcks.share1tmp)

	// h_0 = u_i * pk_0 + e0
	pcks.gaussianSampler.Read(pcks.tmp)
	pcks.ringQP.NTT(pcks.tmp, pcks.tmp)
	pcks.ringQP.Add(pcks.share0tmp, pcks.tmp, pcks.share0tmp)
	// h_1 = u_i * pk_1 + e1
	pcks.gaussianSampler.Read(pcks.tmp)
	pcks.ringQP.NTT(pcks.tmp, pcks.tmp)
	pcks.ringQP.Add(pcks.share1tmp, pcks.tmp, pcks.share1tmp)

	// h_0 = (u_i * pk_0 + e0)/P
	pcks.baseconverter.ModDownNTTPQ(el.Level(), pcks.share0tmp, shareOut.Value[0])

	// h_1 = (u_i * pk_1 + e1)/P
	// Cound be moved to the keyswitch part of the protocol, but the second element of the shares will be larger.
	pcks.baseconverter.ModDownNTTPQ(el.Level(), pcks.share1tmp, shareOut.Value[1])

	var c1NTT *ring.Poly
	if el.IsNTT {
		c1NTT = el.Value[1]
	} else {
		pcks.ringQ.NTT(el.Value[1], pcks.tmp)
		c1NTT = pcks.tmp
	}

	// h_0 = s_i*c_1 + (u_i * pk_0 + e0)/P
	pcks.ringQ.MulCoeffsMontgomeryAndAddLvl(el.Level(), c1NTT, sk.Value, shareOut.Value[0])

	if !el.IsNTT { // if element was not in NTT, put it back in non-ntt form
		pcks.ringQ.InvNTT(shareOut.Value[0], shareOut.Value[0])
		pcks.ringQ.InvNTT(shareOut.Value[1], shareOut.Value[1])
	}
}

// AggregateShares is the second part of the first and unique round of the PCKSProtocol protocol. Each party uppon receiving the j-1 elements from the
// other parties computes :
//
// [ctx[0] + sum(s_i * ctx[0] + u_i * pk[0] + e_0i), sum(u_i * pk[1] + e_1i)]
func (pcks *PCKSProtocol) AggregateShares(share1, share2, shareOut *PCKSShare) {
	level1, level2 := len(share1.Value[0].Coeffs)-1, len(share2.Value[0].Coeffs)-1
	if level1 != level2 {
		panic("cannot aggreate two shares at different levels.")
	}
	pcks.ringQ.AddLvl(level1, share1.Value[0], share2.Value[0], shareOut.Value[0])
	pcks.ringQ.AddLvl(level1, share1.Value[1], share2.Value[1], shareOut.Value[1])
}

// KeySwitch performs the actual keyswitching operation on a ciphertext ct and put the result in ctOut
func (pcks *PCKSProtocol) KeySwitch(combined *PCKSShare, ct, ctOut rlwe.Ciphertext) {
	el, elOut := ct.RLWEElement(), ctOut.RLWEElement()
	pcks.ringQ.AddLvl(el.Level(), el.Value[0], combined.Value[0], elOut.Value[0])
	pcks.ringQ.CopyLvl(el.Level(), combined.Value[1], elOut.Value[1])
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
