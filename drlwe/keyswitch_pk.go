package drlwe

import (
	"github.com/tuneinsight/lattigo/v3/ring"
	"github.com/tuneinsight/lattigo/v3/rlwe"
	"github.com/tuneinsight/lattigo/v3/rlwe/ringqp"
	"github.com/tuneinsight/lattigo/v3/utils"
)

// PublicKeySwitchingProtocol is an interface describing the local steps of a generic RLWE PCKS protocol.
type PublicKeySwitchingProtocol interface {
	AllocateShare(levelQ int) *PCKSShare
	GenShare(skInput *rlwe.SecretKey, pkOutput *rlwe.PublicKey, c1 *ring.Poly, shareOut *PCKSShare)
	AggregateShare(share1, share2, shareOut *PCKSShare)
	KeySwitch(ctIn *rlwe.Ciphertext, combined *PCKSShare, ctOut *rlwe.Ciphertext)
}

// PCKSShare represents a party's share in the PCKS protocol.
type PCKSShare struct {
	Value [2]*ring.Poly
}

// PCKSProtocol is the structure storing the parameters for the collective public key-switching.
type PCKSProtocol struct {
	params        rlwe.Parameters
	sigmaSmudging float64

	tmpQP ringqp.Poly
	tmpP  [2]*ring.Poly

	basisExtender             *ring.BasisExtender
	gaussianSampler           *ring.GaussianSampler
	ternarySamplerMontgomeryQ *ring.TernarySampler
}

// ShallowCopy creates a shallow copy of PCKSProtocol in which all the read-only data-structures are
// shared with the receiver and the temporary buffers are reallocated. The receiver and the returned
// PCKSProtocol can be used concurrently.
func (pcks *PCKSProtocol) ShallowCopy() *PCKSProtocol {
	prng, err := utils.NewPRNG()
	if err != nil {
		panic(err)
	}

	params := pcks.params

	var tmpP [2]*ring.Poly
	if params.RingP() != nil {
		tmpP = [2]*ring.Poly{params.RingP().NewPoly(), params.RingP().NewPoly()}
	}

	return &PCKSProtocol{
		params:                    params,
		sigmaSmudging:             pcks.sigmaSmudging,
		tmpQP:                     params.RingQP().NewPoly(),
		tmpP:                      tmpP,
		basisExtender:             pcks.basisExtender.ShallowCopy(),
		gaussianSampler:           ring.NewGaussianSampler(prng, params.RingQ(), pcks.sigmaSmudging, int(6*pcks.sigmaSmudging)),
		ternarySamplerMontgomeryQ: ring.NewTernarySamplerWithHammingWeight(prng, params.RingQ(), params.HammingWeight(), false),
	}
}

// NewPCKSProtocol creates a new PCKSProtocol object and will be used to re-encrypt a ciphertext ctx encrypted under a secret-shared key among j parties under a new
// collective public-key.
func NewPCKSProtocol(params rlwe.Parameters, sigmaSmudging float64) (pcks *PCKSProtocol) {
	pcks = new(PCKSProtocol)
	pcks.params = params
	pcks.sigmaSmudging = sigmaSmudging

	pcks.tmpQP = params.RingQP().NewPoly()

	if params.RingP() != nil {
		pcks.basisExtender = ring.NewBasisExtender(params.RingQ(), params.RingP())
		pcks.tmpP = [2]*ring.Poly{params.RingP().NewPoly(), params.RingP().NewPoly()}
	}

	prng, err := utils.NewPRNG()
	if err != nil {
		panic(err)
	}
	pcks.gaussianSampler = ring.NewGaussianSampler(prng, params.RingQ(), sigmaSmudging, int(6*sigmaSmudging))
	pcks.ternarySamplerMontgomeryQ = ring.NewTernarySamplerWithHammingWeight(prng, params.RingQ(), params.HammingWeight(), false)

	return pcks
}

// AllocateShare allocates the shares of the PCKS protocol.
func (pcks *PCKSProtocol) AllocateShare(levelQ int) (s *PCKSShare) {
	return &PCKSShare{[2]*ring.Poly{pcks.params.RingQ().NewPolyLvl(levelQ), pcks.params.RingQ().NewPolyLvl(levelQ)}}
}

// GenShare is the first part of the unique round of the PCKSProtocol protocol. Each party computes the following :
//
// [s_i * ct[1] + (u_i * pk[0] + e_0i)/P, (u_i * pk[1] + e_1i)/P]
//
// and broadcasts the result to the other j-1 parties.
// ct1 is the degree 1 element of the rlwe.Ciphertext to keyswitch, i.e. ct1 = rlwe.Ciphertext.Value[1].
// NTT flag for ct1 is expected to be set correctly.
func (pcks *PCKSProtocol) GenShare(sk *rlwe.SecretKey, pk *rlwe.PublicKey, ct1 *ring.Poly, shareOut *PCKSShare) {

	ringQ := pcks.params.RingQ()
	ringP := pcks.params.RingP()
	ringQP := pcks.params.RingQP()

	levelQ := utils.MinInt(shareOut.Value[0].Level(), ct1.Level())
	var levelP int
	if ringP != nil {
		levelP = len(ringP.Modulus) - 1
	}

	// samples MForm(u_i) in Q and P separately
	pcks.ternarySamplerMontgomeryQ.ReadLvl(levelQ, pcks.tmpQP.Q)

	if ringP != nil {
		ringQP.ExtendBasisSmallNormAndCenter(pcks.tmpQP.Q, levelP, nil, pcks.tmpQP.P)
	}

	ringQP.MFormLvl(levelQ, levelP, pcks.tmpQP, pcks.tmpQP)
	ringQP.NTTLvl(levelQ, levelP, pcks.tmpQP, pcks.tmpQP)

	shareOutQP0 := ringqp.Poly{Q: shareOut.Value[0], P: pcks.tmpP[0]}
	shareOutQP1 := ringqp.Poly{Q: shareOut.Value[1], P: pcks.tmpP[1]}

	// h_0 = u_i * pk_0
	// h_1 = u_i * pk_1
	ringQP.MulCoeffsMontgomeryLvl(levelQ, levelP, pcks.tmpQP, pk.Value[0], shareOutQP0)
	ringQP.MulCoeffsMontgomeryLvl(levelQ, levelP, pcks.tmpQP, pk.Value[1], shareOutQP1)

	ringQP.InvNTTLvl(levelQ, levelP, shareOutQP0, shareOutQP0)
	ringQP.InvNTTLvl(levelQ, levelP, shareOutQP1, shareOutQP1)

	// h_0 = u_i * pk_0
	pcks.gaussianSampler.ReadLvl(levelQ, pcks.tmpQP.Q)
	if ringP != nil {
		ringQP.ExtendBasisSmallNormAndCenter(pcks.tmpQP.Q, levelP, nil, pcks.tmpQP.P)
	}

	ringQP.AddLvl(levelQ, levelP, shareOutQP0, pcks.tmpQP, shareOutQP0)

	// h_1 = u_i * pk_1 + e1
	pcks.gaussianSampler.ReadLvl(levelQ, pcks.tmpQP.Q)
	if ringP != nil {
		ringQP.ExtendBasisSmallNormAndCenter(pcks.tmpQP.Q, levelP, nil, pcks.tmpQP.P)
	}

	ringQP.AddLvl(levelQ, levelP, shareOutQP1, pcks.tmpQP, shareOutQP1)

	if ringP != nil {
		// h_0 = (u_i * pk_0 + e0)/P
		pcks.basisExtender.ModDownQPtoQ(levelQ, levelP, shareOutQP0.Q, shareOutQP0.P, shareOutQP0.Q)

		// h_1 = (u_i * pk_1 + e1)/P
		pcks.basisExtender.ModDownQPtoQ(levelQ, levelP, shareOutQP1.Q, shareOutQP1.P, shareOutQP1.Q)
	}

	// h_0 = s_i*c_1 + (u_i * pk_0 + e0)/P
	if ct1.IsNTT {
		ringQ.NTTLvl(levelQ, shareOut.Value[0], shareOut.Value[0])
		ringQ.NTTLvl(levelQ, shareOut.Value[1], shareOut.Value[1])
		ringQ.MulCoeffsMontgomeryAndAddLvl(levelQ, ct1, sk.Value.Q, shareOut.Value[0])
	} else {
		// tmp = s_i*c_1
		ringQ.NTTLazyLvl(levelQ, ct1, pcks.tmpQP.Q)
		ringQ.MulCoeffsMontgomeryConstantLvl(levelQ, pcks.tmpQP.Q, sk.Value.Q, pcks.tmpQP.Q)
		ringQ.InvNTTLvl(levelQ, pcks.tmpQP.Q, pcks.tmpQP.Q)

		// h_0 = s_i*c_1 + (u_i * pk_0 + e0)/P
		ringQ.AddLvl(levelQ, shareOut.Value[0], pcks.tmpQP.Q, shareOut.Value[0])
	}
}

// AggregateShare is the second part of the first and unique round of the PCKSProtocol protocol. Each party uppon receiving the j-1 elements from the
// other parties computes :
//
// [ctx[0] + sum(s_i * ctx[0] + u_i * pk[0] + e_0i), sum(u_i * pk[1] + e_1i)]
func (pcks *PCKSProtocol) AggregateShare(share1, share2, shareOut *PCKSShare) {
	levelQ1, levelQ2 := len(share1.Value[0].Coeffs)-1, len(share2.Value[1].Coeffs)-1
	if levelQ1 != levelQ2 {
		panic("cannot aggreate two shares at different levelQs.")
	}
	pcks.params.RingQ().AddLvl(levelQ1, share1.Value[0], share2.Value[0], shareOut.Value[0])
	pcks.params.RingQ().AddLvl(levelQ1, share1.Value[1], share2.Value[1], shareOut.Value[1])

}

// KeySwitch performs the actual keyswitching operation on a ciphertext ct and put the result in ctOut
func (pcks *PCKSProtocol) KeySwitch(ctIn *rlwe.Ciphertext, combined *PCKSShare, ctOut *rlwe.Ciphertext) {
	level := utils.MinInt(ctIn.Level(), ctOut.Level())
	pcks.params.RingQ().AddLvl(level, ctIn.Value[0], combined.Value[0], ctOut.Value[0])
	ring.CopyValuesLvl(level, combined.Value[1], ctOut.Value[1])
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
