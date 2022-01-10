package drlwe

import (
	"errors"

	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/rlwe"
	"github.com/ldsec/lattigo/v2/utils"
)

// RotationKeyGenerator is an interface for the local operation in the generation of rotation keys
type RotationKeyGenerator interface {
	AllocateShares() (rtgShare *RTGShare)
	GenShare(sk *rlwe.SecretKey, galEl uint64, crp RTGCRP, shareOut *RTGShare)
	Aggregate(share1, share2, shareOut *RTGShare)
	GenRotationKey(share *RTGShare, crp RTGCRP, rotKey *rlwe.SwitchingKey)
}

// RTGShare is represent a Party's share in the RTG protocol
type RTGShare struct {
	Value []rlwe.PolyQP
}

// RTGCRP is a type for common reference polynomials in the RTG protocol.
type RTGCRP []rlwe.PolyQP

// RTGProtocol is the structure storing the parameters for the collective rotation-keys generation.
type RTGProtocol struct {
	params           rlwe.Parameters
	tmpPoly0         rlwe.PolyQP
	tmpPoly1         rlwe.PolyQP
	gaussianSamplerQ *ring.GaussianSampler
}

// ShallowCopy creates a shallow copy of RTGProtocol in which all the read-only data-structures are
// shared with the receiver and the temporary buffers are reallocated. The receiver and the returned
// RTGProtocol can be used concurrently.
func (rtg *RTGProtocol) ShallowCopy() *RTGProtocol {
	prng, err := utils.NewPRNG()
	if err != nil {
		panic(err)
	}

	params := rtg.params

	return &RTGProtocol{
		params:           rtg.params,
		tmpPoly0:         params.RingQP().NewPoly(),
		tmpPoly1:         params.RingQP().NewPoly(),
		gaussianSamplerQ: ring.NewGaussianSampler(prng, params.RingQ(), params.Sigma(), int(6*params.Sigma())),
	}
}

// NewRTGProtocol creates a RTGProtocol instance
func NewRTGProtocol(params rlwe.Parameters) *RTGProtocol {
	rtg := new(RTGProtocol)
	rtg.params = params

	prng, err := utils.NewPRNG()
	if err != nil {
		panic(err)
	}
	rtg.gaussianSamplerQ = ring.NewGaussianSampler(prng, params.RingQ(), params.Sigma(), int(6*params.Sigma()))
	rtg.tmpPoly0 = params.RingQP().NewPoly()
	rtg.tmpPoly1 = params.RingQP().NewPoly()
	return rtg
}

// AllocateShares allocates a party's share in the RTG protocol
func (rtg *RTGProtocol) AllocateShares() (rtgShare *RTGShare) {
	rtgShare = new(RTGShare)
	rtgShare.Value = make([]rlwe.PolyQP, rtg.params.Beta())
	for i := range rtgShare.Value {
		rtgShare.Value[i] = rtg.params.RingQP().NewPoly()
	}
	return
}

// SampleCRP samples a common random polynomial to be used in the RTG protocol from the provided
// common reference string.
func (rtg *RTGProtocol) SampleCRP(crs CRS) RTGCRP {
	crp := make([]rlwe.PolyQP, rtg.params.Beta())
	us := rlwe.NewUniformSamplerQP(rtg.params, crs, rtg.params.RingQP())
	for i := range crp {
		crp[i] = rtg.params.RingQP().NewPoly()
		us.Read(&crp[i])
	}
	return RTGCRP(crp)
}

// GenShare generates a party's share in the RTG protocol
func (rtg *RTGProtocol) GenShare(sk *rlwe.SecretKey, galEl uint64, crp RTGCRP, shareOut *RTGShare) {

	ringQ := rtg.params.RingQ()
	ringP := rtg.params.RingP()
	ringQP := rtg.params.RingQP()
	levelQ := rtg.params.QCount() - 1
	levelP := rtg.params.PCount() - 1

	galElInv := ring.ModExp(galEl, ringQ.NthRoot-1, ringQ.NthRoot)

	ringQ.PermuteNTT(sk.Value.Q, galElInv, rtg.tmpPoly1.Q)
	ringP.PermuteNTT(sk.Value.P, galElInv, rtg.tmpPoly1.P)

	ringQ.MulScalarBigint(sk.Value.Q, ringP.ModulusBigint, rtg.tmpPoly0.Q)

	var index int

	for i := 0; i < rtg.params.Beta(); i++ {

		// e
		rtg.gaussianSamplerQ.Read(shareOut.Value[i].Q)
		ringQP.ExtendBasisSmallNormAndCenter(shareOut.Value[i].Q, levelP, nil, shareOut.Value[i].P)
		ringQP.NTTLazyLvl(levelQ, levelP, shareOut.Value[i], shareOut.Value[i])
		ringQP.MFormLvl(levelQ, levelP, shareOut.Value[i], shareOut.Value[i])

		// a is the CRP

		// e + sk_in * (qiBarre*qiStar) * 2^w
		// (qiBarre*qiStar)%qi = 1, else 0
		for j := 0; j < rtg.params.PCount(); j++ {

			index = i*rtg.params.PCount() + j

			// Handles the case where nb pj does not divides nb qi
			if index >= rtg.params.QCount() {
				break
			}

			qi := ringQ.Modulus[index]
			tmp0 := rtg.tmpPoly0.Q.Coeffs[index]
			tmp1 := shareOut.Value[i].Q.Coeffs[index]

			for w := 0; w < ringQ.N; w++ {
				tmp1[w] = ring.CRed(tmp1[w]+tmp0[w], qi)
			}
		}

		// sk_in * (qiBarre*qiStar) * 2^w - a*sk + e
		ringQP.MulCoeffsMontgomeryAndSubLvl(levelQ, levelP, crp[i], rtg.tmpPoly1, shareOut.Value[i])
	}
}

// Aggregate aggregates two shares in the Rotation Key Generation protocol
func (rtg *RTGProtocol) Aggregate(share1, share2, shareOut *RTGShare) {
	ringQP, levelQ, levelP := rtg.params.RingQP(), rtg.params.QCount()-1, rtg.params.PCount()-1
	for i := 0; i < rtg.params.Beta(); i++ {
		ringQP.AddLvl(levelQ, levelP, share1.Value[i], share2.Value[i], shareOut.Value[i])
	}
}

// GenRotationKey finalizes the RTG protocol and populates the input RotationKey with the computed collective SwitchingKey.
func (rtg *RTGProtocol) GenRotationKey(share *RTGShare, crp RTGCRP, rotKey *rlwe.SwitchingKey) {
	for i := 0; i < rtg.params.Beta(); i++ {
		rotKey.Value[i][0].CopyValues(share.Value[i])
		rotKey.Value[i][1].CopyValues(crp[i])
	}
}

// MarshalBinary encode the target element on a slice of byte.
func (share *RTGShare) MarshalBinary() (data []byte, err error) {
	data = make([]byte, 1+share.Value[0].GetDataLen(true)*len(share.Value))
	if len(share.Value) > 0xFF {
		return []byte{}, errors.New("RKGShare : uint8 overflow on length")
	}
	data[0] = uint8(len(share.Value))
	ptr := 1
	var inc int
	for _, val := range share.Value {
		if inc, err = val.WriteTo(data[ptr:]); err != nil {
			return []byte{}, err
		}
		ptr += inc
	}

	return data, nil
}

// UnmarshalBinary decodes a slice of bytes on the target element.
func (share *RTGShare) UnmarshalBinary(data []byte) (err error) {
	share.Value = make([]rlwe.PolyQP, data[0])
	ptr := 1
	var inc int
	for i := range share.Value {
		if inc, err = share.Value[i].DecodePolyNew(data[ptr:]); err != nil {
			return err
		}
		ptr += inc
	}

	return nil
}
