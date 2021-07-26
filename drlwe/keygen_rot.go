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
	GenShare(sk *rlwe.SecretKey, galEl uint64, crp []*ring.Poly, shareOut *RTGShare)
	Aggregate(share1, share2, shareOut *RTGShare)
	GenRotationKey(share *RTGShare, crp []*ring.Poly, rotKey *rlwe.SwitchingKey)
}

// RTGShare is represent a Party's share in the RTG protocol
type RTGShare struct {
	Value [][2]*ring.Poly
}

// RTGProtocol is the structure storing the parameters for the collective rotation-keys generation.
type RTGProtocol struct { // TODO rename GaloisKeyGen ?
	params           rlwe.Parameters
	ringQ            *ring.Ring
	ringP            *ring.Ring
	tmpPoly0         [2]*ring.Poly
	tmpPoly1         [2]*ring.Poly
	gaussianSamplerQ *ring.GaussianSampler
}

// NewRTGProtocol creates a RTGProtocol instance
func NewRTGProtocol(params rlwe.Parameters) *RTGProtocol {
	rtg := new(RTGProtocol)
	rtg.params = params
	rtg.ringQ = params.RingQ()
	rtg.ringP = params.RingP()

	prng, err := utils.NewPRNG()
	if err != nil {
		panic(err)
	}
	rtg.gaussianSamplerQ = ring.NewGaussianSampler(prng, rtg.ringQ, params.Sigma(), int(6*params.Sigma()))
	rtg.tmpPoly0 = [2]*ring.Poly{rtg.ringQ.NewPoly(), rtg.ringP.NewPoly()}
	rtg.tmpPoly1 = [2]*ring.Poly{rtg.ringQ.NewPoly(), rtg.ringP.NewPoly()}
	return rtg
}

// AllocateShares allocates a party's share in the RTG protocol
func (rtg *RTGProtocol) AllocateShares() (rtgShare *RTGShare) {
	rtgShare = new(RTGShare)
	rtgShare.Value = make([][2]*ring.Poly, rtg.params.Beta())
	for i := range rtgShare.Value {
		rtgShare.Value[i][0] = rtg.ringQ.NewPoly()
		rtgShare.Value[i][1] = rtg.ringP.NewPoly()
	}
	return
}

// GenShare generates a party's share in the RTG protocol
func (rtg *RTGProtocol) GenShare(sk *rlwe.SecretKey, galEl uint64, crp [][2]*ring.Poly, shareOut *RTGShare) {

	ringQ := rtg.ringQ
	ringP := rtg.ringP

	twoN := ringQ.N << 2
	galElInv := ring.ModExp(galEl, int(twoN-1), uint64(twoN))

	ring.PermuteNTT(sk.Value[0], galElInv, rtg.tmpPoly1[0])
	ring.PermuteNTT(sk.Value[1], galElInv, rtg.tmpPoly1[1])

	ringQ.MulScalarBigint(sk.Value[0], ringP.ModulusBigint, rtg.tmpPoly0[0])

	var index int

	for i := 0; i < rtg.params.Beta(); i++ {

		// e
		rtg.gaussianSamplerQ.Read(shareOut.Value[i][0])
		extendBasisSmallNormAndCenter(ringQ, ringP, shareOut.Value[i][0], shareOut.Value[i][1])
		ringQ.NTTLazy(shareOut.Value[i][0], shareOut.Value[i][0])
		ringP.NTTLazy(shareOut.Value[i][1], shareOut.Value[i][1])
		ringQ.MForm(shareOut.Value[i][0], shareOut.Value[i][0])
		ringP.MForm(shareOut.Value[i][1], shareOut.Value[i][1])

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
			tmp0 := rtg.tmpPoly0[0].Coeffs[index]
			tmp1 := shareOut.Value[i][0].Coeffs[index]

			for w := 0; w < ringQ.N; w++ {
				tmp1[w] = ring.CRed(tmp1[w]+tmp0[w], qi)
			}
		}

		// sk_in * (qiBarre*qiStar) * 2^w - a*sk + e
		rtg.ringQ.MulCoeffsMontgomeryAndSub(crp[i][0], rtg.tmpPoly1[0], shareOut.Value[i][0])
		rtg.ringP.MulCoeffsMontgomeryAndSub(crp[i][1], rtg.tmpPoly1[1], shareOut.Value[i][1])
	}

	return
}

// Aggregate aggregates two shares in the Rotation Key Generation protocol
func (rtg *RTGProtocol) Aggregate(share1, share2, shareOut *RTGShare) {
	for i := 0; i < rtg.params.Beta(); i++ {
		rtg.ringQ.Add(share1.Value[i][0], share2.Value[i][0], shareOut.Value[i][0])
		rtg.ringP.Add(share1.Value[i][1], share2.Value[i][1], shareOut.Value[i][1])
	}
}

// GenRotationKey finalizes the RTG protocol and populates the input RotationKey with the computed collective SwitchingKey.
func (rtg *RTGProtocol) GenRotationKey(share *RTGShare, crp [][2]*ring.Poly, rotKey *rlwe.SwitchingKey) {
	for i := 0; i < rtg.params.Beta(); i++ {
		ring.CopyValues(share.Value[i][0], rotKey.Value[i][0][0])
		ring.CopyValues(share.Value[i][1], rotKey.Value[i][0][1])
		ring.CopyValues(crp[i][0], rotKey.Value[i][1][0])
		ring.CopyValues(crp[i][1], rotKey.Value[i][1][1])
	}
}

// MarshalBinary encode the target element on a slice of byte.
func (share *RTGShare) MarshalBinary() (data []byte, err error) {
	data = make([]byte, 1+(share.Value[0][0].GetDataLen(true)+share.Value[0][1].GetDataLen(true))*len(share.Value))
	if len(share.Value) > 0xFF {
		return []byte{}, errors.New("RKGShare : uint8 overflow on length")
	}
	data[0] = uint8(len(share.Value))
	ptr := 1
	var inc int
	for _, val := range share.Value {
		if inc, err = val[0].WriteTo(data[ptr:]); err != nil {
			return []byte{}, err
		}
		ptr += inc

		if inc, err = val[1].WriteTo(data[ptr:]); err != nil {
			return []byte{}, err
		}
		ptr += inc
	}

	return data, nil
}

// UnmarshalBinary decodes a slice of bytes on the target element.
func (share *RTGShare) UnmarshalBinary(data []byte) (err error) {
	share.Value = make([][2]*ring.Poly, data[0])
	ptr := 1
	var inc int
	for i := range share.Value {
		share.Value[i][0] = new(ring.Poly)
		if inc, err = share.Value[i][0].DecodePolyNew(data[ptr:]); err != nil {
			return err
		}
		ptr += inc

		share.Value[i][1] = new(ring.Poly)
		if inc, err = share.Value[i][1].DecodePolyNew(data[ptr:]); err != nil {
			return err
		}
		ptr += inc
	}

	return nil
}
