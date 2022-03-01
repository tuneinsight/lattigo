package drlwe

import (
	"errors"
	"math/big"

	"github.com/tuneinsight/lattigo/v3/ring"
	"github.com/tuneinsight/lattigo/v3/rlwe"
	"github.com/tuneinsight/lattigo/v3/utils"
)

// RotationKeyGenerator is an interface for the local operation in the generation of rotation keys.
type RotationKeyGenerator interface {
	AllocateShare(logBase2 int) (rtgShare *RTGShare)
	GenShare(sk *rlwe.SecretKey, galEl uint64, crp RTGCRP, shareOut *RTGShare)
	AggregateShare(share1, share2, shareOut *RTGShare)
	GenRotationKey(share *RTGShare, crp RTGCRP, rotKey *rlwe.SwitchingKey)
}

// RTGShare is represent a Party's share in the RTG protocol.
type RTGShare struct {
	Value [][]rlwe.PolyQP
}

// RTGCRP is a type for common reference polynomials in the RTG protocol.
type RTGCRP [][]rlwe.PolyQP

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

// NewRTGProtocol creates a RTGProtocol instance.
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

// AllocateShare allocates a party's share in the RTG protocol.
func (rtg *RTGProtocol) AllocateShare(logBase2 int) (rtgShare *RTGShare) {
	rtgShare = new(RTGShare)

	decompRNS := rtg.params.DecompRNS()
	decompBIT := rtg.params.DecompBIT(logBase2)

	rtgShare.Value = make([][]rlwe.PolyQP, decompRNS)

	for i := 0; i < decompRNS; i++ {
		rtgShare.Value[i] = make([]rlwe.PolyQP, decompBIT)
		for j := 0; j < decompBIT; j++ {
			rtgShare.Value[i][j] = rtg.params.RingQP().NewPoly()
		}
	}
	return
}

// SampleCRP samples a common random polynomial to be used in the RTG protocol from the provided
// common reference string.
func (rtg *RTGProtocol) SampleCRP(crs CRS, logBase2 int) RTGCRP {

	decompRNS := rtg.params.DecompRNS()
	decompBIT := rtg.params.DecompBIT(logBase2)

	crp := make([][]rlwe.PolyQP, decompRNS)
	us := rlwe.NewUniformSamplerQP(rtg.params, crs)
	for i := 0; i < decompRNS; i++ {
		crp[i] = make([]rlwe.PolyQP, decompBIT)
		for j := 0; j < decompBIT; j++ {
			crp[i][j] = rtg.params.RingQP().NewPoly()
			us.Read(&crp[i][j])
		}
	}
	return RTGCRP(crp)
}

// GenShare generates a party's share in the RTG protocol.
func (rtg *RTGProtocol) GenShare(sk *rlwe.SecretKey, galEl uint64, crp RTGCRP, shareOut *RTGShare) {

	ringQ := rtg.params.RingQ()
	ringP := rtg.params.RingP()
	ringQP := rtg.params.RingQP()
	levelQ := sk.Value.Q.Level()

	galElInv := ring.ModExp(galEl, ringQ.NthRoot-1, ringQ.NthRoot)

	ringQ.PermuteNTT(sk.Value.Q, galElInv, rtg.tmpPoly1.Q)

<<<<<<< 83ae36f5f9908381fe0d957ce0daa4f037d38e6f
	ringQ.MulScalarBigint(sk.Value.Q, ringP.ModulusAtLevel[levelP], rtg.tmpPoly0.Q)
=======
<<<<<<< dev_bfv_poly
	ringQ.MulScalarBigint(sk.Value.Q, ringP.ModulusBigint[levelP], rtg.tmpPoly0.Q)
=======
	hasModulusP := sk.Value.P != nil
>>>>>>> First step for adding bit-decomp
>>>>>>> First step for adding bit-decomp

	var levelP int
	if hasModulusP {

		ringP.PermuteNTT(sk.Value.P, galElInv, rtg.tmpPoly1.P)

		levelP = sk.Value.P.Level()

		var pBigInt *big.Int
		if levelP == rtg.params.PCount()-1 {
			pBigInt = rtg.params.RingP().ModulusBigint
		} else {
			P := rtg.params.RingP().Modulus
			pBigInt = new(big.Int).SetUint64(P[0])
			for i := 1; i < levelP+1; i++ {
				pBigInt.Mul(pBigInt, ring.NewUint(P[i]))
			}
		}

		// Computes P * sk
		ringQ.MulScalarBigintLvl(levelQ, sk.Value.Q, pBigInt, rtg.tmpPoly0.Q)
	} else {
		ring.CopyLvl(levelQ, sk.Value.Q, rtg.tmpPoly0.Q)
	}

	RNSDecomp := len(shareOut.Value)
	BITDecomp := len(shareOut.Value[0])

	var index int
	for j := 0; j < BITDecomp; j++ {
		for i := 0; i < RNSDecomp; i++ {

			// e
			rtg.gaussianSamplerQ.Read(shareOut.Value[i][j].Q)

			if hasModulusP {
				ringQP.ExtendBasisSmallNormAndCenter(shareOut.Value[i][j].Q, levelP, nil, shareOut.Value[i][j].P)
			}

			ringQP.NTTLazyLvl(levelQ, levelP, shareOut.Value[i][j], shareOut.Value[i][j])
			ringQP.MFormLvl(levelQ, levelP, shareOut.Value[i][j], shareOut.Value[i][j])

			// a is the CRP

			// e + sk_in * (qiBarre*qiStar) * 2^w
			// (qiBarre*qiStar)%qi = 1, else 0
			for k := 0; k < levelP+1; k++ {

				index = i*(levelP+1) + k

				// Handles the case where nb pj does not divides nb qi
				if index >= levelQ+1 {
					break
				}

				qi := ringQ.Modulus[index]
				tmp0 := rtg.tmpPoly0.Q.Coeffs[index]
				tmp1 := shareOut.Value[i][j].Q.Coeffs[index]

				for w := 0; w < ringQ.N; w++ {
					tmp1[w] = ring.CRed(tmp1[w]+tmp0[w], qi)
				}
			}

			// sk_in * (qiBarre*qiStar) * 2^w - a*sk + e
			ringQP.MulCoeffsMontgomeryAndSubLvl(levelQ, levelP, crp[i][j], rtg.tmpPoly1, shareOut.Value[i][j])
		}
	}
}

// AggregateShare aggregates two share in the Rotation Key Generation protocol.
func (rtg *RTGProtocol) AggregateShare(share1, share2, shareOut *RTGShare) {
	ringQP := rtg.params.RingQP()
	levelQ := share1.Value[0][0].Q.Level()

	var levelP int
	if share1.Value[0][0].P != nil {
		levelP = share1.Value[0][0].P.Level()
	}

	RNSDecomp := len(shareOut.Value)
	BITDecomp := len(shareOut.Value[0])
	for i := 0; i < RNSDecomp; i++ {
		for j := 0; j < BITDecomp; j++ {
			ringQP.AddLvl(levelQ, levelP, share1.Value[i][j], share2.Value[i][j], shareOut.Value[i][j])
		}
	}
}

// GenRotationKey finalizes the RTG protocol and populates the input RotationKey with the computed collective SwitchingKey.
func (rtg *RTGProtocol) GenRotationKey(share *RTGShare, crp RTGCRP, rotKey *rlwe.SwitchingKey) {
	RNSDecomp := len(share.Value)
	BITDecomp := len(share.Value[0])
	for i := 0; i < RNSDecomp; i++ {
		for j := 0; j < BITDecomp; j++ {
			rotKey.Value[i][j][0].CopyValues(share.Value[i][j])
			rotKey.Value[i][j][1].CopyValues(crp[i][j])
		}
	}
}

// MarshalBinary encode the target element on a slice of byte.
func (share *RTGShare) MarshalBinary() (data []byte, err error) {
	data = make([]byte, 2+share.Value[0][0].GetDataLen(true)*len(share.Value)*len(share.Value[0]))
	if len(share.Value) > 0xFF {
		return []byte{}, errors.New("RKGShare : uint8 overflow on length")
	}
	data[0] = uint8(len(share.Value))
	data[1] = uint8(len(share.Value[0]))
	ptr := 2
	var inc int
	for i := range share.Value {
		for _, el := range share.Value[i] {
			if inc, err = el.WriteTo(data[ptr:]); err != nil {
				return []byte{}, err
			}
			ptr += inc
		}
	}

	return data, nil
}

// UnmarshalBinary decodes a slice of bytes on the target element.
func (share *RTGShare) UnmarshalBinary(data []byte) (err error) {
	share.Value = make([][]rlwe.PolyQP, data[0])
	ptr := 2
	var inc int
	for i := range share.Value {
		share.Value[i] = make([]rlwe.PolyQP, data[1])
		for j := range share.Value[i] {
			if inc, err = share.Value[i][j].DecodePolyNew(data[ptr:]); err != nil {
				return err
			}
			ptr += inc
		}
	}

	return nil
}
