package drlwe

import (
	"errors"

	"github.com/tuneinsight/lattigo/v3/ring"
	"github.com/tuneinsight/lattigo/v3/rlwe"
	"github.com/tuneinsight/lattigo/v3/rlwe/ringqp"
	"github.com/tuneinsight/lattigo/v3/utils"
)

// RTGShare is represent a Party's share in the RTG protocol.
type RTGShare struct {
	Value [][]ringqp.Poly
}

// RTGCRP is a type for common reference polynomials in the RTG protocol.
type RTGCRP [][]ringqp.Poly

// RTGProtocol is the structure storing the parameters for the collective rotation-keys generation.
type RTGProtocol struct {
	params           rlwe.Parameters
	buff             [2]ringqp.Poly
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
		buff:             [2]ringqp.Poly{params.RingQP().NewPoly(), params.RingQP().NewPoly()},
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
	rtg.buff = [2]ringqp.Poly{params.RingQP().NewPoly(), params.RingQP().NewPoly()}
	return rtg
}

// AllocateShare allocates a party's share in the RTG protocol.
func (rtg *RTGProtocol) AllocateShare() (rtgShare *RTGShare) {
	rtgShare = new(RTGShare)

	params := rtg.params
	decompRNS := rtg.params.DecompRNS(params.QCount()-1, params.PCount()-1)
	decompPw2 := rtg.params.DecompPw2(params.QCount()-1, params.PCount()-1)

	rtgShare.Value = make([][]ringqp.Poly, decompRNS)

	for i := 0; i < decompRNS; i++ {
		rtgShare.Value[i] = make([]ringqp.Poly, decompPw2)
		for j := 0; j < decompPw2; j++ {
			rtgShare.Value[i][j] = rtg.params.RingQP().NewPoly()
		}
	}
	return
}

// SampleCRP samples a common random polynomial to be used in the RTG protocol from the provided
// common reference string.
func (rtg *RTGProtocol) SampleCRP(crs CRS) RTGCRP {

	params := rtg.params
	decompRNS := rtg.params.DecompRNS(params.QCount()-1, params.PCount()-1)
	decompPw2 := rtg.params.DecompPw2(params.QCount()-1, params.PCount()-1)

	crp := make([][]ringqp.Poly, decompRNS)
	us := ringqp.NewUniformSampler(crs, *params.RingQP())
	for i := 0; i < decompRNS; i++ {
		crp[i] = make([]ringqp.Poly, decompPw2)
		for j := 0; j < decompPw2; j++ {
			crp[i][j] = rtg.params.RingQP().NewPoly()
			us.Read(crp[i][j])
		}
	}
	return RTGCRP(crp)
}

// GenShare generates a party's share in the RTG protocol.
func (rtg *RTGProtocol) GenShare(sk *rlwe.SecretKey, galEl uint64, crp RTGCRP, shareOut *RTGShare) {

	ringQ := rtg.params.RingQ()
	ringQP := rtg.params.RingQP()

	levelQ := sk.Value.LevelQ()
	levelP := sk.Value.LevelP()

	hasModulusP := levelP > -1

	galElInv := ring.ModExp(galEl, ringQ.NthRoot-1, ringQ.NthRoot)

	ringQ.PermuteNTT(sk.Value.Q, galElInv, rtg.buff[1].Q)

	if hasModulusP {
		ringQ.PermuteNTT(sk.Value.P, galElInv, rtg.buff[1].P)
		ringQ.MulScalarBigint(sk.Value.Q, ringQP.RingP.ModulusAtLevel[levelP], rtg.buff[0].Q)
	} else {
		levelP = 0
		ring.CopyLvl(levelQ, sk.Value.Q, rtg.buff[0].Q)
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
				tmp0 := rtg.buff[0].Q.Coeffs[index]
				tmp1 := shareOut.Value[i][j].Q.Coeffs[index]

				for w := 0; w < ringQ.N; w++ {
					tmp1[w] = ring.CRed(tmp1[w]+tmp0[w], qi)
				}
			}

			// sk_in * (qiBarre*qiStar) * 2^w - a*sk + e
			ringQP.MulCoeffsMontgomeryAndSubLvl(levelQ, levelP, crp[i][j], rtg.buff[1], shareOut.Value[i][j])
		}

		ringQ.MulScalar(rtg.buff[0].Q, 1<<rtg.params.Pow2Base(), rtg.buff[0].Q)
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
			rotKey.Value[i][j].Value[0].CopyValues(share.Value[i][j])
			rotKey.Value[i][j].Value[1].CopyValues(crp[i][j])
		}
	}
}

// MarshalBinary encode the target element on a slice of byte.
func (share *RTGShare) MarshalBinary() (data []byte, err error) {
	data = make([]byte, 2+share.Value[0][0].GetDataLen64(true)*len(share.Value)*len(share.Value[0]))
	if len(share.Value) > 0xFF {
		return []byte{}, errors.New("RKGShare : uint8 overflow on length")
	}
	data[0] = uint8(len(share.Value))
	data[1] = uint8(len(share.Value[0]))
	ptr := 2
	var inc int
	for i := range share.Value {
		for _, el := range share.Value[i] {
			if inc, err = el.WriteTo64(data[ptr:]); err != nil {
				return []byte{}, err
			}
			ptr += inc
		}
	}

	return data, nil
}

// UnmarshalBinary decodes a slice of bytes on the target element.
func (share *RTGShare) UnmarshalBinary(data []byte) (err error) {
	share.Value = make([][]ringqp.Poly, data[0])
	ptr := 2
	var inc int
	for i := range share.Value {
		share.Value[i] = make([]ringqp.Poly, data[1])
		for j := range share.Value[i] {
			if inc, err = share.Value[i][j].DecodePoly64(data[ptr:]); err != nil {
				return err
			}
			ptr += inc
		}
	}

	return nil
}
