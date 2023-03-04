package drlwe

import (
	"encoding/binary"
	"errors"
	"fmt"

	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/rlwe"
	"github.com/tuneinsight/lattigo/v4/rlwe/ringqp"
	"github.com/tuneinsight/lattigo/v4/utils"
)

// GKGShare is represent a Party's share in the GaloisKey Generation protocol.
type GKGShare struct {
	GaloisElement uint64
	Value         [][]ringqp.Poly
}

// GKGCRP is a type for common reference polynomials in the GaloisKey Generation protocol.
type GKGCRP [][]ringqp.Poly

// GKGProtocol is the structure storing the parameters for the collective GaloisKeys generation.
type GKGProtocol struct {
	params           rlwe.Parameters
	buff             [2]ringqp.Poly
	gaussianSamplerQ *ring.GaussianSampler
}

// ShallowCopy creates a shallow copy of GKGProtocol in which all the read-only data-structures are
// shared with the receiver and the temporary buffers are reallocated. The receiver and the returned
// GKGProtocol can be used concurrently.
func (gkg *GKGProtocol) ShallowCopy() *GKGProtocol {
	prng, err := utils.NewPRNG()
	if err != nil {
		panic(err)
	}

	params := gkg.params

	return &GKGProtocol{
		params:           gkg.params,
		buff:             [2]ringqp.Poly{params.RingQP().NewPoly(), params.RingQP().NewPoly()},
		gaussianSamplerQ: ring.NewGaussianSampler(prng, params.RingQ(), params.Sigma(), int(6*params.Sigma())),
	}
}

// NewGKGProtocol creates a GKGProtocol instance.
func NewGKGProtocol(params rlwe.Parameters) (gkg *GKGProtocol) {
	gkg = new(GKGProtocol)
	gkg.params = params

	prng, err := utils.NewPRNG()
	if err != nil {
		panic(err)
	}
	gkg.gaussianSamplerQ = ring.NewGaussianSampler(prng, params.RingQ(), params.Sigma(), int(6*params.Sigma()))
	gkg.buff = [2]ringqp.Poly{params.RingQP().NewPoly(), params.RingQP().NewPoly()}
	return
}

// AllocateShare allocates a party's share in the GaloisKey Generation.
func (gkg *GKGProtocol) AllocateShare() (gkgShare *GKGShare) {
	gkgShare = new(GKGShare)

	params := gkg.params
	decompRNS := gkg.params.DecompRNS(params.MaxLevelQ(), params.MaxLevelP())
	decompPw2 := gkg.params.DecompPw2(params.MaxLevelQ(), params.MaxLevelP())

	gkgShare.Value = make([][]ringqp.Poly, decompRNS)

	for i := 0; i < decompRNS; i++ {
		gkgShare.Value[i] = make([]ringqp.Poly, decompPw2)
		for j := 0; j < decompPw2; j++ {
			gkgShare.Value[i][j] = gkg.params.RingQP().NewPoly()
		}
	}
	return
}

// SampleCRP samples a common random polynomial to be used in the GaloisKey Generation from the provided
// common reference string.
func (gkg *GKGProtocol) SampleCRP(crs CRS) GKGCRP {

	params := gkg.params
	decompRNS := gkg.params.DecompRNS(params.MaxLevelQ(), params.MaxLevelP())
	decompPw2 := gkg.params.DecompPw2(params.MaxLevelQ(), params.MaxLevelP())

	crp := make([][]ringqp.Poly, decompRNS)
	us := ringqp.NewUniformSampler(crs, *params.RingQP())
	for i := 0; i < decompRNS; i++ {
		crp[i] = make([]ringqp.Poly, decompPw2)
		for j := 0; j < decompPw2; j++ {
			crp[i][j] = gkg.params.RingQP().NewPoly()
			us.Read(crp[i][j])
		}
	}
	return GKGCRP(crp)
}

// GenShare generates a party's share in the GaloisKey Generation.
func (gkg *GKGProtocol) GenShare(sk *rlwe.SecretKey, galEl uint64, crp GKGCRP, shareOut *GKGShare) {

	ringQ := gkg.params.RingQ()
	ringQP := gkg.params.RingQP()

	levelQ := sk.LevelQ()
	levelP := sk.LevelP()

	galElInv := ring.ModExp(galEl, ringQ.NthRoot()-1, ringQ.NthRoot())

	// Important
	shareOut.GaloisElement = galEl

	ringQ.AutomorphismNTT(sk.Value.Q, galElInv, gkg.buff[1].Q)

	var hasModulusP bool

	if levelP > -1 {
		hasModulusP = true
		gkg.params.RingP().AutomorphismNTT(sk.Value.P, galElInv, gkg.buff[1].P)
		ringQ.MulScalarBigint(sk.Value.Q, ringQP.RingP.ModulusAtLevel[levelP], gkg.buff[0].Q)
	} else {
		levelP = 0
		ring.CopyLvl(levelQ, sk.Value.Q, gkg.buff[0].Q)
	}

	RNSDecomp := len(shareOut.Value)
	BITDecomp := len(shareOut.Value[0])
	N := ringQ.N()

	var index int
	for j := 0; j < BITDecomp; j++ {
		for i := 0; i < RNSDecomp; i++ {

			// e
			gkg.gaussianSamplerQ.Read(shareOut.Value[i][j].Q)

			if hasModulusP {
				ringQP.ExtendBasisSmallNormAndCenter(shareOut.Value[i][j].Q, levelP, nil, shareOut.Value[i][j].P)
			}

			ringQP.NTTLazy(shareOut.Value[i][j], shareOut.Value[i][j])
			ringQP.MForm(shareOut.Value[i][j], shareOut.Value[i][j])

			// a is the CRP

			// e + sk_in * (qiBarre*qiStar) * 2^w
			// (qiBarre*qiStar)%qi = 1, else 0
			for k := 0; k < levelP+1; k++ {

				index = i*(levelP+1) + k

				// Handles the case where nb pj does not divides nb qi
				if index >= levelQ+1 {
					break
				}

				qi := ringQ.SubRings[index].Modulus
				tmp0 := gkg.buff[0].Q.Coeffs[index]
				tmp1 := shareOut.Value[i][j].Q.Coeffs[index]

				for w := 0; w < N; w++ {
					tmp1[w] = ring.CRed(tmp1[w]+tmp0[w], qi)
				}
			}

			// sk_in * (qiBarre*qiStar) * 2^w - a*sk + e
			ringQP.MulCoeffsMontgomeryThenSub(crp[i][j], gkg.buff[1], shareOut.Value[i][j])
		}

		ringQ.MulScalar(gkg.buff[0].Q, 1<<gkg.params.Pow2Base(), gkg.buff[0].Q)
	}
}

// AggregateShares aggregates two share in the Rotation Key Generation protocol.
func (gkg *GKGProtocol) AggregateShares(share1, share2, shareOut *GKGShare) {

	if share1.GaloisElement != share2.GaloisElement {
		panic(fmt.Sprintf("cannot aggregate: GKGShares do not share the same GaloisElement: %d != %d", share1.GaloisElement, share2.GaloisElement))
	}

	shareOut.GaloisElement = share1.GaloisElement

	levelQ := share1.Value[0][0].Q.Level()

	var levelP int
	if share1.Value[0][0].P != nil {
		levelP = share1.Value[0][0].P.Level()
	}

	ringQP := gkg.params.RingQP().AtLevel(levelQ, levelP)

	RNSDecomp := len(shareOut.Value)
	BITDecomp := len(shareOut.Value[0])
	for i := 0; i < RNSDecomp; i++ {
		for j := 0; j < BITDecomp; j++ {
			ringQP.Add(share1.Value[i][j], share2.Value[i][j], shareOut.Value[i][j])
		}
	}
}

// GenGaloisKey finalizes the GaloisKey Generation and populates the input GaloisKey with the computed collective GaloisKey.
func (gkg *GKGProtocol) GenGaloisKey(share *GKGShare, crp GKGCRP, gk *rlwe.GaloisKey) {
	RNSDecomp := len(share.Value)
	BITDecomp := len(share.Value[0])
	for i := 0; i < RNSDecomp; i++ {
		for j := 0; j < BITDecomp; j++ {
			gk.Value[i][j].Value[0].Copy(share.Value[i][j])
			gk.Value[i][j].Value[1].Copy(crp[i][j])
		}
	}

	gk.GaloisElement = share.GaloisElement
	gk.NthRoot = gkg.params.RingQ().NthRoot()
}

// MarshalBinary encode the target element on a slice of byte.
func (share *GKGShare) MarshalBinary() (data []byte, err error) {
	data = make([]byte, 2+8+share.Value[0][0].MarshalBinarySize64()*len(share.Value)*len(share.Value[0]))
	if len(share.Value) > 0xFF {
		return []byte{}, errors.New("RKGShare: uint8 overflow on length")
	}
	data[0] = uint8(len(share.Value))
	data[1] = uint8(len(share.Value[0]))
	binary.LittleEndian.PutUint64(data[2:], share.GaloisElement)
	ptr := 10
	var inc int
	for i := range share.Value {
		for _, el := range share.Value[i] {
			if inc, err = el.Encode64(data[ptr:]); err != nil {
				return []byte{}, err
			}
			ptr += inc
		}
	}

	return data, nil
}

// UnmarshalBinary decodes a slice of bytes on the target element.
func (share *GKGShare) UnmarshalBinary(data []byte) (err error) {
	share.Value = make([][]ringqp.Poly, data[0])
	share.GaloisElement = binary.LittleEndian.Uint64(data[2:])
	ptr := 10
	var inc int
	for i := range share.Value {
		share.Value[i] = make([]ringqp.Poly, data[1])
		for j := range share.Value[i] {
			if inc, err = share.Value[i][j].Decode64(data[ptr:]); err != nil {
				return err
			}
			ptr += inc
		}
	}

	return nil
}
