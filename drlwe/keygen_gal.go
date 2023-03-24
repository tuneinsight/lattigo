package drlwe

import (
	"encoding/binary"
	"fmt"

	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/rlwe"
	"github.com/tuneinsight/lattigo/v4/rlwe/ringqp"
	"github.com/tuneinsight/lattigo/v4/utils/sampling"
)

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
	prng, err := sampling.NewPRNG()
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

	prng, err := sampling.NewPRNG()
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

// GKGShare is represent a Party's share in the GaloisKey Generation protocol.
type GKGShare struct {
	GaloisElement uint64
	Value         [][]ringqp.Poly
}

// BinarySize returns the size in bytes that the object once marshalled into a binary form.
func (share *GKGShare) BinarySize() int {
	return 10 + share.Value[0][0].BinarySize()*len(share.Value)*len(share.Value[0])
}

// MarshalBinary encodes the object into a binary form on a newly allocated slice of bytes.
func (share *GKGShare) MarshalBinary() (data []byte, err error) {
	data = make([]byte, share.BinarySize())
	_, err = share.MarshalBinaryInPlace(data)
	return
}

// MarshalBinaryInPlace encodes the object into a binary form on a preallocated slice of bytes
// and returns the number of bytes written.
func (share *GKGShare) MarshalBinaryInPlace(data []byte) (ptr int, err error) {

	if len(share.Value) > 0xFF {
		return ptr, fmt.Errorf("uint8 overflow on length")
	}

	data[ptr] = uint8(len(share.Value))
	ptr++
	data[ptr] = uint8(len(share.Value[0]))
	ptr++

	binary.LittleEndian.PutUint64(data[ptr:ptr+8], share.GaloisElement)
	ptr += 8

	var inc int
	for i := range share.Value {
		for _, el := range share.Value[i] {
			if inc, err = el.Read(data[ptr:]); err != nil {
				return
			}
			ptr += inc
		}
	}

	return
}

// UnmarshalBinary decodes a slice of bytes generated by MarshalBinary
// or MarshalBinaryInPlace on the object.
func (share *GKGShare) UnmarshalBinary(data []byte) (err error) {
	_, err = share.UnmarshalBinaryInPlace(data)
	return
}

// UnmarshalBinaryInPlace decodes a slice of bytes generated by MarshalBinary or
// MarshalBinaryInPlace on the object and returns the number of bytes read.
func (share *GKGShare) UnmarshalBinaryInPlace(data []byte) (ptr int, err error) {

	RNS := int(data[0])
	BIT := int(data[1])

	if share.Value == nil || len(share.Value) != RNS {
		share.Value = make([][]ringqp.Poly, RNS)
	}

	share.GaloisElement = binary.LittleEndian.Uint64(data[2:10])
	ptr = 10
	var inc int
	for i := range share.Value {

		if share.Value[i] == nil {
			share.Value[i] = make([]ringqp.Poly, BIT)
		}

		for j := range share.Value[i] {
			if inc, err = share.Value[i][j].Write(data[ptr:]); err != nil {
				return
			}
			ptr += inc
		}
	}

	return
}
