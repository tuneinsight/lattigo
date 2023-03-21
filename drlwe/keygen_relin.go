package drlwe

import (
	"fmt"

	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/rlwe"
	"github.com/tuneinsight/lattigo/v4/rlwe/ringqp"
	"github.com/tuneinsight/lattigo/v4/utils/sampling"
)

// RKGProtocol is the structure storing the parameters and and precomputations for the collective relinearization key generation protocol.
type RKGProtocol struct {
	params rlwe.Parameters

	gaussianSamplerQ *ring.GaussianSampler
	ternarySamplerQ  *ring.TernarySampler // sampling in Montgomery form

	tmpPoly1 ringqp.Poly
	tmpPoly2 ringqp.Poly
}

// ShallowCopy creates a shallow copy of RKGProtocol in which all the read-only data-structures are
// shared with the receiver and the temporary buffers are reallocated. The receiver and the returned
// RKGProtocol can be used concurrently.
func (ekg *RKGProtocol) ShallowCopy() *RKGProtocol {
	var err error
	prng, err := sampling.NewPRNG()
	if err != nil {
		panic(err)
	}

	params := ekg.params

	return &RKGProtocol{
		params:           ekg.params,
		gaussianSamplerQ: ring.NewGaussianSampler(prng, params.RingQ(), params.Sigma(), int(6*params.Sigma())),
		ternarySamplerQ:  ring.NewTernarySamplerWithHammingWeight(prng, params.RingQ(), params.HammingWeight(), false),
		tmpPoly1:         params.RingQP().NewPoly(),
		tmpPoly2:         params.RingQP().NewPoly(),
	}
}

// RKGCRP is a type for common reference polynomials in the RKG protocol.
type RKGCRP [][]ringqp.Poly

// NewRKGProtocol creates a new RKG protocol struct.
func NewRKGProtocol(params rlwe.Parameters) *RKGProtocol {
	rkg := new(RKGProtocol)
	rkg.params = params

	var err error
	prng, err := sampling.NewPRNG()
	if err != nil {
		panic(err)
	}

	rkg.gaussianSamplerQ = ring.NewGaussianSampler(prng, params.RingQ(), params.Sigma(), int(6*params.Sigma()))
	rkg.ternarySamplerQ = ring.NewTernarySamplerWithHammingWeight(prng, params.RingQ(), params.HammingWeight(), false)
	rkg.tmpPoly1 = params.RingQP().NewPoly()
	rkg.tmpPoly2 = params.RingQP().NewPoly()
	return rkg
}

// AllocateShare allocates the share of the EKG protocol.
func (ekg *RKGProtocol) AllocateShare() (ephSk *rlwe.SecretKey, r1 *RKGShare, r2 *RKGShare) {
	params := ekg.params
	ephSk = rlwe.NewSecretKey(params)
	r1, r2 = new(RKGShare), new(RKGShare)

	decompRNS := params.DecompRNS(params.MaxLevelQ(), params.MaxLevelP())
	decompPw2 := params.DecompPw2(params.MaxLevelQ(), params.MaxLevelP())

	r1.Value = make([][][2]ringqp.Poly, decompRNS)
	r2.Value = make([][][2]ringqp.Poly, decompRNS)

	for i := 0; i < decompRNS; i++ {
		r1.Value[i] = make([][2]ringqp.Poly, decompPw2)
		r2.Value[i] = make([][2]ringqp.Poly, decompPw2)
		for j := 0; j < decompPw2; j++ {
			r1.Value[i][j][0] = ekg.params.RingQP().NewPoly()
			r1.Value[i][j][1] = ekg.params.RingQP().NewPoly()
			r2.Value[i][j][0] = ekg.params.RingQP().NewPoly()
			r2.Value[i][j][1] = ekg.params.RingQP().NewPoly()
		}
	}
	return
}

// SampleCRP samples a common random polynomial to be used in the RKG protocol from the provided
// common reference string.
func (ekg *RKGProtocol) SampleCRP(crs CRS) RKGCRP {
	params := ekg.params
	decompRNS := params.DecompRNS(params.MaxLevelQ(), params.MaxLevelP())
	decompPw2 := params.DecompPw2(params.MaxLevelQ(), params.MaxLevelP())

	crp := make([][]ringqp.Poly, decompRNS)
	us := ringqp.NewUniformSampler(crs, *params.RingQP())
	for i := range crp {
		crp[i] = make([]ringqp.Poly, decompPw2)
		for j := range crp[i] {
			crp[i][j] = params.RingQP().NewPoly()
			us.Read(crp[i][j])
		}
	}
	return RKGCRP(crp)
}

// GenShareRoundOne is the first of three rounds of the RKGProtocol protocol. Each party generates a pseudo encryption of
// its secret share of the key s_i under its ephemeral key u_i : [-u_i*a + s_i*w + e_i] and broadcasts it to the other
// j-1 parties.
//
// round1 = [-u_i * a + s_i * P + e_0i, s_i* a + e_i1]
func (ekg *RKGProtocol) GenShareRoundOne(sk *rlwe.SecretKey, crp RKGCRP, ephSkOut *rlwe.SecretKey, shareOut *RKGShare) {
	// Given a base decomposition w_i (here the CRT decomposition)
	// computes [-u*a_i + P*s_i + e_i, s_i * a + e_i]
	// where a_i = crp_i

	levelQ := sk.LevelQ()
	levelP := sk.LevelP()

	ringQP := ekg.params.RingQP().AtLevel(levelQ, levelP)
	ringQ := ringQP.RingQ

	hasModulusP := levelP > -1

	if hasModulusP {
		// Computes P * sk
		ringQ.MulScalarBigint(sk.Value.Q, ringQP.RingP.ModulusAtLevel[levelP], ekg.tmpPoly1.Q)
	} else {
		levelP = 0
		ring.CopyLvl(levelQ, sk.Value.Q, ekg.tmpPoly1.Q)
	}

	ringQ.IMForm(ekg.tmpPoly1.Q, ekg.tmpPoly1.Q)

	// u
	ekg.ternarySamplerQ.Read(ephSkOut.Value.Q)
	if hasModulusP {
		ringQP.ExtendBasisSmallNormAndCenter(ephSkOut.Value.Q, levelP, nil, ephSkOut.Value.P)
	}
	ringQP.NTT(ephSkOut.Value, ephSkOut.Value)
	ringQP.MForm(ephSkOut.Value, ephSkOut.Value)

	RNSDecomp := len(shareOut.Value)
	BITDecomp := len(shareOut.Value[0])

	N := ringQ.N()

	var index int
	for j := 0; j < BITDecomp; j++ {
		for i := 0; i < RNSDecomp; i++ {
			// h = e
			ekg.gaussianSamplerQ.Read(shareOut.Value[i][j][0].Q)

			if hasModulusP {
				ringQP.ExtendBasisSmallNormAndCenter(shareOut.Value[i][j][0].Q, levelP, nil, shareOut.Value[i][j][0].P)
			}

			ringQP.NTT(shareOut.Value[i][j][0], shareOut.Value[i][j][0])

			// h = sk*CrtBaseDecompQi + e
			for k := 0; k < levelP+1; k++ {

				index = i*(levelP+1) + k

				// Handles the case where nb pj does not divides nb qi
				if index >= levelQ+1 {
					break
				}

				qi := ringQ.SubRings[index].Modulus
				skP := ekg.tmpPoly1.Q.Coeffs[index]
				h := shareOut.Value[i][j][0].Q.Coeffs[index]

				for w := 0; w < N; w++ {
					h[w] = ring.CRed(h[w]+skP[w], qi)
				}
			}

			// h = sk*CrtBaseDecompQi + -u*a + e
			ringQP.MulCoeffsMontgomeryThenSub(ephSkOut.Value, crp[i][j], shareOut.Value[i][j][0])

			// Second Element
			// e_2i
			ekg.gaussianSamplerQ.Read(shareOut.Value[i][j][1].Q)

			if hasModulusP {
				ringQP.ExtendBasisSmallNormAndCenter(shareOut.Value[i][j][1].Q, levelP, nil, shareOut.Value[i][j][1].P)
			}

			ringQP.NTT(shareOut.Value[i][j][1], shareOut.Value[i][j][1])
			// s*a + e_2i
			ringQP.MulCoeffsMontgomeryThenAdd(sk.Value, crp[i][j], shareOut.Value[i][j][1])
		}

		ringQ.MulScalar(ekg.tmpPoly1.Q, 1<<ekg.params.Pow2Base(), ekg.tmpPoly1.Q)
	}
}

// GenShareRoundTwo is the second of three rounds of the RKGProtocol protocol. Upon receiving the j-1 shares, each party computes :
//
// round1 = sum([-u_i * a + s_i * P + e_0i, s_i* a + e_i1])
//
//	= [u * a + s * P + e0, s * a + e1]
//
// round2 = [s_i * round1[0] + e_i2, (u_i - s_i) * round1[1] + e_i3]
//
//	= [s_i * {u * a + s * P + e0} + e_i2, (u_i - s_i) * {s * a + e1} + e_i3]
//
// and broadcasts both values to the other j-1 parties.
func (ekg *RKGProtocol) GenShareRoundTwo(ephSk, sk *rlwe.SecretKey, round1 *RKGShare, shareOut *RKGShare) {

	levelQ := sk.LevelQ()
	levelP := sk.LevelP()

	ringQP := ekg.params.RingQP().AtLevel(levelQ, levelP)

	// (u_i - s_i)
	ringQP.Sub(ephSk.Value, sk.Value, ekg.tmpPoly1)

	RNSDecomp := len(shareOut.Value)
	BITDecomp := len(shareOut.Value[0])

	// Each sample is of the form [-u*a_i + s*w_i + e_i]
	// So for each element of the base decomposition w_i:
	for i := 0; i < RNSDecomp; i++ {
		for j := 0; j < BITDecomp; j++ {

			// Computes [(sum samples)*sk + e_1i, sk*a + e_2i]

			// (AggregateShareRoundTwo samples) * sk
			ringQP.MulCoeffsMontgomeryLazy(round1.Value[i][j][0], sk.Value, shareOut.Value[i][j][0])

			// (AggregateShareRoundTwo samples) * sk + e_1i
			ekg.gaussianSamplerQ.Read(ekg.tmpPoly2.Q)

			if levelP > -1 {
				ringQP.ExtendBasisSmallNormAndCenter(ekg.tmpPoly2.Q, levelP, nil, ekg.tmpPoly2.P)
			}

			ringQP.NTT(ekg.tmpPoly2, ekg.tmpPoly2)
			ringQP.Add(shareOut.Value[i][j][0], ekg.tmpPoly2, shareOut.Value[i][j][0])

			// second part
			// (u_i - s_i) * (sum [x][s*a_i + e_2i]) + e3i
			ekg.gaussianSamplerQ.Read(shareOut.Value[i][j][1].Q)

			if levelP > -1 {
				ringQP.ExtendBasisSmallNormAndCenter(shareOut.Value[i][j][1].Q, levelP, nil, shareOut.Value[i][j][1].P)
			}

			ringQP.NTT(shareOut.Value[i][j][1], shareOut.Value[i][j][1])
			ringQP.MulCoeffsMontgomeryThenAdd(ekg.tmpPoly1, round1.Value[i][j][1], shareOut.Value[i][j][1])
		}
	}
}

// AggregateShares combines two RKG shares into a single one.
func (ekg *RKGProtocol) AggregateShares(share1, share2, shareOut *RKGShare) {

	levelQ := share1.Value[0][0][0].Q.Level()

	var levelP int
	if share1.Value[0][0][0].P != nil {
		levelP = share1.Value[0][0][0].P.Level()
	}

	ringQP := ekg.params.RingQP().AtLevel(levelQ, levelP)

	RNSDecomp := len(shareOut.Value)
	BITDecomp := len(shareOut.Value[0])
	for i := 0; i < RNSDecomp; i++ {
		for j := 0; j < BITDecomp; j++ {
			ringQP.Add(share1.Value[i][j][0], share2.Value[i][j][0], shareOut.Value[i][j][0])
			ringQP.Add(share1.Value[i][j][1], share2.Value[i][j][1], shareOut.Value[i][j][1])
		}
	}
}

// GenRelinearizationKey computes the generated RLK from the public shares and write the result in evalKeyOut.
//
// round1 = [u * a + s * P + e0, s * a + e1]
//
// round2 = sum([s_i * {u * a + s * P + e0} + e_i2, (u_i - s_i) * {s * a + e1} + e_i3])
//
//	= [-sua + P*s^2 + s*e0 + e2, sua + ue1 - s^2a -s*e1 + e3]
//
// [round2[0] + round2[1], round1[1]] = [- s^2a - s*e1 + P*s^2 + s*e0 + u*e1 + e2 + e3, s * a + e1]
//
//	= [s * b + P * s^2 + s*e0 + u*e1 + e2 + e3, b]
func (ekg *RKGProtocol) GenRelinearizationKey(round1 *RKGShare, round2 *RKGShare, evalKeyOut *rlwe.RelinearizationKey) {

	levelQ := round1.Value[0][0][0].Q.Level()

	var levelP int
	if round1.Value[0][0][0].P != nil {
		levelP = round1.Value[0][0][0].P.Level()
	}

	ringQP := ekg.params.RingQP().AtLevel(levelQ, levelP)

	RNSDecomp := len(round1.Value)
	BITDecomp := len(round1.Value[0])
	for i := 0; i < RNSDecomp; i++ {
		for j := 0; j < BITDecomp; j++ {
			ringQP.Add(round2.Value[i][j][0], round2.Value[i][j][1], evalKeyOut.Value[i][j].Value[0])
			evalKeyOut.Value[i][j].Value[1].Copy(round1.Value[i][j][1])
			ringQP.MForm(evalKeyOut.Value[i][j].Value[0], evalKeyOut.Value[i][j].Value[0])
			ringQP.MForm(evalKeyOut.Value[i][j].Value[1], evalKeyOut.Value[i][j].Value[1])
		}
	}
}

// RKGShare is a share in the RKG protocol.
type RKGShare struct {
	Value [][][2]ringqp.Poly
}

// MarshalBinarySize returns the size in bytes that the object once marshalled into a binary form.
func (share *RKGShare) MarshalBinarySize() int {
	return 2 + 2*share.Value[0][0][0].MarshalBinarySize()*len(share.Value)*len(share.Value[0])
}

// MarshalBinary encodes the object into a binary form on a newly allocated slice of bytes.
func (share *RKGShare) MarshalBinary() (data []byte, err error) {
	data = make([]byte, share.MarshalBinarySize())
	_, err = share.MarshalBinaryInPlace(data)
	return
}

// MarshalBinaryInPlace encodes the object into a binary form on a preallocated slice of bytes
// and returns the number of bytes written.
func (share *RKGShare) MarshalBinaryInPlace(data []byte) (ptr int, err error) {

	if len(share.Value) > 0xFF {
		return ptr, fmt.Errorf("uint8 overflow on length")
	}

	if len(share.Value[0]) > 0xFF {
		return ptr, fmt.Errorf("uint8 overflow on length")
	}

	data[ptr] = uint8(len(share.Value))
	ptr++
	data[ptr] = uint8(len(share.Value[0]))
	ptr++

	var inc int
	for i := range share.Value {
		for _, el := range share.Value[i] {

			if inc, err = el[0].Read(data[ptr:]); err != nil {
				return
			}
			ptr += inc

			if inc, err = el[1].Read(data[ptr:]); err != nil {
				return
			}
			ptr += inc
		}
	}

	return
}

// UnmarshalBinary decodes a slice of bytes generated by MarshalBinary
// or MarshalBinaryInPlace on the object.
func (share *RKGShare) UnmarshalBinary(data []byte) (err error) {
	_, err = share.UnmarshalBinaryInPlace(data)
	return
}

// UnmarshalBinaryInPlace decodes a slice of bytes generated by MarshalBinary or
// MarshalBinaryInPlace on the object and returns the number of bytes read.
func (share *RKGShare) UnmarshalBinaryInPlace(data []byte) (ptr int, err error) {

	RNS := int(data[0])
	BIT := int(data[1])

	if share.Value == nil || len(share.Value) != RNS {
		share.Value = make([][][2]ringqp.Poly, RNS)
	}

	ptr = 2
	var inc int
	for i := range share.Value {

		if share.Value[i] == nil || len(share.Value[i]) != BIT {
			share.Value[i] = make([][2]ringqp.Poly, BIT)
		}

		for j := range share.Value[i] {

			if inc, err = share.Value[i][j][0].Write(data[ptr:]); err != nil {
				return
			}
			ptr += inc

			if inc, err = share.Value[i][j][1].Write(data[ptr:]); err != nil {
				return
			}
			ptr += inc
		}
	}

	return
}
