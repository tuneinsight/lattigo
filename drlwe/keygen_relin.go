package drlwe

import (
	"io"

	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/rlwe"
	"github.com/tuneinsight/lattigo/v4/rlwe/ringqp"
	"github.com/tuneinsight/lattigo/v4/utils/sampling"
	"github.com/tuneinsight/lattigo/v4/utils/structs"
)

// RKGProtocol is the structure storing the parameters and and precomputations for the collective relinearization key generation protocol.
type RKGProtocol struct {
	params rlwe.Parameters

	gaussianSamplerQ ring.Sampler
	ternarySamplerQ  ring.Sampler

	buf [2]*ringqp.Poly
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
		buf:              [2]*ringqp.Poly{params.RingQP().NewPoly(), params.RingQP().NewPoly()},
		gaussianSamplerQ: ekg.params.Xe().NewSampler(prng, ekg.params.RingQ(), false),
		ternarySamplerQ:  ekg.params.Xs().NewSampler(prng, ekg.params.RingQ(), false),
	}
}

// RKGCRP is a type for common reference polynomials in the RKG protocol.
type RKGCRP struct {
	Value structs.Matrix[ringqp.Poly]
}

// NewRKGProtocol creates a new RKG protocol struct.
func NewRKGProtocol(params rlwe.Parameters) *RKGProtocol {
	rkg := new(RKGProtocol)
	rkg.params = params

	var err error
	prng, err := sampling.NewPRNG()
	if err != nil {
		panic(err)
	}

	rkg.gaussianSamplerQ = params.Xe().NewSampler(prng, params.RingQ(), false)
	rkg.ternarySamplerQ = params.Xs().NewSampler(prng, params.RingQ(), false)
	rkg.buf = [2]*ringqp.Poly{params.RingQP().NewPoly(), params.RingQP().NewPoly()}
	return rkg
}

// SampleCRP samples a common random polynomial to be used in the RKG protocol from the provided
// common reference string.
func (ekg *RKGProtocol) SampleCRP(crs CRS) RKGCRP {
	params := ekg.params
	decompRNS := params.DecompRNS(params.MaxLevelQ(), params.MaxLevelP())
	decompPw2 := params.DecompPw2(params.MaxLevelQ(), params.MaxLevelP())

	m := make([][]*ringqp.Poly, decompRNS)
	for i := range m {
		vec := make([]*ringqp.Poly, decompPw2)
		for j := range vec {
			vec[j] = ringqp.NewPoly(params.N(), params.MaxLevelQ(), params.MaxLevelP())
		}
		m[i] = vec
	}

	us := ringqp.NewUniformSampler(crs, *params.RingQP())

	for _, v := range m {
		for _, p := range v {
			us.Read(p)
		}
	}

	return RKGCRP{Value: structs.Matrix[ringqp.Poly](m)}
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
		ringQ.MulScalarBigint(sk.Value.Q, ringQP.RingP.ModulusAtLevel[levelP], ekg.buf[0].Q)
	} else {
		levelP = 0
		ring.CopyLvl(levelQ, sk.Value.Q, ekg.buf[0].Q)
	}

	ringQ.IMForm(ekg.buf[0].Q, ekg.buf[0].Q)

	// u
	ekg.ternarySamplerQ.Read(ephSkOut.Value.Q)
	if hasModulusP {
		ringQP.ExtendBasisSmallNormAndCenter(ephSkOut.Value.Q, levelP, nil, ephSkOut.Value.P)
	}
	ringQP.NTT(&ephSkOut.Value, &ephSkOut.Value)
	ringQP.MForm(&ephSkOut.Value, &ephSkOut.Value)

	c := crp.Value

	RNSDecomp := len(shareOut.Value)
	BITDecomp := len(shareOut.Value[0])

	N := ringQ.N()

	var index int
	for j := 0; j < BITDecomp; j++ {
		for i := 0; i < RNSDecomp; i++ {
			// h = e
			ekg.gaussianSamplerQ.Read(shareOut.Value[i][j].Value[0].Q)

			if hasModulusP {
				ringQP.ExtendBasisSmallNormAndCenter(shareOut.Value[i][j].Value[0].Q, levelP, nil, shareOut.Value[i][j].Value[0].P)
			}

			ringQP.NTT(shareOut.Value[i][j].Value[0], shareOut.Value[i][j].Value[0])

			// h = sk*CrtBaseDecompQi + e
			for k := 0; k < levelP+1; k++ {

				index = i*(levelP+1) + k

				// Handles the case where nb pj does not divides nb qi
				if index >= levelQ+1 {
					break
				}

				qi := ringQ.SubRings[index].Modulus
				skP := ekg.buf[0].Q.Coeffs[index]
				h := shareOut.Value[i][j].Value[0].Q.Coeffs[index]

				for w := 0; w < N; w++ {
					h[w] = ring.CRed(h[w]+skP[w], qi)
				}
			}

			// h = sk*CrtBaseDecompQi + -u*a + e
			ringQP.MulCoeffsMontgomeryThenSub(&ephSkOut.Value, c[i][j], shareOut.Value[i][j].Value[0])

			// Second Element
			// e_2i
			ekg.gaussianSamplerQ.Read(shareOut.Value[i][j].Value[1].Q)

			if hasModulusP {
				ringQP.ExtendBasisSmallNormAndCenter(shareOut.Value[i][j].Value[1].Q, levelP, nil, shareOut.Value[i][j].Value[1].P)
			}

			ringQP.NTT(shareOut.Value[i][j].Value[1], shareOut.Value[i][j].Value[1])
			// s*a + e_2i
			ringQP.MulCoeffsMontgomeryThenAdd(&sk.Value, c[i][j], shareOut.Value[i][j].Value[1])
		}

		ringQ.MulScalar(ekg.buf[0].Q, 1<<ekg.params.Pow2Base(), ekg.buf[0].Q)
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
	ringQP.Sub(&ephSk.Value, &sk.Value, ekg.buf[0])

	RNSDecomp := len(shareOut.Value)
	BITDecomp := len(shareOut.Value[0])

	// Each sample is of the form [-u*a_i + s*w_i + e_i]
	// So for each element of the base decomposition w_i:
	for i := 0; i < RNSDecomp; i++ {
		for j := 0; j < BITDecomp; j++ {

			// Computes [(sum samples)*sk + e_1i, sk*a + e_2i]

			// (AggregateShareRoundTwo samples) * sk
			ringQP.MulCoeffsMontgomeryLazy(round1.Value[i][j].Value[0], &sk.Value, shareOut.Value[i][j].Value[0])

			// (AggregateShareRoundTwo samples) * sk + e_1i
			ekg.gaussianSamplerQ.Read(ekg.buf[1].Q)

			if levelP > -1 {
				ringQP.ExtendBasisSmallNormAndCenter(ekg.buf[1].Q, levelP, nil, ekg.buf[1].P)
			}

			ringQP.NTT(ekg.buf[1], ekg.buf[1])
			ringQP.Add(shareOut.Value[i][j].Value[0], ekg.buf[1], shareOut.Value[i][j].Value[0])

			// second part
			// (u_i - s_i) * (sum [x][s*a_i + e_2i]) + e3i
			ekg.gaussianSamplerQ.Read(shareOut.Value[i][j].Value[1].Q)

			if levelP > -1 {
				ringQP.ExtendBasisSmallNormAndCenter(shareOut.Value[i][j].Value[1].Q, levelP, nil, shareOut.Value[i][j].Value[1].P)
			}

			ringQP.NTT(shareOut.Value[i][j].Value[1], shareOut.Value[i][j].Value[1])
			ringQP.MulCoeffsMontgomeryThenAdd(ekg.buf[0], round1.Value[i][j].Value[1], shareOut.Value[i][j].Value[1])
		}
	}
}

// AggregateShares combines two RKG shares into a single one.
func (ekg *RKGProtocol) AggregateShares(share1, share2, shareOut *RKGShare) {

	levelQ := share1.Value[0][0].LevelQ()
	levelP := share1.Value[0][0].LevelP()

	ringQP := ekg.params.RingQP().AtLevel(levelQ, levelP)

	RNSDecomp := len(shareOut.Value)
	BITDecomp := len(shareOut.Value[0])
	for i := 0; i < RNSDecomp; i++ {
		for j := 0; j < BITDecomp; j++ {
			ringQP.Add(share1.Value[i][j].Value[0], share2.Value[i][j].Value[0], shareOut.Value[i][j].Value[0])
			ringQP.Add(share1.Value[i][j].Value[1], share2.Value[i][j].Value[1], shareOut.Value[i][j].Value[1])
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

	levelQ := round1.Value[0][0].LevelQ()
	levelP := round1.Value[0][0].LevelP()

	ringQP := ekg.params.RingQP().AtLevel(levelQ, levelP)

	RNSDecomp := len(round1.Value)
	BITDecomp := len(round1.Value[0])
	for i := 0; i < RNSDecomp; i++ {
		for j := 0; j < BITDecomp; j++ {
			ringQP.Add(round2.Value[i][j].Value[0], round2.Value[i][j].Value[1], evalKeyOut.Value[i][j].Value[0])
			evalKeyOut.Value[i][j].Value[1].Copy(round1.Value[i][j].Value[1])
			ringQP.MForm(evalKeyOut.Value[i][j].Value[0], evalKeyOut.Value[i][j].Value[0])
			ringQP.MForm(evalKeyOut.Value[i][j].Value[1], evalKeyOut.Value[i][j].Value[1])
		}
	}
}

// RKGShare is a share in the RKG protocol.
type RKGShare struct {
	rlwe.GadgetCiphertext
}

// AllocateShare allocates the share of the EKG protocol.
func (ekg *RKGProtocol) AllocateShare() (ephSk *rlwe.SecretKey, r1 *RKGShare, r2 *RKGShare) {
	params := ekg.params
	ephSk = rlwe.NewSecretKey(params)

	decompRNS := params.DecompRNS(params.MaxLevelQ(), params.MaxLevelP())
	decompPw2 := params.DecompPw2(params.MaxLevelQ(), params.MaxLevelP())

	r1 = &RKGShare{GadgetCiphertext: *rlwe.NewGadgetCiphertext(params, params.MaxLevelQ(), params.MaxLevelP(), decompRNS, decompPw2)}
	r2 = &RKGShare{GadgetCiphertext: *rlwe.NewGadgetCiphertext(params, params.MaxLevelQ(), params.MaxLevelP(), decompRNS, decompPw2)}

	return
}

// BinarySize returns the size in bytes of the object
// when encoded using Encode.
func (share *RKGShare) BinarySize() int {
	return share.GadgetCiphertext.BinarySize()
}

// MarshalBinary encodes the object into a binary form on a newly allocated slice of bytes.
func (share *RKGShare) MarshalBinary() (data []byte, err error) {
	return share.GadgetCiphertext.MarshalBinary()
}

// Encode encodes the object into a binary form on a preallocated slice of bytes
// and returns the number of bytes written.
func (share *RKGShare) Encode(data []byte) (n int, err error) {
	return share.GadgetCiphertext.Encode(data)
}

// WriteTo writes the object on an io.Writer.
// To ensure optimal efficiency and minimal allocations, the user is encouraged
// to provide a struct implementing the interface buffer.Writer, which defines
// a subset of the method of the bufio.Writer.
// If w is not compliant to the buffer.Writer interface, it will be wrapped in
// a new bufio.Writer.
// For additional information, see lattigo/utils/buffer/writer.go.
func (share *RKGShare) WriteTo(w io.Writer) (n int64, err error) {
	return share.GadgetCiphertext.WriteTo(w)
}

// UnmarshalBinary decodes a slice of bytes generated by
// MarshalBinary or WriteTo on the object.
func (share *RKGShare) UnmarshalBinary(data []byte) (err error) {
	return share.GadgetCiphertext.UnmarshalBinary(data)
}

// Decode decodes a slice of bytes generated by Encode
// on the object and returns the number of bytes read.
func (share *RKGShare) Decode(data []byte) (n int, err error) {
	return share.GadgetCiphertext.Decode(data)
}

// ReadFrom reads on the object from an io.Writer.
// To ensure optimal efficiency and minimal allocations, the user is encouraged
// to provide a struct implementing the interface buffer.Reader, which defines
// a subset of the method of the bufio.Reader.
// If r is not compliant to the buffer.Reader interface, it will be wrapped in
// a new bufio.Reader.
// For additional information, see lattigo/utils/buffer/reader.go.
func (share *RKGShare) ReadFrom(r io.Reader) (n int64, err error) {
	return share.GadgetCiphertext.ReadFrom(r)
}
