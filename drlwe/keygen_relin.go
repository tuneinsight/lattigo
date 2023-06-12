package drlwe

import (
	"io"

	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/rlwe"
	"github.com/tuneinsight/lattigo/v4/rlwe/ringqp"
	"github.com/tuneinsight/lattigo/v4/utils/sampling"
	"github.com/tuneinsight/lattigo/v4/utils/structs"
)

// RelinKeyGenProtocol is the structure storing the parameters and and precomputations for the collective relinearization key generation protocol.
type RelinKeyGenProtocol struct {
	params rlwe.Parameters

	gaussianSamplerQ ring.Sampler
	ternarySamplerQ  ring.Sampler

	buf [2]*ringqp.Poly
}

// RelinKeyGenShare is a share in the RelinKeyGen protocol.
type RelinKeyGenShare struct {
	rlwe.GadgetCiphertext
}

// RelinKeyGenCRP is a type for common reference polynomials in the RelinKeyGen protocol.
type RelinKeyGenCRP struct {
	Value structs.Matrix[ringqp.Poly]
}

// ShallowCopy creates a shallow copy of RelinKeyGenProtocol in which all the read-only data-structures are
// shared with the receiver and the temporary buffers are reallocated. The receiver and the returned
// RelinKeyGenProtocol can be used concurrently.
func (ekg *RelinKeyGenProtocol) ShallowCopy() *RelinKeyGenProtocol {
	var err error
	prng, err := sampling.NewPRNG()
	if err != nil {
		panic(err)
	}

	params := ekg.params

	return &RelinKeyGenProtocol{
		params:           ekg.params,
		buf:              [2]*ringqp.Poly{params.RingQP().NewPoly(), params.RingQP().NewPoly()},
		gaussianSamplerQ: ring.NewSampler(prng, ekg.params.RingQ(), ekg.params.Xe(), false),
		ternarySamplerQ:  ring.NewSampler(prng, ekg.params.RingQ(), ekg.params.Xs(), false),
	}
}

// NewRelinKeyGenProtocol creates a new RelinKeyGen protocol struct.
func NewRelinKeyGenProtocol(params rlwe.Parameters) *RelinKeyGenProtocol {
	rkg := new(RelinKeyGenProtocol)
	rkg.params = params

	var err error
	prng, err := sampling.NewPRNG()
	if err != nil {
		panic(err)
	}

	rkg.gaussianSamplerQ = ring.NewSampler(prng, params.RingQ(), params.Xe(), false)
	rkg.ternarySamplerQ = ring.NewSampler(prng, params.RingQ(), params.Xs(), false)
	rkg.buf = [2]*ringqp.Poly{params.RingQP().NewPoly(), params.RingQP().NewPoly()}
	return rkg
}

// SampleCRP samples a common random polynomial to be used in the RelinKeyGen protocol from the provided
// common reference string.
func (ekg *RelinKeyGenProtocol) SampleCRP(crs CRS) RelinKeyGenCRP {
	params := ekg.params
	decompRNS := params.DecompRNS(params.MaxLevelQ(), params.MaxLevelP())
	decompPw2 := params.DecompPw2(params.MaxLevelQ(), params.MaxLevelP())

	m := make([][]ringqp.Poly, decompRNS)
	for i := range m {
		vec := make([]ringqp.Poly, decompPw2)
		for j := range vec {
			vec[j] = *ringqp.NewPoly(params.N(), params.MaxLevelQ(), params.MaxLevelP())
		}
		m[i] = vec
	}

	us := ringqp.NewUniformSampler(crs, *params.RingQP())

	for _, v := range m {
		for _, p := range v {
			us.Read(&p)
		}
	}

	return RelinKeyGenCRP{Value: structs.Matrix[ringqp.Poly](m)}
}

// GenShareRoundOne is the first of three rounds of the RelinKeyGenProtocol protocol. Each party generates a pseudo encryption of
// its secret share of the key s_i under its ephemeral key u_i : [-u_i*a + s_i*w + e_i] and broadcasts it to the other
// j-1 parties.
//
// round1 = [-u_i * a + s_i * P + e_0i, s_i* a + e_i1]
func (ekg *RelinKeyGenProtocol) GenShareRoundOne(sk *rlwe.SecretKey, crp RelinKeyGenCRP, ephSkOut *rlwe.SecretKey, shareOut *RelinKeyGenShare) {
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

			ringQP.NTT(&shareOut.Value[i][j].Value[0], &shareOut.Value[i][j].Value[0])

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
			ringQP.MulCoeffsMontgomeryThenSub(&ephSkOut.Value, &c[i][j], &shareOut.Value[i][j].Value[0])

			// Second Element
			// e_2i
			ekg.gaussianSamplerQ.Read(shareOut.Value[i][j].Value[1].Q)

			if hasModulusP {
				ringQP.ExtendBasisSmallNormAndCenter(shareOut.Value[i][j].Value[1].Q, levelP, nil, shareOut.Value[i][j].Value[1].P)
			}

			ringQP.NTT(&shareOut.Value[i][j].Value[1], &shareOut.Value[i][j].Value[1])
			// s*a + e_2i
			ringQP.MulCoeffsMontgomeryThenAdd(&sk.Value, &c[i][j], &shareOut.Value[i][j].Value[1])
		}

		ringQ.MulScalar(ekg.buf[0].Q, 1<<ekg.params.Pow2Base(), ekg.buf[0].Q)
	}
}

// GenShareRoundTwo is the second of three rounds of the RelinKeyGenProtocol protocol. Upon receiving the j-1 shares, each party computes :
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
func (ekg *RelinKeyGenProtocol) GenShareRoundTwo(ephSk, sk *rlwe.SecretKey, round1 *RelinKeyGenShare, shareOut *RelinKeyGenShare) {

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
			ringQP.MulCoeffsMontgomeryLazy(&round1.Value[i][j].Value[0], &sk.Value, &shareOut.Value[i][j].Value[0])

			// (AggregateShareRoundTwo samples) * sk + e_1i
			ekg.gaussianSamplerQ.Read(ekg.buf[1].Q)

			if levelP > -1 {
				ringQP.ExtendBasisSmallNormAndCenter(ekg.buf[1].Q, levelP, nil, ekg.buf[1].P)
			}

			ringQP.NTT(ekg.buf[1], ekg.buf[1])
			ringQP.Add(&shareOut.Value[i][j].Value[0], ekg.buf[1], &shareOut.Value[i][j].Value[0])

			// second part
			// (u_i - s_i) * (sum [x][s*a_i + e_2i]) + e3i
			ekg.gaussianSamplerQ.Read(shareOut.Value[i][j].Value[1].Q)

			if levelP > -1 {
				ringQP.ExtendBasisSmallNormAndCenter(shareOut.Value[i][j].Value[1].Q, levelP, nil, shareOut.Value[i][j].Value[1].P)
			}

			ringQP.NTT(&shareOut.Value[i][j].Value[1], &shareOut.Value[i][j].Value[1])
			ringQP.MulCoeffsMontgomeryThenAdd(ekg.buf[0], &round1.Value[i][j].Value[1], &shareOut.Value[i][j].Value[1])
		}
	}
}

// AggregateShares combines two RelinKeyGen shares into a single one.
func (ekg *RelinKeyGenProtocol) AggregateShares(share1, share2, shareOut *RelinKeyGenShare) {

	levelQ := share1.Value[0][0].LevelQ()
	levelP := share1.Value[0][0].LevelP()

	ringQP := ekg.params.RingQP().AtLevel(levelQ, levelP)

	RNSDecomp := len(shareOut.Value)
	BITDecomp := len(shareOut.Value[0])
	for i := 0; i < RNSDecomp; i++ {
		for j := 0; j < BITDecomp; j++ {
			ringQP.Add(&share1.Value[i][j].Value[0], &share2.Value[i][j].Value[0], &shareOut.Value[i][j].Value[0])
			ringQP.Add(&share1.Value[i][j].Value[1], &share2.Value[i][j].Value[1], &shareOut.Value[i][j].Value[1])
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
func (ekg *RelinKeyGenProtocol) GenRelinearizationKey(round1 *RelinKeyGenShare, round2 *RelinKeyGenShare, evalKeyOut *rlwe.RelinearizationKey) {

	levelQ := round1.Value[0][0].LevelQ()
	levelP := round1.Value[0][0].LevelP()

	ringQP := ekg.params.RingQP().AtLevel(levelQ, levelP)

	RNSDecomp := len(round1.Value)
	BITDecomp := len(round1.Value[0])
	for i := 0; i < RNSDecomp; i++ {
		for j := 0; j < BITDecomp; j++ {
			ringQP.Add(&round2.Value[i][j].Value[0], &round2.Value[i][j].Value[1], &evalKeyOut.Value[i][j].Value[0])
			evalKeyOut.Value[i][j].Value[1].Copy(&round1.Value[i][j].Value[1])
			ringQP.MForm(&evalKeyOut.Value[i][j].Value[0], &evalKeyOut.Value[i][j].Value[0])
			ringQP.MForm(&evalKeyOut.Value[i][j].Value[1], &evalKeyOut.Value[i][j].Value[1])
		}
	}
}

// AllocateShare allocates the share of the EKG protocol.
func (ekg *RelinKeyGenProtocol) AllocateShare() (ephSk *rlwe.SecretKey, r1 *RelinKeyGenShare, r2 *RelinKeyGenShare) {
	params := ekg.params
	ephSk = rlwe.NewSecretKey(params)

	decompRNS := params.DecompRNS(params.MaxLevelQ(), params.MaxLevelP())
	decompPw2 := params.DecompPw2(params.MaxLevelQ(), params.MaxLevelP())

	r1 = &RelinKeyGenShare{GadgetCiphertext: *rlwe.NewGadgetCiphertext(params, params.MaxLevelQ(), params.MaxLevelP(), decompRNS, decompPw2)}
	r2 = &RelinKeyGenShare{GadgetCiphertext: *rlwe.NewGadgetCiphertext(params, params.MaxLevelQ(), params.MaxLevelP(), decompRNS, decompPw2)}

	return
}

// BinarySize returns the serialized size of the object in bytes.
func (share *RelinKeyGenShare) BinarySize() int {
	return share.GadgetCiphertext.BinarySize()
}

// WriteTo writes the object on an io.Writer. It implements the io.WriterTo
// interface, and will write exactly object.BinarySize() bytes on w.
//
// Unless w implements the buffer.Writer interface (see lattigo/utils/buffer/writer.go),
// it will be wrapped into a bufio.Writer. Since this requires allocations, it
// is preferable to pass a buffer.Writer directly:
//
//   - When writing multiple times to a io.Writer, it is preferable to first wrap the
//     io.Writer in a pre-allocated bufio.Writer.
//   - When writing to a pre-allocated var b []byte, it is preferable to pass
//     buffer.NewBuffer(b) as w (see lattigo/utils/buffer/buffer.go).
func (share *RelinKeyGenShare) WriteTo(w io.Writer) (n int64, err error) {
	return share.GadgetCiphertext.WriteTo(w)
}

// ReadFrom reads on the object from an io.Writer. It implements the
// io.ReaderFrom interface.
//
// Unless r implements the buffer.Reader interface (see see lattigo/utils/buffer/reader.go),
// it will be wrapped into a bufio.Reader. Since this requires allocation, it
// is preferable to pass a buffer.Reader directly:
//
//   - When reading multiple values from a io.Reader, it is preferable to first
//     first wrap io.Reader in a pre-allocated bufio.Reader.
//   - When reading from a var b []byte, it is preferable to pass a buffer.NewBuffer(b)
//     as w (see lattigo/utils/buffer/buffer.go).
func (share *RelinKeyGenShare) ReadFrom(r io.Reader) (n int64, err error) {
	return share.GadgetCiphertext.ReadFrom(r)
}

// MarshalBinary encodes the object into a binary form on a newly allocated slice of bytes.
func (share *RelinKeyGenShare) MarshalBinary() (data []byte, err error) {
	return share.GadgetCiphertext.MarshalBinary()
}

// UnmarshalBinary decodes a slice of bytes generated by
// MarshalBinary or WriteTo on the object.
func (share *RelinKeyGenShare) UnmarshalBinary(data []byte) (err error) {
	return share.GadgetCiphertext.UnmarshalBinary(data)
}
