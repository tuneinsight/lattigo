package multiparty

import (
	"io"
	"slices"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/ring/ringqp"
	"github.com/tuneinsight/lattigo/v6/utils/sampling"
	"github.com/tuneinsight/lattigo/v6/utils/structs"
)

// RelinearizationKeyGenProtocol is the structure storing the parameters and and precomputations for the collective relinearization key generation protocol.
type RelinearizationKeyGenProtocol struct {
	params rlwe.Parameters

	gaussianSamplerQ ring.Sampler
	ternarySamplerQ  ring.Sampler

	buf [2]ringqp.Poly
}

// RelinearizationKeyGenShare is a share in the RelinearizationKeyGen protocol.
type RelinearizationKeyGenShare struct {
	rlwe.GadgetCiphertext
}

// RelinearizationKeyGenCRP is a type for common reference polynomials in the RelinearizationKeyGen protocol.
type RelinearizationKeyGenCRP struct {
	Value structs.Matrix[ringqp.Poly]
}

// ShallowCopy creates a shallow copy of [RelinearizationKeyGenProtocol] in which all the read-only data-structures are
// shared with the receiver and the temporary buffers are reallocated. The receiver and the returned
// [RelinearizationKeyGenProtocol] can be used concurrently.
func (ekg *RelinearizationKeyGenProtocol) ShallowCopy() RelinearizationKeyGenProtocol {
	var err error
	prng, err := sampling.NewPRNG()

	// Sanity check, this error should not happen.
	if err != nil {
		panic(err)
	}

	params := ekg.params

	Xe, err := ring.NewSampler(prng, ekg.params.RingQ(), ekg.params.Xe(), false)

	// Sanity check, this error should not happen.
	if err != nil {
		panic(err)
	}

	Xs, err := ring.NewSampler(prng, ekg.params.RingQ(), ekg.params.Xs(), false)

	// Sanity check, this error should not happen.
	if err != nil {
		panic(err)
	}

	return RelinearizationKeyGenProtocol{
		params:           ekg.params,
		buf:              [2]ringqp.Poly{params.RingQP().NewPoly(), params.RingQP().NewPoly()},
		gaussianSamplerQ: Xe,
		ternarySamplerQ:  Xs,
	}
}

// NewRelinearizationKeyGenProtocol creates a new RelinearizationKeyGen protocol struct.
func NewRelinearizationKeyGenProtocol(params rlwe.ParameterProvider) RelinearizationKeyGenProtocol {
	rkg := RelinearizationKeyGenProtocol{}
	rkg.params = *params.GetRLWEParameters()

	var err error
	prng, err := sampling.NewPRNG()

	// Sanity check, this error should not happen.
	if err != nil {
		panic(err)
	}

	rkg.gaussianSamplerQ, err = ring.NewSampler(prng, rkg.params.RingQ(), rkg.params.Xe(), false)

	// Sanity check, this error should not happen.
	if err != nil {
		panic(err)
	}

	rkg.ternarySamplerQ, err = ring.NewSampler(prng, rkg.params.RingQ(), rkg.params.Xs(), false)

	// Sanity check, this error should not happen.
	if err != nil {
		panic(err)
	}

	rkg.buf = [2]ringqp.Poly{rkg.params.RingQP().NewPoly(), rkg.params.RingQP().NewPoly()}
	return rkg
}

// SampleCRP samples a common random polynomial to be used in the RelinearizationKeyGen protocol from the provided
// common reference string.
func (ekg RelinearizationKeyGenProtocol) SampleCRP(crs CRS, evkParams ...rlwe.EvaluationKeyParameters) RelinearizationKeyGenCRP {
	params := ekg.params

	levelQ, levelP, BaseTwoDecomposition, _ := rlwe.ResolveEvaluationKeyParameters(ekg.params, evkParams)

	BaseRNSDecompositionVectorSize := params.BaseRNSDecompositionVectorSize(levelQ, levelP)
	BaseTwoDecompositionVectorSize := params.BaseTwoDecompositionVectorSize(levelQ, levelP, BaseTwoDecomposition)

	us := ringqp.NewUniformSampler(crs, params.RingQP().AtLevel(levelQ, levelP))

	m := make([][]ringqp.Poly, BaseRNSDecompositionVectorSize)
	for i := range m {
		vec := make([]ringqp.Poly, BaseTwoDecompositionVectorSize[i])
		for j := range vec {
			vec[j] = us.ReadNew()
		}
		m[i] = vec
	}

	return RelinearizationKeyGenCRP{Value: structs.Matrix[ringqp.Poly](m)}
}

// GenShareRoundOne is the first of three rounds of the [RelinearizationKeyGenProtocol] protocol. Each party generates a pseudo encryption of
// its secret share of the key s_i under its ephemeral key u_i : [-u_i*a + s_i*w + e_i] and broadcasts it to the other
// j-1 parties.
//
// round1 = [-u_i * a + s_i * P + e_0i, s_i* a + e_i1]
func (ekg RelinearizationKeyGenProtocol) GenShareRoundOne(sk *rlwe.SecretKey, crp RelinearizationKeyGenCRP, ephSkOut *rlwe.SecretKey, shareOut *RelinearizationKeyGenShare) {
	// Given a base decomposition w_i (here the CRT decomposition)
	// computes [-u*a_i + P*s_i + e_i, s_i * a + e_i]
	// where a_i = crp_i

	levelQ := shareOut.LevelQ()
	levelP := shareOut.LevelP()

	ringQP := ekg.params.RingQP().AtLevel(levelQ, levelP)
	ringQ := ringQP.RingQ

	hasModulusP := levelP > -1

	if hasModulusP {
		// Computes P * sk
		ringQ.MulScalarBigint(sk.Value.Q, ringQP.RingP.ModulusAtLevel[levelP], ekg.buf[0].Q)
	} else {
		levelP = 0
		ekg.buf[0].Q.CopyLvl(levelQ, sk.Value.Q)
	}

	ringQ.IMForm(ekg.buf[0].Q, ekg.buf[0].Q)

	// u
	ekg.ternarySamplerQ.Read(ephSkOut.Value.Q)
	if hasModulusP {
		ringQP.ExtendBasisSmallNormAndCenter(ephSkOut.Value.Q, levelP, ephSkOut.Value.Q, ephSkOut.Value.P)
	}
	ringQP.NTT(ephSkOut.Value, ephSkOut.Value)
	ringQP.MForm(ephSkOut.Value, ephSkOut.Value)

	c := crp.Value

	BaseRNSDecompositionVectorSize := shareOut.BaseRNSDecompositionVectorSize()
	BaseTwoDecompositionVectorSize := shareOut.BaseTwoDecompositionVectorSize()

	N := ringQ.N()

	sampler := ekg.gaussianSamplerQ.AtLevel(levelQ)

	var index int
	for j := 0; j < slices.Max(BaseTwoDecompositionVectorSize); j++ {
		for i := 0; i < BaseRNSDecompositionVectorSize; i++ {

			if j < BaseTwoDecompositionVectorSize[i] {
				// h = e
				sampler.Read(shareOut.Value[i][j][0].Q)

				if hasModulusP {
					ringQP.ExtendBasisSmallNormAndCenter(shareOut.Value[i][j][0].Q, levelP, shareOut.Value[i][j][0].Q, shareOut.Value[i][j][0].P)
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
					skP := ekg.buf[0].Q.Coeffs[index]
					h := shareOut.Value[i][j][0].Q.Coeffs[index]

					for w := 0; w < N; w++ {
						h[w] = ring.CRed(h[w]+skP[w], qi)
					}
				}

				// h = sk*CrtBaseDecompQi + -u*a + e
				ringQP.MulCoeffsMontgomeryThenSub(ephSkOut.Value, c[i][j], shareOut.Value[i][j][0])

				// Second Element
				// e_2i
				sampler.Read(shareOut.Value[i][j][1].Q)

				if hasModulusP {
					ringQP.ExtendBasisSmallNormAndCenter(shareOut.Value[i][j][1].Q, levelP, shareOut.Value[i][j][1].Q, shareOut.Value[i][j][1].P)
				}

				ringQP.NTT(shareOut.Value[i][j][1], shareOut.Value[i][j][1])
				// s*a + e_2i
				ringQP.MulCoeffsMontgomeryThenAdd(sk.Value, c[i][j], shareOut.Value[i][j][1])
			}
		}

		ringQ.MulScalar(ekg.buf[0].Q, 1<<shareOut.BaseTwoDecomposition, ekg.buf[0].Q)
	}
}

// GenShareRoundTwo is the second of three rounds of the [RelinearizationKeyGenProtocol] protocol. Upon receiving the j-1 shares, each party computes :
//
//   - round1 = sum([-u_i * a + s_i * P + e_0i, s_i* a + e_i1]) = [-ua + sP + e0, sa + e1]
//
//   - round2 = [s_i * round1[0] + (u_i - s_i) * round1[1] + e_i2] = [s_i * {-ua + s * P + e0} + (u_i - s_i) * {sa + e1} + e_i2]
//
// and broadcasts both values to the other j-1 parties.
func (ekg RelinearizationKeyGenProtocol) GenShareRoundTwo(ephSk, sk *rlwe.SecretKey, round1 RelinearizationKeyGenShare, shareOut *RelinearizationKeyGenShare) {

	levelQ := shareOut.LevelQ()
	levelP := shareOut.LevelP()
	BaseRNSDecompositionVectorSize := shareOut.BaseRNSDecompositionVectorSize()
	BaseTwoDecompositionVectorSize := shareOut.BaseTwoDecompositionVectorSize()

	ringQP := ekg.params.RingQP().AtLevel(levelQ, levelP)

	// (u_i - s_i)
	ringQP.Sub(ephSk.Value, sk.Value, ekg.buf[0])

	sampler := ekg.gaussianSamplerQ.AtLevel(levelQ)

	// Each sample is of the form [-u*a_i + s*w_i + e_i]
	// So for each element of the base decomposition w_i:
	for i := 0; i < BaseRNSDecompositionVectorSize; i++ {
		for j := 0; j < BaseTwoDecompositionVectorSize[i]; j++ {

			// Computes [(sum samples)*sk + e_1i, sk*a + e_2i]

			// (AggregateShareRoundTwo samples) * sk
			ringQP.MulCoeffsMontgomeryLazy(round1.Value[i][j][0], sk.Value, shareOut.Value[i][j][0])

			// (AggregateShareRoundTwo samples) * sk + e_1i
			sampler.Read(ekg.buf[1].Q)

			if levelP > -1 {
				ringQP.ExtendBasisSmallNormAndCenter(ekg.buf[1].Q, levelP, ekg.buf[1].Q, ekg.buf[1].P)
			}

			ringQP.NTT(ekg.buf[1], ekg.buf[1])
			ringQP.Add(shareOut.Value[i][j][0], ekg.buf[1], shareOut.Value[i][j][0])

			// second part
			// (AggRound1Samples[0])*sk + (u_i - s_i) * (AggRound1Samples[1]) + e_1
			ringQP.MulCoeffsMontgomeryThenAdd(ekg.buf[0], round1.Value[i][j][1], shareOut.Value[i][j][0])
		}
	}
}

// AggregateShares combines two RelinearizationKeyGen shares into a single one.
func (ekg RelinearizationKeyGenProtocol) AggregateShares(share1, share2 RelinearizationKeyGenShare, shareOut *RelinearizationKeyGenShare) {

	levelQ := share1.LevelQ()
	levelP := share1.LevelP()
	BaseRNSDecompositionVectorSize := share1.BaseRNSDecompositionVectorSize()
	BaseTwoDecompositionVectorSize := share1.BaseTwoDecompositionVectorSize()

	ringQP := ekg.params.RingQP().AtLevel(levelQ, levelP)

	for i := 0; i < BaseRNSDecompositionVectorSize; i++ {
		for j := 0; j < BaseTwoDecompositionVectorSize[i]; j++ {
			// deg(round 1 shares) = 1, deg(round 2 shares) = 0
			for k := 0; k <= share1.Degree(); k++ {
				ringQP.Add(share1.Value[i][j][k], share2.Value[i][j][k], shareOut.Value[i][j][k])
			}
		}
	}
}

// GenRelinearizationKey computes the generated RLK from the public shares and write the result in evalKeyOut.
//
//   - round1 = [-ua + sP + e0, sa + e1]
//   - round2 = sum([s_i * {-ua + sP + e0} + (u_i - s_i) * {sa + e1} + e_i2]) = [-sua + Ps^2 + se0 + e2, sua + ue1 - s^2a -se1]
//   - [round2[0] + round2[1], round1[1]] = [-{s^2a + se1} + Ps^2 + {se0 + ue1 + e2}, sa + e1] = [sb + Ps^2 + e, b]
func (ekg RelinearizationKeyGenProtocol) GenRelinearizationKey(round1 RelinearizationKeyGenShare, round2 RelinearizationKeyGenShare, evalKeyOut *rlwe.RelinearizationKey) {

	levelQ := round1.LevelQ()
	levelP := round1.LevelP()
	BaseRNSDecompositionVectorSize := round1.BaseRNSDecompositionVectorSize()
	BaseTwoDecompositionVectorSize := round1.BaseTwoDecompositionVectorSize()

	ringQP := ekg.params.RingQP().AtLevel(levelQ, levelP)

	for i := 0; i < BaseRNSDecompositionVectorSize; i++ {
		for j := 0; j < BaseTwoDecompositionVectorSize[i]; j++ {
			ringQP.MForm(round2.Value[i][j][0], evalKeyOut.Value[i][j][0])
			ringQP.MForm(round1.Value[i][j][1], evalKeyOut.Value[i][j][1])
		}
	}
}

// AllocateShare allocates the share of the EKG protocol.
// To satisfy the correctness of the multi-party protocol, linearization keys shares cannot be allocated in the compressed format.
func (ekg RelinearizationKeyGenProtocol) AllocateShare(evkParams ...rlwe.EvaluationKeyParameters) (ephSk *rlwe.SecretKey, r1 RelinearizationKeyGenShare, r2 RelinearizationKeyGenShare) {
	params := ekg.params
	ephSk = rlwe.NewSecretKey(params)

	levelQ, levelP, BaseTwoDecomposition, _ := rlwe.ResolveEvaluationKeyParameters(ekg.params, evkParams)

	r1 = RelinearizationKeyGenShare{GadgetCiphertext: *rlwe.NewGadgetCiphertext(params, 1, levelQ, levelP, BaseTwoDecomposition)}
	r2 = RelinearizationKeyGenShare{GadgetCiphertext: *rlwe.NewGadgetCiphertext(params, 0, levelQ, levelP, BaseTwoDecomposition)}

	return
}

// BinarySize returns the serialized size of the object in bytes.
func (share RelinearizationKeyGenShare) BinarySize() int {
	return share.GadgetCiphertext.BinarySize()
}

// WriteTo writes the object on an [io.Writer]. It implements the [io.WriterTo]
// interface, and will write exactly object.BinarySize() bytes on w.
//
// Unless w implements the [buffer.Writer] interface (see lattigo/utils/buffer/writer.go),
// it will be wrapped into a [bufio.Writer]. Since this requires allocations, it
// is preferable to pass a [buffer.Writer] directly:
//
//   - When writing multiple times to a io.Writer, it is preferable to first wrap the
//     io.Writer in a pre-allocated [bufio.Writer].
//   - When writing to a pre-allocated var b []byte, it is preferable to pass
//     buffer.NewBuffer(b) as w (see lattigo/utils/buffer/buffer.go).
func (share RelinearizationKeyGenShare) WriteTo(w io.Writer) (n int64, err error) {
	return share.GadgetCiphertext.WriteTo(w)
}

// ReadFrom reads on the object from an [io.Writer]. It implements the
// [io.ReaderFrom] interface.
//
// Unless r implements the [buffer.Reader] interface (see see lattigo/utils/buffer/reader.go),
// it will be wrapped into a [bufio.Reader]. Since this requires allocation, it
// is preferable to pass a [buffer.Reader] directly:
//
//   - When reading multiple values from a io.Reader, it is preferable to first
//     first wrap [io.Reader] in a pre-allocated bufio.Reader.
//   - When reading from a var b []byte, it is preferable to pass a buffer.NewBuffer(b)
//     as w (see lattigo/utils/buffer/buffer.go).
func (share *RelinearizationKeyGenShare) ReadFrom(r io.Reader) (n int64, err error) {
	return share.GadgetCiphertext.ReadFrom(r)
}

// MarshalBinary encodes the object into a binary form on a newly allocated slice of bytes.
func (share RelinearizationKeyGenShare) MarshalBinary() (data []byte, err error) {
	return share.GadgetCiphertext.MarshalBinary()
}

// UnmarshalBinary decodes a slice of bytes generated by
// [RelinearizationKeyGenShare.MarshalBinary] or [RelinearizationKeyGenShare.WriteTo] on the object.
func (share *RelinearizationKeyGenShare) UnmarshalBinary(data []byte) (err error) {
	return share.GadgetCiphertext.UnmarshalBinary(data)
}
