package drlwe

import (
	"io"

	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/rlwe"
	"github.com/tuneinsight/lattigo/v4/rlwe/ringqp"
	"github.com/tuneinsight/lattigo/v4/utils"
	"github.com/tuneinsight/lattigo/v4/utils/sampling"
	"github.com/tuneinsight/lattigo/v4/utils/structs"
)

// EvaluationKeyGenCRP is a type for common reference polynomials in the EvaluationKey Generation protocol.
type EvaluationKeyGenCRP struct {
	Value structs.Matrix[ringqp.Poly]
}

// EvaluationKeyGenProtocol is the structure storing the parameters for the collective EvaluationKey generation.
type EvaluationKeyGenProtocol struct {
	params           rlwe.Parameters
	buff             [2]ringqp.Poly
	gaussianSamplerQ ring.Sampler
}

// ShallowCopy creates a shallow copy of EvaluationKeyGenProtocol in which all the read-only data-structures are
// shared with the receiver and the temporary buffers are reallocated. The receiver and the returned
// EvaluationKeyGenProtocol can be used concurrently.
func (evkg EvaluationKeyGenProtocol) ShallowCopy() EvaluationKeyGenProtocol {
	prng, err := sampling.NewPRNG()
	if err != nil {
		panic(err)
	}

	params := evkg.params

	return EvaluationKeyGenProtocol{
		params:           evkg.params,
		buff:             [2]ringqp.Poly{params.RingQP().NewPoly(), params.RingQP().NewPoly()},
		gaussianSamplerQ: ring.NewSampler(prng, evkg.params.RingQ(), evkg.params.Xe(), false),
	}
}

// NewEvaluationKeyGenProtocol creates a EvaluationKeyGenProtocol instance.
func NewEvaluationKeyGenProtocol(params rlwe.Parameters) (evkg EvaluationKeyGenProtocol) {

	prng, err := sampling.NewPRNG()
	if err != nil {
		panic(err)
	}

	return EvaluationKeyGenProtocol{
		params:           params,
		gaussianSamplerQ: ring.NewSampler(prng, params.RingQ(), params.Xe(), false),
		buff:             [2]ringqp.Poly{params.RingQP().NewPoly(), params.RingQP().NewPoly()},
	}
}

// AllocateShare allocates a party's share in the EvaluationKey Generation.
func (evkg EvaluationKeyGenProtocol) AllocateShare() EvaluationKeyGenShare {
	params := evkg.params
	decompRNS := params.DecompRNS(params.MaxLevelQ(), params.MaxLevelP())
	decompPw2 := params.DecompPw2(params.MaxLevelQ(), params.MaxLevelP())

	p := make([][]ringqp.Poly, decompRNS)
	for i := range p {
		vec := make([]ringqp.Poly, decompPw2)
		for j := range vec {
			vec[j] = ringqp.NewPoly(params.N(), params.MaxLevelQ(), params.MaxLevelP())
		}
		p[i] = vec
	}

	return EvaluationKeyGenShare{Value: structs.Matrix[ringqp.Poly](p)}
}

// SampleCRP samples a common random polynomial to be used in the EvaluationKey Generation from the provided
// common reference string.
func (evkg EvaluationKeyGenProtocol) SampleCRP(crs CRS) EvaluationKeyGenCRP {

	params := evkg.params
	decompRNS := params.DecompRNS(params.MaxLevelQ(), params.MaxLevelP())
	decompPw2 := params.DecompPw2(params.MaxLevelQ(), params.MaxLevelP())

	m := make([][]ringqp.Poly, decompRNS)
	for i := range m {
		vec := make([]ringqp.Poly, decompPw2)
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

	return EvaluationKeyGenCRP{Value: structs.Matrix[ringqp.Poly](m)}
}

// GenShare generates a party's share in the EvaluationKey Generation.
func (evkg EvaluationKeyGenProtocol) GenShare(skIn, skOut *rlwe.SecretKey, crp EvaluationKeyGenCRP, shareOut *EvaluationKeyGenShare) {

	ringQ := evkg.params.RingQ()
	ringQP := evkg.params.RingQP()

	levelQ := utils.Min(skIn.LevelQ(), skOut.LevelQ())
	levelP := utils.Min(skIn.LevelP(), skOut.LevelP())

	var hasModulusP bool

	if levelP > -1 {
		hasModulusP = true
		ringQ.MulScalarBigint(skIn.Value.Q, ringQP.RingP.ModulusAtLevel[levelP], evkg.buff[0].Q)
	} else {
		levelP = 0
		ring.CopyLvl(levelQ, skIn.Value.Q, evkg.buff[0].Q)
	}

	m := shareOut.Value
	c := crp.Value

	RNSDecomp := len(m)
	BITDecomp := len(m[0])

	N := ringQ.N()

	var index int
	for j := 0; j < BITDecomp; j++ {
		for i := 0; i < RNSDecomp; i++ {

			// e
			evkg.gaussianSamplerQ.Read(m[i][j].Q)

			if hasModulusP {
				ringQP.ExtendBasisSmallNormAndCenter(m[i][j].Q, levelP, m[i][j].Q, m[i][j].P)
			}

			ringQP.NTTLazy(m[i][j], m[i][j])
			ringQP.MForm(m[i][j], m[i][j])

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
				tmp0 := evkg.buff[0].Q.Coeffs[index]
				tmp1 := m[i][j].Q.Coeffs[index]

				for w := 0; w < N; w++ {
					tmp1[w] = ring.CRed(tmp1[w]+tmp0[w], qi)
				}
			}

			// sk_in * (qiBarre*qiStar) * 2^w - a*sk + e
			ringQP.MulCoeffsMontgomeryThenSub(c[i][j], skOut.Value, m[i][j])
		}

		ringQ.MulScalar(evkg.buff[0].Q, 1<<evkg.params.Pow2Base(), evkg.buff[0].Q)
	}
}

// AggregateShares computes share3 = share1 + share2.
func (evkg EvaluationKeyGenProtocol) AggregateShares(share1, share2 EvaluationKeyGenShare, share3 *EvaluationKeyGenShare) {

	m1 := share1.Value
	m2 := share2.Value
	m3 := share3.Value

	levelQ := m1[0][0].Q.Level()
	levelP := m1[0][0].P.Level()

	ringQP := evkg.params.RingQP().AtLevel(levelQ, levelP)

	RNSDecomp := len(m1)
	BITDecomp := len(m1[0])
	for i := 0; i < RNSDecomp; i++ {
		for j := 0; j < BITDecomp; j++ {
			ringQP.Add(m1[i][j], m2[i][j], m3[i][j])
		}
	}
}

// GenEvaluationKey finalizes the EvaluationKey Generation and populates the input Evaluationkey with the computed collective EvaluationKey.
func (evkg EvaluationKeyGenProtocol) GenEvaluationKey(share EvaluationKeyGenShare, crp EvaluationKeyGenCRP, evk *rlwe.EvaluationKey) {

	m := share.Value
	p := crp.Value

	RNSDecomp := len(m)
	BITDecomp := len(m[0])
	for i := 0; i < RNSDecomp; i++ {
		for j := 0; j < BITDecomp; j++ {
			evk.Value[i][j][0].Copy(m[i][j])
			evk.Value[i][j][1].Copy(p[i][j])
		}
	}
}

// EvaluationKeyGenShare is represent a Party's share in the EvaluationKey Generation protocol.
type EvaluationKeyGenShare struct {
	Value structs.Matrix[ringqp.Poly]
}

// LevelQ returns the level of the ciphertext modulus of the target share.
func (share EvaluationKeyGenShare) LevelQ() int {
	return share.Value[0][0].LevelQ()
}

// LevelP returns the level of the auxiliary switching key modulus of the target share.
func (share EvaluationKeyGenShare) LevelP() int {
	return share.Value[0][0].LevelP()
}

// BinarySize returns the serialized size of the object in bytes.
func (share EvaluationKeyGenShare) BinarySize() int {
	return share.Value.BinarySize()
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
func (share EvaluationKeyGenShare) WriteTo(w io.Writer) (n int64, err error) {
	return share.Value.WriteTo(w)
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
func (share *EvaluationKeyGenShare) ReadFrom(r io.Reader) (n int64, err error) {
	return share.Value.ReadFrom(r)
}

// MarshalBinary encodes the object into a binary form on a newly allocated slice of bytes.
func (share EvaluationKeyGenShare) MarshalBinary() (p []byte, err error) {
	return share.Value.MarshalBinary()
}

// UnmarshalBinary decodes a slice of bytes generated by
// MarshalBinary or WriteTo on the object.
func (share *EvaluationKeyGenShare) UnmarshalBinary(p []byte) (err error) {
	return share.Value.UnmarshalBinary(p)
}
