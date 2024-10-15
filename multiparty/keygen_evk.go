package multiparty

import (
	"fmt"
	"io"
	"slices"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/ring/ringqp"
	"github.com/tuneinsight/lattigo/v6/utils"
	"github.com/tuneinsight/lattigo/v6/utils/sampling"
	"github.com/tuneinsight/lattigo/v6/utils/structs"
)

// EvaluationKeyGenProtocol is the structure storing the parameters for the collective EvaluationKey generation.
type EvaluationKeyGenProtocol struct {
	params           rlwe.Parameters
	buff             [2]ringqp.Poly
	gaussianSamplerQ ring.Sampler
}

// ShallowCopy creates a shallow copy of [EvaluationKeyGenProtocol] in which all the read-only data-structures are
// shared with the receiver and the temporary buffers are reallocated. The receiver and the returned
// [EvaluationKeyGenProtocol] can be used concurrently.
func (evkg EvaluationKeyGenProtocol) ShallowCopy() EvaluationKeyGenProtocol {
	prng, err := sampling.NewPRNG()

	// Sanity check, this error should not happen.
	if err != nil {
		panic(err)
	}

	params := evkg.params

	Xe, err := ring.NewSampler(prng, evkg.params.RingQ(), evkg.params.Xe(), false)

	// Sanity check, this error should not happen.
	if err != nil {
		panic(err)
	}

	return EvaluationKeyGenProtocol{
		params:           evkg.params,
		buff:             [2]ringqp.Poly{params.RingQP().NewPoly(), params.RingQP().NewPoly()},
		gaussianSamplerQ: Xe,
	}
}

// NewEvaluationKeyGenProtocol creates a [EvaluationKeyGenProtocol] instance.
func NewEvaluationKeyGenProtocol(params rlwe.ParameterProvider) (evkg EvaluationKeyGenProtocol) {

	prng, err := sampling.NewPRNG()

	// Sanity check, this error should not happen.
	if err != nil {
		panic(err)
	}

	pRLWE := *params.GetRLWEParameters()

	Xe, err := ring.NewSampler(prng, pRLWE.RingQ(), pRLWE.Xe(), false)

	// Sanity check, this error should not happen.
	if err != nil {
		panic(err)
	}

	return EvaluationKeyGenProtocol{
		params:           pRLWE,
		gaussianSamplerQ: Xe,
		buff:             [2]ringqp.Poly{pRLWE.RingQP().NewPoly(), pRLWE.RingQP().NewPoly()},
	}
}

// AllocateShare allocates a party's share in the EvaluationKey Generation.
func (evkg EvaluationKeyGenProtocol) AllocateShare(evkParams ...rlwe.EvaluationKeyParameters) EvaluationKeyGenShare {
	levelQ, levelP, BaseTwoDecomposition, _ := rlwe.ResolveEvaluationKeyParameters(evkg.params, evkParams)
	return evkg.allocateShare(levelQ, levelP, BaseTwoDecomposition)
}

func (evkg EvaluationKeyGenProtocol) allocateShare(levelQ, levelP, BaseTwoDecomposition int) EvaluationKeyGenShare {
	return EvaluationKeyGenShare{*rlwe.NewGadgetCiphertext(evkg.params, 0, levelQ, levelP, BaseTwoDecomposition)}
}

// SampleCRP samples a common random polynomial to be used in the EvaluationKey Generation from the provided
// common reference string.
func (evkg EvaluationKeyGenProtocol) SampleCRP(crs CRS, evkParams ...rlwe.EvaluationKeyParameters) EvaluationKeyGenCRP {
	levelQ, levelP, BaseTwoDecomposition, _ := rlwe.ResolveEvaluationKeyParameters(evkg.params, evkParams)
	return evkg.sampleCRP(crs, levelQ, levelP, BaseTwoDecomposition)
}

func (evkg EvaluationKeyGenProtocol) sampleCRP(crs CRS, levelQ, levelP, BaseTwoDecomposition int) EvaluationKeyGenCRP {

	params := evkg.params

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

	return EvaluationKeyGenCRP{Value: structs.Matrix[ringqp.Poly](m)}
}

// GenShare generates a party's share in the EvaluationKey Generation.
func (evkg EvaluationKeyGenProtocol) GenShare(skIn, skOut *rlwe.SecretKey, crp EvaluationKeyGenCRP, shareOut *EvaluationKeyGenShare) (err error) {

	levelQ := shareOut.LevelQ()
	levelP := shareOut.LevelP()

	if levelQ > utils.Min(skIn.LevelQ(), skOut.LevelQ()) {
		return fmt.Errorf("cannot GenShare: min(skIn, skOut) LevelQ < shareOut LevelQ")
	}

	if shareOut.LevelP() != levelP {
		return fmt.Errorf("cannot GenShare: min(skIn, skOut) LevelP != shareOut LevelP")
	}

	if shareOut.BaseRNSDecompositionVectorSize() != crp.BaseRNSDecompositionVectorSize() {
		return fmt.Errorf("cannot GenShare: crp.BaseRNSDecompositionVectorSize() != shareOut.BaseRNSDecompositionVectorSize()")
	}

	if !slices.Equal(shareOut.BaseTwoDecompositionVectorSize(), crp.BaseTwoDecompositionVectorSize()) {
		return fmt.Errorf("cannot GenShare: crp.BaseTwoDecompositionVectorSize() != shareOut.BaseTwoDecompositionVectorSize()")
	}

	ringQP := evkg.params.RingQP().AtLevel(levelQ, levelP)
	ringQ := ringQP.RingQ

	var hasModulusP bool

	if levelP > -1 {
		hasModulusP = true
		ringQ.MulScalarBigint(skIn.Value.Q, ringQP.RingP.ModulusAtLevel[levelP], evkg.buff[0].Q)
	} else {
		levelP = 0
		evkg.buff[0].Q.CopyLvl(levelQ, skIn.Value.Q)
	}

	m := shareOut.Value
	c := crp.Value

	N := ringQ.N()

	sampler := evkg.gaussianSamplerQ.AtLevel(levelQ)

	BaseTwoDecompositionVectorSize := shareOut.BaseTwoDecompositionVectorSize()
	BaseRNSDecompositionVectorSize := shareOut.BaseRNSDecompositionVectorSize()

	var index int

	for j := 0; j < slices.Max(BaseTwoDecompositionVectorSize); j++ {

		for i := 0; i < BaseRNSDecompositionVectorSize; i++ {

			if j < BaseTwoDecompositionVectorSize[i] {

				mij := m[i][j][0]

				// e
				sampler.Read(mij.Q)

				if hasModulusP {
					ringQP.ExtendBasisSmallNormAndCenter(mij.Q, levelP, mij.Q, mij.P)
				}

				ringQP.NTTLazy(mij, mij)
				ringQP.MForm(mij, mij)

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
					tmp1 := mij.Q.Coeffs[index]

					for w := 0; w < N; w++ {
						tmp1[w] = ring.CRed(tmp1[w]+tmp0[w], qi)
					}
				}

				// sk_in * (qiBarre*qiStar) * 2^w - a*sk + e
				ringQP.MulCoeffsMontgomeryThenSub(c[i][j], skOut.Value, mij)
			}
		}

		ringQ.MulScalar(evkg.buff[0].Q, 1<<shareOut.BaseTwoDecomposition, evkg.buff[0].Q)
	}

	return
}

// AggregateShares computes share3 = share1 + share2.
func (evkg EvaluationKeyGenProtocol) AggregateShares(share1, share2 EvaluationKeyGenShare, share3 *EvaluationKeyGenShare) (err error) {

	if share1.LevelQ() != share2.LevelQ() || share1.LevelQ() != share3.LevelQ() {
		return fmt.Errorf("cannot AggregateShares: share LevelQ do not match")
	}

	if share1.LevelP() != share2.LevelP() || share1.LevelP() != share3.LevelP() {
		return fmt.Errorf("cannot AggregateShares: share LevelP do not match")
	}

	m1 := share1.Value
	m2 := share2.Value
	m3 := share3.Value

	levelQ := share1.LevelQ()
	levelP := share2.LevelP()

	ringQP := evkg.params.RingQP().AtLevel(levelQ, levelP)

	BaseRNSDecompositionVectorSize := share1.BaseRNSDecompositionVectorSize()
	BaseTwoDecompositionVectorSize := share1.BaseTwoDecompositionVectorSize()

	for i := 0; i < BaseRNSDecompositionVectorSize; i++ {
		for j := 0; j < BaseTwoDecompositionVectorSize[i]; j++ {
			ringQP.Add(m1[i][j][0], m2[i][j][0], m3[i][j][0])
		}
	}

	return
}

// GenEvaluationKey finalizes the EvaluationKey Generation and populates the input Evaluationkey with the computed collective EvaluationKey.
func (evkg EvaluationKeyGenProtocol) GenEvaluationKey(share EvaluationKeyGenShare, crp EvaluationKeyGenCRP, evk *rlwe.EvaluationKey) (err error) {

	if share.LevelQ() != evk.LevelQ() {
		return fmt.Errorf("cannot GenEvaluationKey: share LevelQ != evk LevelQ")
	}

	if share.LevelP() != evk.LevelP() {
		return fmt.Errorf("cannot GenEvaluationKey: share LevelP != evk LevelP")
	}

	m := share.Value
	p := crp.Value

	BaseRNSDecompositionVectorSize := len(m)
	BaseTwoDecompositionVectorSize := len(m[0])
	for i := 0; i < BaseRNSDecompositionVectorSize; i++ {
		for j := 0; j < BaseTwoDecompositionVectorSize; j++ {
			evk.Value[i][j][0].Copy(m[i][j][0])
			evk.Value[i][j][1].Copy(p[i][j])
		}
	}

	return
}

// EvaluationKeyGenCRP is a type for common reference polynomials in the EvaluationKey Generation protocol.
type EvaluationKeyGenCRP struct {
	Value structs.Matrix[ringqp.Poly]
}

// LevelQ returns the level of the ciphertext modulus of the target share.
func (crp EvaluationKeyGenCRP) LevelQ() int {
	return crp.Value[0][0].LevelQ()
}

// LevelP returns the level of the auxiliary switching key modulus of the target share.
func (crp EvaluationKeyGenCRP) LevelP() int {
	return crp.Value[0][0].LevelP()
}

// BaseTwoDecompositionVectorSize returns the number of element in the Power of two decomposition basis for each prime of Q.
func (crp EvaluationKeyGenCRP) BaseTwoDecompositionVectorSize() (base []int) {
	base = make([]int, len(crp.Value))
	for i := range crp.Value {
		base[i] = len(crp.Value[i])
	}
	return
}

// BaseRNSDecompositionVectorSize returns the number of element in the RNS decomposition basis: Ceil(lenQi / lenPi)
func (crp EvaluationKeyGenCRP) BaseRNSDecompositionVectorSize() int {
	return len(crp.Value)
}

// EvaluationKeyGenShare is represent a Party's share in the EvaluationKey Generation protocol.
type EvaluationKeyGenShare struct {
	rlwe.GadgetCiphertext
}

// BinarySize returns the serialized size of the object in bytes.
func (share EvaluationKeyGenShare) BinarySize() int {
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
func (share EvaluationKeyGenShare) WriteTo(w io.Writer) (n int64, err error) {
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
func (share *EvaluationKeyGenShare) ReadFrom(r io.Reader) (n int64, err error) {
	return share.GadgetCiphertext.ReadFrom(r)
}

// MarshalBinary encodes the object into a binary form on a newly allocated slice of bytes.
func (share EvaluationKeyGenShare) MarshalBinary() (p []byte, err error) {
	return share.GadgetCiphertext.MarshalBinary()
}

// UnmarshalBinary decodes a slice of bytes generated by
// [EvaluationKeyGenShare.MarshalBinary] or [EvaluationKeyGenShare.WriteTo] on the object.
func (share *EvaluationKeyGenShare) UnmarshalBinary(p []byte) (err error) {
	return share.GadgetCiphertext.UnmarshalBinary(p)
}
