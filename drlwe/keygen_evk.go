package drlwe

import (
	"fmt"
	"io"

	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/rlwe"
	"github.com/tuneinsight/lattigo/v4/rlwe/ringqp"
	"github.com/tuneinsight/lattigo/v4/utils"
	"github.com/tuneinsight/lattigo/v4/utils/sampling"
	"github.com/tuneinsight/lattigo/v4/utils/structs"
)

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
func (evkg EvaluationKeyGenProtocol) AllocateShare(levelQ, levelP, BaseTwoDecomposition int) EvaluationKeyGenShare {
	return EvaluationKeyGenShare{*rlwe.NewGadgetCiphertext(evkg.params, 0, levelQ, levelP, BaseTwoDecomposition)}
}

// SampleCRP samples a common random polynomial to be used in the EvaluationKey Generation from the provided
// common reference string.
func (evkg EvaluationKeyGenProtocol) SampleCRP(crs CRS, levelQ, levelP, BaseTwoDecomposition int) EvaluationKeyGenCRP {

	params := evkg.params
	decompRNS := params.DecompRNS(levelQ, levelP)
	decompPw2 := params.DecompPw2(levelQ, levelP, BaseTwoDecomposition)

	m := make([][]ringqp.Poly, decompRNS)
	for i := range m {
		vec := make([]ringqp.Poly, decompPw2)
		for j := range vec {
			vec[j] = ringqp.NewPoly(params.N(), levelQ, levelP)
		}
		m[i] = vec
	}

	us := ringqp.NewUniformSampler(crs, params.RingQP().AtLevel(levelQ, levelP))

	for _, v := range m {
		for _, p := range v {
			us.Read(p)
		}
	}

	return EvaluationKeyGenCRP{Value: structs.Matrix[ringqp.Poly](m)}
}

// GenShare generates a party's share in the EvaluationKey Generation.
func (evkg EvaluationKeyGenProtocol) GenShare(skIn, skOut *rlwe.SecretKey, crp EvaluationKeyGenCRP, shareOut *EvaluationKeyGenShare) {

	levelQ := shareOut.LevelQ()
	levelP := shareOut.LevelP()

	if levelQ > utils.Min(skIn.LevelQ(), skOut.LevelQ()) {
		panic(fmt.Errorf("cannot GenShare: min(skIn, skOut) LevelQ < shareOut LevelQ"))
	}

	if shareOut.LevelP() != levelP {
		panic(fmt.Errorf("cannot GenShare: min(skIn, skOut) LevelP != shareOut LevelP"))
	}

	if shareOut.DecompRNS() != crp.DecompRNS() {
		panic(fmt.Errorf("cannot GenSahre: crp.DecompRNS() != shareOut.DecompRNS()"))
	}

	if shareOut.DecompPw2() != crp.DecompPw2() {
		panic(fmt.Errorf("cannot GenSahre: crp.DecompPw2() != shareOut.DecompPw2()"))
	}

	ringQP := evkg.params.RingQP().AtLevel(levelQ, levelP)
	ringQ := ringQP.RingQ

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

	N := ringQ.N()

	sampler := evkg.gaussianSamplerQ.AtLevel(levelQ)

	var index int
	for j := 0; j < shareOut.DecompPw2(); j++ {
		for i := 0; i < shareOut.DecompRNS(); i++ {

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

		ringQ.MulScalar(evkg.buff[0].Q, 1<<shareOut.BaseTwoDecomposition, evkg.buff[0].Q)
	}
}

// AggregateShares computes share3 = share1 + share2.
func (evkg EvaluationKeyGenProtocol) AggregateShares(share1, share2 EvaluationKeyGenShare, share3 *EvaluationKeyGenShare) {

	if share1.LevelQ() != share2.LevelQ() || share1.LevelQ() != share3.LevelQ() {
		panic(fmt.Errorf("cannot AggregateShares: share LevelQ do not match"))
	}

	if share1.LevelP() != share2.LevelP() || share1.LevelP() != share3.LevelP() {
		panic(fmt.Errorf("cannot AggregateShares: share LevelP do not match"))
	}

	m1 := share1.Value
	m2 := share2.Value
	m3 := share3.Value

	levelQ := share1.LevelQ()
	levelP := share2.LevelP()

	ringQP := evkg.params.RingQP().AtLevel(levelQ, levelP)

	DecompRNS := share1.DecompRNS()
	DecompPw2 := share1.DecompPw2()

	for i := 0; i < DecompRNS; i++ {
		for j := 0; j < DecompPw2; j++ {
			ringQP.Add(m1[i][j][0], m2[i][j][0], m3[i][j][0])
		}
	}
}

// GenEvaluationKey finalizes the EvaluationKey Generation and populates the input Evaluationkey with the computed collective EvaluationKey.
func (evkg EvaluationKeyGenProtocol) GenEvaluationKey(share EvaluationKeyGenShare, crp EvaluationKeyGenCRP, evk *rlwe.EvaluationKey) {

	if share.LevelQ() != evk.LevelQ() {
		panic(fmt.Errorf("cannot GenEvaluationKey: share LevelQ != evk LevelQ"))
	}

	if share.LevelP() != evk.LevelP() {
		panic(fmt.Errorf("cannot GenEvaluationKey: share LevelP != evk LevelP"))
	}

	m := share.Value
	p := crp.Value

	DecompRNS := len(m)
	DecompPw2 := len(m[0])
	for i := 0; i < DecompRNS; i++ {
		for j := 0; j < DecompPw2; j++ {
			evk.Value[i][j][0].Copy(m[i][j][0])
			evk.Value[i][j][1].Copy(p[i][j])
		}
	}
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

// DecompPw2 returns ceil(p.MaxBitQ(levelQ, levelP)/DecompPw2).
func (crp EvaluationKeyGenCRP) DecompPw2() int {
	return len(crp.Value[0])
}

// DecompRNS returns the number of element in the RNS decomposition basis: Ceil(lenQi / lenPi)
func (crp EvaluationKeyGenCRP) DecompRNS() int {
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
func (share *EvaluationKeyGenShare) ReadFrom(r io.Reader) (n int64, err error) {
	return share.GadgetCiphertext.ReadFrom(r)
}

// MarshalBinary encodes the object into a binary form on a newly allocated slice of bytes.
func (share EvaluationKeyGenShare) MarshalBinary() (p []byte, err error) {
	return share.GadgetCiphertext.MarshalBinary()
}

// UnmarshalBinary decodes a slice of bytes generated by
// MarshalBinary or WriteTo on the object.
func (share *EvaluationKeyGenShare) UnmarshalBinary(p []byte) (err error) {
	return share.GadgetCiphertext.UnmarshalBinary(p)
}
