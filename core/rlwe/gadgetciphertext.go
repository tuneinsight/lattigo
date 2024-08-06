package rlwe

import (
	"bufio"
	"fmt"
	"io"
	"slices"

	"github.com/google/go-cmp/cmp"

	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/ring/ringqp"
	"github.com/tuneinsight/lattigo/v6/utils/buffer"
	"github.com/tuneinsight/lattigo/v6/utils/structs"
)

// GadgetCiphertext is a struct for storing an encrypted
// plaintext times the gadget power matrix.
type GadgetCiphertext struct {
	BaseTwoDecomposition int
	Value                structs.Matrix[VectorQP]
}

// NewGadgetCiphertext returns a new [GadgetCiphertext] key with pre-allocated zero-value.
// [GadgetCiphertext] is always in the NTT domain.
// A [GadgetCiphertext] is created by default at degree 1 with the the maximum levelQ and levelP and with no base 2 decomposition.
// Give the optional GadgetCiphertextParameters struct to create a GadgetCiphertext with at a specific degree, levelQ, levelP and/or base 2 decomposition.
func NewGadgetCiphertext(params ParameterProvider, Degree, LevelQ, LevelP, BaseTwoDecomposition int) *GadgetCiphertext {

	p := params.GetRLWEParameters()

	BaseRNSDecompositionVectorSize := p.BaseRNSDecompositionVectorSize(LevelQ, LevelP)
	BaseTwoDecompositionVectorSize := p.BaseTwoDecompositionVectorSize(LevelQ, LevelP, BaseTwoDecomposition)

	m := make(structs.Matrix[VectorQP], BaseRNSDecompositionVectorSize)
	for i := 0; i < BaseRNSDecompositionVectorSize; i++ {
		m[i] = make([]VectorQP, BaseTwoDecompositionVectorSize[i])
		for j := range m[i] {
			m[i][j] = NewVectorQP(params, Degree+1, LevelQ, LevelP)
		}
	}

	return &GadgetCiphertext{BaseTwoDecomposition: BaseTwoDecomposition, Value: m}
}

// Degree returns the degree of the target ciphertext.
func (ct GadgetCiphertext) Degree() int {
	return len(ct.Value[0][0]) - 1
}

// LevelQ returns the level of the modulus Q of the target ciphertext.
func (ct GadgetCiphertext) LevelQ() int {
	return ct.Value[0][0][0].LevelQ()
}

// LevelP returns the level of the modulus P of the target ciphertext.
func (ct GadgetCiphertext) LevelP() int {
	return ct.Value[0][0][0].LevelP()
}

// BaseRNSDecompositionVectorSize returns the number of element in the RNS decomposition basis.
func (ct GadgetCiphertext) BaseRNSDecompositionVectorSize() int {
	return len(ct.Value)
}

// BaseTwoDecompositionVectorSize returns the number of element in the power of two decomposition basis for each prime of Q.
func (ct GadgetCiphertext) BaseTwoDecompositionVectorSize() (base []int) {
	base = make([]int, len(ct.Value))
	for i := range ct.Value {
		base[i] = len(ct.Value[i])
	}
	return
}

// Equal checks two ciphertexts for equality.
func (ct GadgetCiphertext) Equal(other *GadgetCiphertext) bool {
	return (ct.BaseTwoDecomposition == other.BaseTwoDecomposition) && cmp.Equal(ct.Value, other.Value)
}

// CopyNew creates a deep copy of the receiver ciphertext and returns it.
func (ct GadgetCiphertext) CopyNew() (ctCopy *GadgetCiphertext) {
	return &GadgetCiphertext{BaseTwoDecomposition: ct.BaseTwoDecomposition, Value: ct.Value.CopyNew()}
}

// BinarySize returns the serialized size of the object in bytes.
func (ct GadgetCiphertext) BinarySize() (dataLen int) {
	return 8 + ct.Value.BinarySize()
}

// WriteTo writes the object on an [io.Writer]. It implements the [io.WriterTo]
// interface, and will write exactly object.BinarySize() bytes on w.
//
// Unless w implements the [buffer.Writer] interface (see lattigo/utils/buffer/writer.go),
// it will be wrapped into a [bufio.Writer]. Since this requires allocations, it
// is preferable to pass a [buffer.Writer] directly:
//
//   - When writing multiple times to a [io.Writer], it is preferable to first wrap the
//     io.Writer in a pre-allocated [bufio.Writer].
//   - When writing to a pre-allocated var b []byte, it is preferable to pass
//     buffer.NewBuffer(b) as w (see lattigo/utils/buffer/buffer.go).
func (ct GadgetCiphertext) WriteTo(w io.Writer) (n int64, err error) {

	switch w := w.(type) {
	case buffer.Writer:

		var inc int64

		if inc, err = buffer.WriteAsUint64(w, ct.BaseTwoDecomposition); err != nil {
			return n + inc, err
		}

		n += inc

		inc, err = ct.Value.WriteTo(w)

		return n + inc, err

	default:
		return ct.WriteTo(bufio.NewWriter(w))
	}
}

// ReadFrom reads on the object from an [io.Writer]. It implements the
// [io.ReaderFrom] interface.
//
// Unless r implements the [buffer.Reader] interface (see see lattigo/utils/buffer/reader.go),
// it will be wrapped into a [bufio.Reader]. Since this requires allocation, it
// is preferable to pass a [buffer.Reader] directly:
//
//   - When reading multiple values from a [io.Reader], it is preferable to first
//     first wrap [io.Reader] in a pre-allocated [bufio.Reader].
//   - When reading from a var b []byte, it is preferable to pass a buffer.NewBuffer(b)
//     as w (see lattigo/utils/buffer/buffer.go).
func (ct *GadgetCiphertext) ReadFrom(r io.Reader) (n int64, err error) {
	switch r := r.(type) {
	case buffer.Reader:

		var inc int64

		if inc, err = buffer.ReadAsUint64(r, &ct.BaseTwoDecomposition); err != nil {
			return n + inc, err
		}

		n += inc

		inc, err = ct.Value.ReadFrom(r)

		return n + inc, err

	default:
		return ct.ReadFrom(bufio.NewReader(r))
	}
}

// MarshalBinary encodes the object into a binary form on a newly allocated slice of bytes.
func (ct GadgetCiphertext) MarshalBinary() (data []byte, err error) {
	buf := buffer.NewBufferSize(ct.BinarySize())
	_, err = ct.WriteTo(buf)
	return buf.Bytes(), err
}

// UnmarshalBinary decodes a slice of bytes generated by
// [GadgetCiphertext.MarshalBinary] or [GadgetCiphertext.WriteTo] on the object.
func (ct *GadgetCiphertext) UnmarshalBinary(p []byte) (err error) {
	_, err = ct.ReadFrom(buffer.NewBuffer(p))
	return
}

// AddPolyTimesGadgetVectorToGadgetCiphertext takes a plaintext polynomial and a list of [GadgetCiphertext] and adds the
// plaintext times the RNS and BIT decomposition to the i-th element of the i-th ciphertexts. This method return
// an error if len(cts) > 2.
func AddPolyTimesGadgetVectorToGadgetCiphertext(pt ring.Poly, cts []GadgetCiphertext, ringQP ringqp.Ring, buff ring.Poly) (err error) {

	levelQ := cts[0].LevelQ()
	levelP := cts[0].LevelP()

	ringQ := ringQP.RingQ.AtLevel(levelQ)

	if len(cts) > 2 {
		return fmt.Errorf("cannot AddPolyTimesGadgetVectorToGadgetCiphertext: len(cts) should be <= 2")
	}

	if levelP != -1 {
		ringQ.MulScalarBigint(pt, ringQP.RingP.AtLevel(levelP).Modulus(), buff) // P * pt
	} else {
		levelP = 0
		buff.CopyLvl(levelQ, pt) // 1 * pt
	}

	BaseRNSDecompositionVectorSize := len(cts[0].Value)

	BaseTwoDecompositionVectorSize := make([]int, len(cts[0].Value))
	for i := range BaseTwoDecompositionVectorSize {
		BaseTwoDecompositionVectorSize[i] = len(cts[0].Value[i])
	}

	N := ringQ.N()

	var index int
	for j := 0; j < slices.Max(BaseTwoDecompositionVectorSize); j++ {

		for i := 0; i < BaseRNSDecompositionVectorSize; i++ {

			if j < BaseTwoDecompositionVectorSize[i] {

				// e + (m * P * w^2j) * (q_star * q_tild) mod QP
				//
				// q_prod = prod(q[i*#Pi+j])
				// q_star = Q/qprod
				// q_tild = q_star^-1 mod q_prod
				//
				// Therefore : (pt * P * w^2j) * (q_star * q_tild) = pt*P*w^2j mod q[i*#Pi+j], else 0
				for k := 0; k < levelP+1; k++ {

					index = i*(levelP+1) + k

					// Handle cases where #pj does not divide #qi
					if index >= levelQ+1 {
						break
					}

					qi := ringQ.SubRings[index].Modulus
					p0tmp := buff.Coeffs[index]

					for u, ct := range cts {
						p1tmp := ct.Value[i][j][u].Q.Coeffs[index]
						for w := 0; w < N; w++ {
							p1tmp[w] = ring.CRed(p1tmp[w]+p0tmp[w], qi)
						}
					}

				}
			}
		}

		// w^2j
		ringQ.MulScalar(buff, 1<<cts[0].BaseTwoDecomposition, buff)
	}

	return
}

// GadgetPlaintext stores a plaintext value times the gadget vector.
type GadgetPlaintext struct {
	Value structs.Vector[ring.Poly]
}

// NewGadgetPlaintext creates a new [GadgetPlaintext] from value, which can be either uint64, int64 or *[ring.Poly].
// Plaintext is returned in the NTT and Montgomery domain.
func NewGadgetPlaintext(params Parameters, value interface{}, levelQ, levelP, baseTwoDecomposition int) (pt *GadgetPlaintext, err error) {

	ringQ := params.RingQP().RingQ.AtLevel(levelQ)

	BaseTwoDecompositionVectorSize := slices.Max(params.BaseTwoDecompositionVectorSize(levelQ, levelP, baseTwoDecomposition))

	pt = new(GadgetPlaintext)
	pt.Value = make([]ring.Poly, BaseTwoDecompositionVectorSize)

	switch el := value.(type) {
	case uint64:
		pt.Value[0] = ringQ.NewPoly()
		for i := 0; i < levelQ+1; i++ {
			pt.Value[0].Coeffs[i][0] = el
		}
	case int64:
		pt.Value[0] = ringQ.NewPoly()
		if el < 0 {
			for i := 0; i < levelQ+1; i++ {
				pt.Value[0].Coeffs[i][0] = ringQ.SubRings[i].Modulus - uint64(-el)
			}
		} else {
			for i := 0; i < levelQ+1; i++ {
				pt.Value[0].Coeffs[i][0] = uint64(el)
			}
		}
	case ring.Poly:
		pt.Value[0] = *el.CopyNew()
	default:
		return nil, fmt.Errorf("cannot NewGadgetPlaintext: unsupported type, must be either int64, uint64 or ring.Poly but is %T", el)
	}

	if levelP > -1 {
		ringQ.MulScalarBigint(pt.Value[0], params.RingP().AtLevel(levelP).Modulus(), pt.Value[0])
	}

	ringQ.NTT(pt.Value[0], pt.Value[0])
	ringQ.MForm(pt.Value[0], pt.Value[0])

	for i := 1; i < len(pt.Value); i++ {

		pt.Value[i] = *pt.Value[0].CopyNew()

		for j := 0; j < i; j++ {
			ringQ.MulScalar(pt.Value[i], 1<<baseTwoDecomposition, pt.Value[i])
		}
	}

	return
}
