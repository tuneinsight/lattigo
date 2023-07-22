package rlwe

import (
	"bufio"
	"fmt"
	"io"

	"github.com/google/go-cmp/cmp"
	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/rlwe/ringqp"
	"github.com/tuneinsight/lattigo/v4/utils/buffer"
	"github.com/tuneinsight/lattigo/v4/utils/sampling"
	"github.com/tuneinsight/lattigo/v4/utils/structs"
)

// OperandInterface is a common interface for Ciphertext and Plaintext types.
type OperandInterface[T ring.Poly | ringqp.Poly] interface {
	El() *Operand[T]
	Degree() int
	Level() int
}

type Operand[T ring.Poly | ringqp.Poly] struct {
	*MetaData
	Value structs.Vector[T]
}

func NewOperandQ(params ParametersInterface, degree, levelQ int) *Operand[ring.Poly] {
	ringQ := params.RingQ().AtLevel(levelQ)

	Value := make([]ring.Poly, degree+1)
	for i := range Value {
		Value[i] = ringQ.NewPoly()
	}

	return &Operand[ring.Poly]{
		Value: Value,
		MetaData: &MetaData{
			CiphertextMetaData: CiphertextMetaData{
				IsNTT: params.NTTFlag(),
			},
		},
	}
}

func NewOperandQP(params ParametersInterface, degree, levelQ, levelP int) *Operand[ringqp.Poly] {
	ringQP := params.RingQP().AtLevel(levelQ, levelP)

	Value := make([]ringqp.Poly, degree+1)
	for i := range Value {
		Value[i] = ringQP.NewPoly()
	}

	return &Operand[ringqp.Poly]{
		Value: Value,
		MetaData: &MetaData{
			CiphertextMetaData: CiphertextMetaData{
				IsNTT: params.NTTFlag(),
			},
		},
	}
}

// NewOperandQAtLevelFromPoly constructs a new Operand at a specific level
// where the message is set to the passed poly. No checks are performed on poly and
// the returned Operand will share its backing array of coefficients.
// Returned Operand's MetaData is nil.
func NewOperandQAtLevelFromPoly(level int, poly []ring.Poly) (*Operand[ring.Poly], error) {
	Value := make([]ring.Poly, len(poly))
	for i := range Value {

		if len(poly[i].Coeffs) < level+1 {
			return nil, fmt.Errorf("cannot NewOperandQAtLevelFromPoly: provided ring.Poly[%d] level is too small", i)
		}

		Value[i].Coeffs = poly[i].Coeffs[:level+1]
		Value[i].Buff = poly[i].Buff[:poly[i].N()*(level+1)]
	}

	return &Operand[ring.Poly]{Value: Value}, nil
}

// Equal performs a deep equal.
func (op Operand[T]) Equal(other *Operand[T]) bool {
	return cmp.Equal(op.MetaData, other.MetaData) && cmp.Equal(op.Value, other.Value)
}

// Degree returns the degree of the target Operand.
func (op Operand[T]) Degree() int {
	return len(op.Value) - 1
}

// Level returns the level of the target Operand.
func (op Operand[T]) Level() int {
	return op.LevelQ()
}

func (op Operand[T]) LevelQ() int {
	switch el := any(op.Value[0]).(type) {
	case ring.Poly:
		return el.Level()
	case ringqp.Poly:
		return el.LevelQ()
	default:
		panic("invalid Operand[type]")
	}
}

func (op Operand[T]) LevelP() int {
	switch el := any(op.Value[0]).(type) {
	case ring.Poly:
		panic("cannot levelP on Operand[ring.Poly]")
	case ringqp.Poly:
		return el.LevelP()
	default:
		panic("invalid Operand[type]")
	}
}

func (op *Operand[T]) El() *Operand[T] {
	return op
}

// Resize resizes the degree of the target element.
// Sets the NTT flag of the added poly equal to the NTT flag
// to the poly at degree zero.
func (op *Operand[T]) Resize(degree, level int) {

	switch op := any(op).(type) {
	case *Operand[ring.Poly]:
		if op.Level() != level {
			for i := range op.Value {
				op.Value[i].Resize(level)
			}
		}

		if op.Degree() > degree {
			op.Value = op.Value[:degree+1]
		} else if op.Degree() < degree {

			for op.Degree() < degree {
				op.Value = append(op.Value, []ring.Poly{ring.NewPoly(op.Value[0].N(), level)}...)
			}
		}
	default:
		panic(fmt.Errorf("can only resize Operand[ring.Poly] but is %T", op))
	}
}

// CopyNew creates a deep copy of the object and returns it.
func (op Operand[T]) CopyNew() *Operand[T] {
	return &Operand[T]{Value: *op.Value.CopyNew(), MetaData: op.MetaData.CopyNew()}
}

// Copy copies the input element and its parameters on the target element.
func (op *Operand[T]) Copy(opCopy *Operand[T]) {

	if op != opCopy {
		switch any(op.Value).(type) {
		case structs.Vector[ring.Poly]:

			op0 := any(op.Value).(structs.Vector[ring.Poly])
			op1 := any(opCopy.Value).(structs.Vector[ring.Poly])

			for i := range opCopy.Value {
				op0[i].Copy(op1[i])
			}

		case structs.Vector[ringqp.Poly]:

			op0 := any(op.Value).(structs.Vector[ringqp.Poly])
			op1 := any(opCopy.Value).(structs.Vector[ringqp.Poly])

			for i := range opCopy.Value {
				op0[i].Copy(op1[i])
			}
		}

		*op.MetaData = *opCopy.MetaData
	}
}

// GetSmallestLargest returns the provided element that has the smallest degree as a first
// returned value and the largest degree as second return value. If the degree match, the
// order is the same as for the input.
func GetSmallestLargest[T ring.Poly | ringqp.Poly](el0, el1 *Operand[T]) (smallest, largest *Operand[T], sameDegree bool) {
	switch {
	case el0.Degree() > el1.Degree():
		return el1, el0, false
	case el0.Degree() < el1.Degree():
		return el0, el1, false
	}
	return el0, el1, true
}

// PopulateElementRandom creates a new rlwe.Element with random coefficients.
func PopulateElementRandom(prng sampling.PRNG, params ParametersInterface, ct *Operand[ring.Poly]) {
	sampler := ring.NewUniformSampler(prng, params.RingQ()).AtLevel(ct.Level())
	for i := range ct.Value {
		sampler.Read(ct.Value[i])
	}
}

// SwitchCiphertextRingDegreeNTT changes the ring degree of ctIn to the one of opOut.
// Maps Y^{N/n} -> X^{N} or X^{N} -> Y^{N/n}.
// If the ring degree of opOut is larger than the one of ctIn, then the ringQ of opOut
// must be provided (otherwise, a nil pointer).
// The ctIn must be in the NTT domain and opOut will be in the NTT domain.
func SwitchCiphertextRingDegreeNTT(ctIn *Operand[ring.Poly], ringQLargeDim *ring.Ring, opOut *Operand[ring.Poly]) {

	NIn, NOut := len(ctIn.Value[0].Coeffs[0]), len(opOut.Value[0].Coeffs[0])

	if NIn > NOut {

		gap := NIn / NOut
		buff := make([]uint64, NIn)
		for i := range opOut.Value {
			for j := range opOut.Value[i].Coeffs {

				tmpIn, tmpOut := ctIn.Value[i].Coeffs[j], opOut.Value[i].Coeffs[j]

				ringQLargeDim.SubRings[j].INTT(tmpIn, buff)

				for w0, w1 := 0, 0; w0 < NOut; w0, w1 = w0+1, w1+gap {
					tmpOut[w0] = buff[w1]
				}

				s := ringQLargeDim.SubRings[j]

				switch ringQLargeDim.Type() {
				case ring.Standard:
					ring.NTTStandard(tmpOut, tmpOut, NOut, s.Modulus, s.MRedConstant, s.BRedConstant, s.RootsForward)
				case ring.ConjugateInvariant:
					ring.NTTConjugateInvariant(tmpOut, tmpOut, NOut, s.Modulus, s.MRedConstant, s.BRedConstant, s.RootsForward)
				}
			}
		}

	} else {
		for i := range opOut.Value {
			ring.MapSmallDimensionToLargerDimensionNTT(ctIn.Value[i], opOut.Value[i])
		}
	}

	*opOut.MetaData = *ctIn.MetaData
}

// SwitchCiphertextRingDegree changes the ring degree of ctIn to the one of opOut.
// Maps Y^{N/n} -> X^{N} or X^{N} -> Y^{N/n}.
// If the ring degree of opOut is larger than the one of ctIn, then the ringQ of ctIn
// must be provided (otherwise, a nil pointer).
func SwitchCiphertextRingDegree(ctIn, opOut *Operand[ring.Poly]) {

	NIn, NOut := len(ctIn.Value[0].Coeffs[0]), len(opOut.Value[0].Coeffs[0])

	gapIn, gapOut := NOut/NIn, 1
	if NIn > NOut {
		gapIn, gapOut = 1, NIn/NOut
	}

	for i := range opOut.Value {
		for j := range opOut.Value[i].Coeffs {
			tmp0, tmp1 := opOut.Value[i].Coeffs[j], ctIn.Value[i].Coeffs[j]
			for w0, w1 := 0, 0; w0 < NOut; w0, w1 = w0+gapIn, w1+gapOut {
				tmp0[w0] = tmp1[w1]
			}
		}
	}

	*opOut.MetaData = *ctIn.MetaData
}

// BinarySize returns the serialized size of the object in bytes.
func (op Operand[T]) BinarySize() (size int) {
	size++
	if op.MetaData != nil {
		size += op.MetaData.BinarySize()
	}

	return size + op.Value.BinarySize()
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
func (op Operand[T]) WriteTo(w io.Writer) (n int64, err error) {

	switch w := w.(type) {
	case buffer.Writer:

		var inc int64

		if op.MetaData != nil {

			if inc, err = buffer.WriteUint8(w, 1); err != nil {
				return n, err
			}

			n += inc

			if inc, err = op.MetaData.WriteTo(w); err != nil {
				return n, err
			}

			n += inc

		} else {
			if inc, err = buffer.WriteUint8(w, 0); err != nil {
				return n, err
			}

			n += inc
		}

		inc, err = op.Value.WriteTo(w)

		return n + inc, err

	default:
		return op.WriteTo(bufio.NewWriter(w))
	}
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
func (op *Operand[T]) ReadFrom(r io.Reader) (n int64, err error) {

	switch r := r.(type) {
	case buffer.Reader:

		if op == nil {
			return 0, fmt.Errorf("cannot ReadFrom: target object is nil")
		}

		var inc int64

		var hasMetaData uint8

		if inc, err = buffer.ReadUint8(r, &hasMetaData); err != nil {
			return n, err
		}

		n += inc

		if hasMetaData == 1 {

			if op.MetaData == nil {
				op.MetaData = &MetaData{}
			}

			if inc, err = op.MetaData.ReadFrom(r); err != nil {
				return n, err
			}

			n += inc
		}

		inc, err = op.Value.ReadFrom(r)

		return n + inc, err

	default:
		return op.ReadFrom(bufio.NewReader(r))
	}
}

// MarshalBinary encodes the object into a binary form on a newly allocated slice of bytes.
func (op Operand[T]) MarshalBinary() (data []byte, err error) {
	buf := buffer.NewBufferSize(op.BinarySize())
	_, err = op.WriteTo(buf)
	return buf.Bytes(), err
}

// UnmarshalBinary decodes a slice of bytes generated by
// MarshalBinary or WriteTo on the object.
func (op *Operand[T]) UnmarshalBinary(p []byte) (err error) {
	_, err = op.ReadFrom(buffer.NewBuffer(p))
	return
}
