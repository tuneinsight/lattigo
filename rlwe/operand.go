package rlwe

import (
	"bytes"
	"fmt"
	"io"

	"github.com/google/go-cmp/cmp"
	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/rlwe/ringqp"
	"github.com/tuneinsight/lattigo/v4/utils/sampling"
	"github.com/tuneinsight/lattigo/v4/utils/structs"
)

// Operand is a common interface for OperandQ types.
type Operand interface {
	El() *OperandQ
	Degree() int
	Level() int
	GetScale() Scale
	SetScale(Scale)
}

type OperandQ struct {
	MetaData
	Value structs.Vector[ring.Poly]
}

func NewOperandQ(params Parameters, degree, levelQ int) *OperandQ {
	ringQ := params.RingQ().AtLevel(levelQ)

	Value := make([]*ring.Poly, degree+1)
	for i := range Value {
		Value[i] = ringQ.NewPoly()
	}

	return &OperandQ{
		Value: Value,
		MetaData: MetaData{
			IsNTT: params.DefaultNTTFlag(),
		},
	}
}

// NewOperandQAtLevelFromPoly constructs a new OperandQ at a specific level
// where the message is set to the passed poly. No checks are performed on poly and
// the returned OperandQ will share its backing array of coefficients.
// Returned OperandQ's MetaData is empty.
func NewOperandQAtLevelFromPoly(level int, poly []*ring.Poly) *OperandQ {
	Value := make([]*ring.Poly, len(poly))
	for i := range Value {

		if len(poly[i].Coeffs) < level+1 {
			panic(fmt.Errorf("cannot NewOperandQAtLevelFromPoly: provided ring.Poly[%d] level is too small", i))
		}

		Value[i] = new(ring.Poly)
		Value[i].Coeffs = poly[i].Coeffs[:level+1]
		Value[i].Buff = poly[i].Buff[:poly[i].N()*(level+1)]
	}

	return &OperandQ{Value: Value}
}

// Equal performs a deep equal.
func (op *OperandQ) Equal(other *OperandQ) bool {
	return cmp.Equal(op.MetaData, other.MetaData) && cmp.Equal(op.Value, other.Value)
}

// Degree returns the degree of the target OperandQ.
func (op *OperandQ) Degree() int {
	return len(op.Value) - 1
}

// Level returns the level of the target OperandQ.
func (op *OperandQ) Level() int {
	return len(op.Value[0].Coeffs) - 1
}

// GetScale gets the scale of the target OperandQ.
func (op *OperandQ) GetScale() Scale {
	return op.Scale
}

// SetScale sets the scale of the target OperandQ.
func (op *OperandQ) SetScale(scale Scale) {
	op.Scale = scale
}

func (op *OperandQ) El() *OperandQ {
	return op
}

// Resize resizes the degree of the target element.
// Sets the NTT flag of the added poly equal to the NTT flag
// to the poly at degree zero.
func (op *OperandQ) Resize(degree, level int) {

	if op.Level() != level {
		for i := range op.Value {
			op.Value[i].Resize(level)
		}
	}

	if op.Degree() > degree {
		op.Value = op.Value[:degree+1]
	} else if op.Degree() < degree {
		for op.Degree() < degree {
			op.Value = append(op.Value, []*ring.Poly{ring.NewPoly(op.Value[0].N(), level)}...)
		}
	}
}

// CopyNew creates a deep copy of the object and returns it.
func (op *OperandQ) CopyNew() *OperandQ {

	Value := make([]*ring.Poly, len(op.Value))

	for i := range Value {
		Value[i] = op.Value[i].CopyNew()
	}

	return &OperandQ{Value: Value, MetaData: op.MetaData}
}

// Copy copies the input element and its parameters on the target element.
func (op *OperandQ) Copy(opCopy *OperandQ) {

	if op != opCopy {
		for i := range opCopy.Value {
			op.Value[i].Copy(opCopy.Value[i])
		}

		op.MetaData = opCopy.MetaData
	}
}

// GetSmallestLargest returns the provided element that has the smallest degree as a first
// returned value and the largest degree as second return value. If the degree match, the
// order is the same as for the input.
func GetSmallestLargest(el0, el1 *OperandQ) (smallest, largest *OperandQ, sameDegree bool) {
	switch {
	case el0.Degree() > el1.Degree():
		return el1, el0, false
	case el0.Degree() < el1.Degree():
		return el0, el1, false
	}
	return el0, el1, true
}

// PopulateElementRandom creates a new rlwe.Element with random coefficients.
func PopulateElementRandom(prng sampling.PRNG, params Parameters, ct *OperandQ) {
	sampler := ring.NewUniformSampler(prng, params.RingQ()).AtLevel(ct.Level())
	for i := range ct.Value {
		sampler.Read(ct.Value[i])
	}
}

// SwitchCiphertextRingDegreeNTT changes the ring degree of ctIn to the one of ctOut.
// Maps Y^{N/n} -> X^{N} or X^{N} -> Y^{N/n}.
// If the ring degree of ctOut is larger than the one of ctIn, then the ringQ of ctOut
// must be provided (otherwise, a nil pointer).
// The ctIn must be in the NTT domain and ctOut will be in the NTT domain.
func SwitchCiphertextRingDegreeNTT(ctIn *OperandQ, ringQLargeDim *ring.Ring, ctOut *OperandQ) {

	NIn, NOut := len(ctIn.Value[0].Coeffs[0]), len(ctOut.Value[0].Coeffs[0])

	if NIn > NOut {

		gap := NIn / NOut
		buff := make([]uint64, NIn)
		for i := range ctOut.Value {
			for j := range ctOut.Value[i].Coeffs {

				tmpIn, tmpOut := ctIn.Value[i].Coeffs[j], ctOut.Value[i].Coeffs[j]

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
		for i := range ctOut.Value {
			ring.MapSmallDimensionToLargerDimensionNTT(ctIn.Value[i], ctOut.Value[i])
		}
	}

	ctOut.MetaData = ctIn.MetaData
}

// SwitchCiphertextRingDegree changes the ring degree of ctIn to the one of ctOut.
// Maps Y^{N/n} -> X^{N} or X^{N} -> Y^{N/n}.
// If the ring degree of ctOut is larger than the one of ctIn, then the ringQ of ctIn
// must be provided (otherwise, a nil pointer).
func SwitchCiphertextRingDegree(ctIn, ctOut *OperandQ) {

	NIn, NOut := len(ctIn.Value[0].Coeffs[0]), len(ctOut.Value[0].Coeffs[0])

	gapIn, gapOut := NOut/NIn, 1
	if NIn > NOut {
		gapIn, gapOut = 1, NIn/NOut
	}

	for i := range ctOut.Value {
		for j := range ctOut.Value[i].Coeffs {
			tmp0, tmp1 := ctOut.Value[i].Coeffs[j], ctIn.Value[i].Coeffs[j]
			for w0, w1 := 0, 0; w0 < NOut; w0, w1 = w0+gapIn, w1+gapOut {
				tmp0[w0] = tmp1[w1]
			}
		}
	}

	ctOut.MetaData = ctIn.MetaData
}

// MarshalBinary encodes the object into a binary form on a newly allocated slice of bytes.
func (op *OperandQ) MarshalBinary() (data []byte, err error) {
	buf := bytes.NewBuffer([]byte{})
	_, err = op.WriteTo(buf)
	return buf.Bytes(), err
}

// UnmarshalBinary decodes a slice of bytes generated by MarshalBinary
// or Read on the objeop.
func (op *OperandQ) UnmarshalBinary(p []byte) (err error) {
	_, err = op.ReadFrom(bytes.NewBuffer(p))
	return
}

// WriteTo writes the object on an io.Writer.
// To ensure optimal efficiency and minimal allocations, the user is encouraged
// to provide a struct implementing the interface buffer.Writer, which defines
// a subset of the method of the bufio.Writer.
// If w is not compliant to the buffer.Writer interface, it will be wrapped in
// a new bufio.Writer.
// For additional information, see lattigo/utils/buffer/writer.go.
func (op *OperandQ) WriteTo(w io.Writer) (n int64, err error) {

	if n, err = op.MetaData.WriteTo(w); err != nil {
		return n, err
	}

	inc, err := op.Value.WriteTo(w)

	return n + inc, err
}

// ReadFrom reads on the object from an io.Writer.
// To ensure optimal efficiency and minimal allocations, the user is encouraged
// to provide a struct implementing the interface buffer.Reader, which defines
// a subset of the method of the bufio.Reader.
// If r is not compliant to the buffer.Reader interface, it will be wrapped in
// a new bufio.Reader.
// For additional information, see lattigo/utils/buffer/reader.go.
func (op *OperandQ) ReadFrom(r io.Reader) (n int64, err error) {

	if op == nil {
		return 0, fmt.Errorf("cannot ReadFrom: target object is nil")
	}

	if n, err = op.MetaData.ReadFrom(r); err != nil {
		return n, err
	}

	inc, err := op.Value.ReadFrom(r)

	return n + inc, err
}

// BinarySize returns the size in bytes of the object
// when encoded using Encode.
func (op *OperandQ) BinarySize() int {
	return op.MetaData.BinarySize() + op.Value.BinarySize()
}

// Encode encodes the object into a binary form on a preallocated slice of bytes
// and returns the number of bytes written.
func (op *OperandQ) Encode(p []byte) (n int, err error) {

	if len(p) < op.BinarySize() {
		return 0, fmt.Errorf("cannot Encode: len(p) is too small")
	}

	if n, err = op.MetaData.Encode(p); err != nil {
		return
	}

	inc, err := op.Value.Encode(p[n:])

	return n + inc, err
}

// Decode decodes a slice of bytes generated by Encode
// on the object and returns the number of bytes read.
func (op *OperandQ) Decode(p []byte) (n int, err error) {

	if n, err = op.MetaData.Decode(p); err != nil {
		return
	}

	inc, err := op.Value.Decode(p[n:])

	return n + inc, err
}

type OperandQP struct {
	MetaData
	Value structs.Vector[ringqp.Poly]
}

func NewOperandQP(params Parameters, degree, levelQ, levelP int) *OperandQP {
	ringQP := params.RingQP().AtLevel(levelQ, levelP)

	Value := make([]*ringqp.Poly, degree+1)
	for i := range Value {
		Value[i] = ringQP.NewPoly()
	}

	return &OperandQP{
		Value: Value,
		MetaData: MetaData{
			IsNTT: params.DefaultNTTFlag(),
		},
	}
}

// GetScale gets the scale of the target OperandQP.
func (op *OperandQP) GetScale() Scale {
	return op.Scale
}

// SetScale sets the scale of the target OperandQP.
func (op *OperandQP) SetScale(scale Scale) {
	op.Scale = scale
}

// Equal performs a deep equal.
func (op *OperandQP) Equal(other *OperandQP) bool {
	return cmp.Equal(op.MetaData, other.MetaData) && cmp.Equal(op.Value, other.Value)
}

// LevelQ returns the level of the modulus Q of the first element of the objeop.
func (op *OperandQP) LevelQ() int {
	return op.Value[0].LevelQ()
}

// LevelP returns the level of the modulus P of the first element of the objeop.
func (op *OperandQP) LevelP() int {
	return op.Value[0].LevelP()
}

// CopyNew creates a deep copy of the object and returns it.
func (op *OperandQP) CopyNew() *OperandQP {

	Value := make([]*ringqp.Poly, len(op.Value))

	for i := range Value {
		Value[i] = op.Value[i].CopyNew()
	}

	return &OperandQP{Value: Value, MetaData: op.MetaData}
}

// MarshalBinary encodes the object into a binary form on a newly allocated slice of bytes.
func (op *OperandQP) MarshalBinary() (data []byte, err error) {
	buf := bytes.NewBuffer([]byte{})
	_, err = op.WriteTo(buf)
	return buf.Bytes(), err
}

// UnmarshalBinary decodes a slice of bytes generated by MarshalBinary
// or Read on the objeop.
func (op *OperandQP) UnmarshalBinary(p []byte) (err error) {
	_, err = op.ReadFrom(bytes.NewBuffer(p))
	return
}

// WriteTo writes the object on an io.Writer.
// To ensure optimal efficiency and minimal allocations, the user is encouraged
// to provide a struct implementing the interface buffer.Writer, which defines
// a subset of the method of the bufio.Writer.
// If w is not compliant to the buffer.Writer interface, it will be wrapped in
// a new bufio.Writer.
// For additional information, see lattigo/utils/buffer/writer.go.
func (op *OperandQP) WriteTo(w io.Writer) (n int64, err error) {

	if n, err = op.MetaData.WriteTo(w); err != nil {
		return n, err
	}

	inc, err := op.Value.WriteTo(w)

	return n + inc, err
}

// ReadFrom reads on the object from an io.Writer.
// To ensure optimal efficiency and minimal allocations, the user is encouraged
// to provide a struct implementing the interface buffer.Reader, which defines
// a subset of the method of the bufio.Reader.
// If r is not compliant to the buffer.Reader interface, it will be wrapped in
// a new bufio.Reader.
// For additional information, see lattigo/utils/buffer/reader.go.
func (op *OperandQP) ReadFrom(r io.Reader) (n int64, err error) {

	if op == nil {
		return 0, fmt.Errorf("cannot ReadFrom: target object is nil")
	}

	if n, err = op.MetaData.ReadFrom(r); err != nil {
		return n, err
	}

	inc, err := op.Value.ReadFrom(r)

	return n + inc, err
}

// BinarySize returns the size in bytes of the object
// when encoded using Encode.
func (op *OperandQP) BinarySize() int {
	return op.MetaData.BinarySize() + op.Value.BinarySize()
}

// Encode encodes the object into a binary form on a preallocated slice of bytes
// and returns the number of bytes written.
func (op *OperandQP) Encode(p []byte) (n int, err error) {

	if len(p) < op.BinarySize() {
		return 0, fmt.Errorf("cannote Encode: len(p) is too small")
	}

	if n, err = op.MetaData.Encode(p); err != nil {
		return
	}

	inc, err := op.Value.Encode(p[n:])

	return n + inc, err
}

// Decode decodes a slice of bytes generated by Encode
// on the object and returns the number of bytes read.
func (op *OperandQP) Decode(p []byte) (n int, err error) {

	if n, err = op.MetaData.Decode(p); err != nil {
		return
	}

	inc, err := op.Value.Decode(p[n:])

	return n + inc, err
}
