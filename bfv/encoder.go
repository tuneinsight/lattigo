// Package bfv implements a RNS-accelerated Fan-Vercauteren version of Brakerski's scale invariant homomorphic encryption scheme. It provides modular arithmetic over the integers.
package bfv

import (
	"fmt"

	"github.com/tuneinsight/lattigo/v3/ring"
	"github.com/tuneinsight/lattigo/v3/utils"
)

// GaloisGen is an integer of order N=2^d modulo M=2N and that spans Z_M with the integer -1.
// The j-th ring automorphism takes the root zeta to zeta^(5j).
const GaloisGen uint64 = 5

// Encoder is an interface for plaintext encoding and decoding operations. It provides methods to embed []uint64 and []int64 types into
// the various plaintext types and the inverse operations. It also provides methodes to convert between the different plaintext types.
// The different plaintext types represent different embeddings of the message in the polynomial space. This relation is illustrated in
// the figure below:
//                   ┌-> Encoder.Encode(.) -----------------------------------------------------┐
// []uint64/[]int64 -┼-> Encoder.EncodeRingT(.) ---> PlaintextRingT -┬-> Encoder.ScaleUp(.) ----┴-> Plaintext
//                   |                                               └-> Encoder.RingTToMul(.) -┬-> PlaintextMul
//                   └-> Encoder.EncodeMul(.) --------------------------------------------------┘
//
// The different plaintext types have different efficiency-related caracteristics that we summarize in the Table below. For more information
// about the different plaintext types, see plaintext.go.
//
// Relative efficiency of operations
//  -------------------------------------------------------------------------
// |                      |  PlaintextRingT  |  Plaintext  | PlaintextMul    |
//  -------------------------------------------------------------------------
// | Encoding/Decoding    |    Faster        |    Slower   |    Slower       |
// | Memory size          |    Smaller       |    Larger   |    Larger       |
// | Ct-Pt Add / Sub      |    Slower        |    Faster   |    N/A          |
// | Ct-Pt Mul            |    Faster        |    Slower   |    Much Faster  |
//  -------------------------------------------------------------------------
//
type Encoder interface {
	Encode(coeffs interface{}, pt *Plaintext)
	EncodeNew(coeffs interface{}, level int) (pt *Plaintext)
	EncodeRingT(coeffs interface{}, pt *PlaintextRingT)
	EncodeRingTNew(coeffs interface{}) (pt *PlaintextRingT)
	EncodeMul(coeffs interface{}, pt *PlaintextMul)
	EncodeMulNew(coeffs interface{}, level int) (pt *PlaintextMul)

	ScaleUp(*PlaintextRingT, *Plaintext)
	ScaleDown(pt *Plaintext, ptRt *PlaintextRingT)
	RingTToMul(ptRt *PlaintextRingT, ptmul *PlaintextMul)
	MulToRingT(pt *PlaintextMul, ptRt *PlaintextRingT)

	DecodeRingT(pt interface{}, ptRt *PlaintextRingT)
	DecodeUint(pt interface{}, coeffs []uint64)
	DecodeInt(pt interface{}, coeffs []int64)
	DecodeUintNew(pt interface{}) (coeffs []uint64)
	DecodeIntNew(pt interface{}) (coeffs []int64)

	ShallowCopy() Encoder
}

// encoder is a structure that stores the parameters to encode values on a plaintext in a SIMD (Single-Instruction Multiple-Data) fashion.
type encoder struct {
	params Parameters

	indexMatrix []uint64
	scaler      Scaler

	tInvModQ []uint64

	tmpPoly *ring.Poly
	tmpPtRt *PlaintextRingT
}

// NewEncoder creates a new encoder from the provided parameters.
func NewEncoder(params Parameters) Encoder {

	ringQ := params.RingQ()
	ringT := params.RingT()

	var m, pos, index1, index2 uint64

	slots := params.N()

	indexMatrix := make([]uint64, slots)

	logN := uint64(params.LogN())

	rowSize := params.N() >> 1
	m = uint64(params.N()) << 1
	pos = 1

	for i := 0; i < rowSize; i++ {

		index1 = (pos - 1) >> 1
		index2 = (m - pos - 1) >> 1

		indexMatrix[i] = utils.BitReverse64(index1, logN)
		indexMatrix[i|rowSize] = utils.BitReverse64(index2, logN)

		pos *= GaloisGen
		pos &= (m - 1)
	}

	return &encoder{
		params:      params,
		indexMatrix: indexMatrix,
		scaler:      NewRNSScaler(ringQ, ringT.Modulus[0]),
		tmpPoly:     ringQ.NewPoly(),
		tmpPtRt:     NewPlaintextRingT(params),
	}
}

// EncodeNew encodes a slice of integers of type []uint64 or []int64 of size at most N on a newly allocated plaintext.
func (ecd *encoder) EncodeNew(values interface{}, level int) (pt *Plaintext) {
	pt = NewPlaintextLvl(ecd.params, level)
	ecd.Encode(values, pt)
	return
}

// Encode encodes a slice of integers of type []uint64 or []int64 of size at most N into a pre-allocated plaintext.
func (ecd *encoder) Encode(values interface{}, pt *Plaintext) {
	ptRt := &PlaintextRingT{pt.Plaintext}

	// Encodes the values in RingT
	ecd.EncodeRingT(values, ptRt)

	// Scales by Q/t
	ecd.ScaleUp(ptRt, pt)
}

// EncodeRingTNew encodes a slice of integers of type []uint64 or []int64 of size at most N into a newly allocated PlaintextRingT.
func (ecd *encoder) EncodeRingTNew(values interface{}) (pt *PlaintextRingT) {
	pt = NewPlaintextRingT(ecd.params)
	ecd.EncodeRingT(values, pt)
	return
}

// EncodeRingT encodes a slice of integers of type []uint64 or []int64 of size at most N into a pre-allocated PlaintextRingT.
// The input values are reduced modulo T before encoding.
func (ecd *encoder) EncodeRingT(values interface{}, ptOut *PlaintextRingT) {

	if len(ptOut.Value.Coeffs[0]) != len(ecd.indexMatrix) {
		panic("cannot EncodeRingT: invalid plaintext to receive encoding: number of coefficients does not match the ring degree")
	}

	pt := ptOut.Value.Coeffs[0]

	ringT := ecd.params.RingT()

	var valLen int
	switch values := values.(type) {
	case []uint64:
		for i, c := range values {
			pt[ecd.indexMatrix[i]] = c
		}
		ringT.Reduce(ptOut.Value, ptOut.Value)
		valLen = len(values)
	case []int64:

		T := ringT.Modulus[0]
		bredparamsT := ringT.BredParams[0]

		var sign, abs uint64
		for i, c := range values {
			sign = uint64(c) >> 63
			abs = ring.BRedAdd(uint64(c*((int64(sign)^1)-int64(sign))), T, bredparamsT)
			pt[ecd.indexMatrix[i]] = sign*(T-abs) | (sign^1)*abs
		}
		valLen = len(values)
	default:
		panic("cannot EncodeRingT: coeffs must be either []uint64 or []int64")
	}

	for i := valLen; i < len(ecd.indexMatrix); i++ {
		pt[ecd.indexMatrix[i]] = 0
	}

	ringT.InvNTT(ptOut.Value, ptOut.Value)
}

// EncodeMulNew encodes a slice of integers of type []uint64 or []int64 of size at most N into a newly allocated PlaintextMul (optimized for ciphertext-plaintext multiplication).
func (ecd *encoder) EncodeMulNew(coeffs interface{}, level int) (pt *PlaintextMul) {
	pt = NewPlaintextMulLvl(ecd.params, level)
	ecd.EncodeMul(coeffs, pt)
	return
}

// EncodeMul encodes a slice of integers of type []uint64 or []int64 of size at most N into a pre-allocated PlaintextMul (optimized for ciphertext-plaintext multiplication).
func (ecd *encoder) EncodeMul(coeffs interface{}, pt *PlaintextMul) {

	ptRt := &PlaintextRingT{pt.Plaintext}

	// Encodes the values in RingT
	ecd.EncodeRingT(coeffs, ptRt)

	// Puts in NTT+Montgomery domains of ringQ
	ecd.RingTToMul(ptRt, pt)
}

// ScaleUp transforms a PlaintextRingT (R_t) into a Plaintext (R_q) by scaling up the coefficient by Q/t.
func (ecd *encoder) ScaleUp(ptRt *PlaintextRingT, pt *Plaintext) {
	ecd.scaler.ScaleUpByQOverTLvl(pt.Level(), ptRt.Value, pt.Value)
}

// ScaleDown transforms a Plaintext (R_q) into a PlaintextRingT (R_t) by scaling down the coefficient by t/Q and rounding.
func (ecd *encoder) ScaleDown(pt *Plaintext, ptRt *PlaintextRingT) {
	ecd.scaler.DivByQOverTRoundedLvl(pt.Level(), pt.Value, ptRt.Value)
}

// RingTToMul transforms a PlaintextRingT into a PlaintextMul by performing the NTT transform
// of R_q and putting the coefficients in Montgomery form.
func (ecd *encoder) RingTToMul(ptRt *PlaintextRingT, ptMul *PlaintextMul) {

	level := ptMul.Level()

	if ptRt.Value != ptMul.Value {
		copy(ptMul.Value.Coeffs[0], ptRt.Value.Coeffs[0])
	}
	for i := 1; i < level+1; i++ {
		copy(ptMul.Value.Coeffs[i], ptRt.Value.Coeffs[0])
	}

	ecd.params.RingQ().NTTLazyLvl(level, ptMul.Value, ptMul.Value)
	ecd.params.RingQ().MFormLvl(level, ptMul.Value, ptMul.Value)
}

// MulToRingT transforms a PlaintextMul into PlaintextRingT by performing the inverse NTT transform of R_q and
// putting the coefficients out of the Montgomery form.
func (ecd *encoder) MulToRingT(pt *PlaintextMul, ptRt *PlaintextRingT) {
	ecd.params.RingQ().InvNTTLvl(0, pt.Value, ptRt.Value)
	ecd.params.RingQ().InvMFormLvl(0, ptRt.Value, ptRt.Value)
}

// DecodeRingT decodes any plaintext type into a PlaintextRingT. It panics if p is not PlaintextRingT, Plaintext or PlaintextMul.
func (ecd *encoder) DecodeRingT(p interface{}, ptRt *PlaintextRingT) {
	switch pt := p.(type) {
	case *Plaintext:
		ecd.ScaleDown(pt, ptRt)
	case *PlaintextMul:
		ecd.MulToRingT(pt, ptRt)
	case *PlaintextRingT:
		ptRt.Copy(pt.Plaintext)
	default:
		panic(fmt.Errorf("cannot DecodeRingT: unsupported plaintext type (%T)", pt))
	}
}

// DecodeUint decodes a any plaintext type and write the coefficients in coeffs. It panics if p is not PlaintextRingT, Plaintext or PlaintextMul.
func (ecd *encoder) DecodeUint(p interface{}, coeffs []uint64) {

	var ptRt *PlaintextRingT
	var isInRingT bool
	if ptRt, isInRingT = p.(*PlaintextRingT); !isInRingT {
		ecd.DecodeRingT(p, ecd.tmpPtRt)
		ptRt = ecd.tmpPtRt
	}

	ecd.params.RingT().NTT(ptRt.Value, ecd.tmpPoly)

	for i := 0; i < ecd.params.RingQ().N; i++ {
		coeffs[i] = ecd.tmpPoly.Coeffs[0][ecd.indexMatrix[i]]
	}
}

// DecodeUintNew decodes any plaintext type and returns the coefficients in a new []uint64.
// It panics if p is not PlaintextRingT, Plaintext or PlaintextMul.
func (ecd *encoder) DecodeUintNew(p interface{}) (coeffs []uint64) {
	coeffs = make([]uint64, ecd.params.RingQ().N)
	ecd.DecodeUint(p, coeffs)
	return
}

// DecodeInt decodes a any plaintext type and writes the coefficients in coeffs. It also decodes the sign
// modulus (by centering the values around the plaintext). It panics if p is not PlaintextRingT, Plaintext or PlaintextMul.
func (ecd *encoder) DecodeInt(p interface{}, coeffs []int64) {

	ecd.DecodeRingT(p, ecd.tmpPtRt)

	ecd.params.RingT().NTT(ecd.tmpPtRt.Value, ecd.tmpPoly)

	modulus := int64(ecd.params.T())
	modulusHalf := modulus >> 1
	var value int64
	for i := 0; i < ecd.params.RingQ().N; i++ {

		value = int64(ecd.tmpPoly.Coeffs[0][ecd.indexMatrix[i]])
		coeffs[i] = value
		if value >= modulusHalf {
			coeffs[i] -= modulus
		}
	}
}

// DecodeIntNew decodes any plaintext type and returns the coefficients in a new []int64. It also decodes the sign
// modulus (by centering the values around the plaintext). It panics if p is not PlaintextRingT, Plaintext or PlaintextMul.
func (ecd *encoder) DecodeIntNew(p interface{}) (coeffs []int64) {
	coeffs = make([]int64, ecd.params.RingQ().N)
	ecd.DecodeInt(p, coeffs)
	return
}

// ShallowCopy creates a shallow copy of Encoder in which all the read-only data-structures are
// shared with the receiver and the temporary buffers are reallocated. The receiver and the returned
// Encoder can be used concurrently.
func (ecd *encoder) ShallowCopy() Encoder {
	return &encoder{
		params:      ecd.params,
		indexMatrix: ecd.indexMatrix,
		scaler:      NewRNSScaler(ecd.params.RingQ(), ecd.params.RingT().Modulus[0]),
		tmpPoly:     ecd.params.RingQ().NewPoly(),
		tmpPtRt:     NewPlaintextRingT(ecd.params),
	}
}
