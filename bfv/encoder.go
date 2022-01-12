// Package bfv implements a RNS-accelerated Fan-Vercauteren version of Brakerski's scale invariant homomorphic encryption scheme. It provides modular arithmetic over the integers.
package bfv

import (
	"fmt"
	"math/big"
	"unsafe"

	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/utils"
)

// GaloisGen is an integer of order N=2^d modulo M=2N and that spans Z_M with the integer -1.
// The j-th ring automorphism takes the root zeta to zeta^(5j).
const GaloisGen uint64 = 5

// Encoder is an interface for plaintext encoding and decoding operations. It provides methods to embed []uint64 and []int64 types into
// the various plaintext types and the inverse operations. It also provides methodes to convert between the different plaintext types.
// The different plaintext types represent different embeddings of the message in the polynomial space. This relation is illustrated in
// The figure below:
//
// []uint64 --- Encoder.EncodeUintRingT(.) -┬-> PlaintextRingT -┬-> Encoder.ScaleUp(.) -----> Plaintext
// []uint64 --- Encoder.EncodeIntRingT(.) --┘                   └-> Encoder.RingTToMul(.) ---> PlaintextMul
//
//
// The different plaintext types have different efficiency-related caracteristics that we summarize in the Table below. For more information
// about the different plaintext types, see plaintext.go.
//
// Relative efficiency of operation
//  -----------------------------------------------------------------------
// |                      |  PlaintextRingT  |  Plaintext  | PlaintextMul  |
//  -----------------------------------------------------------------------
// | Encoding/Decoding    |    Faster      |    Slower   |    Slower       |
// | Memory size          |    Smaller     |    Larger   |    Larger       |
// | Ct-Pt Add / Sub      |    Slower      |    Faster   |    N/A          |
// | Ct-Pt Mul            |    Faster      |    Slower   |    Much Faster  |
//  -----------------------------------------------------------------------
//
type Encoder interface {
	EncodeUint(coeffs []uint64, pt *Plaintext)
	EncodeUintRingT(coeffs []uint64, pt *PlaintextRingT)
	EncodeUintMul(coeffs []uint64, pt *PlaintextMul)
	EncodeInt(coeffs []int64, pt *Plaintext)
	EncodeIntRingT(coeffs []int64, pt *PlaintextRingT)
	EncodeIntMul(coeffs []int64, pt *PlaintextMul)

	ScaleUp(*PlaintextRingT, *Plaintext)
	ScaleDown(pt *Plaintext, ptRt *PlaintextRingT)
	RingTToMul(ptRt *PlaintextRingT, ptmul *PlaintextMul)
	MulToRingT(pt *PlaintextMul, ptRt *PlaintextRingT)

	DecodeRingT(pt interface{}, ptRt *PlaintextRingT)
	DecodeUint(pt interface{}, coeffs []uint64)
	DecodeInt(pt interface{}, coeffs []int64)
	DecodeUintNew(pt interface{}) (coeffs []uint64)
	DecodeIntNew(pt interface{}) (coeffs []int64)
}

// Encoder is a structure that stores the parameters to encode values on a plaintext in a SIMD (Single-Instruction Multiple-Data) fashion.
type encoder struct {
	params Parameters

	indexMatrix []uint64
	scaler      ring.Scaler

	rescaleParams []uint64

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

	rescaleParams := make([]uint64, len(ringQ.Modulus))
	for i, qi := range ringQ.Modulus {
		rescaleParams[i] = ring.MForm(ring.ModExp(params.T(), qi-2, qi), qi, ringQ.BredParams[i])
	}

	return &encoder{
		params:        params,
		indexMatrix:   indexMatrix,
		scaler:        ring.NewRNSScaler(ringQ, ringT),
		rescaleParams: rescaleParams,
		tmpPoly:       ringT.NewPoly(),
		tmpPtRt:       NewPlaintextRingT(params),
	}
}

// EncodeUintRingT encodes a slice of uint64 into a Plaintext in R_t
func (encoder *encoder) EncodeUintRingT(coeffs []uint64, p *PlaintextRingT) {
	if len(coeffs) > len(encoder.indexMatrix) {
		panic("invalid input to encode: number of coefficients must be smaller or equal to the ring degree")
	}

	if len(p.Value.Coeffs[0]) != len(encoder.indexMatrix) {
		panic("invalid plaintext to receive encoding: number of coefficients does not match the ring degree")
	}

	for i := 0; i < len(coeffs); i++ {
		p.Value.Coeffs[0][encoder.indexMatrix[i]] = coeffs[i]
	}

	for i := len(coeffs); i < len(encoder.indexMatrix); i++ {
		p.Value.Coeffs[0][encoder.indexMatrix[i]] = 0
	}

	encoder.params.RingT().InvNTT(p.Value, p.Value)
}

// EncodeUint encodes an uint64 slice of size at most N on a plaintext.
func (encoder *encoder) EncodeUint(coeffs []uint64, p *Plaintext) {
	ptRt := &PlaintextRingT{p.Plaintext}

	// Encodes the values in RingT
	encoder.EncodeUintRingT(coeffs, ptRt)

	// Scales by Q/t
	encoder.ScaleUp(ptRt, p)
}

func (encoder *encoder) EncodeUintMul(coeffs []uint64, p *PlaintextMul) {

	ptRt := &PlaintextRingT{p.Plaintext}

	// Encodes the values in RingT
	encoder.EncodeUintRingT(coeffs, ptRt)

	// Puts in NTT+Montgomery domains of ringQ
	encoder.RingTToMul(ptRt, p)
}

// EncodeInt encodes an int64 slice of size at most N on a plaintext. It also encodes the sign of the given integer (as its inverse modulo the plaintext modulus).
// The sign will correctly decode as long as the absolute value of the coefficient does not exceed half of the plaintext modulus.
func (encoder *encoder) EncodeIntRingT(coeffs []int64, p *PlaintextRingT) {

	if len(coeffs) > len(encoder.indexMatrix) {
		panic("invalid input to encode: number of coefficients must be smaller or equal to the ring degree")
	}

	if len(p.Value.Coeffs[0]) != len(encoder.indexMatrix) {
		panic("invalid plaintext to receive encoding: number of coefficients does not match the ring degree")
	}

	for i := 0; i < len(coeffs); i++ {

		if coeffs[i] < 0 {
			p.Value.Coeffs[0][encoder.indexMatrix[i]] = uint64(int64(encoder.params.T()) + coeffs[i])
		} else {
			p.Value.Coeffs[0][encoder.indexMatrix[i]] = uint64(coeffs[i])
		}
	}

	for i := len(coeffs); i < len(encoder.indexMatrix); i++ {
		p.Value.Coeffs[0][encoder.indexMatrix[i]] = 0
	}

	encoder.params.RingT().InvNTTLazy(p.Value, p.Value)
}

func (encoder *encoder) EncodeInt(coeffs []int64, p *Plaintext) {
	ptRt := &PlaintextRingT{p.Plaintext}

	// Encodes the values in RingT
	encoder.EncodeIntRingT(coeffs, ptRt)

	// Scales by Q/t
	encoder.ScaleUp(ptRt, p)
}

func (encoder *encoder) EncodeIntMul(coeffs []int64, p *PlaintextMul) {
	ptRt := &PlaintextRingT{p.Plaintext}

	// Encodes the values in RingT
	encoder.EncodeIntRingT(coeffs, ptRt)

	// Puts in NTT+Montgomery domains of ringQ
	encoder.RingTToMul(ptRt, p)
}

// ScaleUp transforms a PlaintextRingT (R_t) into a Plaintext (R_q) by scaling up the coefficient by Q/t.
func (encoder *encoder) ScaleUp(ptRt *PlaintextRingT, pt *Plaintext) {
	encoder.scaleUp(encoder.params.RingQ(), encoder.params.RingT(), encoder.tmpPoly.Coeffs[0], ptRt.Value, pt.Value)
}

// takes m mod T and returns round((m*Q)/T) mod Q
func (encoder *encoder) scaleUp(ringQ, ringT *ring.Ring, tmp []uint64, pIn, pOut *ring.Poly) {

	qModTmontgomery := ring.MForm(new(big.Int).Mod(ringQ.ModulusBigint, ringT.ModulusBigint).Uint64(), ringT.Modulus[0], ringT.BredParams[0])

	t := ringT.Modulus[0]
	tHalf := t >> 1
	tInv := ringT.MredParams[0]

	// (x * Q + T/2) mod T
	for i := 0; i < ringQ.N; i = i + 8 {
		x := (*[8]uint64)(unsafe.Pointer(&pIn.Coeffs[0][i]))
		z := (*[8]uint64)(unsafe.Pointer(&tmp[i]))

		z[0] = ring.CRed(ring.MRed(x[0], qModTmontgomery, t, tInv)+tHalf, t)
		z[1] = ring.CRed(ring.MRed(x[1], qModTmontgomery, t, tInv)+tHalf, t)
		z[2] = ring.CRed(ring.MRed(x[2], qModTmontgomery, t, tInv)+tHalf, t)
		z[3] = ring.CRed(ring.MRed(x[3], qModTmontgomery, t, tInv)+tHalf, t)
		z[4] = ring.CRed(ring.MRed(x[4], qModTmontgomery, t, tInv)+tHalf, t)
		z[5] = ring.CRed(ring.MRed(x[5], qModTmontgomery, t, tInv)+tHalf, t)
		z[6] = ring.CRed(ring.MRed(x[6], qModTmontgomery, t, tInv)+tHalf, t)
		z[7] = ring.CRed(ring.MRed(x[7], qModTmontgomery, t, tInv)+tHalf, t)
	}

	// (x * T^-1 - T/2) mod Qi
	for i := 0; i < len(pOut.Coeffs); i++ {
		p0tmp := tmp
		p1tmp := pOut.Coeffs[i]
		qi := ringQ.Modulus[i]
		bredParams := ringQ.BredParams[i]
		mredParams := ringQ.MredParams[i]
		rescaleParams := qi - encoder.rescaleParams[i]

		tHalfNegQi := qi - ring.BRedAdd(tHalf, qi, bredParams)

		for j := 0; j < ringQ.N; j = j + 8 {

			x := (*[8]uint64)(unsafe.Pointer(&p0tmp[j]))
			z := (*[8]uint64)(unsafe.Pointer(&p1tmp[j]))

			z[0] = ring.MRed(x[0]+tHalfNegQi, rescaleParams, qi, mredParams)
			z[1] = ring.MRed(x[1]+tHalfNegQi, rescaleParams, qi, mredParams)
			z[2] = ring.MRed(x[2]+tHalfNegQi, rescaleParams, qi, mredParams)
			z[3] = ring.MRed(x[3]+tHalfNegQi, rescaleParams, qi, mredParams)
			z[4] = ring.MRed(x[4]+tHalfNegQi, rescaleParams, qi, mredParams)
			z[5] = ring.MRed(x[5]+tHalfNegQi, rescaleParams, qi, mredParams)
			z[6] = ring.MRed(x[6]+tHalfNegQi, rescaleParams, qi, mredParams)
			z[7] = ring.MRed(x[7]+tHalfNegQi, rescaleParams, qi, mredParams)
		}
	}
}

// ScaleDown transforms a Plaintext (R_q) into a PlaintextRingT (R_t) by scaling down the coefficient by t/Q and rounding.
func (encoder *encoder) ScaleDown(pt *Plaintext, ptRt *PlaintextRingT) {
	encoder.scaler.DivByQOverTRounded(pt.Value, ptRt.Value)
}

// RingTToMul transforms a PlaintextRingT into a PlaintextMul by operating the NTT transform
// of R_q and putting the coefficients in Montgomery form.
func (encoder *encoder) RingTToMul(ptRt *PlaintextRingT, ptMul *PlaintextMul) {
	if ptRt.Value != ptMul.Value {
		copy(ptMul.Value.Coeffs[0], ptRt.Value.Coeffs[0])
	}
	for i := 1; i < len(encoder.params.RingQ().Modulus); i++ {
		copy(ptMul.Value.Coeffs[i], ptRt.Value.Coeffs[0])
	}

	encoder.params.RingQ().NTTLazy(ptMul.Value, ptMul.Value)
	encoder.params.RingQ().MForm(ptMul.Value, ptMul.Value)
}

// MulToRingT transforms a PlaintextMul into PlaintextRingT by operating the inverse NTT transform of R_q and
// putting the coefficients out of the Montgomery form.
func (encoder *encoder) MulToRingT(pt *PlaintextMul, ptRt *PlaintextRingT) {
	encoder.params.RingQ().InvNTTLvl(0, pt.Value, ptRt.Value)
	encoder.params.RingQ().InvMFormLvl(0, ptRt.Value, ptRt.Value)
}

// DecodeRingT decodes any plaintext type into a PlaintextRingT. It panics if p is not PlaintextRingT, Plaintext or PlaintextMul.
func (encoder *encoder) DecodeRingT(p interface{}, ptRt *PlaintextRingT) {
	switch pt := p.(type) {
	case *Plaintext:
		encoder.ScaleDown(pt, ptRt)
	case *PlaintextMul:
		encoder.MulToRingT(pt, ptRt)
	case *PlaintextRingT:
		ptRt.Copy(pt.Plaintext)
	default:
		panic(fmt.Errorf("unsupported plaintext type (%T)", pt))
	}
}

// DecodeUint decodes a any plaintext type and write the coefficients in coeffs. It panics if p is not PlaintextRingT, Plaintext or PlaintextMul.
func (encoder *encoder) DecodeUint(p interface{}, coeffs []uint64) {

	var ptRt *PlaintextRingT
	var isInRingT bool
	if ptRt, isInRingT = p.(*PlaintextRingT); !isInRingT {
		encoder.DecodeRingT(p, encoder.tmpPtRt)
		ptRt = encoder.tmpPtRt
	}

	encoder.params.RingT().NTT(ptRt.Value, encoder.tmpPoly)

	for i := 0; i < encoder.params.RingQ().N; i++ {
		coeffs[i] = encoder.tmpPoly.Coeffs[0][encoder.indexMatrix[i]]
	}
}

// DecodeUintNew decodes any plaintext type and returns the coefficients in a new []uint64.
// It panics if p is not PlaintextRingT, Plaintext or PlaintextMul.
func (encoder *encoder) DecodeUintNew(p interface{}) (coeffs []uint64) {
	coeffs = make([]uint64, encoder.params.RingQ().N)
	encoder.DecodeUint(p, coeffs)
	return
}

// DecodeInt decodes a any plaintext type and write the coefficients in coeffs. It also decodes the sign
// modulus (by centering the values around the plaintext). It panics if p is not PlaintextRingT, Plaintext or PlaintextMul.
func (encoder *encoder) DecodeInt(p interface{}, coeffs []int64) {

	encoder.DecodeRingT(p, encoder.tmpPtRt)

	encoder.params.RingT().NTT(encoder.tmpPtRt.Value, encoder.tmpPoly)

	modulus := int64(encoder.params.T())
	modulusHalf := modulus >> 1
	var value int64
	for i := 0; i < encoder.params.RingQ().N; i++ {

		value = int64(encoder.tmpPoly.Coeffs[0][encoder.indexMatrix[i]])
		coeffs[i] = value
		if value >= modulusHalf {
			coeffs[i] -= modulus
		}
	}
}

// DecodeIntNew decodes any plaintext type and returns the coefficients in a new []int64. It also decodes the sign
// modulus (by centering the values around the plaintext). It panics if p is not PlaintextRingT, Plaintext or PlaintextMul.
func (encoder *encoder) DecodeIntNew(p interface{}) (coeffs []int64) {
	coeffs = make([]int64, encoder.params.RingQ().N)
	encoder.DecodeInt(p, coeffs)
	return
}
