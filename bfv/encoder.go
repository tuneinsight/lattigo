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

	ShallowCopy() Encoder
}

// Encoder is a structure that stores the parameters to encode values on a plaintext in a SIMD (Single-Instruction Multiple-Data) fashion.
type encoder struct {
	params Parameters

	indexMatrix []uint64
	scaler      ring.Scaler

	tInvModQ []uint64

	tmpPoly *ring.Poly
	tmpPtRt *PlaintextRingT

	tdividesQ bool
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

	var rescaleParams []uint64
	var scaler ring.Scaler
	var tdividesQ bool
	if !utils.IsInSliceUint64(params.T(), params.Q()) {
		rescaleParams = make([]uint64, len(ringQ.Modulus))
		for i, qi := range ringQ.Modulus {
			rescaleParams[i] = ring.MForm(ring.ModExp(params.T(), qi-2, qi), qi, ringQ.BredParams[i])
		}
		scaler = ring.NewRNSScaler(ringQ, ringT)
	} else {
		tdividesQ = true
	}

	return &encoder{
		params:      params,
		indexMatrix: indexMatrix,
		scaler:      scaler,
		tInvModQ:    rescaleParams,
		tmpPoly:     ringQ.NewPoly(),
		tmpPtRt:     NewPlaintextRingT(params),
		tdividesQ:   tdividesQ,
	}
}

// EncodeUint encodes an uint64 slice of size at most N on a plaintext.
func (ecd *encoder) EncodeUint(coeffs []uint64, p *Plaintext) {
	ptRt := &PlaintextRingT{p.Plaintext}

	// Encodes the values in RingT
	ecd.EncodeUintRingT(coeffs, ptRt)

	// Scales by Q/t
	ecd.ScaleUp(ptRt, p)
}

// EncodeUintRingT encodes a slice of uint64 into a Plaintext in R_t
func (ecd *encoder) EncodeUintRingT(coeffs []uint64, p *PlaintextRingT) {
	if len(coeffs) > len(ecd.indexMatrix) {
		panic("invalid input to encode: number of coefficients must be smaller or equal to the ring degree")
	}

	if len(p.Value.Coeffs[0]) != len(ecd.indexMatrix) {
		panic("invalid plaintext to receive encoding: number of coefficients does not match the ring degree")
	}

	for i := 0; i < len(coeffs); i++ {
		p.Value.Coeffs[0][ecd.indexMatrix[i]] = coeffs[i]
	}

	for i := len(coeffs); i < len(ecd.indexMatrix); i++ {
		p.Value.Coeffs[0][ecd.indexMatrix[i]] = 0
	}

	ecd.params.RingT().InvNTT(p.Value, p.Value)
}

// EncodeUintMul encodes an uint64 slice of size at most N on a PlaintextRingT (R_t) optimized for ciphertext-plaintext multiplication.
func (ecd *encoder) EncodeUintMul(coeffs []uint64, p *PlaintextMul) {

	ptRt := &PlaintextRingT{p.Plaintext}

	// Encodes the values in RingT
	ecd.EncodeUintRingT(coeffs, ptRt)

	// Puts in NTT+Montgomery domains of ringQ
	ecd.RingTToMul(ptRt, p)
}

// EncodeIntRingT encodes an int64 slice of size at most N on a plaintext. It also encodes the sign of the given integer (as its inverse modulo the plaintext modulus).
// The sign will correctly decode as long as the absolute value of the coefficient does not exceed half of the plaintext modulus.
func (ecd *encoder) EncodeIntRingT(coeffs []int64, p *PlaintextRingT) {

	if len(coeffs) > len(ecd.indexMatrix) {
		panic("invalid input to encode: number of coefficients must be smaller or equal to the ring degree")
	}

	if len(p.Value.Coeffs[0]) != len(ecd.indexMatrix) {
		panic("invalid plaintext to receive encoding: number of coefficients does not match the ring degree")
	}

	for i := 0; i < len(coeffs); i++ {

		if coeffs[i] < 0 {
			p.Value.Coeffs[0][ecd.indexMatrix[i]] = uint64(int64(ecd.params.T()) + coeffs[i])
		} else {
			p.Value.Coeffs[0][ecd.indexMatrix[i]] = uint64(coeffs[i])
		}
	}

	for i := len(coeffs); i < len(ecd.indexMatrix); i++ {
		p.Value.Coeffs[0][ecd.indexMatrix[i]] = 0
	}

	ecd.params.RingT().InvNTTLazy(p.Value, p.Value)
}

// EncodeInt encodes an int64 slice of size at most N on a plaintext.
func (ecd *encoder) EncodeInt(coeffs []int64, p *Plaintext) {
	ptRt := &PlaintextRingT{p.Plaintext}

	// Encodes the values in RingT
	ecd.EncodeIntRingT(coeffs, ptRt)

	// Scales by Q/t
	ecd.ScaleUp(ptRt, p)
}

// EncodeIntMul encodes an int64 slice of size at most N on a PlaintextRingT (R_t) optimized for ciphertext-plaintext multiplication.
func (ecd *encoder) EncodeIntMul(coeffs []int64, p *PlaintextMul) {
	ptRt := &PlaintextRingT{p.Plaintext}

	// Encodes the values in RingT
	ecd.EncodeIntRingT(coeffs, ptRt)

	// Puts in NTT+Montgomery domains of ringQ
	ecd.RingTToMul(ptRt, p)
}

// ScaleUp transforms a PlaintextRingT (R_t) into a Plaintext (R_q) by scaling up the coefficient by Q/t.
func (ecd *encoder) ScaleUp(ptRt *PlaintextRingT, pt *Plaintext) {
	if ecd.tdividesQ {
		ScaleUpTDividesQVec(ecd.params.RingQ(), ptRt.Value, pt.Value)
	} else {
		ScaleUpTCoprimeWithQVec(ecd.params.RingQ(), ecd.params.RingT(), ecd.tInvModQ, ecd.tmpPoly.Coeffs[0], ptRt.Value, pt.Value)
	}
}

// ScaleDown transforms a Plaintext (R_q) into a PlaintextRingT (R_t) by scaling down the coefficient by t/Q and rounding.
func (ecd *encoder) ScaleDown(pt *Plaintext, ptRt *PlaintextRingT) {
	if ecd.tdividesQ {
		ecd.params.RingQ().DivRoundByLastModulusManyLvl(pt.Level(), pt.Level(), pt.Value, ecd.tmpPoly, ptRt.Value)
	} else {
		ecd.scaler.DivByQOverTRoundedLvl(pt.Level(), pt.Value, ptRt.Value)
	}
}

// RingTToMul transforms a PlaintextRingT into a PlaintextMul by operating the NTT transform
// of R_q and putting the coefficients in Montgomery form.
func (ecd *encoder) RingTToMul(ptRt *PlaintextRingT, ptMul *PlaintextMul) {
	if ptRt.Value != ptMul.Value {
		copy(ptMul.Value.Coeffs[0], ptRt.Value.Coeffs[0])
	}
	for i := 1; i < len(ecd.params.RingQ().Modulus); i++ {
		copy(ptMul.Value.Coeffs[i], ptRt.Value.Coeffs[0])
	}

	ecd.params.RingQ().NTTLazy(ptMul.Value, ptMul.Value)
	ecd.params.RingQ().MForm(ptMul.Value, ptMul.Value)
}

// MulToRingT transforms a PlaintextMul into PlaintextRingT by operating the inverse NTT transform of R_q and
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
		panic(fmt.Errorf("unsupported plaintext type (%T)", pt))
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

// DecodeInt decodes a any plaintext type and write the coefficients in coeffs. It also decodes the sign
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

	var rnsScaler ring.Scaler
	if !ecd.tdividesQ {
		rnsScaler = ring.NewRNSScaler(ecd.params.RingQ(), ecd.params.RingT())
	}

	return &encoder{
		params:      ecd.params,
		indexMatrix: ecd.indexMatrix,
		scaler:      rnsScaler,
		tInvModQ:    ecd.tInvModQ,
		tmpPoly:     ecd.params.RingQ().NewPoly(),
		tmpPtRt:     NewPlaintextRingT(ecd.params),
		tdividesQ:   ecd.tdividesQ,
	}
}
