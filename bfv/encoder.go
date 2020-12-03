// Package bfv implements a RNS-accelerated Fan-Vercauteren version of Brakerski's scale invariant homomorphic encryption scheme. It provides modular arithmetic over the integers.
package bfv

import (
	"fmt"
	"math/big"
	"unsafe"

	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/utils"
)

// GaloisGen is an integer of order N/2 modulo M and that spans Z_M with the integer -1.
// The j-th ring automorphism takes the root zeta to zeta^(5j).
const GaloisGen uint64 = 5

// Encoder is an interface implementing the encoder.
type Encoder interface {
	EncodeUint(coeffs []uint64, plaintext *Plaintext)
	EncodeUintZt(coeffs []uint64, plaintext *PlaintextZT)
	EncodeUintMul(coeffs []uint64, plaintext *PlaintextMul)
	EncodeInt(coeffs []int64, plaintext *Plaintext)
	EncodeIntZt(coeffs []int64, plaintext *PlaintextZT)
	EncodeIntMul(coeffs []int64, plaintext *PlaintextMul)
	ZtToZq(*PlaintextZT, *Plaintext)
	ZqToZt(pt *Plaintext, ptZt *PlaintextZT)
	MulToZt(pt *PlaintextMul, ptZt *PlaintextZT)
	DecodeZt(p interface{}, ptZt *PlaintextZT)
	DecodeUint(plaintext interface{}, coeffs []uint64)
	DecodeInt(plaintext interface{}, coeffs []int64)
	DecodeUintNew(plaintext interface{}) (coeffs []uint64)
	DecodeIntNew(plaintext interface{}) (coeffs []int64)
}

// Encoder is a structure that stores the parameters to encode values on a plaintext in a SIMD (Single-Instruction Multiple-Data) fashion.
type encoder struct {
	params *Parameters

	ringQ *ring.Ring
	ringT *ring.Ring

	indexMatrix []uint64
	scaler      ring.Scaler
	deltaMont   []uint64

	tmpPoly *ring.Poly
	tmpPtZt *PlaintextZT
}

// NewEncoder creates a new encoder from the provided parameters.
func NewEncoder(params *Parameters) Encoder {

	var ringQ, ringT *ring.Ring
	var err error

	if ringQ, err = ring.NewRing(params.N(), params.qi); err != nil {
		panic(err)
	}

	if ringT, err = ring.NewRing(params.N(), []uint64{params.t}); err != nil {
		panic(err)
	}

	var m, pos, index1, index2 uint64

	slots := params.N()

	indexMatrix := make([]uint64, slots)

	logN := params.LogN()

	rowSize := params.N() >> 1
	m = (params.N() << 1)
	pos = 1

	for i := uint64(0); i < rowSize; i++ {

		index1 = (pos - 1) >> 1
		index2 = (m - pos - 1) >> 1

		indexMatrix[i] = utils.BitReverse64(index1, logN)
		indexMatrix[i|rowSize] = utils.BitReverse64(index2, logN)

		pos *= GaloisGen
		pos &= (m - 1)
	}

	return &encoder{
		params:      params.Copy(),
		ringQ:       ringQ,
		ringT:       ringT,
		indexMatrix: indexMatrix,
		deltaMont:   GenLiftParams(ringQ, params.t),
		scaler:      ring.NewRNSScaler(params.t, ringQ),
		tmpPoly:     ringT.NewPoly(),
		tmpPtZt:     NewPlaintextZT(params),
	}
}

// GenLiftParams generates the lifting parameters.
func GenLiftParams(ringQ *ring.Ring, t uint64) (deltaMont []uint64) {

	delta := new(big.Int).Quo(ringQ.ModulusBigint, ring.NewUint(t))

	deltaMont = make([]uint64, len(ringQ.Modulus))

	tmp := new(big.Int)
	bredParams := ringQ.GetBredParams()
	for i, Qi := range ringQ.Modulus {
		deltaMont[i] = tmp.Mod(delta, ring.NewUint(Qi)).Uint64()
		deltaMont[i] = ring.MForm(deltaMont[i], Qi, bredParams[i])
	}

	return
}

func (encoder *encoder) EncodeUintZt(coeffs []uint64, p *PlaintextZT) {
	if len(coeffs) > len(encoder.indexMatrix) {
		panic("invalid input to encode: number of coefficients must be smaller or equal to the ring degree")
	}

	if len(p.value.Coeffs[0]) != len(encoder.indexMatrix) {
		panic("invalid plaintext to receive encoding: number of coefficients does not match the ring degree")
	}

	for i := 0; i < len(coeffs); i++ {
		p.value.Coeffs[0][encoder.indexMatrix[i]] = coeffs[i]
	}

	for i := len(coeffs); i < len(encoder.indexMatrix); i++ {
		p.value.Coeffs[0][encoder.indexMatrix[i]] = 0
	}

	encoder.ringT.InvNTT(p.value, p.value)
}

// EncodeUint encodes an uint64 slice of size at most N on a plaintext.
func (encoder *encoder) EncodeUint(coeffs []uint64, p *Plaintext) {
	ptZt := &PlaintextZT{p.Element, p.Element.value[0]}

	// Encodes the values in Zt
	encoder.EncodeUintZt(coeffs, ptZt)

	// Scales by Q/t
	encoder.ZtToZq(ptZt, p)
}

func (encoder *encoder) EncodeUintMul(coeffs []uint64, p *PlaintextMul) {

	// Encodes the values in Zt
	encoder.EncodeUintZt(coeffs, &PlaintextZT{p.Element, p.Element.value[0]})

	// Puts in NTT+Montgomerry domains of ringQ
	encoder.nttMontZQ(p)
}

// EncodeInt encodes an int64 slice of size at most N on a plaintext. It also encodes the sign of the given integer (as its inverse modulo the plaintext modulus).
// The sign will correctly decode as long as the absolute value of the coefficient does not exceed half of the plaintext modulus.
func (encoder *encoder) EncodeIntZt(coeffs []int64, p *PlaintextZT) {

	if len(coeffs) > len(encoder.indexMatrix) {
		panic("invalid input to encode: number of coefficients must be smaller or equal to the ring degree")
	}

	if len(p.value.Coeffs[0]) != len(encoder.indexMatrix) {
		panic("invalid plaintext to receive encoding: number of coefficients does not match the ring degree")
	}

	for i := 0; i < len(coeffs); i++ {

		if coeffs[i] < 0 {
			p.value.Coeffs[0][encoder.indexMatrix[i]] = uint64(int64(encoder.params.t) + coeffs[i])
		} else {
			p.value.Coeffs[0][encoder.indexMatrix[i]] = uint64(coeffs[i])
		}
	}

	for i := len(coeffs); i < len(encoder.indexMatrix); i++ {
		p.value.Coeffs[0][encoder.indexMatrix[i]] = 0
	}

	encoder.ringT.InvNTTLazy(p.value, p.value)
}

func (encoder *encoder) EncodeInt(coeffs []int64, p *Plaintext) {
	ptZt := &PlaintextZT{p.Element, p.value}

	// Encodes the values in Zt
	encoder.EncodeIntZt(coeffs, ptZt)

	// Scales by Q/t
	encoder.ZtToZq(ptZt, p)
}

func (encoder *encoder) EncodeIntMul(coeffs []int64, p *PlaintextMul) {
	ptZt := &PlaintextZT{p.Element, p.value}

	// Encodes the values in Zt
	encoder.EncodeIntZt(coeffs, ptZt)

	// Puts in NTT+Montgomerry domains of ringQ
	encoder.nttMontZQ(p)
}

func (encoder *encoder) ZtToZq(ptzt *PlaintextZT, pol *Plaintext) {

	for i := len(encoder.ringQ.Modulus) - 1; i >= 0; i-- {
		tmp1 := pol.value.Coeffs[i]
		tmp2 := ptzt.value.Coeffs[0]
		d := encoder.deltaMont[i]
		qi := encoder.ringQ.Modulus[i]
		mredParams := encoder.ringQ.GetMredParams()[i]

		for j := uint64(0); j < encoder.ringQ.N; j = j + 8 {

			x := (*[8]uint64)(unsafe.Pointer(&tmp2[j]))
			z := (*[8]uint64)(unsafe.Pointer(&tmp1[j]))

			z[0] = ring.MRed(x[0], d, qi, mredParams)
			z[1] = ring.MRed(x[1], d, qi, mredParams)
			z[2] = ring.MRed(x[2], d, qi, mredParams)
			z[3] = ring.MRed(x[3], d, qi, mredParams)
			z[4] = ring.MRed(x[4], d, qi, mredParams)
			z[5] = ring.MRed(x[5], d, qi, mredParams)
			z[6] = ring.MRed(x[6], d, qi, mredParams)
			z[7] = ring.MRed(x[7], d, qi, mredParams)
		}
	}
}

func (encoder *encoder) nttMontZQ(p *PlaintextMul) {
	ringQ := encoder.ringQ
	for i := 1; i < len(ringQ.Modulus); i++ {
		copy(p.value.Coeffs[i], p.value.Coeffs[0])
	}

	ringQ.NTTLazy(p.value, p.value)
	ringQ.MForm(p.value, p.value)
}

func (encoder *encoder) ZqToZt(pt *Plaintext, ptZt *PlaintextZT) {
	encoder.scaler.DivByQOverTRounded(pt.value, ptZt.value)
}

func (encoder *encoder) MulToZt(pt *PlaintextMul, ptZt *PlaintextZT) {
	encoder.ringQ.InvNTTLvl(0, pt.value, ptZt.value)
	encoder.ringQ.InvMFormLvl(0, ptZt.value, ptZt.value)
}

func (encoder *encoder) DecodeZt(p interface{}, ptZt *PlaintextZT) {
	switch pt := p.(type) {
	case *Plaintext:
		encoder.ZqToZt(pt, ptZt)
	case *PlaintextMul:
		encoder.MulToZt(pt, ptZt)
	case *PlaintextZT:
		ptZt.Copy(pt.Element)
	default:
		panic(fmt.Errorf("unsuported plaintext type (%T)", pt))
	}
}

// DecodeUint decodes a batched plaintext and returns the coefficients in a uint64 slice.
func (encoder *encoder) DecodeUint(p interface{}, coeffs []uint64) {

	encoder.DecodeZt(p, encoder.tmpPtZt)

	encoder.ringT.NTT(encoder.tmpPtZt.value, encoder.tmpPoly)

	for i := uint64(0); i < encoder.ringQ.N; i++ {
		coeffs[i] = encoder.tmpPoly.Coeffs[0][encoder.indexMatrix[i]]
	}
}

func (encoder *encoder) DecodeUintNew(p interface{}) (coeffs []uint64) {
	coeffs = make([]uint64, encoder.ringQ.N)
	encoder.DecodeUint(p, coeffs)
	return
}

// DecodeInt decodes a batched plaintext and returns the coefficients in an int64 slice. It also decodes the sign (by centering the values around the plaintext
// modulus).
func (encoder *encoder) DecodeInt(p interface{}, coeffs []int64) {

	encoder.DecodeZt(p, encoder.tmpPtZt)

	encoder.ringT.NTT(encoder.tmpPtZt.value, encoder.tmpPoly)

	modulus := int64(encoder.params.t)
	var value int64
	for i := uint64(0); i < encoder.ringQ.N; i++ {

		value = int64(encoder.tmpPoly.Coeffs[0][encoder.indexMatrix[i]])
		coeffs[i] = value
		if value > modulus>>1 {
			coeffs[i] -= modulus
		}
	}
}

func (encoder *encoder) DecodeIntNew(p interface{}) (coeffs []int64) {
	coeffs = make([]int64, encoder.ringQ.N)
	encoder.DecodeInt(p, coeffs)
	return
}
