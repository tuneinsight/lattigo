// Package bfv implements a RNS-accelerated Fan-Vercauteren version of Brakerski's scale invariant homomorphic encryption scheme. It provides modular arithmetic over the integers.
package bfv

import (
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
	EncodeAndLiftUint(coeffs []uint64, plaintext *Plaintext)
	EncodeInt(coeffs []int64, plaintext *Plaintext)
	EncodeAndLiftInt(coeffs []int64, plaintext *Plaintext)
	DecodeUint(plaintext *Plaintext) (coeffs []uint64)
	DecodeInt(plaintext *Plaintext) (coeffs []int64)
}

// Encoder is a structure that stores the parameters to encode values on a plaintext in a SIMD (Single-Instruction Multiple-Data) fashion.
type encoder struct {
	params *Parameters

	ringQ *ring.Ring
	ringT *ring.Ring

	indexMatrix []uint64
	scaler      ring.Scaler
	polypool    *ring.Poly
	deltaMont   []uint64
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
		polypool:    ringT.NewPoly(),
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

// EncodeUint encodes an uint64 slice of size at most N on a plaintext.
func (encoder *encoder) EncodeUint(coeffs []uint64, plaintext *Plaintext) {

	if len(coeffs) > len(encoder.indexMatrix) {
		panic("invalid input to encode: number of coefficients must be smaller or equal to the ring degree")
	}

	if len(plaintext.value.Coeffs[0]) != len(encoder.indexMatrix) {
		panic("invalid plaintext to receive encoding: number of coefficients does not match the ring degree")
	}

	for i := 0; i < len(coeffs); i++ {
		plaintext.value.Coeffs[0][encoder.indexMatrix[i]] = coeffs[i]
	}

	for i := len(coeffs); i < len(encoder.indexMatrix); i++ {
		plaintext.value.Coeffs[0][encoder.indexMatrix[i]] = 0
	}

	encoder.ringT.InvNTT(plaintext.value, plaintext.value)
}

// EncodeUint encodes an uint64 slice of size at most N on a plaintext.
func (encoder *encoder) EncodeAndLiftUint(coeffs []uint64, plaintext *Plaintext) {
	encoder.EncodeUint(coeffs, plaintext)
	encoder.TtoQ(plaintext)
}

// EncodeInt encodes an int64 slice of size at most N on a plaintext. It also encodes the sign of the given integer (as its inverse modulo the plaintext modulus).
// The sign will correctly decode as long as the absolute value of the coefficient does not exceed half of the plaintext modulus.
func (encoder *encoder) EncodeInt(coeffs []int64, plaintext *Plaintext) {

	if len(coeffs) > len(encoder.indexMatrix) {
		panic("invalid input to encode: number of coefficients must be smaller or equal to the ring degree")
	}

	if len(plaintext.value.Coeffs[0]) != len(encoder.indexMatrix) {
		panic("invalid plaintext to receive encoding: number of coefficients does not match the ring degree")
	}

	for i := 0; i < len(coeffs); i++ {

		if coeffs[i] < 0 {
			plaintext.value.Coeffs[0][encoder.indexMatrix[i]] = uint64(int64(encoder.params.t) + coeffs[i])
		} else {
			plaintext.value.Coeffs[0][encoder.indexMatrix[i]] = uint64(coeffs[i])
		}
	}

	for i := len(coeffs); i < len(encoder.indexMatrix); i++ {
		plaintext.value.Coeffs[0][encoder.indexMatrix[i]] = 0
	}

	encoder.ringT.InvNTTLazy(plaintext.value, plaintext.value)
}

func (encoder *encoder) EncodeAndLiftInt(coeffs []int64, plaintext *Plaintext) {
	encoder.EncodeInt(coeffs, plaintext)
	encoder.TtoQ(plaintext)
}

func (encoder *encoder) TtoQ(p *Plaintext) {

	ringQ := encoder.ringQ

	if !p.inZQ {
		additionalCoeffs := make([][]uint64, len(ringQ.Modulus)-1)
		for i := 0; i < len(ringQ.Modulus)-1; i++ {
			additionalCoeffs[i] = make([]uint64, ringQ.N)
		}
		p.value.Coeffs = append(p.value.Coeffs, additionalCoeffs...)
	}

	for i := len(ringQ.Modulus) - 1; i >= 0; i-- {
		tmp1 := p.value.Coeffs[i]
		tmp2 := p.value.Coeffs[0]
		deltaMont := encoder.deltaMont[i]
		qi := ringQ.Modulus[i]
		bredParams := ringQ.GetMredParams()[i]

		for j := uint64(0); j < ringQ.N; j = j + 8 {

			x := (*[8]uint64)(unsafe.Pointer(&tmp2[j]))
			z := (*[8]uint64)(unsafe.Pointer(&tmp1[j]))

			z[0] = ring.MRed(x[0], deltaMont, qi, bredParams)
			z[1] = ring.MRed(x[1], deltaMont, qi, bredParams)
			z[2] = ring.MRed(x[2], deltaMont, qi, bredParams)
			z[3] = ring.MRed(x[3], deltaMont, qi, bredParams)
			z[4] = ring.MRed(x[4], deltaMont, qi, bredParams)
			z[5] = ring.MRed(x[5], deltaMont, qi, bredParams)
			z[6] = ring.MRed(x[6], deltaMont, qi, bredParams)
			z[7] = ring.MRed(x[7], deltaMont, qi, bredParams)
		}
	}

	p.inZQ = true
}

// DecodeUint decodes a batched plaintext and returns the coefficients in a uint64 slice.
func (encoder *encoder) DecodeUint(plaintext *Plaintext) (coeffs []uint64) {

	if plaintext.inZQ {
		encoder.scaler.DivByQOverTRounded(plaintext.value, encoder.polypool)
		encoder.ringT.NTT(encoder.polypool, encoder.polypool)
	} else {
		encoder.ringT.NTT(plaintext.value, encoder.polypool)
	}

	coeffs = make([]uint64, encoder.ringQ.N)

	for i := uint64(0); i < encoder.ringQ.N; i++ {
		coeffs[i] = encoder.polypool.Coeffs[0][encoder.indexMatrix[i]]
	}

	return

}

// DecodeInt decodes a batched plaintext and returns the coefficients in an int64 slice. It also decodes the sign (by centering the values around the plaintext
// modulus).
func (encoder *encoder) DecodeInt(plaintext *Plaintext) (coeffs []int64) {

	var value int64

	if plaintext.inZQ {
		encoder.scaler.DivByQOverTRounded(plaintext.value, encoder.polypool)
		encoder.ringT.NTT(encoder.polypool, encoder.polypool)
	} else {
		encoder.ringT.NTT(plaintext.value, encoder.polypool)
	}

	coeffs = make([]int64, encoder.ringQ.N)

	modulus := int64(encoder.params.t)

	for i := uint64(0); i < encoder.ringQ.N; i++ {

		value = int64(encoder.polypool.Coeffs[0][encoder.indexMatrix[i]])

		coeffs[i] = value

		if value > modulus>>1 {
			coeffs[i] -= modulus
		}
	}

	return coeffs
}
