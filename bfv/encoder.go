package bfv

import (
	"github.com/ldsec/lattigo/ring"
	"github.com/ldsec/lattigo/utils"
	"math/bits"
)

// Encoder is an interface implementing the encoder.
type Encoder interface {
	EncodeUint(coeffs []uint64, plaintext *Plaintext)
	EncodeInt(coeffs []int64, plaintext *Plaintext)
	DecodeUint(plaintext *Plaintext) (coeffs []uint64)
	DecodeInt(plaintext *Plaintext) (coeffs []int64)
}

// Encoder is a structure that stores the parameters to encode values on a plaintext in a SIMD (Single-Instruction Multiple-Data) fashion.
type encoder struct {
	params       *Parameters
	bfvContext   *bfvContext
	indexMatrix  []uint64
	simplescaler *ring.SimpleScaler
	polypool     *ring.Poly
	deltaMont    []uint64
}

// NewEncoder creates a new encoder from the provided parameters.
func NewEncoder(params *Parameters) Encoder {

	if !params.isValid {
		panic("cannot NewEncoder: params not valid (check if they were generated properly)")
	}

	bfvContext := newBFVContext(params)

	var m, pos, index1, index2 uint64

	slots := bfvContext.n

	indexMatrix := make([]uint64, slots)

	logN := uint64(bits.Len64(bfvContext.n) - 1)

	rowSize := bfvContext.n >> 1
	m = (bfvContext.n << 1)
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
		params:       params.Copy(),
		bfvContext:   bfvContext,
		indexMatrix:  indexMatrix,
		deltaMont:    GenLiftParams(bfvContext.contextQ, params.T),
		simplescaler: ring.NewSimpleScaler(params.T, bfvContext.contextQ),
		polypool:     bfvContext.contextT.NewPoly(),
	}
}

// EncodeUint encodes an uint64 slice of size at most N on a plaintext.
func (encoder *encoder) EncodeUint(coeffs []uint64, plaintext *Plaintext) {

	if len(coeffs) > len(encoder.indexMatrix) {
		panic("cannot EncodeUint: invalid input to encode (number of coefficients must be smaller or equal to the context)")
	}

	if len(plaintext.value.Coeffs[0]) != len(encoder.indexMatrix) {
		panic("cannot EncodeUint: invalid plaintext to receive encoding (number of coefficients does not match the context of the encoder)")
	}

	for i := 0; i < len(coeffs); i++ {
		plaintext.value.Coeffs[0][encoder.indexMatrix[i]] = coeffs[i]
	}

	for i := len(coeffs); i < len(encoder.indexMatrix); i++ {
		plaintext.value.Coeffs[0][encoder.indexMatrix[i]] = 0
	}

	encoder.encodePlaintext(plaintext)

}

// EncodeInt encodes an int64 slice of size at most N on a plaintext. It also encodes the sign of the given integer (as its inverse modulo the plaintext modulus).
// The sign will correctly decode as long as the absolute value of the coefficient does not exceed half of the plaintext modulus.
func (encoder *encoder) EncodeInt(coeffs []int64, plaintext *Plaintext) {

	if len(coeffs) > len(encoder.indexMatrix) {
		panic("cannot EncodeInt: invalid input to encode (number of coefficients must be smaller or equal to the context)")
	}

	if len(plaintext.value.Coeffs[0]) != len(encoder.indexMatrix) {
		panic("cannot EncodeInt: invalid plaintext to receive encoding (number of coefficients does not match the context of the encoder)")
	}

	for i := 0; i < len(coeffs); i++ {

		if coeffs[i] < 0 {
			plaintext.value.Coeffs[0][encoder.indexMatrix[i]] = uint64(int64(encoder.params.T) + coeffs[i])
		} else {
			plaintext.value.Coeffs[0][encoder.indexMatrix[i]] = uint64(coeffs[i])
		}
	}

	for i := len(coeffs); i < len(encoder.indexMatrix); i++ {
		plaintext.value.Coeffs[0][encoder.indexMatrix[i]] = 0
	}

	encoder.encodePlaintext(plaintext)
}

func (encoder *encoder) encodePlaintext(p *Plaintext) {

	encoder.bfvContext.contextT.InvNTT(p.value, p.value)

	ringContext := encoder.bfvContext.contextQ

	for i := len(ringContext.Modulus) - 1; i >= 0; i-- {
		tmp1 := p.value.Coeffs[i]
		tmp2 := p.value.Coeffs[0]
		deltaMont := encoder.deltaMont[i]
		qi := ringContext.Modulus[i]
		bredParams := ringContext.GetMredParams()[i]
		for j := uint64(0); j < ringContext.N; j++ {
			tmp1[j] = ring.MRed(tmp2[j], deltaMont, qi, bredParams)
		}
	}
}

// DecodeUint decodes a batched plaintext and returns the coefficients in a uint64 slice.
func (encoder *encoder) DecodeUint(plaintext *Plaintext) (coeffs []uint64) {

	encoder.simplescaler.Scale(plaintext.value, encoder.polypool)

	encoder.bfvContext.contextT.NTT(encoder.polypool, encoder.polypool)

	coeffs = make([]uint64, encoder.bfvContext.n)

	for i := uint64(0); i < encoder.bfvContext.n; i++ {
		coeffs[i] = encoder.polypool.Coeffs[0][encoder.indexMatrix[i]]
	}

	return

}

// DecodeInt decodes a batched plaintext and returns the coefficients in an int64 slice. It also decodes the sign (by centering the values around the plaintext
// modulus).
func (encoder *encoder) DecodeInt(plaintext *Plaintext) (coeffs []int64) {

	var value int64

	encoder.simplescaler.Scale(plaintext.value, encoder.polypool)

	encoder.bfvContext.contextT.NTT(encoder.polypool, encoder.polypool)

	coeffs = make([]int64, encoder.bfvContext.n)

	modulus := int64(encoder.params.T)

	for i := uint64(0); i < encoder.bfvContext.n; i++ {

		value = int64(encoder.polypool.Coeffs[0][encoder.indexMatrix[i]])

		coeffs[i] = value

		if value > modulus>>1 {
			coeffs[i] -= modulus
		}
	}

	return coeffs
}
