package bfv

import (
	"github.com/ldsec/lattigo/ring"
	"github.com/ldsec/lattigo/utils"
	"math/big"
	"math/bits"
)

type EncoderContext struct {
	// Polynomial degree
	n uint64

	// Plaintext Modulus
	t uint64

	// floor(Q/T) mod each Qi in montgomeryform
	deltaMont []uint64
	delta     []uint64

	// Polynomial contexts
	contextT *ring.Context
	contextQ *ring.Context

	// Galois elements used to permute the batched plaintext in the encrypted domain
	gen uint64
}

func NewEncoderContext(params *Parameters) *EncoderContext {
	n := params.N
	t := params.T

	contextT := ring.NewContext()
	// Plaintext NTT Parameters
	// We do not check for an error since the plaintext NTT is optional
	// it will still compute the other relevant parameters
	contextT.SetParameters(n, []uint64{t})
	if err := contextT.GenNTTParams(); err != nil {
		panic(err)
	}

	contextQ := ring.NewContext()
	contextQ.SetParameters(n, params.Qi)
	if err := contextQ.GenNTTParams(); err != nil {
		panic(err)
	}

	delta0 := new(big.Int).Quo(contextQ.ModulusBigint, ring.NewUint(t))
	tmpBig := new(big.Int)

	deltaMont := make([]uint64, len(params.Qi))
	delta := make([]uint64, len(params.Qi))

	for i, Qi := range params.Qi {
		delta[i] = tmpBig.Mod(delta0, ring.NewUint(Qi)).Uint64()
		deltaMont[i] = ring.MForm(delta[i], Qi, contextQ.GetBredParams()[i])
	}

	return &EncoderContext{
		n:         n,
		t:         t,
		deltaMont: deltaMont,
		delta:     delta,

		contextT: contextT,
		contextQ: contextQ,
		gen:      GaloisGen,
	}
}

// Encoder is a structure storing the parameters encode values on a plaintext in a SIMD fashion.
type Encoder struct {
	indexMatrix  []uint64
	context      *EncoderContext
	simplescaler *ring.SimpleScaler
	polypool     *ring.Poly
}

// NewEncoder creates a new encoder from the provided parameters
func NewEncoder(params *Parameters) (encoder *Encoder) {
	context := NewEncoderContext(params)

	if !context.contextT.AllowsNTT() {
		panic("cannot create batch encoder : plaintext modulus does not allow NTT")
	}

	var m, gen, pos, index1, index2 uint64

	encoder = new(Encoder)

	encoder.context = context

	slots := context.n

	encoder.indexMatrix = make([]uint64, slots)

	logN := uint64(bits.Len64(context.n) - 1)

	rowSize := context.n >> 1
	m = (context.n << 1)
	gen = context.gen
	pos = 1

	for i := uint64(0); i < rowSize; i++ {

		index1 = (pos - 1) >> 1
		index2 = (m - pos - 1) >> 1

		encoder.indexMatrix[i] = utils.BitReverse64(index1, logN)
		encoder.indexMatrix[i|rowSize] = utils.BitReverse64(index2, logN)

		pos *= gen
		pos &= (m - 1)
	}

	encoder.simplescaler = ring.NewSimpleScaler(context.t, context.contextQ)
	encoder.polypool = context.contextT.NewPoly()

	return encoder
}

// EncodeUint encodes an uint64 slice of size at most N on a plaintext.
func (encoder *Encoder) EncodeUint(coeffs []uint64, plaintext *Plaintext) {

	if len(coeffs) > len(encoder.indexMatrix) {
		panic("invalid input to encode (number of coefficients must be smaller or equal to the context)")
	}

	if len(plaintext.value.Coeffs[0]) != len(encoder.indexMatrix) {
		panic("invalid plaintext to receive encoding (number of coefficients does not match the context of the encoder")
	}

	for i := 0; i < len(coeffs); i++ {
		plaintext.value.Coeffs[0][encoder.indexMatrix[i]] = coeffs[i]
	}

	for i := len(coeffs); i < len(encoder.indexMatrix); i++ {
		plaintext.value.Coeffs[0][encoder.indexMatrix[i]] = 0
	}

	plaintext.InvNTTPlainModulus(encoder.context)

	plaintext.Lift(encoder.context)
}

// EncodeInt encodes an int64 slice of size at most N on a plaintext. Also encodes the sign of the given integer (as its inverse modulo the plaintext modulus).
// The sign will correctly decode as long as the absolute value of the coefficient does not exceed half of the plaintext modulus.
func (encoder *Encoder) EncodeInt(coeffs []int64, plaintext *Plaintext) {

	if len(coeffs) > len(encoder.indexMatrix) {
		panic("invalid input to encode (number of coefficients must be smaller or equal to the context)")
	}

	if len(plaintext.value.Coeffs[0]) != len(encoder.indexMatrix) {
		panic("invalid plaintext to receive encoding (number of coefficients does not match the context of the encoder)")
	}

	for i := 0; i < len(coeffs); i++ {

		if coeffs[i] < 0 {
			plaintext.value.Coeffs[0][encoder.indexMatrix[i]] = uint64(int64(encoder.context.t) + coeffs[i])
		} else {
			plaintext.value.Coeffs[0][encoder.indexMatrix[i]] = uint64(coeffs[i])
		}
	}

	for i := len(coeffs); i < len(encoder.indexMatrix); i++ {
		plaintext.value.Coeffs[0][encoder.indexMatrix[i]] = 0
	}

	plaintext.InvNTTPlainModulus(encoder.context)
	plaintext.Lift(encoder.context)
}

// DecodeUint decodes a batched plaintext and returns the coefficients in a uint64 slice.
func (encoder *Encoder) DecodeUint(plaintext *Plaintext) (coeffs []uint64) {

	encoder.simplescaler.Scale(plaintext.value, encoder.polypool)

	encoder.context.contextT.NTT(encoder.polypool, encoder.polypool)

	coeffs = make([]uint64, encoder.context.n)

	for i := uint64(0); i < encoder.context.n; i++ {
		coeffs[i] = encoder.polypool.Coeffs[0][encoder.indexMatrix[i]]
	}

	return

}

// DecodeInt decodes a batched plaintext and returns the coefficients in an int64 slice. Also decodes the sign (by centering the values around the plaintext
// modulus).
func (encoder *Encoder) DecodeInt(plaintext *Plaintext) (coeffs []int64) {

	var value int64

	encoder.simplescaler.Scale(plaintext.value, encoder.polypool)

	encoder.context.contextT.NTT(encoder.polypool, encoder.polypool)

	coeffs = make([]int64, encoder.context.n)

	modulus := int64(encoder.context.t)

	for i := uint64(0); i < encoder.context.n; i++ {

		value = int64(encoder.polypool.Coeffs[0][encoder.indexMatrix[i]])

		coeffs[i] = value

		if value > modulus>>1 {
			coeffs[i] -= modulus
		}
	}

	return coeffs
}
