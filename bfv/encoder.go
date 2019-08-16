package bfv

import (
	"errors"
	"github.com/lca1/lattigo-private/ring"
	"math/bits"
)

type BatchEncoder struct {
	indexMatrix  []uint64
	bfvcontext   *BfvContext
	simplescaler *ring.SimpleScaler
}

func (bfvcontext *BfvContext) NewBatchEncoder() *BatchEncoder {

	var m, gen, pos, index1, index2 uint64

	batchencoder := new(BatchEncoder)

	batchencoder.bfvcontext = bfvcontext

	slots := bfvcontext.n

	batchencoder.indexMatrix = make([]uint64, slots)

	logN := uint64(bits.Len64(bfvcontext.n) - 1)

	row_size := bfvcontext.n >> 1
	m = (bfvcontext.n << 1)
	gen = bfvcontext.gen
	pos = 1

	for i := uint64(0); i < row_size; i++ {

		index1 = (pos - 1) >> 1
		index2 = (m - pos - 1) >> 1

		batchencoder.indexMatrix[i] = bitReverse64(index1, logN)
		batchencoder.indexMatrix[i|row_size] = bitReverse64(index2, logN)

		pos *= gen
		pos &= (m - 1)
	}

	batchencoder.simplescaler, _ = ring.NewSimpleScaler(bfvcontext.t, bfvcontext.contextQ)

	return batchencoder
}

func (batchencoder *BatchEncoder) EncodeUint(coeffs []uint64, plaintext *Plaintext) error {

	if len(coeffs) > len(batchencoder.indexMatrix) {
		return errors.New("error : invalid input to encode (number of coefficients must be smaller or equal to the context)")
	}

	if len(plaintext.value[0].Coeffs[0]) != len(batchencoder.indexMatrix) {
		return errors.New("error : invalid plaintext to receive encoding (number of coefficients does not match the context of the encoder")
	}

	for i := 0; i < len(coeffs); i++ {
		plaintext.value[0].Coeffs[0][batchencoder.indexMatrix[i]] = coeffs[i]
	}

	for i := len(coeffs); i < len(batchencoder.indexMatrix); i++ {
		plaintext.value[0].Coeffs[0][batchencoder.indexMatrix[i]] = 0
	}

	if err := plaintext.EMB(); err != nil {
		return err
	}

	plaintext.Lift()

	return nil
}

func (batchencoder *BatchEncoder) EncodeInt(coeffs []int64, plaintext *Plaintext) error {

	if len(coeffs) > len(batchencoder.indexMatrix) {
		return errors.New("error : invalid input to encode (number of coefficients must be smaller or equal to the context)")
	}

	if len(plaintext.value[0].Coeffs[0]) != len(batchencoder.indexMatrix) {
		return errors.New("error : invalid plaintext to receive encoding (number of coefficients does not match the context of the encoder")
	}

	for i := 0; i < len(coeffs); i++ {

		if coeffs[i] < 0 {
			plaintext.value[0].Coeffs[0][batchencoder.indexMatrix[i]] = uint64(int64(plaintext.bfvcontext.t) + coeffs[i])
		} else {
			plaintext.value[0].Coeffs[0][batchencoder.indexMatrix[i]] = uint64(coeffs[i])
		}
	}

	for i := len(coeffs); i < len(batchencoder.indexMatrix); i++ {
		plaintext.value[0].Coeffs[0][batchencoder.indexMatrix[i]] = 0
	}

	if err := plaintext.EMB(); err != nil {
		return err
	}

	plaintext.Lift()

	return nil
}

func (batchencoder *BatchEncoder) DecodeUint(plaintext *Plaintext) ([]uint64, error) {

	if len(plaintext.value[0].Coeffs[0]) != len(batchencoder.indexMatrix) {
		return nil, errors.New("error : invalid plaintext to decode (number of coefficients does not match the context of the encoder")
	}

	tmp := batchencoder.bfvcontext.contextQ.NewPoly()

	batchencoder.simplescaler.Scale(plaintext.value[0], tmp)

	batchencoder.bfvcontext.contextT.NTT(tmp, tmp)

	coeffs := make([]uint64, batchencoder.bfvcontext.n)

	for i := uint64(0); i < batchencoder.bfvcontext.n; i++ {
		coeffs[i] = tmp.Coeffs[0][batchencoder.indexMatrix[i]]
	}

	return coeffs, nil

}

func (batchencoder *BatchEncoder) DecodeInt(plaintext *Plaintext) ([]int64, error) {

	var value int64

	if len(plaintext.value[0].Coeffs[0]) != len(batchencoder.indexMatrix) {
		return nil, errors.New("error : invalid plaintext to decode (number of coefficients does not match the context of the encoder")
	}

	tmp := batchencoder.bfvcontext.contextQ.NewPoly()

	batchencoder.simplescaler.Scale(plaintext.value[0], tmp)

	batchencoder.bfvcontext.contextT.NTT(tmp, tmp)

	coeffs := make([]int64, batchencoder.bfvcontext.n)

	modulus := int64(batchencoder.bfvcontext.t)

	for i := uint64(0); i < batchencoder.bfvcontext.n; i++ {

		value = int64(tmp.Coeffs[0][batchencoder.indexMatrix[i]])

		coeffs[i] = value

		if value > modulus>>1 {
			coeffs[i] -= modulus
		}
	}

	return coeffs, nil

}

type IntEncoder struct {
	base         int64
	simplescaler *ring.SimpleScaler
}

func (bfvcontext *BfvContext) NewIntEncoder(base int64) *IntEncoder {
	encoder := new(IntEncoder)
	encoder.base = base
	encoder.simplescaler, _ = ring.NewSimpleScaler(bfvcontext.t, bfvcontext.contextQ)
	return encoder
}

func (encoder *IntEncoder) Encode(msg int64, plaintext *Plaintext) {
	plaintext.value[0].Coeffs[0] = intEncode(msg, encoder.base, int64(plaintext.bfvcontext.t), plaintext.value[0].Coeffs[0])
	plaintext.Lift()
}

func (encoder *IntEncoder) Decode(plaintext *Plaintext) int64 {
	tmp := plaintext.bfvcontext.contextQ.NewPoly()
	encoder.simplescaler.Scale(plaintext.value[0], tmp)
	return intDecode(tmp.Coeffs[0], encoder.base, int64(plaintext.bfvcontext.t))
}

// Encodes an integer on a ring F given a base W.
// Decomposes the integer in base W, then set the coefficients
// of F accordingly.
// One has to be carefule about two things :
// 		- Make sure that base > T/2 (or it might not decode properly)
//		- Make sure that log_base(msg) <= N (or it will not fit in the ring)
func intEncode(msg, base, modulus int64, coeffs []uint64) []uint64 {

	//if base > plaintext.Context.Qi[0]>>1{
	//	return errors.New("Encoding base > T/2")
	//}
	//if msg != 0{
	//	msgLen := uint64(math.Round(math.Log(math.Abs(float64(msg)))/math.Log(float64(base))))
	//}
	//if msgLen > plaintext.Context.N{
	//	return errors.New("Messages does not fit on the ring with the given base")
	//}

	is_negative := msg < 0

	if is_negative {
		for i := 0; msg != 0; i++ {
			coeffs[i] = uint64(modulus + (msg % base))
			msg /= base
		}
	} else {
		for i := 0; msg != 0; i++ {
			coeffs[i] = uint64(msg % base)
			msg /= base
		}
	}

	return coeffs
}

// Decodes a integer from a ring F given a base W
// This is equivalent to evaluating F with the base W
// This evaluation is done with Horner's method (O(N))
func intDecode(coeffs []uint64, base, modulus int64) (msg int64) {

	msg = 0
	modulus_half := uint64(modulus >> 1)
	for i := len(coeffs) - 1; i >= 0; i-- {
		msg *= base
		msg += int64(coeffs[i])
		if coeffs[i] > modulus_half {
			msg -= modulus
		}
	}

	return msg
}
