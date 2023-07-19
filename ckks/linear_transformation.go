package ckks

import (
	"fmt"
	"math/big"

	"github.com/tuneinsight/lattigo/v4/he"
	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/rlwe"
	"github.com/tuneinsight/lattigo/v4/rlwe/ringqp"
	"github.com/tuneinsight/lattigo/v4/utils"
	"github.com/tuneinsight/lattigo/v4/utils/bignum"
)

// NewLinearTransformation allocates a new LinearTransformation with zero values according to the parameters specified by the LinearTranfromationParameters.
func NewLinearTransformation[T float64 | complex128 | *big.Float | *bignum.Complex](params rlwe.ParametersInterface, lt he.LinearTranfromationParameters[T]) he.LinearTransformation {
	return he.NewLinearTransformation(params, lt)
}

// EncodeLinearTransformation encodes a linear transformation on a pre-allocated linear transformation.
// The method will return an error if the non-zero diagonals between the pre-allocated linear transformation and the parameters of the linear transformation to encode do not match.
func EncodeLinearTransformation[T float64 | complex128 | *big.Float | *bignum.Complex](allocated he.LinearTransformation, params he.LinearTranfromationParameters[T], ecd *Encoder) (err error) {
	return he.EncodeLinearTransformation[T](allocated, params, &encoder[T, ringqp.Poly]{ecd})
}

// TraceNew maps X -> sum((-1)^i * X^{i*n+1}) for 0 <= i < N and returns the result on a new ciphertext.
// For log(n) = logSlots.
func (eval Evaluator) TraceNew(ctIn *rlwe.Ciphertext, logSlots int) (opOut *rlwe.Ciphertext, err error) {
	opOut = NewCiphertext(eval.parameters, 1, ctIn.Level())
	return opOut, eval.Trace(ctIn, logSlots, opOut)
}

// Average returns the average of vectors of batchSize elements.
// The operation assumes that ctIn encrypts SlotCount/'batchSize' sub-vectors of size 'batchSize'.
// It then replaces all values of those sub-vectors by the component-wise average between all the sub-vectors.
// Example for batchSize=4 and slots=8: [{a, b, c, d}, {e, f, g, h}] -> [0.5*{a+e, b+f, c+g, d+h}, 0.5*{a+e, b+f, c+g, d+h}]
// Operation requires log2(SlotCout/'batchSize') rotations.
// Required rotation keys can be generated with 'RotationsForInnerSumLog(batchSize, SlotCount/batchSize)”
func (eval Evaluator) Average(ctIn *rlwe.Ciphertext, logBatchSize int, opOut *rlwe.Ciphertext) (err error) {

	if ctIn.Degree() != 1 || opOut.Degree() != 1 {
		return fmt.Errorf("cannot Average: ctIn.Degree() != 1 or opOut.Degree() != 1")
	}

	if logBatchSize > ctIn.PlaintextLogDimensions.Cols {
		return fmt.Errorf("cannot Average: batchSize must be smaller or equal to the number of slots")
	}

	ringQ := eval.parameters.RingQ()

	level := utils.Min(ctIn.Level(), opOut.Level())

	n := 1 << (ctIn.PlaintextLogDimensions.Cols - logBatchSize)

	// pre-multiplication by n^-1
	for i, s := range ringQ.SubRings[:level+1] {

		invN := ring.ModExp(uint64(n), s.Modulus-2, s.Modulus)
		invN = ring.MForm(invN, s.Modulus, s.BRedConstant)

		s.MulScalarMontgomery(ctIn.Value[0].Coeffs[i], invN, opOut.Value[0].Coeffs[i])
		s.MulScalarMontgomery(ctIn.Value[1].Coeffs[i], invN, opOut.Value[1].Coeffs[i])
	}

	return eval.InnerSum(opOut, 1<<logBatchSize, n, opOut)
}
