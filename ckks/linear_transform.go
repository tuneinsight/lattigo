package ckks

import (
	"math/big"

	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/rlwe"
	"github.com/tuneinsight/lattigo/v4/rlwe/ringqp"
	"github.com/tuneinsight/lattigo/v4/utils"
	"github.com/tuneinsight/lattigo/v4/utils/bignum"
)

// NewLinearTransform allocates a new LinearTransform with zero plaintexts at the specified level.
// If LogBSGSRatio < 0, the LinearTransform is set to not use the BSGS approach.
func NewLinearTransform(params Parameters, nonZeroDiags []int, level int, scale rlwe.Scale, LogSlots, LogBSGSRatio int) rlwe.LinearTransform {
	return rlwe.NewLinearTransform(params, nonZeroDiags, level, scale, [2]int{0, LogSlots}, LogBSGSRatio)
}

func EncodeLinearTransform[T float64 | complex128 | *big.Float | *bignum.Complex](LT rlwe.LinearTransform, diagonals map[int][]T, ecd *Encoder) (err error) {
	return rlwe.EncodeLinearTransform[T](LT, diagonals, &encoder[T, ringqp.Poly]{ecd})
}

func GenLinearTransform[T float64 | complex128 | *big.Float | *bignum.Complex](diagonals map[int][]T, ecd *Encoder, level int, scale rlwe.Scale, LogSlots, LogBSGSRatio int) (LT rlwe.LinearTransform, err error) {
	return rlwe.GenLinearTransform[T](diagonals, &encoder[T, ringqp.Poly]{ecd}, level, scale, [2]int{0, LogSlots}, LogBSGSRatio)
}

// TraceNew maps X -> sum((-1)^i * X^{i*n+1}) for 0 <= i < N and returns the result on a new ciphertext.
// For log(n) = logSlots.
func (eval Evaluator) TraceNew(ctIn *rlwe.Ciphertext, logSlots int) (opOut *rlwe.Ciphertext) {
	opOut = NewCiphertext(eval.parameters, 1, ctIn.Level())
	eval.Trace(ctIn, logSlots, opOut)
	return
}

// Average returns the average of vectors of batchSize elements.
// The operation assumes that ctIn encrypts SlotCount/'batchSize' sub-vectors of size 'batchSize'.
// It then replaces all values of those sub-vectors by the component-wise average between all the sub-vectors.
// Example for batchSize=4 and slots=8: [{a, b, c, d}, {e, f, g, h}] -> [0.5*{a+e, b+f, c+g, d+h}, 0.5*{a+e, b+f, c+g, d+h}]
// Operation requires log2(SlotCout/'batchSize') rotations.
// Required rotation keys can be generated with 'RotationsForInnerSumLog(batchSize, SlotCount/batchSize)”
func (eval Evaluator) Average(ctIn *rlwe.Ciphertext, logBatchSize int, opOut *rlwe.Ciphertext) {

	if ctIn.Degree() != 1 || opOut.Degree() != 1 {
		panic("ctIn.Degree() != 1 or opOut.Degree() != 1")
	}

	if logBatchSize > ctIn.PlaintextLogDimensions[1] {
		panic("cannot Average: batchSize must be smaller or equal to the number of slots")
	}

	ringQ := eval.parameters.RingQ()

	level := utils.Min(ctIn.Level(), opOut.Level())

	n := 1 << (ctIn.PlaintextLogDimensions[1] - logBatchSize)

	// pre-multiplication by n^-1
	for i, s := range ringQ.SubRings[:level+1] {

		invN := ring.ModExp(uint64(n), s.Modulus-2, s.Modulus)
		invN = ring.MForm(invN, s.Modulus, s.BRedConstant)

		s.MulScalarMontgomery(ctIn.Value[0].Coeffs[i], invN, opOut.Value[0].Coeffs[i])
		s.MulScalarMontgomery(ctIn.Value[1].Coeffs[i], invN, opOut.Value[1].Coeffs[i])
	}

	eval.InnerSum(opOut, 1<<logBatchSize, n, opOut)
}
