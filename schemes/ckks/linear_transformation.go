package ckks

import (
	"fmt"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/utils"
)

// TraceNew maps X -> sum((-1)^i * X^{i*n+1}) for 0 <= i < N and returns the result on a new ciphertext.
// For log(n) = logSlots.
func (eval Evaluator) TraceNew(ctIn *rlwe.Ciphertext, logSlots int) (opOut *rlwe.Ciphertext, err error) {
	opOut = NewCiphertext(*eval.GetParameters(), 1, ctIn.Level())
	return opOut, eval.Trace(ctIn, logSlots, opOut)
}

// Average returns the average of vectors of batchSize elements.
// The operation assumes that ctIn encrypts SlotCount/batchSize sub-vectors of size batchSize.
// It then replaces all values of those sub-vectors by the component-wise average between all the sub-vectors.
// Example for batchSize=4 and slots=8: [{a, b, c, d}, {e, f, g, h}] -> [0.5*{a+e, b+f, c+g, d+h}, 0.5*{a+e, b+f, c+g, d+h}]
// Operation requires log2(SlotCout/batchSize) rotations.
// Required rotation keys can be generated with RotationsForInnerSumLog(batchSize, SlotCount/batchSize).
func (eval Evaluator) Average(ctIn *rlwe.Ciphertext, logBatchSize int, opOut *rlwe.Ciphertext) (err error) {

	if ctIn.Degree() != 1 || opOut.Degree() != 1 {
		return fmt.Errorf("cannot Average: ctIn.Degree() != 1 or opOut.Degree() != 1")
	}

	if logBatchSize > ctIn.LogDimensions.Cols {
		return fmt.Errorf("cannot Average: batchSize must be smaller or equal to the number of slots")
	}

	ringQ := eval.GetParameters().RingQ()

	level := utils.Min(ctIn.Level(), opOut.Level())

	n := 1 << (ctIn.LogDimensions.Cols - logBatchSize)

	// pre-multiplication by n^-1
	for i, s := range ringQ.SubRings[:level+1] {

		invN := ring.ModExp(uint64(n), s.Modulus-2, s.Modulus)
		invN = ring.MForm(invN, s.Modulus, s.BRedConstant)

		s.MulScalarMontgomery(ctIn.Value[0].Coeffs[i], invN, opOut.Value[0].Coeffs[i])
		s.MulScalarMontgomery(ctIn.Value[1].Coeffs[i], invN, opOut.Value[1].Coeffs[i])
	}

	return eval.InnerSum(opOut, 1<<logBatchSize, n, opOut)
}
