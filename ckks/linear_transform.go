package ckks

import (
	"fmt"
	"math/big"

	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/rlwe"
	"github.com/tuneinsight/lattigo/v4/rlwe/ringqp"
	"github.com/tuneinsight/lattigo/v4/utils"
	"github.com/tuneinsight/lattigo/v4/utils/bignum"
)

type LinearTransformEncoder[T float64 | complex128 | *big.Float | *bignum.Complex] struct {
	*Encoder
	diagonals map[int][]T
	values    []T
}

func NewLinearTransformEncoder[T float64 | complex128 | *big.Float | *bignum.Complex](ecd *Encoder, diagonals map[int][]T) LinearTransformEncoder[T] {
	return LinearTransformEncoder[T]{
		Encoder:   ecd,
		diagonals: diagonals,
		values:    make([]T, ecd.Parameters().MaxSlots()),
	}
}

func (l LinearTransformEncoder[_]) Parameters() rlwe.Parameters {
	return l.Encoder.Parameters().Parameters
}

// Diagonals returns the list of non-zero diagonals.
func (l LinearTransformEncoder[_]) NonZeroDiagonals() []int {
	return utils.GetKeys(l.diagonals)
}

func (l LinearTransformEncoder[_]) EncodeLinearTransformDiagonalNaive(i int, scale rlwe.Scale, LogSlots int, output ringqp.Poly) (err error) {

	if diag, ok := l.diagonals[i]; ok {
		return l.Embed(diag, LogSlots, scale, true, output)
	}

	return fmt.Errorf("cannot EncodeLinearTransformDiagonalNaive: diagonal [%d] doesn't exist", i)
}

func (l LinearTransformEncoder[_]) EncodeLinearTransformDiagonalBSGS(i, rot int, scale rlwe.Scale, logSlots int, output ringqp.Poly) (err error) {

	ecd := l.Encoder
	slots := 1 << logSlots
	values := l.values

	// manages inputs that have rotation between 0 and slots-1 or between -slots/2 and slots/2-1
	v, ok := l.diagonals[i]
	if !ok {
		v = l.diagonals[i-slots]
	}

	copyRotInterface(values[:slots], v, rot)

	return ecd.Embed(values[:slots], logSlots, scale, true, output)
}

func copyRotInterface[T any](a, b []T, rot int) {
	n := len(a)

	if len(b) >= rot {
		copy(a[:n-rot], b[rot:])
		copy(a[n-rot:], b[:rot])
	} else {
		copy(a[n-rot:], b)
	}
}

// TraceNew maps X -> sum((-1)^i * X^{i*n+1}) for 0 <= i < N and returns the result on a new ciphertext.
// For log(n) = logSlots.
func (eval *Evaluator) TraceNew(ctIn *rlwe.Ciphertext, logSlots int) (ctOut *rlwe.Ciphertext) {
	ctOut = NewCiphertext(eval.params, 1, ctIn.Level())
	eval.Trace(ctIn, logSlots, ctOut)
	return
}

// Average returns the average of vectors of batchSize elements.
// The operation assumes that ctIn encrypts SlotCount/'batchSize' sub-vectors of size 'batchSize'.
// It then replaces all values of those sub-vectors by the component-wise average between all the sub-vectors.
// Example for batchSize=4 and slots=8: [{a, b, c, d}, {e, f, g, h}] -> [0.5*{a+e, b+f, c+g, d+h}, 0.5*{a+e, b+f, c+g, d+h}]
// Operation requires log2(SlotCout/'batchSize') rotations.
// Required rotation keys can be generated with 'RotationsForInnerSumLog(batchSize, SlotCount/batchSize)”
func (eval *Evaluator) Average(ctIn *rlwe.Ciphertext, logBatchSize int, ctOut *rlwe.Ciphertext) {

	if ctIn.Degree() != 1 || ctOut.Degree() != 1 {
		panic("ctIn.Degree() != 1 or ctOut.Degree() != 1")
	}

	if logBatchSize > ctIn.LogSlots {
		panic("cannot Average: batchSize must be smaller or equal to the number of slots")
	}

	ringQ := eval.params.RingQ()

	level := utils.Min(ctIn.Level(), ctOut.Level())

	n := 1 << (ctIn.LogSlots - logBatchSize)

	// pre-multiplication by n^-1
	for i, s := range ringQ.SubRings[:level+1] {

		invN := ring.ModExp(uint64(n), s.Modulus-2, s.Modulus)
		invN = ring.MForm(invN, s.Modulus, s.BRedConstant)

		s.MulScalarMontgomery(ctIn.Value[0].Coeffs[i], invN, ctOut.Value[0].Coeffs[i])
		s.MulScalarMontgomery(ctIn.Value[1].Coeffs[i], invN, ctOut.Value[1].Coeffs[i])
	}

	eval.InnerSum(ctOut, 1<<logBatchSize, n, ctOut)
}
