// Package lintrans implements homomorphic linear transformations for the BFV/BGV schemes.
package lintrans

import (
	"github.com/tuneinsight/lattigo/v6/circuits/common/lintrans"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/schemes"
	"github.com/tuneinsight/lattigo/v6/schemes/bgv"
	"github.com/tuneinsight/lattigo/v6/utils"
)

// Diagonals is a wrapper of [lintrans.Diagonals].
type Diagonals[T bgv.Integer] lintrans.Diagonals[T]

// DiagonalsIndexList returns the list of the non-zero diagonals of the square matrix.
// A non zero diagonals is a diagonal with a least one non-zero element.
func (m Diagonals[T]) DiagonalsIndexList() (indexes []int) {
	return lintrans.Diagonals[T](m).DiagonalsIndexList()
}

// Evaluate evaluates the linear transformation on the provided vector.
// add: c = a + b
// muladd: c = c + a * b
func (m Diagonals[T]) Evaluate(vector []T, newVec func(size int) []T, add func(a, b, c []T), muladd func(a, b, c []T)) (res []T) {

	slots := len(vector) >> 1 // 2 x n/2 matrix

	keys := utils.GetKeys(m)

	N1 := lintrans.FindBestBSGSRatio(keys, slots, 1)

	index, _, _ := lintrans.BSGSIndex(keys, slots, N1)

	res = newVec(2 * slots)

	for j := range index {

		rot := -j & (slots - 1)

		tmp := newVec(2 * slots)

		for _, i := range index[j] {

			v, ok := m[j+i]
			if !ok {
				v = m[j+i-slots]
			}

			muladd(utils.RotateSlice(vector[:slots], i), utils.RotateSlice(v[:slots], rot), tmp[:slots])
			muladd(utils.RotateSlice(vector[slots:], i), utils.RotateSlice(v[slots:], rot), tmp[slots:])
		}

		add(res[:slots], utils.RotateSlice(tmp[:slots], j), res[:slots])
		add(res[slots:], utils.RotateSlice(tmp[slots:], j), res[slots:])
	}

	return
}

// PermutationMapping is a struct storing
// a mapping: From -> To and a scaling value.
type PermutationMapping[T bgv.Integer] struct {
	From    int
	To      int
	Scaling T
}

// Permutation is a struct that defines generic permutations
// over 2 x n/2 matrices, treating each row independently.
//
// For example, given the 2x2 matrix [[a, b] [c, d]] that would
// be mapped to [[1b, 2a] [3c, 4d]] then the Permutation would
// contain the following map:
// {{0:{1, 2}, 1:{0: 1}}, {0:{0, 3}: 1:{1, 4}}}
type Permutation[T bgv.Integer] [2][]PermutationMapping[T]

// GetDiagonals returns the non-zero diagonals of the matrix
// representation of the permutation, which can be used to
// instantiate [Parameters].
func (p Permutation[T]) GetDiagonals(logSlots int) Diagonals[T] {

	slots := 1 << (logSlots - 1) // matrix of 2 x 2^{logSlots}

	diagonals := map[int][]T{}

	for i := range p {

		offset := i * slots

		for _, pm := range p[i] {

			From := pm.From
			To := pm.To
			Scaling := pm.Scaling

			diagIndex := (slots + From - To) & (slots - 1)

			if _, ok := diagonals[diagIndex]; !ok {
				diagonals[diagIndex] = make([]T, 2*slots)
			}

			diagonals[diagIndex][To+offset] = Scaling
		}
	}

	return Diagonals[T](diagonals)
}

// Parameters is a wrapper of [lintrans.Parameters].
type Parameters lintrans.Parameters

// LinearTransformation is a wrapper of [lintrans.Parameters].
type LinearTransformation lintrans.LinearTransformation

// GaloisElements returns the list of Galois elements required to evaluate the linear transformation.
func (lt LinearTransformation) GaloisElements(params rlwe.ParameterProvider) []uint64 {
	return lintrans.LinearTransformation(lt).GaloisElements(params)
}

// NewLinearTransformation instantiates a new [LinearTransformation] and is a wrapper of [lintrans.LinearTransformation].
func NewLinearTransformation(params rlwe.ParameterProvider, lt Parameters) LinearTransformation {
	return LinearTransformation(lintrans.NewLinearTransformation(params, lintrans.Parameters(lt)))
}

// Encode is a method used to encode a [LinearTransformation] and a wrapper of [lintrans.Encode].
func Encode[T bgv.Integer](ecd schemes.Encoder, diagonals Diagonals[T], allocated LinearTransformation) (err error) {
	return lintrans.Encode(
		ecd,
		lintrans.Diagonals[T](diagonals),
		lintrans.LinearTransformation(allocated))
}

// Evaluator is a struct for evaluating linear transformations on [rlwe.Ciphertexts].
// All fields of this struct are public, enabling custom instantiations.
type Evaluator struct {
	lintrans.Evaluator
}

// NewEvaluator instantiates a new [Evaluator] from a [schemes.Evaluator].
// The default [bgv.Evaluator] is compliant to the [schemes.Evaluator] interface.
func NewEvaluator(eval schemes.Evaluator) (linTransEval *Evaluator) {
	return &Evaluator{
		lintrans.Evaluator{
			Evaluator: eval,
		},
	}
}

// EvaluateNew takes as input a ciphertext ctIn and a linear transformation M and evaluate and returns opOut: M(ctIn).
func (eval Evaluator) EvaluateNew(ctIn *rlwe.Ciphertext, linearTransformation LinearTransformation) (opOut *rlwe.Ciphertext, err error) {
	ops, err := eval.EvaluateManyNew(ctIn, []LinearTransformation{linearTransformation})
	if err != nil {
		return nil, err
	}
	return ops[0], nil
}

// Evaluate takes as input a ciphertext ctIn, a linear transformation M and evaluates opOut: M(ctIn).
func (eval Evaluator) Evaluate(ctIn *rlwe.Ciphertext, linearTransformation LinearTransformation, opOut *rlwe.Ciphertext) (err error) {
	return eval.Evaluator.EvaluateMany(ctIn, []lintrans.LinearTransformation{lintrans.LinearTransformation(linearTransformation)}, []*rlwe.Ciphertext{opOut})
}

// EvaluateManyNew takes as input a ciphertext ctIn and a list of linear transformations [M0, M1, M2, ...] and returns opOut:[M0(ctIn), M1(ctIn), M2(ctInt), ...].
func (eval Evaluator) EvaluateManyNew(ctIn *rlwe.Ciphertext, linearTransformations []LinearTransformation) (opOut []*rlwe.Ciphertext, err error) {
	params := eval.GetRLWEParameters()
	opOut = make([]*rlwe.Ciphertext, len(linearTransformations))
	for i := range opOut {
		opOut[i] = rlwe.NewCiphertext(params, 1, linearTransformations[i].LevelQ)
	}
	return opOut, eval.EvaluateMany(ctIn, linearTransformations, opOut)
}

// EvaluateMany takes as input a ciphertext ctIn, a list of linear transformations [M0, M1, M2, ...] and a list of pre-allocated receiver opOut
// and evaluates opOut: [M0(ctIn), M1(ctIn), M2(ctIn), ...]
func (eval Evaluator) EvaluateMany(ctIn *rlwe.Ciphertext, linearTransformations []LinearTransformation, opOut []*rlwe.Ciphertext) (err error) {
	circuitLTs := make([]lintrans.LinearTransformation, len(linearTransformations))
	for i := range circuitLTs {
		circuitLTs[i] = lintrans.LinearTransformation(linearTransformations[i])
	}
	return eval.Evaluator.EvaluateMany(ctIn, circuitLTs, opOut)
}

// EvaluateSequentialNew takes as input a ciphertext ctIn and a list of linear transformations [M0, M1, M2, ...] and returns opOut:...M2(M1(M0(ctIn))
func (eval Evaluator) EvaluateSequentialNew(ctIn *rlwe.Ciphertext, linearTransformations []LinearTransformation) (opOut *rlwe.Ciphertext, err error) {
	opOut = rlwe.NewCiphertext(eval.GetRLWEParameters(), 1, linearTransformations[0].LevelQ)
	return opOut, eval.EvaluateSequential(ctIn, linearTransformations, opOut)
}

// EvaluateSequential takes as input a ciphertext ctIn and a list of linear transformations [M0, M1, M2, ...] and returns opOut:...M2(M1(M0(ctIn))
func (eval Evaluator) EvaluateSequential(ctIn *rlwe.Ciphertext, linearTransformations []LinearTransformation, opOut *rlwe.Ciphertext) (err error) {
	circuitLTs := make([]lintrans.LinearTransformation, len(linearTransformations))
	for i := range circuitLTs {
		circuitLTs[i] = lintrans.LinearTransformation(linearTransformations[i])
	}
	return eval.Evaluator.EvaluateSequential(ctIn, circuitLTs, opOut)
}
