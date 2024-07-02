package ltfloat

import (
	"github.com/tuneinsight/lattigo/v5/circuits/linear_transformation"
	"github.com/tuneinsight/lattigo/v5/core/rlwe"
	"github.com/tuneinsight/lattigo/v5/he"
	"github.com/tuneinsight/lattigo/v5/schemes"
	"github.com/tuneinsight/lattigo/v5/schemes/ckks"
	"github.com/tuneinsight/lattigo/v5/utils"
)

// Diagonals is a wrapper of [he.Diagonals].
// See [he.Diagonals] for the documentation.
type Diagonals[T ckks.Float] he.Diagonals[T]

// DiagonalsIndexList returns the list of the non-zero diagonals of the square matrix.
// A non zero diagonals is a diagonal with a least one non-zero element.
func (m Diagonals[T]) DiagonalsIndexList() (indexes []int) {
	return he.Diagonals[T](m).DiagonalsIndexList()
}

// Evaluate evaluates the linear transformation on the provided vector.
// add: c = a + b
// muladd: c = c + a * b
func (m Diagonals[T]) Evaluate(vector []T, newVec func(size int) []T, add func(a, b, c []T), muladd func(a, b, c []T)) (res []T) {

	slots := len(vector)

	keys := utils.GetKeys(m)

	N1 := he.FindBestBSGSRatio(keys, slots, 1)

	index, _, _ := he.BSGSIndex(keys, slots, N1)

	res = newVec(slots)

	for j := range index {

		rot := -j & (slots - 1)

		tmp := newVec(slots)

		for _, i := range index[j] {

			v, ok := m[j+i]
			if !ok {
				v = m[j+i-slots]
			}

			muladd(utils.RotateSlice(vector, i), utils.RotateSlice(v, rot), tmp)
		}

		add(res, utils.RotateSlice(tmp, j), res)
	}

	return
}

// PermutationMapping is a struct storing
// a mapping: From -> To and a scaling value.
type PermutationMapping[T ckks.Float] struct {
	From    int
	To      int
	Scaling T
}

// Permutation is a struct that defines generic permutations
// over vectors.
//
// For example, given the vector [a, b, c, d] that would be
// mapped to the vector [1b, 2c, 3d, 4a] then the Permutation
// would contain the following map:
// {0: {3, 4}, 1:{0, 1}, 2:{1, 2}, 3:{2, 3}}
type Permutation[T ckks.Float] []PermutationMapping[T]

// GetDiagonals returns the non-zero diagonals of the matrix
// representation of the permutation, which can be used to
// instantiate [LinearTransformationParameters].
func (p Permutation[T]) GetDiagonals(logSlots int) Diagonals[T] {

	slots := 1 << logSlots

	diagonals := map[int][]T{}

	for _, pm := range p {

		From := pm.From
		To := pm.To
		Scaling := pm.Scaling

		diagIndex := (slots + From - To) & (slots - 1)

		if _, ok := diagonals[diagIndex]; !ok {
			diagonals[diagIndex] = make([]T, slots)
		}

		diagonals[diagIndex][To] = Scaling
	}

	return Diagonals[T](diagonals)
}

// LinearTransformationParameters is a wrapper of he.LinearTransformationParameters.
// See he.LinearTransformationParameters for the documentation.
type LinearTransformationParameters linear_transformation.LinearTransformationParameters

// LinearTransformation is a wrapper of [he.LinearTransformation].
// See [he.LinearTransformation] for the documentation.
type LinearTransformation linear_transformation.LinearTransformation

// GaloisElements returns the list of Galois elements required to evaluate the linear transformation.
func (lt LinearTransformation) GaloisElements(params rlwe.ParameterProvider) []uint64 {
	return linear_transformation.LinearTransformation(lt).GaloisElements(params)
}

// NewLinearTransformation instantiates a new LinearTransformation and is a wrapper of [he.LinearTransformation].
// See [he.LinearTransformation] for the documentation.
func NewLinearTransformation(params rlwe.ParameterProvider, lt LinearTransformationParameters) LinearTransformation {
	return LinearTransformation(linear_transformation.NewLinearTransformation(params, linear_transformation.LinearTransformationParameters(lt)))
}

// EncodeLinearTransformation is a method used to encode EncodeLinearTransformation and a wrapper of [he.EncodeLinearTransformation].
// See [he.EncodeLinearTransformation] for the documentation.
func EncodeLinearTransformation[T ckks.Float](ecd schemes.Encoder, diagonals Diagonals[T], allocated LinearTransformation) (err error) {
	return linear_transformation.EncodeLinearTransformation(
		ecd,
		linear_transformation.Diagonals[T](diagonals),
		linear_transformation.LinearTransformation(allocated))
}

// // GaloisElementsForLinearTransformation returns the list of Galois elements required to evaluate the linear transformation.
func GaloisElementsForLinearTransformation(params rlwe.ParameterProvider, lt LinearTransformationParameters) (galEls []uint64) {
	return linear_transformation.GaloisElementsForLinearTransformation(params, lt.DiagonalsIndexList, 1<<lt.LogDimensions.Cols, lt.LogBabyStepGianStepRatio)
}

// LinearTransformationEvaluator is an evaluator providing an API to evaluate linear transformations on [rlwe.Ciphertexts].
// All fields of this struct are public, enabling custom instantiations.
type LinearTransformationEvaluator struct {
	linear_transformation.LinearTransformationEvaluator
	//he.EvaluatorForLinearTransformation
	//he.EvaluatorForDiagonalMatrix
}

// NewLinearTransformationEvaluator instantiates a new [LinearTransformationEvaluator] from a circuit.EvaluatorForLinearTransformation.
// The default [hefloat.Evaluator] is compliant to the [he.EvaluatorForLinearTransformation] interface.
// This method is allocation free.
func NewLinearTransformationEvaluator(eval schemes.Evaluator) (linTransEval *LinearTransformationEvaluator) {
	return &LinearTransformationEvaluator{
		linear_transformation.LinearTransformationEvaluator{
			Evaluator: eval,
		},
	}
}

// EvaluateNew takes as input a ciphertext ctIn and a linear transformation M and evaluate and returns opOut: M(ctIn).
func (eval LinearTransformationEvaluator) EvaluateNew(ctIn *rlwe.Ciphertext, linearTransformation LinearTransformation) (opOut *rlwe.Ciphertext, err error) {
	ops, err := eval.EvaluateManyNew(ctIn, []LinearTransformation{linearTransformation})
	if err != nil {
		return nil, err
	}
	return ops[0], nil
}

// Evaluate takes as input a ciphertext ctIn, a linear transformation M and evaluates opOut: M(ctIn).
func (eval LinearTransformationEvaluator) Evaluate(ctIn *rlwe.Ciphertext, linearTransformation LinearTransformation, opOut *rlwe.Ciphertext) (err error) {
	return eval.EvaluateLinearTransformationsMany(ctIn, []linear_transformation.LinearTransformation{linear_transformation.LinearTransformation(linearTransformation)}, []*rlwe.Ciphertext{opOut})
}

// EvaluateManyNew takes as input a ciphertext ctIn and a list of linear transformations [M0, M1, M2, ...] and returns opOut:[M0(ctIn), M1(ctIn), M2(ctInt), ...].
func (eval LinearTransformationEvaluator) EvaluateManyNew(ctIn *rlwe.Ciphertext, linearTransformations []LinearTransformation) (opOut []*rlwe.Ciphertext, err error) {
	params := eval.GetRLWEParameters()
	opOut = make([]*rlwe.Ciphertext, len(linearTransformations))
	for i := range opOut {
		opOut[i] = rlwe.NewCiphertext(params, 1, linearTransformations[i].LevelQ)
	}
	return opOut, eval.EvaluateMany(ctIn, linearTransformations, opOut)
}

// EvaluateMany takes as input a ciphertext ctIn, a list of linear transformations [M0, M1, M2, ...] and a list of pre-allocated receiver opOut
// and evaluates opOut: [M0(ctIn), M1(ctIn), M2(ctIn), ...]
func (eval LinearTransformationEvaluator) EvaluateMany(ctIn *rlwe.Ciphertext, linearTransformations []LinearTransformation, opOut []*rlwe.Ciphertext) (err error) {
	circuitLTs := make([]linear_transformation.LinearTransformation, len(linearTransformations))
	for i := range circuitLTs {
		circuitLTs[i] = linear_transformation.LinearTransformation(linearTransformations[i])
	}
	return eval.EvaluateLinearTransformationsMany(ctIn, circuitLTs, opOut)
}

// EvaluateSequentialNew takes as input a ciphertext ctIn and a list of linear transformations [M0, M1, M2, ...] and returns opOut:...M2(M1(M0(ctIn))
func (eval LinearTransformationEvaluator) EvaluateSequentialNew(ctIn *rlwe.Ciphertext, linearTransformations []LinearTransformation) (opOut *rlwe.Ciphertext, err error) {
	opOut = rlwe.NewCiphertext(eval.GetRLWEParameters(), 1, linearTransformations[0].LevelQ)
	return opOut, eval.EvaluateSequential(ctIn, linearTransformations, opOut)
}

// EvaluateSequential takes as input a ciphertext ctIn and a list of linear transformations [M0, M1, M2, ...] and returns opOut:...M2(M1(M0(ctIn))
func (eval LinearTransformationEvaluator) EvaluateSequential(ctIn *rlwe.Ciphertext, linearTransformations []LinearTransformation, opOut *rlwe.Ciphertext) (err error) {
	circuitLTs := make([]linear_transformation.LinearTransformation, len(linearTransformations))
	for i := range circuitLTs {
		circuitLTs[i] = linear_transformation.LinearTransformation(linearTransformations[i])
	}
	return eval.EvaluateLinearTranformationSequential(ctIn, circuitLTs, opOut)
}
