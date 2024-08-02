package heint

import (
	"github.com/tuneinsight/lattigo/v5/core/rlwe"
	"github.com/tuneinsight/lattigo/v5/he"
	"github.com/tuneinsight/lattigo/v5/ring"
	"github.com/tuneinsight/lattigo/v5/ring/ringqp"
	"github.com/tuneinsight/lattigo/v5/utils"
)

type intEncoder[T Integer, U ring.Poly | ringqp.Poly | *rlwe.Plaintext] struct {
	*Encoder
}

func (e intEncoder[T, U]) Encode(values []T, metadata *rlwe.MetaData, output U) (err error) {
	return e.Embed(values, false, metadata, output)
}

// Diagonals is a wrapper of [he.Diagonals].
// See [he.Diagonals] for the documentation.
type Diagonals[T Integer] he.Diagonals[T]

// DiagonalsIndexList returns the list of the non-zero diagonals of the square matrix.
// A non zero diagonals is a diagonal with a least one non-zero element.
func (m Diagonals[T]) DiagonalsIndexList() (indexes []int) {
	return he.Diagonals[T](m).DiagonalsIndexList()
}

// Evaluate evaluates the linear transformation on the provided vector.
// add: c = a + b
// muladd: c = c + a * b
func (m Diagonals[T]) Evaluate(vector []T, newVec func(size int) []T, add func(a, b, c []T), muladd func(a, b, c []T)) (res []T) {

	slots := len(vector) >> 1 // 2 x n/2 matrix

	keys := utils.GetKeys(m)

	N1 := he.FindBestBSGSRatio(keys, slots, 1)

	index, _, _ := he.BSGSIndex(keys, slots, N1)

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
type PermutationMapping[T Integer] struct {
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
type Permutation[T Integer] [2][]PermutationMapping[T]

// GetDiagonals returns the non-zero diagonals of the matrix
// representation of the permutation, which can be used to
// instantiate [LinearTransformationParameters].
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

// LinearTransformationParameters is a wrapper of [he.LinearTransformationParameters].
// See [he.LinearTransformationParameters] for the documentation.
type LinearTransformationParameters he.LinearTransformationParameters

// LinearTransformation is a wrapper of [he.LinearTransformation].
// See [he.LinearTransformation] for the documentation.
type LinearTransformation he.LinearTransformation

// GaloisElements returns the list of Galois elements required to evaluate the linear transformation.
func (lt LinearTransformation) GaloisElements(params rlwe.ParameterProvider) []uint64 {
	return he.LinearTransformation(lt).GaloisElements(params)
}

// NewLinearTransformation instantiates a new [LinearTransformation] and is a wrapper of [he.LinearTransformation].
// See he.LinearTransformation for the documentation.
func NewLinearTransformation(params rlwe.ParameterProvider, lt LinearTransformationParameters) LinearTransformation {
	return LinearTransformation(he.NewLinearTransformation(params, he.LinearTransformationParameters(lt)))
}

// EncodeLinearTransformation is a method used to encode EncodeLinearTransformation and a wrapper of [he.EncodeLinearTransformation].
// See [he.EncodeLinearTransformation] for the documentation.
func EncodeLinearTransformation[T Integer](ecd *Encoder, diagonals Diagonals[T], allocated LinearTransformation) (err error) {
	return he.EncodeLinearTransformation[T](
		&intEncoder[T, ringqp.Poly]{ecd},
		he.Diagonals[T](diagonals),
		he.LinearTransformation(allocated))
}

// GaloisElementsForLinearTransformation returns the list of Galois elements required to evaluate the linear transformation.
func GaloisElementsForLinearTransformation(params rlwe.ParameterProvider, lt LinearTransformationParameters) (galEls []uint64) {
	return he.GaloisElementsForLinearTransformation(params, lt.DiagonalsIndexList, 1<<lt.LogDimensions.Cols, lt.LogBabyStepGianStepRatio)
}

// LinearTransformationEvaluator is a struct for evaluating linear transformations on [rlwe.Ciphertexts].
// All fields of this struct are public, enabling custom instantiations.
type LinearTransformationEvaluator struct {
	he.EvaluatorForLinearTransformation
	he.EvaluatorForDiagonalMatrix
}

// NewLinearTransformationEvaluator instantiates a new [LinearTransformationEvaluator] from a circuit.EvaluatorForLinearTransformation.
// The default [heint.Evaluator] is compliant to the circuit.EvaluatorForLinearTransformation interface.
func NewLinearTransformationEvaluator(eval he.EvaluatorForLinearTransformation) (linTransEval *LinearTransformationEvaluator) {
	return &LinearTransformationEvaluator{
		EvaluatorForLinearTransformation: eval,
		EvaluatorForDiagonalMatrix:       &defaultDiagonalMatrixEvaluator{eval},
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
	return he.EvaluateLinearTransformationsMany(eval.EvaluatorForLinearTransformation, eval.EvaluatorForDiagonalMatrix, ctIn, []he.LinearTransformation{he.LinearTransformation(linearTransformation)}, []*rlwe.Ciphertext{opOut})
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
	circuitLTs := make([]he.LinearTransformation, len(linearTransformations))
	for i := range circuitLTs {
		circuitLTs[i] = he.LinearTransformation(linearTransformations[i])
	}
	return he.EvaluateLinearTransformationsMany(eval.EvaluatorForLinearTransformation, eval.EvaluatorForDiagonalMatrix, ctIn, circuitLTs, opOut)
}

// EvaluateSequentialNew takes as input a ciphertext ctIn and a list of linear transformations [M0, M1, M2, ...] and returns opOut:...M2(M1(M0(ctIn))
func (eval LinearTransformationEvaluator) EvaluateSequentialNew(ctIn *rlwe.Ciphertext, linearTransformations []LinearTransformation) (opOut *rlwe.Ciphertext, err error) {
	opOut = rlwe.NewCiphertext(eval.GetRLWEParameters(), 1, linearTransformations[0].LevelQ)
	return opOut, eval.EvaluateSequential(ctIn, linearTransformations, opOut)
}

// EvaluateSequential takes as input a ciphertext ctIn and a list of linear transformations [M0, M1, M2, ...] and returns opOut:...M2(M1(M0(ctIn))
func (eval LinearTransformationEvaluator) EvaluateSequential(ctIn *rlwe.Ciphertext, linearTransformations []LinearTransformation, opOut *rlwe.Ciphertext) (err error) {
	circuitLTs := make([]he.LinearTransformation, len(linearTransformations))
	for i := range circuitLTs {
		circuitLTs[i] = he.LinearTransformation(linearTransformations[i])
	}
	return he.EvaluateLinearTransformationSequential(eval.EvaluatorForLinearTransformation, eval.EvaluatorForDiagonalMatrix, ctIn, circuitLTs, opOut)
}

// defaultDiagonalMatrixEvaluator is a struct implementing the interface [he.EvaluatorForDiagonalMatrix].
type defaultDiagonalMatrixEvaluator struct {
	he.EvaluatorForLinearTransformation
}

// Decompose applies the RNS decomposition on ct[1] at the given level and stores the result in BuffDecompQP.
func (eval defaultDiagonalMatrixEvaluator) Decompose(levelQ, levelP int, ct *rlwe.Ciphertext, BuffDecompQP []ringqp.Poly) {
	eval.DecomposeNTT(levelQ, levelP, levelP+1, ct.Value[1], ct.IsNTT, BuffDecompQP)
}

// GetPreRotatedCiphertextForDiagonalMatrixMultiplication populates ctPreRot with the pre-rotated ciphertext for the rotations rots and deletes rotated ciphertexts that are not in rots.
func (eval defaultDiagonalMatrixEvaluator) GetPreRotatedCiphertextForDiagonalMatrixMultiplication(levelQ, levelP int, ctIn *rlwe.Ciphertext, BuffDecompQP []ringqp.Poly, rots []int, ctPreRot map[int]*rlwe.Element[ringqp.Poly]) (err error) {
	return he.GetPreRotatedCiphertextForDiagonalMatrixMultiplication(levelQ, levelP, eval, ctIn, BuffDecompQP, rots, ctPreRot)
}

// MultiplyByDiagMatrix multiplies the Ciphertext ctIn by the plaintext matrix and returns the result on the Ciphertext
// opOut. Memory buffers for the decomposed ciphertext BuffDecompQP, BuffDecompQP must be provided, those are list of poly of ringQ and ringP
// respectively, each of size params.Beta().
// The naive approach is used (single hoisting and no baby-step giant-step), which is faster than MultiplyByDiagMatrixBSGS
// for matrix of only a few non-zero diagonals but uses more keys.
func (eval defaultDiagonalMatrixEvaluator) MultiplyByDiagMatrix(ctIn *rlwe.Ciphertext, matrix he.LinearTransformation, BuffDecompQP []ringqp.Poly, opOut *rlwe.Ciphertext) (err error) {
	return he.MultiplyByDiagMatrix(eval.EvaluatorForLinearTransformation, ctIn, matrix, BuffDecompQP, opOut)
}

// MultiplyByDiagMatrixBSGS multiplies the Ciphertext ctIn by the plaintext matrix and returns the result on the Ciphertext opOut.
// ctInPreRotated can be obtained with GetPreRotatedCiphertextForDiagonalMatrixMultiplication.
// The BSGS approach is used (double hoisting with baby-step giant-step), which is faster than MultiplyByDiagMatrix
// for matrix with more than a few non-zero diagonals and uses significantly less keys.
func (eval defaultDiagonalMatrixEvaluator) MultiplyByDiagMatrixBSGS(ctIn *rlwe.Ciphertext, matrix he.LinearTransformation, ctPreRot map[int]*rlwe.Element[ringqp.Poly], opOut *rlwe.Ciphertext) (err error) {
	return he.MultiplyByDiagMatrixBSGS(eval.EvaluatorForLinearTransformation, ctIn, matrix, ctPreRot, opOut)
}
