package float

import (
	"github.com/tuneinsight/lattigo/v4/circuits"
	"github.com/tuneinsight/lattigo/v4/ckks"
	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/rlwe"
	"github.com/tuneinsight/lattigo/v4/rlwe/ringqp"
)

type floatEncoder[T Float, U ring.Poly | ringqp.Poly | *rlwe.Plaintext] struct {
	*ckks.Encoder
}

func (e floatEncoder[T, U]) Encode(values []T, metadata *rlwe.MetaData, output U) (err error) {
	return e.Encoder.Embed(values, metadata, output)
}

// Diagonals is a wrapper of circuits.Diagonals.
// See circuits.Diagonals for the documentation.
type Diagonals[T Float] circuits.Diagonals[T]

// DiagonalsIndexList returns the list of the non-zero diagonals of the square matrix.
// A non zero diagonals is a diagonal with a least one non-zero element.
func (m Diagonals[T]) DiagonalsIndexList() (indexes []int) {
	return circuits.Diagonals[T](m).DiagonalsIndexList()
}

// LinearTransformationParameters is a wrapper of circuits.LinearTransformationParameters.
// See circuits.LinearTransformationParameters for the documentation.
type LinearTransformationParameters circuits.LinearTransformationParameters

// LinearTransformation is a wrapper of circuits.LinearTransformation.
// See circuits.LinearTransformation for the documentation.
type LinearTransformation circuits.LinearTransformation

// GaloisElements returns the list of Galois elements required to evaluate the linear transformation.
func (lt LinearTransformation) GaloisElements(params rlwe.ParameterProvider) []uint64 {
	return circuits.LinearTransformation(lt).GaloisElements(params)
}

// NewLinearTransformation instantiates a new LinearTransformation and is a wrapper of circuits.LinearTransformation.
// See circuits.LinearTransformation for the documentation.
func NewLinearTransformation(params rlwe.ParameterProvider, lt LinearTransformationParameters) LinearTransformation {
	return LinearTransformation(circuits.NewLinearTransformation(params, circuits.LinearTransformationParameters(lt)))
}

// EncodeLinearTransformation is a method used to encode EncodeLinearTransformation and a wrapper of circuits.EncodeLinearTransformation.
// See circuits.EncodeLinearTransformation for the documentation.
func EncodeLinearTransformation[T Float](ecd *ckks.Encoder, diagonals Diagonals[T], allocated LinearTransformation) (err error) {
	return circuits.EncodeLinearTransformation[T](
		&floatEncoder[T, ringqp.Poly]{ecd},
		circuits.Diagonals[T](diagonals),
		circuits.LinearTransformation(allocated))
}

// GaloisElementsForLinearTransformation returns the list of Galois elements required to evaluate the linear transformation.
func GaloisElementsForLinearTransformation(params rlwe.ParameterProvider, lt LinearTransformationParameters) (galEls []uint64) {
	return circuits.GaloisElementsForLinearTransformation(params, lt.DiagonalsIndexList, 1<<lt.LogDimensions.Cols, lt.LogBabyStepGianStepRatio)
}

// LinearTransformationEvaluator is an evaluator providing an API to evaluate linear transformations on rlwe.Ciphertexts.
// All fields of this struct are public, enabling custom instantiations.
type LinearTransformationEvaluator struct {
	circuits.EvaluatorForLinearTransformation
	circuits.EvaluatorForDiagonalMatrix
}

// NewLinearTransformationEvaluator instantiates a new LinearTransformationEvaluator from a circuit.EvaluatorForLinearTransformation.
// The default ckks.Evaluator is compliant to the circuits.EvaluatorForLinearTransformation interface.
// This method is allocation free.
func NewLinearTransformationEvaluator(eval circuits.EvaluatorForLinearTransformation) (linTransEval *LinearTransformationEvaluator) {
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
	return circuits.EvaluateLinearTransformationsMany(eval.EvaluatorForLinearTransformation, eval.EvaluatorForDiagonalMatrix, ctIn, []circuits.LinearTransformation{circuits.LinearTransformation(linearTransformation)}, []*rlwe.Ciphertext{opOut})
}

// EvaluateManyNew takes as input a ciphertext ctIn and a list of linear transformations [M0, M1, M2, ...] and returns opOut:[M0(ctIn), M1(ctIn), M2(ctInt), ...].
func (eval LinearTransformationEvaluator) EvaluateManyNew(ctIn *rlwe.Ciphertext, linearTransformations []LinearTransformation) (opOut []*rlwe.Ciphertext, err error) {
	params := eval.GetRLWEParameters()
	opOut = make([]*rlwe.Ciphertext, len(linearTransformations))
	for i := range opOut {
		opOut[i] = rlwe.NewCiphertext(params, 1, linearTransformations[i].Level)
	}
	return opOut, eval.EvaluateMany(ctIn, linearTransformations, opOut)
}

// EvaluateMany takes as input a ciphertext ctIn, a list of linear transformations [M0, M1, M2, ...] and a list of pre-allocated receiver opOut
// and evaluates opOut: [M0(ctIn), M1(ctIn), M2(ctIn), ...]
func (eval LinearTransformationEvaluator) EvaluateMany(ctIn *rlwe.Ciphertext, linearTransformations []LinearTransformation, opOut []*rlwe.Ciphertext) (err error) {
	circuitLTs := make([]circuits.LinearTransformation, len(linearTransformations))
	for i := range circuitLTs {
		circuitLTs[i] = circuits.LinearTransformation(linearTransformations[i])
	}
	return circuits.EvaluateLinearTransformationsMany(eval.EvaluatorForLinearTransformation, eval.EvaluatorForDiagonalMatrix, ctIn, circuitLTs, opOut)
}

// EvaluateSequentialNew takes as input a ciphertext ctIn and a list of linear transformations [M0, M1, M2, ...] and returns opOut:...M2(M1(M0(ctIn))
func (eval LinearTransformationEvaluator) EvaluateSequentialNew(ctIn *rlwe.Ciphertext, linearTransformations []LinearTransformation) (opOut *rlwe.Ciphertext, err error) {
	opOut = rlwe.NewCiphertext(eval.GetRLWEParameters(), 1, linearTransformations[0].Level)
	return opOut, eval.EvaluateSequential(ctIn, linearTransformations, opOut)
}

// EvaluateSequential takes as input a ciphertext ctIn and a list of linear transformations [M0, M1, M2, ...] and returns opOut:...M2(M1(M0(ctIn))
func (eval LinearTransformationEvaluator) EvaluateSequential(ctIn *rlwe.Ciphertext, linearTransformations []LinearTransformation, opOut *rlwe.Ciphertext) (err error) {
	circuitLTs := make([]circuits.LinearTransformation, len(linearTransformations))
	for i := range circuitLTs {
		circuitLTs[i] = circuits.LinearTransformation(linearTransformations[i])
	}
	return circuits.EvaluateLinearTranformationSequential(eval.EvaluatorForLinearTransformation, eval.EvaluatorForDiagonalMatrix, ctIn, circuitLTs, opOut)
}

// defaultDiagonalMatrixEvaluator is a struct implementing the interface circuits.EvaluatorForDiagonalMatrix.
type defaultDiagonalMatrixEvaluator struct {
	circuits.EvaluatorForLinearTransformation
}

// Decompose applies the RNS decomposition on ct[1] at the given level and stores the result in BuffDecompQP.
func (eval defaultDiagonalMatrixEvaluator) Decompose(level int, ct *rlwe.Ciphertext, BuffDecompQP []ringqp.Poly) {
	params := eval.GetRLWEParameters()
	eval.DecomposeNTT(level, params.MaxLevelP(), params.PCount(), ct.Value[1], ct.IsNTT, BuffDecompQP)
}

// GetPreRotatedCiphertextForDiagonalMatrixMultiplication populates ctPreRot with the pre-rotated ciphertext for the rotations rots and deletes rotated ciphertexts that are not in rots.
func (eval defaultDiagonalMatrixEvaluator) GetPreRotatedCiphertextForDiagonalMatrixMultiplication(levelQ int, ctIn *rlwe.Ciphertext, BuffDecompQP []ringqp.Poly, rots []int, ctPreRot map[int]*rlwe.Element[ringqp.Poly]) (err error) {
	return circuits.GetPreRotatedCiphertextForDiagonalMatrixMultiplication(levelQ, eval, ctIn, BuffDecompQP, rots, ctPreRot)
}

// MultiplyByDiagMatrix multiplies the Ciphertext "ctIn" by the plaintext matrix "matrix" and returns the result on the Ciphertext
// "opOut". Memory buffers for the decomposed ciphertext BuffDecompQP, BuffDecompQP must be provided, those are list of poly of ringQ and ringP
// respectively, each of size params.Beta().
// The naive approach is used (single hoisting and no baby-step giant-step), which is faster than MultiplyByDiagMatrixBSGS
// for matrix of only a few non-zero diagonals but uses more keys.
func (eval defaultDiagonalMatrixEvaluator) MultiplyByDiagMatrix(ctIn *rlwe.Ciphertext, matrix circuits.LinearTransformation, BuffDecompQP []ringqp.Poly, opOut *rlwe.Ciphertext) (err error) {
	return circuits.MultiplyByDiagMatrix(eval.EvaluatorForLinearTransformation, ctIn, matrix, BuffDecompQP, opOut)
}

// MultiplyByDiagMatrixBSGS multiplies the Ciphertext "ctIn" by the plaintext matrix "matrix" and returns the result on the Ciphertext "opOut".
// ctInPreRotated can be obtained with GetPreRotatedCiphertextForDiagonalMatrixMultiplication.
// The BSGS approach is used (double hoisting with baby-step giant-step), which is faster than MultiplyByDiagMatrix
// for matrix with more than a few non-zero diagonals and uses significantly less keys.
func (eval defaultDiagonalMatrixEvaluator) MultiplyByDiagMatrixBSGS(ctIn *rlwe.Ciphertext, matrix circuits.LinearTransformation, ctPreRot map[int]*rlwe.Element[ringqp.Poly], opOut *rlwe.Ciphertext) (err error) {
	return circuits.MultiplyByDiagMatrixBSGS(eval.EvaluatorForLinearTransformation, ctIn, matrix, ctPreRot, opOut)
}
