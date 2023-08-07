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

func (e *floatEncoder[T, U]) Encode(values []T, metadata *rlwe.MetaData, output U) (err error) {
	return e.Encoder.Embed(values, metadata, output)
}

type Diagonals[T Float] circuits.Diagonals[T]

func (m Diagonals[T]) DiagonalsIndexList() (indexes []int) {
	return circuits.Diagonals[T](m).DiagonalsIndexList()
}

type LinearTransformationParameters circuits.LinearTransformationParameters

type LinearTransformation circuits.LinearTransformation

// NewLinearTransformationEvaluator instantiates a new LinearTransformationEvaluator from an EvaluatorForLinearTransformation.
// The method is allocation free if the underlying EvaluatorForLinearTransformation returns a non-nil
// *rlwe.EvaluatorBuffers.
func NewLinearTransformationEvaluator(eval circuits.EvaluatorForLinearTransformation) (linTransEval *LinearTransformationEvaluator) {
	return &LinearTransformationEvaluator{*circuits.NewLinearTransformationEvaluator(eval)}
}

func NewLinearTransformation(params rlwe.ParameterProvider, lt LinearTransformationParameters) LinearTransformation {
	return LinearTransformation(circuits.NewLinearTransformation(params, circuits.LinearTransformationParameters(lt)))
}

func EncodeLinearTransformation[T Float](params LinearTransformationParameters, ecd *ckks.Encoder, diagonals Diagonals[T], allocated LinearTransformation) (err error) {
	return circuits.EncodeLinearTransformation[T](
		circuits.LinearTransformationParameters(params),
		&floatEncoder[T, ringqp.Poly]{ecd},
		circuits.Diagonals[T](diagonals),
		circuits.LinearTransformation(allocated))
}

// GaloisElementsForLinearTransformation returns the list of Galois elements required to evaluate the linear transformation.
func GaloisElementsForLinearTransformation(params rlwe.ParameterProvider, lt LinearTransformationParameters) (galEls []uint64) {
	return circuits.GaloisElementsForLinearTransformation(params, lt.DiagonalsIndexList, 1<<lt.LogDimensions.Cols, lt.LogBabyStepGianStepRatio)
}

type LinearTransformationEvaluator struct {
	circuits.LinearTransformationEvaluator
}

// LinearTransformationsNew takes as input a ciphertext ctIn and a list of linear transformations [M0, M1, M2, ...] and returns opOut:[M0(ctIn), M1(ctIn), M2(ctInt), ...].
func (eval LinearTransformationEvaluator) LinearTransformationsNew(ctIn *rlwe.Ciphertext, linearTransformations []LinearTransformation) (opOut []*rlwe.Ciphertext, err error) {
	circuitLTs := make([]circuits.LinearTransformation, len(linearTransformations))
	for i := range circuitLTs {
		circuitLTs[i] = circuits.LinearTransformation(linearTransformations[i])
	}
	return eval.LinearTransformationEvaluator.LinearTransformationsNew(ctIn, circuitLTs)
}

// LinearTransformationNew takes as input a ciphertext ctIn and a linear transformation M and evaluate and returns opOut: M(ctIn).
func (eval LinearTransformationEvaluator) LinearTransformationNew(ctIn *rlwe.Ciphertext, linearTransformation LinearTransformation) (opOut *rlwe.Ciphertext, err error) {
	cts, err := eval.LinearTransformationsNew(ctIn, []LinearTransformation{linearTransformation})
	return cts[0], err
}

// LinearTransformation takes as input a ciphertext ctIn, a linear transformation M and evaluates opOut: M(ctIn).
func (eval LinearTransformationEvaluator) LinearTransformation(ctIn *rlwe.Ciphertext, linearTransformation LinearTransformation, opOut *rlwe.Ciphertext) (err error) {
	return eval.LinearTransformations(ctIn, []LinearTransformation{linearTransformation}, []*rlwe.Ciphertext{opOut})
}

// LinearTransformations takes as input a ciphertext ctIn, a list of linear transformations [M0, M1, M2, ...] and a list of pre-allocated receiver opOut
// and evaluates opOut: [M0(ctIn), M1(ctIn), M2(ctIn), ...]
func (eval LinearTransformationEvaluator) LinearTransformations(ctIn *rlwe.Ciphertext, linearTransformations []LinearTransformation, opOut []*rlwe.Ciphertext) (err error) {
	circuitLTs := make([]circuits.LinearTransformation, len(linearTransformations))
	for i := range circuitLTs {
		circuitLTs[i] = circuits.LinearTransformation(linearTransformations[i])
	}
	return eval.LinearTransformationEvaluator.LinearTransformations(ctIn, circuitLTs, opOut)
}

// MultiplyByDiagMatrix multiplies the Ciphertext "ctIn" by the plaintext matrix "matrix" and returns the result on the Ciphertext
// "opOut". Memory buffers for the decomposed ciphertext BuffDecompQP, BuffDecompQP must be provided, those are list of poly of ringQ and ringP
// respectively, each of size params.Beta().
// The naive approach is used (single hoisting and no baby-step giant-step), which is faster than MultiplyByDiagMatrixBSGS
// for matrix of only a few non-zero diagonals but uses more keys.
func (eval LinearTransformationEvaluator) MultiplyByDiagMatrix(ctIn *rlwe.Ciphertext, matrix LinearTransformation, BuffDecompQP []ringqp.Poly, opOut *rlwe.Ciphertext) (err error) {
	return eval.LinearTransformationEvaluator.MultiplyByDiagMatrix(ctIn, circuits.LinearTransformation(matrix), BuffDecompQP, opOut)
}

// MultiplyByDiagMatrixBSGS multiplies the Ciphertext "ctIn" by the plaintext matrix "matrix" and returns the result on the Ciphertext
// "opOut". Memory buffers for the decomposed Ciphertext BuffDecompQP, BuffDecompQP must be provided, those are list of poly of ringQ and ringP
// respectively, each of size params.Beta().
// The BSGS approach is used (double hoisting with baby-step giant-step), which is faster than MultiplyByDiagMatrix
// for matrix with more than a few non-zero diagonals and uses significantly less keys.
func (eval LinearTransformationEvaluator) MultiplyByDiagMatrixBSGS(ctIn *rlwe.Ciphertext, matrix LinearTransformation, BuffDecompQP []ringqp.Poly, opOut *rlwe.Ciphertext) (err error) {
	return eval.LinearTransformationEvaluator.MultiplyByDiagMatrixBSGS(ctIn, circuits.LinearTransformation(matrix), BuffDecompQP, opOut)
}
