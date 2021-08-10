package advanced

import (
	"github.com/ldsec/lattigo/v2/ckks"
	"github.com/ldsec/lattigo/v2/rlwe"
)

// Evaluator is an interface embeding the ckks.Evaluator interface with
// additional advanced arithmetic features.
type Evaluator interface {
	ckks.Evaluator
	// Homomorphic Encoding
	CoeffsToSlots(ctIn *ckks.Ciphertext, ctsMatrices EncodingMatrix) (ctReal, ctImag *ckks.Ciphertext)
	SlotsToCoeffs(ctReal, ctImag *ckks.Ciphertext, stcMatrices EncodingMatrix) (ctOut *ckks.Ciphertext)
	// Homomorphic Modular Reduction
	EvalMod(ctIn *ckks.Ciphertext, evalModPoly EvalModPoly) (ctOut *ckks.Ciphertext)
}

type evaluator struct {
	ckks.Evaluator
	params ckks.Parameters
}

// NewEvaluator creates a new Evaluator.
func NewEvaluator(params ckks.Parameters, evaluationKey rlwe.EvaluationKey) Evaluator {
	return &evaluator{ckks.NewEvaluator(params, evaluationKey), params}
}

// CoeffsToSlots applies the homomorphic encoding.
func (eval *evaluator) CoeffsToSlots(ctIn *ckks.Ciphertext, ctsMatrices EncodingMatrix) (ctReal, ctImag *ckks.Ciphertext) {

	var zV, zVconj *ckks.Ciphertext

	zV = eval.dft(ctIn, ctsMatrices.Matrices)

	zVconj = eval.ConjugateNew(zV)

	// The real part is stored in ct0
	ctReal = eval.AddNew(zV, zVconj)

	// The imaginary part is stored in ct1
	ctImag = eval.SubNew(zV, zVconj)

	eval.DivByi(ctImag, ctImag)

	// If repacking, then ct0 and ct1 right n/2 slots are zero.
	if eval.params.LogSlots() < eval.params.LogN()-1 {
		eval.Rotate(ctImag, eval.params.Slots(), ctImag)
		eval.Add(ctReal, ctImag, ctReal)
		return ctReal, nil
	}

	zV = nil
	zVconj = nil

	return ctReal, ctImag
}

// SlotsToCoeffs applies the homomorphic decoding.
func (eval *evaluator) SlotsToCoeffs(ctReal, ctImag *ckks.Ciphertext, stcMatrices EncodingMatrix) (ctOut *ckks.Ciphertext) {

	// If full packing, the repacking can be done directly using ct0 and ct1.
	if ctImag != nil {
		eval.MultByi(ctImag, ctImag)
		eval.Add(ctReal, ctImag, ctReal)
	}

	ctImag = nil

	return eval.dft(ctReal, stcMatrices.Matrices)
}

func (eval *evaluator) dft(vec *ckks.Ciphertext, plainVectors []ckks.PtDiagMatrix) *ckks.Ciphertext {

	// Sequentially multiplies w with the provided dft matrices.
	for _, plainVector := range plainVectors {
		scale := vec.Scale
		vec = eval.LinearTransform(vec, plainVector)[0]
		if err := eval.Rescale(vec, scale, vec); err != nil {
			panic(err)
		}
	}

	return vec
}
