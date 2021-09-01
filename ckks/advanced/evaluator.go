package advanced

import (
	"github.com/ldsec/lattigo/v2/ckks"
	"github.com/ldsec/lattigo/v2/rlwe"
	"math"
)

// Evaluator is an interface embeding the ckks.Evaluator interface with
// additional advanced arithmetic features.
type Evaluator interface {
	ckks.Evaluator
	// Homomorphic Encoding
	CoeffsToSlotsNew(ctIn *ckks.Ciphertext, ctsMatrices EncodingMatrix) (ctReal, ctImag *ckks.Ciphertext)
	SlotsToCoeffsNew(ctReal, ctImag *ckks.Ciphertext, stcMatrices EncodingMatrix) (ctOut *ckks.Ciphertext)
	// Homomorphic Modular Reduction
	EvalModNew(ctIn *ckks.Ciphertext, evalModPoly EvalModPoly) (ctOut *ckks.Ciphertext)
}

type evaluator struct {
	ckks.Evaluator
	params ckks.Parameters
}

// NewEvaluator creates a new Evaluator.
func NewEvaluator(params ckks.Parameters, evaluationKey rlwe.EvaluationKey) Evaluator {
	return &evaluator{ckks.NewEvaluator(params, evaluationKey), params}
}

// CoeffsToSlotsNew applies the homomorphic encoding.
func (eval *evaluator) CoeffsToSlotsNew(ctIn *ckks.Ciphertext, ctsMatrices EncodingMatrix) (ctReal, ctImag *ckks.Ciphertext) {

	var zV *ckks.Ciphertext

	zV = eval.dft(ctIn, ctsMatrices.matrices)

	ctReal = eval.ConjugateNew(zV)

	// Imaginary part
	ctImag = eval.SubNew(zV, ctReal)

	// Real part
	eval.Add(ctReal, zV, ctReal)

	eval.DivByi(ctImag, ctImag)

	// If repacking, then ct0 and ct1 right n/2 slots are zero.
	if eval.params.LogSlots() < eval.params.LogN()-1 {
		eval.Rotate(ctImag, eval.params.Slots(), ctImag)
		eval.Add(ctReal, ctImag, ctReal)
		return ctReal, nil
	}

	zV = nil

	return ctReal, ctImag
}

// SlotsToCoeffsNew applies the homomorphic decoding.
func (eval *evaluator) SlotsToCoeffsNew(ctReal, ctImag *ckks.Ciphertext, stcMatrices EncodingMatrix) (ctOut *ckks.Ciphertext) {

	// If full packing, the repacking can be done directly using ct0 and ct1.
	if ctImag != nil {
		eval.MultByi(ctImag, ctImag)
		eval.Add(ctReal, ctImag, ctReal)
	}

	ctImag = nil

	return eval.dft(ctReal, stcMatrices.matrices)
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

// EvalModNew does :
//
//	1) Delta * (Q/Delta * I(X) + m(X)) (Delta = scaling factor, I(X) integer poly, m(X) message)
//	2) Delta * (I(X) + Delta/Q * m(X)) (divide by Q/Delta)
//	3) Delta * (Delta/Q * m(X)) (x mod 1)
//	4) Delta * (m(X)) (multiply back by Q/Delta)
//
// Since Q is not a power of two, but Delta is, then does an approximate division by the closest
// power of two to Q instead. Hence, it assumes that the input plaintext is already scaled by
// the correcting factor Q/2^{round(log(Q))}.
//
// !! Assumes that the input is normalized by 1/K for K the range of the approximation.
//
// Scaling back error correction by 2^{round(log(Q))}/Q afterward is included in the polynomial
func (eval *evaluator) EvalModNew(ct *ckks.Ciphertext, evalModPoly EvalModPoly) *ckks.Ciphertext {

	// Stores default scales
	prevScaleCt := ct.Scale

	// Normalize the modular reduction to mod by 1 (division by Q)
	ct.Scale = evalModPoly.scalingFactor

	var err error

	// Compute the scales that the ciphertext should have before the double angle
	// formula such that after it it has the scale it had before the polynomial
	// evaluation
	targetScale := ct.Scale
	for i := 0; i < evalModPoly.doubleAngle; i++ {
		targetScale = math.Sqrt(targetScale * eval.params.QiFloat64(evalModPoly.levelStart-evalModPoly.sinePoly.Depth()-evalModPoly.doubleAngle+i+1))
	}

	// Division by 1/2^r and change of variable for the Chebysehev evaluation
	if evalModPoly.sineType == Cos1 || evalModPoly.sineType == Cos2 {
		eval.AddConst(ct, -0.5/(complex(evalModPoly.scFac, 0)*(evalModPoly.sinePoly.B-evalModPoly.sinePoly.A)), ct)
	}

	// Chebyshev evaluation
	if ct, err = eval.EvaluatePoly(ct, evalModPoly.sinePoly, targetScale); err != nil {
		panic(err)
	}

	// Double angle
	sqrt2pi := evalModPoly.sqrt2Pi
	for i := 0; i < evalModPoly.doubleAngle; i++ {
		sqrt2pi *= sqrt2pi
		eval.MulRelin(ct, ct, ct)
		eval.Add(ct, ct, ct)
		eval.AddConst(ct, -sqrt2pi, ct)
		if err := eval.Rescale(ct, targetScale, ct); err != nil {
			panic(err)
		}
	}

	// ArcSine
	if evalModPoly.arcSinePoly != nil {
		if ct, err = eval.EvaluatePoly(ct, evalModPoly.arcSinePoly, ct.Scale); err != nil {
			panic(err)
		}
	}

	// Multiplies back by q
	ct.Scale = prevScaleCt
	return ct
}
