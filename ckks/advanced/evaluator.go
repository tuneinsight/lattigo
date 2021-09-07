package advanced

import (
	"math"

	"github.com/ldsec/lattigo/v2/ckks"
	"github.com/ldsec/lattigo/v2/rlwe"
)

// Evaluator is an interface embeding the ckks.Evaluator interface with
// additional advanced arithmetic features.
type Evaluator interface {

	// =======================================
	// === Original ckks.Evaluator methods ===
	// =======================================

	Add(op0, op1 ckks.Operand, ctOut *ckks.Ciphertext)
	AddNoMod(op0, op1 ckks.Operand, ctOut *ckks.Ciphertext)
	AddNew(op0, op1 ckks.Operand) (ctOut *ckks.Ciphertext)
	AddNoModNew(op0, op1 ckks.Operand) (ctOut *ckks.Ciphertext)
	Sub(op0, op1 ckks.Operand, ctOut *ckks.Ciphertext)
	SubNoMod(op0, op1 ckks.Operand, ctOut *ckks.Ciphertext)
	SubNew(op0, op1 ckks.Operand) (ctOut *ckks.Ciphertext)
	SubNoModNew(op0, op1 ckks.Operand) (ctOut *ckks.Ciphertext)
	Neg(ctIn *ckks.Ciphertext, ctOut *ckks.Ciphertext)
	NegNew(ctIn *ckks.Ciphertext) (ctOut *ckks.Ciphertext)
	AddConstNew(ctIn *ckks.Ciphertext, constant interface{}) (ctOut *ckks.Ciphertext)
	AddConst(ctIn *ckks.Ciphertext, constant interface{}, ctOut *ckks.Ciphertext)
	MultByConstNew(ctIn *ckks.Ciphertext, constant interface{}) (ctOut *ckks.Ciphertext)
	MultByConst(ctIn *ckks.Ciphertext, constant interface{}, ctOut *ckks.Ciphertext)
	MultByGaussianInteger(ctIn *ckks.Ciphertext, cReal, cImag interface{}, ctOut *ckks.Ciphertext)
	MultByConstAndAdd(ctIn *ckks.Ciphertext, constant interface{}, ctOut *ckks.Ciphertext)
	MultByGaussianIntegerAndAdd(ctIn *ckks.Ciphertext, cReal, cImag interface{}, ctOut *ckks.Ciphertext)
	MultByiNew(ctIn *ckks.Ciphertext) (ctOut *ckks.Ciphertext)
	MultByi(ctIn *ckks.Ciphertext, ctOut *ckks.Ciphertext)
	DivByiNew(ctIn *ckks.Ciphertext) (ctOut *ckks.Ciphertext)
	DivByi(ctIn *ckks.Ciphertext, ctOut *ckks.Ciphertext)
	ConjugateNew(ctIn *ckks.Ciphertext) (ctOut *ckks.Ciphertext)
	Conjugate(ctIn *ckks.Ciphertext, ctOut *ckks.Ciphertext)
	Mul(op0, op1 ckks.Operand, ctOut *ckks.Ciphertext)
	MulNew(op0, op1 ckks.Operand) (ctOut *ckks.Ciphertext)
	MulRelin(op0, op1 ckks.Operand, ctOut *ckks.Ciphertext)
	MulRelinNew(op0, op1 ckks.Operand) (ctOut *ckks.Ciphertext)
	RotateNew(ctIn *ckks.Ciphertext, k int) (ctOut *ckks.Ciphertext)
	Rotate(ctIn *ckks.Ciphertext, k int, ctOut *ckks.Ciphertext)
	RotateHoistedNew(ctIn *ckks.Ciphertext, rotations []int) (ctOut map[int]*ckks.Ciphertext)
	RotateHoisted(ctIn *ckks.Ciphertext, rotations []int, ctOut map[int]*ckks.Ciphertext)
	MulByPow2New(ctIn *ckks.Ciphertext, pow2 int) (ctOut *ckks.Ciphertext)
	MulByPow2(ctIn *ckks.Ciphertext, pow2 int, ctOut *ckks.Ciphertext)
	PowerOf2(ctIn *ckks.Ciphertext, logPow2 int, ctOut *ckks.Ciphertext)
	Power(ctIn *ckks.Ciphertext, degree int, ctOut *ckks.Ciphertext)
	PowerNew(ctIn *ckks.Ciphertext, degree int) (ctOut *ckks.Ciphertext)
	EvaluatePoly(ctIn *ckks.Ciphertext, pol *ckks.Polynomial, targetScale float64) (ctOut *ckks.Ciphertext, err error)
	InverseNew(ctIn *ckks.Ciphertext, steps int) (ctOut *ckks.Ciphertext)
	LinearTransformNew(ctIn *ckks.Ciphertext, linearTransform interface{}) (ctOut []*ckks.Ciphertext)
	LinearTransform(ctIn *ckks.Ciphertext, linearTransform interface{}, ctOut []*ckks.Ciphertext)
	MultiplyByDiagMatrix(ctIn *ckks.Ciphertext, matrix ckks.PtDiagMatrix, c2DecompQP []rlwe.PolyQP, ctOut *ckks.Ciphertext)
	MultiplyByDiagMatrixBSGS(ctIn *ckks.Ciphertext, matrix ckks.PtDiagMatrix, c2DecompQP []rlwe.PolyQP, ctOut *ckks.Ciphertext)
	InnerSumLog(ctIn *ckks.Ciphertext, batch, n int, ctOut *ckks.Ciphertext)
	InnerSum(ctIn *ckks.Ciphertext, batch, n int, ctOut *ckks.Ciphertext)
	ReplicateLog(ctIn *ckks.Ciphertext, batch, n int, ctOut *ckks.Ciphertext)
	Replicate(ctIn *ckks.Ciphertext, batch, n int, ctOut *ckks.Ciphertext)
	TraceNew(ctIn *ckks.Ciphertext, logSlotsStart, logSlotsEnd int) *ckks.Ciphertext
	Trace(ctIn *ckks.Ciphertext, logSlotsStart, logSlotsEnd int, ctOut *ckks.Ciphertext)
	SwitchKeysNew(ctIn *ckks.Ciphertext, switchingKey *rlwe.SwitchingKey) (ctOut *ckks.Ciphertext)
	SwitchKeys(ctIn *ckks.Ciphertext, switchingKey *rlwe.SwitchingKey, ctOut *ckks.Ciphertext)
	RelinearizeNew(ctIn *ckks.Ciphertext) (ctOut *ckks.Ciphertext)
	Relinearize(ctIn *ckks.Ciphertext, ctOut *ckks.Ciphertext)
	ScaleUpNew(ctIn *ckks.Ciphertext, scale float64) (ctOut *ckks.Ciphertext)
	ScaleUp(ctIn *ckks.Ciphertext, scale float64, ctOut *ckks.Ciphertext)
	SetScale(ctIn *ckks.Ciphertext, scale float64)
	Rescale(ctIn *ckks.Ciphertext, minScale float64, ctOut *ckks.Ciphertext) (err error)
	DropLevelNew(ctIn *ckks.Ciphertext, levels int) (ctOut *ckks.Ciphertext)
	DropLevel(ctIn *ckks.Ciphertext, levels int)
	ReduceNew(ctIn *ckks.Ciphertext) (ctOut *ckks.Ciphertext)
	Reduce(ctIn *ckks.Ciphertext, ctOut *ckks.Ciphertext) error

	// ======================================
	// === advanced.Evaluator new methods ===
	// ======================================

	CoeffsToSlotsNew(ctIn *ckks.Ciphertext, ctsMatrices EncodingMatrix) (ctReal, ctImag *ckks.Ciphertext)
	SlotsToCoeffsNew(ctReal, ctImag *ckks.Ciphertext, stcMatrices EncodingMatrix) (ctOut *ckks.Ciphertext)
	EvalModNew(ctIn *ckks.Ciphertext, evalModPoly EvalModPoly) (ctOut *ckks.Ciphertext)

	// =================================================
	// === original ckks.Evaluator redefined methods ===
	// =================================================

	ShallowCopy() Evaluator
	WithKey(evakey rlwe.EvaluationKey) Evaluator
}

type evaluator struct {
	ckks.Evaluator
	params ckks.Parameters
}

// NewEvaluator creates a new Evaluator.
func NewEvaluator(params ckks.Parameters, evaluationKey rlwe.EvaluationKey) Evaluator {
	return &evaluator{ckks.NewEvaluator(params, evaluationKey), params}
}

// ShallowCopy creates a shallow copy of this evaluator in which all the read-only data-structures are
// shared with the receiver and the temporary buffers are reallocated. The receiver and the returned
// Evaluators can be used concurrently.
func (eval *evaluator) ShallowCopy() Evaluator {
	return &evaluator{eval.Evaluator.ShallowCopy(), eval.params}
}

// WithKey creates a shallow copy of the receiver Evaluator for which the new EvaluationKey is evaluationKey
// and where the temporary buffers are shared. The receiver and the returned Evaluators cannot be used concurrently.
func (eval *evaluator) WithKey(evaluationKey rlwe.EvaluationKey) Evaluator {
	return &evaluator{eval.Evaluator.WithKey(evaluationKey), eval.params}
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
		vec = eval.LinearTransformNew(vec, plainVector)[0]
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
