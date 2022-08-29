// Package advanced implements advanced operations for the CKKS scheme.
package advanced

import (
	"math"

	"github.com/tuneinsight/lattigo/v3/ckks"
	"github.com/tuneinsight/lattigo/v3/ring"
	"github.com/tuneinsight/lattigo/v3/rlwe"
	"github.com/tuneinsight/lattigo/v3/rlwe/ringqp"
)

// Evaluator is an interface embeding the ckks.Evaluator interface with
// additional advanced arithmetic features.
type Evaluator interface {

	// =======================================
	// === Original ckks.Evaluator methods ===
	// =======================================

	Add(ctIn *ckks.Ciphertext, op1 ckks.Operand, ctOut *ckks.Ciphertext)
	AddNoMod(ctIn *ckks.Ciphertext, op1 ckks.Operand, ctOut *ckks.Ciphertext)
	AddNew(ctIn *ckks.Ciphertext, op1 ckks.Operand) (ctOut *ckks.Ciphertext)
	AddNoModNew(ctIn *ckks.Ciphertext, op1 ckks.Operand) (ctOut *ckks.Ciphertext)
	Sub(ctIn *ckks.Ciphertext, op1 ckks.Operand, ctOut *ckks.Ciphertext)
	SubNoMod(ctIn *ckks.Ciphertext, op1 ckks.Operand, ctOut *ckks.Ciphertext)
	SubNew(ctIn *ckks.Ciphertext, op1 ckks.Operand) (ctOut *ckks.Ciphertext)
	SubNoModNew(ctIn *ckks.Ciphertext, op1 ckks.Operand) (ctOut *ckks.Ciphertext)
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
	Mul(ctIn *ckks.Ciphertext, op1 ckks.Operand, ctOut *ckks.Ciphertext)
	MulNew(ctIn *ckks.Ciphertext, op1 ckks.Operand) (ctOut *ckks.Ciphertext)
	MulRelin(ctIn *ckks.Ciphertext, op1 ckks.Operand, ctOut *ckks.Ciphertext)
	MulRelinNew(ctIn *ckks.Ciphertext, op1 ckks.Operand) (ctOut *ckks.Ciphertext)
	RotateNew(ctIn *ckks.Ciphertext, k int) (ctOut *ckks.Ciphertext)
	Rotate(ctIn *ckks.Ciphertext, k int, ctOut *ckks.Ciphertext)
	RotateHoistedNew(ctIn *ckks.Ciphertext, rotations []int) (ctOut map[int]*ckks.Ciphertext)
	RotateHoisted(ctIn *ckks.Ciphertext, rotations []int, ctOut map[int]*ckks.Ciphertext)
	MulByPow2New(ctIn *ckks.Ciphertext, pow2 int) (ctOut *ckks.Ciphertext)
	MulByPow2(ctIn *ckks.Ciphertext, pow2 int, ctOut *ckks.Ciphertext)
	PowerOf2(ctIn *ckks.Ciphertext, logPow2 int, ctOut *ckks.Ciphertext)
	Power(ctIn *ckks.Ciphertext, degree int, ctOut *ckks.Ciphertext)
	PowerNew(ctIn *ckks.Ciphertext, degree int) (ctOut *ckks.Ciphertext)
	EvaluatePoly(input interface{}, pol *ckks.Polynomial, targetScale float64) (ctOut *ckks.Ciphertext, err error)
	EvaluatePolyVector(input interface{}, pols []*ckks.Polynomial, encoder ckks.Encoder, slotIndex map[int][]int, targetScale float64) (ctOut *ckks.Ciphertext, err error)
	InverseNew(ctIn *ckks.Ciphertext, steps int) (ctOut *ckks.Ciphertext)
	LinearTransformNew(ctIn *ckks.Ciphertext, linearTransform interface{}) (ctOut []*ckks.Ciphertext)
	LinearTransform(ctIn *ckks.Ciphertext, linearTransform interface{}, ctOut []*ckks.Ciphertext)
	MultiplyByDiagMatrix(ctIn *ckks.Ciphertext, matrix ckks.LinearTransform, c2DecompQP []ringqp.Poly, ctOut *ckks.Ciphertext)
	MultiplyByDiagMatrixBSGS(ctIn *ckks.Ciphertext, matrix ckks.LinearTransform, c2DecompQP []ringqp.Poly, ctOut *ckks.Ciphertext)
	InnerSumLog(ctIn *ckks.Ciphertext, batch, n int, ctOut *ckks.Ciphertext)
	InnerSum(ctIn *ckks.Ciphertext, batch, n int, ctOut *ckks.Ciphertext)
	ReplicateLog(ctIn *ckks.Ciphertext, batch, n int, ctOut *ckks.Ciphertext)
	Replicate(ctIn *ckks.Ciphertext, batch, n int, ctOut *ckks.Ciphertext)
	TraceNew(ctIn *ckks.Ciphertext, logSlots int) *ckks.Ciphertext
	Trace(ctIn *ckks.Ciphertext, logSlots int, ctOut *ckks.Ciphertext)
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
	CoeffsToSlots(ctIn *ckks.Ciphertext, ctsMatrices EncodingMatrix, ctReal, ctImag *ckks.Ciphertext)
	SlotsToCoeffsNew(ctReal, ctImag *ckks.Ciphertext, stcMatrices EncodingMatrix) (ctOut *ckks.Ciphertext)
	SlotsToCoeffs(ctReal, ctImag *ckks.Ciphertext, stcMatrices EncodingMatrix, ctOut *ckks.Ciphertext)
	EvalModNew(ctIn *ckks.Ciphertext, evalModPoly EvalModPoly) (ctOut *ckks.Ciphertext)

	// =================================================
	// === original ckks.Evaluator redefined methods ===
	// =================================================

	GetRLWEEvaluator() *rlwe.Evaluator
	BuffQ() [3]*ring.Poly
	BuffCt() *ckks.Ciphertext
	ShallowCopy() Evaluator
	WithKey(rlwe.EvaluationKey) Evaluator
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

// CoeffsToSlotsNew applies the homomorphic encoding and returns the result on new ciphertexts.
// Homomorphically encodes a complex vector vReal + i*vImag.
// If the packing is sparse (n < N/2), then returns ctReal = Ecd(vReal || vImag) and ctImag = nil.
// If the packing is dense (n == N/2), then returns ctReal = Ecd(vReal) and ctImag = Ecd(vImag).
func (eval *evaluator) CoeffsToSlotsNew(ctIn *ckks.Ciphertext, ctsMatrices EncodingMatrix) (ctReal, ctImag *ckks.Ciphertext) {
	ctReal = ckks.NewCiphertext(eval.params, 1, ctsMatrices.LevelStart, 0)

	if eval.params.LogSlots() == eval.params.LogN()-1 {
		ctImag = ckks.NewCiphertext(eval.params, 1, ctsMatrices.LevelStart, 0)
	}

	eval.CoeffsToSlots(ctIn, ctsMatrices, ctReal, ctImag)
	return
}

// CoeffsToSlots applies the homomorphic encoding and returns the results on the provided ciphertexts.
// Homomorphically encodes a complex vector vReal + i*vImag of size n on a real vector of size 2n.
// If the packing is sparse (n < N/2), then returns ctReal = Ecd(vReal || vImag) and ctImag = nil.
// If the packing is dense (n == N/2), then returns ctReal = Ecd(vReal) and ctImag = Ecd(vImag).
func (eval *evaluator) CoeffsToSlots(ctIn *ckks.Ciphertext, ctsMatrices EncodingMatrix, ctReal, ctImag *ckks.Ciphertext) {

	if ctsMatrices.RepackImag2Real {

		zV := ctIn.CopyNew()

		eval.dft(ctIn, ctsMatrices.matrices, zV)

		eval.Conjugate(zV, ctReal)

		var tmp *ckks.Ciphertext
		if ctImag != nil {
			tmp = ctImag
		} else {
			tmp = ckks.NewCiphertextAtLevelFromPoly(ctReal.Level(), [2]*ring.Poly{eval.BuffCt().Value[0], eval.BuffCt().Value[1]})
		}

		// Imag part
		eval.Sub(zV, ctReal, tmp)
		eval.DivByi(tmp, tmp)

		// Real part
		eval.Add(ctReal, zV, ctReal)

		// If repacking, then ct0 and ct1 right n/2 slots are zero.
		if eval.params.LogSlots() < eval.params.LogN()-1 {
			eval.Rotate(tmp, eval.params.Slots(), tmp)
			eval.Add(ctReal, tmp, ctReal)
		}

		zV = nil
	} else {
		eval.dft(ctIn, ctsMatrices.matrices, ctReal)
	}
}

// SlotsToCoeffsNew applies the homomorphic decoding and returns the result on a new ciphertext.
// Homomorphically decodes a real vector of size 2n on a complex vector vReal + i*vImag of size n.
// If the packing is sparse (n < N/2) then ctReal = Ecd(vReal || vImag) and ctImag = nil.
// If the packing is dense (n == N/2), then ctReal = Ecd(vReal) and ctImag = Ecd(vImag).
func (eval *evaluator) SlotsToCoeffsNew(ctReal, ctImag *ckks.Ciphertext, stcMatrices EncodingMatrix) (ctOut *ckks.Ciphertext) {

	if ctReal.Level() < stcMatrices.LevelStart || (ctImag != nil && ctImag.Level() < stcMatrices.LevelStart) {
		panic("ctReal.Level() or ctImag.Level() < EncodingMatrix.LevelStart")
	}

	ctOut = ckks.NewCiphertext(eval.params, 1, stcMatrices.LevelStart, ctReal.Scale())
	eval.SlotsToCoeffs(ctReal, ctImag, stcMatrices, ctOut)
	return

}

// SlotsToCoeffsNew applies the homomorphic decoding and returns the result on the provided ciphertext.
// Homomorphically decodes a real vector of size 2n on a complex vector vReal + i*vImag of size n.
// If the packing is sparse (n < N/2) then ctReal = Ecd(vReal || vImag) and ctImag = nil.
// If the packing is dense (n == N/2), then ctReal = Ecd(vReal) and ctImag = Ecd(vImag).
func (eval *evaluator) SlotsToCoeffs(ctReal, ctImag *ckks.Ciphertext, stcMatrices EncodingMatrix, ctOut *ckks.Ciphertext) {
	// If full packing, the repacking can be done directly using ct0 and ct1.
	if ctImag != nil {
		eval.MultByi(ctImag, ctOut)
		eval.Add(ctOut, ctReal, ctOut)
		eval.dft(ctOut, stcMatrices.matrices, ctOut)
	} else {
		eval.dft(ctReal, stcMatrices.matrices, ctOut)
	}
}

func (eval *evaluator) dft(ctIn *ckks.Ciphertext, plainVectors []ckks.LinearTransform, ctOut *ckks.Ciphertext) {
	// Sequentially multiplies w with the provided dft matrices.
	scale := ctIn.Scale()
	var in, out *ckks.Ciphertext
	for i, plainVector := range plainVectors {
		in, out = ctOut, ctOut
		if i == 0 {
			in, out = ctIn, ctOut
		}
		eval.LinearTransform(in, plainVector, []*ckks.Ciphertext{out})
		if err := eval.Rescale(out, scale, out); err != nil {
			panic(err)
		}
	}
}

// EvalModNew applies a homomorphic mod Q on a vector scaled by Delta, scaled down to mod 1 :
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

	if ct.Level() < evalModPoly.LevelStart() {
		panic("ct.Level() < evalModPoly.LevelStart")
	}

	if ct.Level() > evalModPoly.LevelStart() {
		eval.DropLevel(ct, ct.Level()-evalModPoly.LevelStart())
	}

	// Stores default scales
	prevScaleCt := ct.Scale()

	// Normalize the modular reduction to mod by 1 (division by Q)
	ct.SetScale(evalModPoly.scalingFactor)

	var err error

	// Compute the scales that the ciphertext should have before the double angle
	// formula such that after it it has the scale it had before the polynomial
	// evaluation
	targetScale := ct.Scale()
	for i := 0; i < evalModPoly.doubleAngle; i++ {
		targetScale = math.Sqrt(targetScale * eval.params.QiFloat64(evalModPoly.levelStart-evalModPoly.sinePoly.Depth()-evalModPoly.doubleAngle+i+1))
	}

	// Division by 1/2^r and change of variable for the Chebysehev evaluation
	if evalModPoly.sineType == Cos1 || evalModPoly.sineType == Cos2 {
		eval.AddConst(ct, -0.5/(evalModPoly.scFac*(evalModPoly.sinePoly.B-evalModPoly.sinePoly.A)), ct)
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
		if ct, err = eval.EvaluatePoly(ct, evalModPoly.arcSinePoly, ct.Scale()); err != nil {
			panic(err)
		}
	}

	// Multiplies back by q
	ct.SetScale(prevScaleCt)
	return ct
}
