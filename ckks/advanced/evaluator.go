// Package advanced implements advanced operations for the CKKS scheme.
package advanced

import (
	"math"

	"github.com/tuneinsight/lattigo/v4/ckks"
	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/rlwe"
	"github.com/tuneinsight/lattigo/v4/rlwe/ringqp"
)

// Evaluator is an interface embedding the ckks.Evaluator interface with
// additional advanced arithmetic features.
type Evaluator interface {

	// =======================================
	// === Original ckks.Evaluator methods ===
	// =======================================

	Add(ctIn *rlwe.Ciphertext, op1 rlwe.Operand, ctOut *rlwe.Ciphertext)
	AddNew(ctIn *rlwe.Ciphertext, op1 rlwe.Operand) (ctOut *rlwe.Ciphertext)
	Sub(ctIn *rlwe.Ciphertext, op1 rlwe.Operand, ctOut *rlwe.Ciphertext)
	SubNew(ctIn *rlwe.Ciphertext, op1 rlwe.Operand) (ctOut *rlwe.Ciphertext)
	Neg(ctIn *rlwe.Ciphertext, ctOut *rlwe.Ciphertext)
	NegNew(ctIn *rlwe.Ciphertext) (ctOut *rlwe.Ciphertext)
	AddConstNew(ctIn *rlwe.Ciphertext, constant interface{}) (ctOut *rlwe.Ciphertext)
	AddConst(ctIn *rlwe.Ciphertext, constant interface{}, ctOut *rlwe.Ciphertext)
	MultByConstNew(ctIn *rlwe.Ciphertext, constant interface{}) (ctOut *rlwe.Ciphertext)
	MultByConst(ctIn *rlwe.Ciphertext, constant interface{}, ctOut *rlwe.Ciphertext)
	MultByConstThenAdd(ctIn *rlwe.Ciphertext, constant interface{}, ctOut *rlwe.Ciphertext)
	ConjugateNew(ctIn *rlwe.Ciphertext) (ctOut *rlwe.Ciphertext)
	Conjugate(ctIn *rlwe.Ciphertext, ctOut *rlwe.Ciphertext)
	Mul(ctIn *rlwe.Ciphertext, op1 rlwe.Operand, ctOut *rlwe.Ciphertext)
	MulNew(ctIn *rlwe.Ciphertext, op1 rlwe.Operand) (ctOut *rlwe.Ciphertext)
	MulRelin(ctIn *rlwe.Ciphertext, op1 rlwe.Operand, ctOut *rlwe.Ciphertext)
	MulRelinNew(ctIn *rlwe.Ciphertext, op1 rlwe.Operand) (ctOut *rlwe.Ciphertext)
	RotateNew(ctIn *rlwe.Ciphertext, k int) (ctOut *rlwe.Ciphertext)
	Rotate(ctIn *rlwe.Ciphertext, k int, ctOut *rlwe.Ciphertext)
	RotateHoistedNew(ctIn *rlwe.Ciphertext, rotations []int) (ctOut map[int]*rlwe.Ciphertext)
	RotateHoisted(ctIn *rlwe.Ciphertext, rotations []int, ctOut map[int]*rlwe.Ciphertext)
	EvaluatePoly(input interface{}, pol *ckks.Polynomial, targetscale rlwe.Scale) (ctOut *rlwe.Ciphertext, err error)
	EvaluatePolyVector(input interface{}, pols []*ckks.Polynomial, encoder ckks.Encoder, slotIndex map[int][]int, targetscale rlwe.Scale) (ctOut *rlwe.Ciphertext, err error)
	InverseNew(ctIn *rlwe.Ciphertext, steps int) (ctOut *rlwe.Ciphertext, err error)
	LinearTransformNew(ctIn *rlwe.Ciphertext, linearTransform interface{}) (ctOut []*rlwe.Ciphertext)
	LinearTransform(ctIn *rlwe.Ciphertext, linearTransform interface{}, ctOut []*rlwe.Ciphertext)
	MultiplyByDiagMatrix(ctIn *rlwe.Ciphertext, matrix ckks.LinearTransform, c2DecompQP []ringqp.Poly, ctOut *rlwe.Ciphertext)
	MultiplyByDiagMatrixBSGS(ctIn *rlwe.Ciphertext, matrix ckks.LinearTransform, c2DecompQP []ringqp.Poly, ctOut *rlwe.Ciphertext)
	InnerSum(ctIn *rlwe.Ciphertext, batch, n int, ctOut *rlwe.Ciphertext)
	Replicate(ctIn *rlwe.Ciphertext, batch, n int, ctOut *rlwe.Ciphertext)
	TraceNew(ctIn *rlwe.Ciphertext, logSlots int) *rlwe.Ciphertext
	Trace(ctIn *rlwe.Ciphertext, logSlots int, ctOut *rlwe.Ciphertext)
	SwitchKeysNew(ctIn *rlwe.Ciphertext, switchingKey *rlwe.SwitchingKey) (ctOut *rlwe.Ciphertext)
	SwitchKeys(ctIn *rlwe.Ciphertext, switchingKey *rlwe.SwitchingKey, ctOut *rlwe.Ciphertext)
	RelinearizeNew(ctIn *rlwe.Ciphertext) (ctOut *rlwe.Ciphertext)
	Relinearize(ctIn *rlwe.Ciphertext, ctOut *rlwe.Ciphertext)
	ScaleUpNew(ctIn *rlwe.Ciphertext, scale rlwe.Scale) (ctOut *rlwe.Ciphertext)
	ScaleUp(ctIn *rlwe.Ciphertext, scale rlwe.Scale, ctOut *rlwe.Ciphertext)
	SetScale(ctIn *rlwe.Ciphertext, scale rlwe.Scale)
	Rescale(ctIn *rlwe.Ciphertext, minscale rlwe.Scale, ctOut *rlwe.Ciphertext) (err error)
	DropLevelNew(ctIn *rlwe.Ciphertext, levels int) (ctOut *rlwe.Ciphertext)
	DropLevel(ctIn *rlwe.Ciphertext, levels int)

	// ======================================
	// === advanced.Evaluator new methods ===
	// ======================================

	CoeffsToSlotsNew(ctIn *rlwe.Ciphertext, ctsMatrices EncodingMatrix) (ctReal, ctImag *rlwe.Ciphertext)
	CoeffsToSlots(ctIn *rlwe.Ciphertext, ctsMatrices EncodingMatrix, ctReal, ctImag *rlwe.Ciphertext)
	SlotsToCoeffsNew(ctReal, ctImag *rlwe.Ciphertext, stcMatrices EncodingMatrix) (ctOut *rlwe.Ciphertext)
	SlotsToCoeffs(ctReal, ctImag *rlwe.Ciphertext, stcMatrices EncodingMatrix, ctOut *rlwe.Ciphertext)
	EvalModNew(ctIn *rlwe.Ciphertext, evalModPoly EvalModPoly) (ctOut *rlwe.Ciphertext)

	// =================================================
	// === original ckks.Evaluator redefined methods ===
	// =================================================

	GetRLWEEvaluator() *rlwe.Evaluator
	BuffQ() [3]*ring.Poly
	BuffCt() *rlwe.Ciphertext
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
func (eval *evaluator) CoeffsToSlotsNew(ctIn *rlwe.Ciphertext, ctsMatrices EncodingMatrix) (ctReal, ctImag *rlwe.Ciphertext) {
	ctReal = ckks.NewCiphertext(eval.params, 1, ctsMatrices.LevelStart)

	if eval.params.LogSlots() == eval.params.LogN()-1 {
		ctImag = ckks.NewCiphertext(eval.params, 1, ctsMatrices.LevelStart)
	}

	eval.CoeffsToSlots(ctIn, ctsMatrices, ctReal, ctImag)
	return
}

// CoeffsToSlots applies the homomorphic encoding and returns the results on the provided ciphertexts.
// Homomorphically encodes a complex vector vReal + i*vImag of size n on a real vector of size 2n.
// If the packing is sparse (n < N/2), then returns ctReal = Ecd(vReal || vImag) and ctImag = nil.
// If the packing is dense (n == N/2), then returns ctReal = Ecd(vReal) and ctImag = Ecd(vImag).
func (eval *evaluator) CoeffsToSlots(ctIn *rlwe.Ciphertext, ctsMatrices EncodingMatrix, ctReal, ctImag *rlwe.Ciphertext) {

	if ctsMatrices.RepackImag2Real {

		zV := ctIn.CopyNew()

		eval.dft(ctIn, ctsMatrices.matrices, zV)

		eval.Conjugate(zV, ctReal)

		var tmp *rlwe.Ciphertext
		if ctImag != nil {
			tmp = ctImag
		} else {
			tmp = rlwe.NewCiphertextAtLevelFromPoly(ctReal.Level(), eval.BuffCt().Value[:2])
			tmp.IsNTT = true
		}

		// Imag part
		eval.Sub(zV, ctReal, tmp)
		eval.MultByConst(tmp, -1i, tmp)

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
func (eval *evaluator) SlotsToCoeffsNew(ctReal, ctImag *rlwe.Ciphertext, stcMatrices EncodingMatrix) (ctOut *rlwe.Ciphertext) {

	if ctReal.Level() < stcMatrices.LevelStart || (ctImag != nil && ctImag.Level() < stcMatrices.LevelStart) {
		panic("ctReal.Level() or ctImag.Level() < EncodingMatrix.LevelStart")
	}

	ctOut = ckks.NewCiphertext(eval.params, 1, stcMatrices.LevelStart)
	eval.SlotsToCoeffs(ctReal, ctImag, stcMatrices, ctOut)
	return

}

// SlotsToCoeffsNew applies the homomorphic decoding and returns the result on the provided ciphertext.
// Homomorphically decodes a real vector of size 2n on a complex vector vReal + i*vImag of size n.
// If the packing is sparse (n < N/2) then ctReal = Ecd(vReal || vImag) and ctImag = nil.
// If the packing is dense (n == N/2), then ctReal = Ecd(vReal) and ctImag = Ecd(vImag).
func (eval *evaluator) SlotsToCoeffs(ctReal, ctImag *rlwe.Ciphertext, stcMatrices EncodingMatrix, ctOut *rlwe.Ciphertext) {
	// If full packing, the repacking can be done directly using ct0 and ct1.
	if ctImag != nil {
		eval.MultByConst(ctImag, 1i, ctOut)
		eval.Add(ctOut, ctReal, ctOut)
		eval.dft(ctOut, stcMatrices.matrices, ctOut)
	} else {
		eval.dft(ctReal, stcMatrices.matrices, ctOut)
	}
}

func (eval *evaluator) dft(ctIn *rlwe.Ciphertext, plainVectors []ckks.LinearTransform, ctOut *rlwe.Ciphertext) {
	// Sequentially multiplies w with the provided dft matrices.
	scale := ctIn.Scale
	var in, out *rlwe.Ciphertext
	for i, plainVector := range plainVectors {
		in, out = ctOut, ctOut
		if i == 0 {
			in, out = ctIn, ctOut
		}
		eval.LinearTransform(in, plainVector, []*rlwe.Ciphertext{out})
		if err := eval.Rescale(out, scale, out); err != nil {
			panic(err)
		}
	}
}

// EvalModNew applies a homomorphic mod Q on a vector scaled by Delta, scaled down to mod 1 :
//
//  1. Delta * (Q/Delta * I(X) + m(X)) (Delta = scaling factor, I(X) integer poly, m(X) message)
//  2. Delta * (I(X) + Delta/Q * m(X)) (divide by Q/Delta)
//  3. Delta * (Delta/Q * m(X)) (x mod 1)
//  4. Delta * (m(X)) (multiply back by Q/Delta)
//
// Since Q is not a power of two, but Delta is, then does an approximate division by the closest
// power of two to Q instead. Hence, it assumes that the input plaintext is already scaled by
// the correcting factor Q/2^{round(log(Q))}.
//
// !! Assumes that the input is normalized by 1/K for K the range of the approximation.
//
// Scaling back error correction by 2^{round(log(Q))}/Q afterward is included in the polynomial
func (eval *evaluator) EvalModNew(ct *rlwe.Ciphertext, evalModPoly EvalModPoly) *rlwe.Ciphertext {

	if ct.Level() < evalModPoly.LevelStart() {
		panic("ct.Level() < evalModPoly.LevelStart")
	}

	if ct.Level() > evalModPoly.LevelStart() {
		eval.DropLevel(ct, ct.Level()-evalModPoly.LevelStart())
	}

	// Stores default scales
	prevScaleCt := ct.Scale

	// Normalize the modular reduction to mod by 1 (division by Q)
	ct.Scale = evalModPoly.scalingFactor

	var err error

	// Compute the scales that the ciphertext should have before the double angle
	// formula such that after it it has the scale it had before the polynomial
	// evaluation

	targetScale := ct.Scale.Float64()
	for i := 0; i < evalModPoly.doubleAngle; i++ {
		targetScale = math.Sqrt(targetScale * eval.params.QiFloat64(evalModPoly.levelStart-evalModPoly.sinePoly.Depth()-evalModPoly.doubleAngle+i+1))
	}

	// Division by 1/2^r and change of variable for the Chebyshev evaluation
	if evalModPoly.sineType == Cos1 || evalModPoly.sineType == Cos2 {
		eval.AddConst(ct, -0.5/(evalModPoly.scFac*(evalModPoly.sinePoly.B-evalModPoly.sinePoly.A)), ct)
	}

	// Chebyshev evaluation
	if ct, err = eval.EvaluatePoly(ct, evalModPoly.sinePoly, rlwe.NewScale(targetScale)); err != nil {
		panic(err)
	}

	// Double angle
	sqrt2pi := evalModPoly.sqrt2Pi
	for i := 0; i < evalModPoly.doubleAngle; i++ {
		sqrt2pi *= sqrt2pi
		eval.MulRelin(ct, ct, ct)
		eval.Add(ct, ct, ct)
		eval.AddConst(ct, -sqrt2pi, ct)
		if err := eval.Rescale(ct, rlwe.NewScale(targetScale), ct); err != nil {
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
