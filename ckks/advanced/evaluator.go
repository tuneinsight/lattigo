// Package advanced implements advanced operations for the CKKS scheme.
package advanced

import (
	"math/big"

	"github.com/tuneinsight/lattigo/v4/ckks"
	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/rlwe"
	"github.com/tuneinsight/lattigo/v4/rlwe/ringqp"
	"github.com/tuneinsight/lattigo/v4/utils/bignum"
)

// Evaluator is an interface embedding the ckks.Evaluator interface with
// additional advanced arithmetic features.
type Evaluator interface {

	// =======================================
	// === Original ckks.Evaluator methods ===
	// =======================================

	// ========================
	// === Basic Arithmetic ===
	// ========================

	// Addition
	Add(op0 *rlwe.Ciphertext, op1 interface{}, op2 *rlwe.Ciphertext)
	AddNew(op0 *rlwe.Ciphertext, op1 interface{}) (op2 *rlwe.Ciphertext)

	// Subtraction
	Sub(op0 *rlwe.Ciphertext, op1 interface{}, op2 *rlwe.Ciphertext)
	SubNew(op0 *rlwe.Ciphertext, op1 interface{}) (op2 *rlwe.Ciphertext)

	// Complex Conjugation
	ConjugateNew(op0 *rlwe.Ciphertext) (ctOut *rlwe.Ciphertext)
	Conjugate(op0 *rlwe.Ciphertext, ctOut *rlwe.Ciphertext)

	// Multiplication
	Mul(op0 *rlwe.Ciphertext, op1 interface{}, ctOut *rlwe.Ciphertext)
	MulNew(op0 *rlwe.Ciphertext, op1 interface{}) (ctOut *rlwe.Ciphertext)
	MulRelin(op0 *rlwe.Ciphertext, op1 rlwe.Operand, ctOut *rlwe.Ciphertext)
	MulRelinNew(op0 *rlwe.Ciphertext, op1 rlwe.Operand) (ctOut *rlwe.Ciphertext)

	MulThenAdd(op0 *rlwe.Ciphertext, op1 interface{}, ctOut *rlwe.Ciphertext)
	MulRelinThenAdd(op0 *rlwe.Ciphertext, op1 rlwe.Operand, ctOut *rlwe.Ciphertext)

	// Slot Rotations
	RotateNew(op0 *rlwe.Ciphertext, k int) (ctOut *rlwe.Ciphertext)
	Rotate(op0 *rlwe.Ciphertext, k int, ctOut *rlwe.Ciphertext)
	RotateHoistedNew(op0 *rlwe.Ciphertext, rotations []int) (ctOut map[int]*rlwe.Ciphertext)
	RotateHoisted(op0 *rlwe.Ciphertext, rotations []int, ctOut map[int]*rlwe.Ciphertext)
	RotateHoistedLazyNew(level int, rotations []int, ct *rlwe.Ciphertext, c2DecompQP []ringqp.Poly) (cOut map[int]rlwe.CiphertextQP)

	// ===========================
	// === Advanced Arithmetic ===
	// ===========================

	// Polynomial evaluation
	EvaluatePoly(input interface{}, pol *bignum.Polynomial, targetScale rlwe.Scale) (ctOut *rlwe.Ciphertext, err error)
	EvaluatePolyVector(input interface{}, pols []*bignum.Polynomial, encoder *ckks.Encoder, slotIndex map[int][]int, targetScale rlwe.Scale) (ctOut *rlwe.Ciphertext, err error)

	// GoldschmidtDivision
	GoldschmidtDivisionNew(ct *rlwe.Ciphertext, minValue, log2Targetprecision float64, btp rlwe.Bootstrapper) (ctInv *rlwe.Ciphertext, err error)

	// Linear Transformations
	LinearTransformNew(op0 *rlwe.Ciphertext, linearTransform interface{}) (ctOut []*rlwe.Ciphertext)
	LinearTransform(op0 *rlwe.Ciphertext, linearTransform interface{}, ctOut []*rlwe.Ciphertext)
	MultiplyByDiagMatrix(op0 *rlwe.Ciphertext, matrix ckks.LinearTransform, c2DecompQP []ringqp.Poly, ctOut *rlwe.Ciphertext)
	MultiplyByDiagMatrixBSGS(op0 *rlwe.Ciphertext, matrix ckks.LinearTransform, c2DecompQP []ringqp.Poly, ctOut *rlwe.Ciphertext)

	// Inner sum
	InnerSum(op0 *rlwe.Ciphertext, batch, n int, ctOut *rlwe.Ciphertext)
	Average(op0 *rlwe.Ciphertext, batch int, ctOut *rlwe.Ciphertext)

	// Replication (inverse of Inner sum)
	Replicate(op0 *rlwe.Ciphertext, batch, n int, ctOut *rlwe.Ciphertext)

	// Trace
	Trace(op0 *rlwe.Ciphertext, logSlots int, ctOut *rlwe.Ciphertext)
	TraceNew(op0 *rlwe.Ciphertext, logSlots int) (ctOut *rlwe.Ciphertext)

	// =============================
	// === Ciphertext Management ===
	// =============================

	// Generic EvaluationKeys
	ApplyEvaluationKeyNew(op0 *rlwe.Ciphertext, evk *rlwe.EvaluationKey) (ctOut *rlwe.Ciphertext)
	ApplyEvaluationKey(op0 *rlwe.Ciphertext, evk *rlwe.EvaluationKey, ctOut *rlwe.Ciphertext)

	// Degree Management
	RelinearizeNew(op0 *rlwe.Ciphertext) (ctOut *rlwe.Ciphertext)
	Relinearize(op0 *rlwe.Ciphertext, ctOut *rlwe.Ciphertext)

	// Scale Management
	ScaleUpNew(op0 *rlwe.Ciphertext, scale rlwe.Scale) (ctOut *rlwe.Ciphertext)
	ScaleUp(op0 *rlwe.Ciphertext, scale rlwe.Scale, ctOut *rlwe.Ciphertext)
	SetScale(op0 *rlwe.Ciphertext, scale rlwe.Scale)
	Rescale(op0 *rlwe.Ciphertext, minScale rlwe.Scale, ctOut *rlwe.Ciphertext) (err error)

	// Level Management
	DropLevelNew(op0 *rlwe.Ciphertext, levels int) (ctOut *rlwe.Ciphertext)
	DropLevel(op0 *rlwe.Ciphertext, levels int)

	// ======================================
	// === advanced.Evaluator new methods ===
	// ======================================

	CoeffsToSlotsNew(ctIn *rlwe.Ciphertext, ctsMatrices HomomorphicDFTMatrix) (ctReal, ctImag *rlwe.Ciphertext)
	CoeffsToSlots(ctIn *rlwe.Ciphertext, ctsMatrices HomomorphicDFTMatrix, ctReal, ctImag *rlwe.Ciphertext)
	SlotsToCoeffsNew(ctReal, ctImag *rlwe.Ciphertext, stcMatrices HomomorphicDFTMatrix) (ctOut *rlwe.Ciphertext)
	SlotsToCoeffs(ctReal, ctImag *rlwe.Ciphertext, stcMatrices HomomorphicDFTMatrix, ctOut *rlwe.Ciphertext)
	EvalModNew(ctIn *rlwe.Ciphertext, evalModPoly EvalModPoly) (ctOut *rlwe.Ciphertext)

	// =================================================
	// === original ckks.Evaluator redefined methods ===
	// =================================================

	CheckBinary(op0, op1, opOut rlwe.Operand, opOutMinDegree int) (degree, level int)
	CheckUnary(op0, opOut rlwe.Operand) (degree, level int)
	GetRLWEEvaluator() *rlwe.Evaluator
	BuffQ() [3]*ring.Poly
	BuffCt() *rlwe.Ciphertext
	ShallowCopy() Evaluator
	WithKey(evk rlwe.EvaluationKeySetInterface) Evaluator
}

type evaluator struct {
	ckks.Evaluator
	params ckks.Parameters
}

// NewEvaluator creates a new Evaluator.
func NewEvaluator(params ckks.Parameters, evk rlwe.EvaluationKeySetInterface) Evaluator {
	return &evaluator{ckks.NewEvaluator(params, evk), params}
}

// ShallowCopy creates a shallow copy of this evaluator in which all the read-only data-structures are
// shared with the receiver and the temporary buffers are reallocated. The receiver and the returned
// Evaluators can be used concurrently.
func (eval *evaluator) ShallowCopy() Evaluator {
	return &evaluator{eval.Evaluator.ShallowCopy(), eval.params}
}

// Parameters returns the ckks.Parameters of the target Evaluator.
func (eval *evaluator) Parameters() ckks.Parameters {
	return eval.params
}

// WithKey creates a shallow copy of the receiver Evaluator for which the new EvaluationKey is evaluationKey
// and where the temporary buffers are shared. The receiver and the returned Evaluators cannot be used concurrently.
func (eval *evaluator) WithKey(evk rlwe.EvaluationKeySetInterface) Evaluator {
	return &evaluator{eval.Evaluator.WithKey(evk), eval.params}
}

// CoeffsToSlotsNew applies the homomorphic encoding and returns the result on new ciphertexts.
// Homomorphically encodes a complex vector vReal + i*vImag.
// If the packing is sparse (n < N/2), then returns ctReal = Ecd(vReal || vImag) and ctImag = nil.
// If the packing is dense (n == N/2), then returns ctReal = Ecd(vReal) and ctImag = Ecd(vImag).
func (eval *evaluator) CoeffsToSlotsNew(ctIn *rlwe.Ciphertext, ctsMatrices HomomorphicDFTMatrix) (ctReal, ctImag *rlwe.Ciphertext) {
	ctReal = ckks.NewCiphertext(eval.params, 1, ctsMatrices.LevelStart)

	if ctsMatrices.LogSlots == eval.params.MaxLogSlots() {
		ctImag = ckks.NewCiphertext(eval.params, 1, ctsMatrices.LevelStart)
	}

	eval.CoeffsToSlots(ctIn, ctsMatrices, ctReal, ctImag)
	return
}

// CoeffsToSlots applies the homomorphic encoding and returns the results on the provided ciphertexts.
// Homomorphically encodes a complex vector vReal + i*vImag of size n on a real vector of size 2n.
// If the packing is sparse (n < N/2), then returns ctReal = Ecd(vReal || vImag) and ctImag = nil.
// If the packing is dense (n == N/2), then returns ctReal = Ecd(vReal) and ctImag = Ecd(vImag).
func (eval *evaluator) CoeffsToSlots(ctIn *rlwe.Ciphertext, ctsMatrices HomomorphicDFTMatrix, ctReal, ctImag *rlwe.Ciphertext) {

	if ctsMatrices.RepackImag2Real {

		zV := ctIn.CopyNew()

		eval.dft(ctIn, ctsMatrices.Matrices, zV)

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
		eval.Mul(tmp, -1i, tmp)

		// Real part
		eval.Add(ctReal, zV, ctReal)

		// If repacking, then ct0 and ct1 right n/2 slots are zero.
		if ctsMatrices.LogSlots < eval.params.MaxLogSlots() {
			eval.Rotate(tmp, ctIn.Slots(), tmp)
			eval.Add(ctReal, tmp, ctReal)
		}

		zV = nil

	} else {
		eval.dft(ctIn, ctsMatrices.Matrices, ctReal)
	}
}

// SlotsToCoeffsNew applies the homomorphic decoding and returns the result on a new ciphertext.
// Homomorphically decodes a real vector of size 2n on a complex vector vReal + i*vImag of size n.
// If the packing is sparse (n < N/2) then ctReal = Ecd(vReal || vImag) and ctImag = nil.
// If the packing is dense (n == N/2), then ctReal = Ecd(vReal) and ctImag = Ecd(vImag).
func (eval *evaluator) SlotsToCoeffsNew(ctReal, ctImag *rlwe.Ciphertext, stcMatrices HomomorphicDFTMatrix) (ctOut *rlwe.Ciphertext) {

	if ctReal.Level() < stcMatrices.LevelStart || (ctImag != nil && ctImag.Level() < stcMatrices.LevelStart) {
		panic("ctReal.Level() or ctImag.Level() < HomomorphicDFTMatrix.LevelStart")
	}

	ctOut = ckks.NewCiphertext(eval.params, 1, stcMatrices.LevelStart)
	eval.SlotsToCoeffs(ctReal, ctImag, stcMatrices, ctOut)
	return

}

// SlotsToCoeffsNew applies the homomorphic decoding and returns the result on the provided ciphertext.
// Homomorphically decodes a real vector of size 2n on a complex vector vReal + i*vImag of size n.
// If the packing is sparse (n < N/2) then ctReal = Ecd(vReal || vImag) and ctImag = nil.
// If the packing is dense (n == N/2), then ctReal = Ecd(vReal) and ctImag = Ecd(vImag).
func (eval *evaluator) SlotsToCoeffs(ctReal, ctImag *rlwe.Ciphertext, stcMatrices HomomorphicDFTMatrix, ctOut *rlwe.Ciphertext) {
	// If full packing, the repacking can be done directly using ct0 and ct1.
	if ctImag != nil {
		eval.Mul(ctImag, 1i, ctOut)
		eval.Add(ctOut, ctReal, ctOut)
		eval.dft(ctOut, stcMatrices.Matrices, ctOut)
	} else {
		eval.dft(ctReal, stcMatrices.Matrices, ctOut)
	}
}

func (eval *evaluator) dft(ctIn *rlwe.Ciphertext, plainVectors []ckks.LinearTransform, ctOut *rlwe.Ciphertext) {

	inputLogSlots := ctIn.LogSlots

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

	// Encoding matrices are a special case of `fractal` linear transform
	// that doesn't change the underlying plaintext polynomial Y = X^{N/n}
	// of the input ciphertext.
	ctOut.LogSlots = inputLogSlots
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
	ct.Scale = evalModPoly.ScalingFactor()

	var err error

	// Compute the scales that the ciphertext should have before the double angle
	// formula such that after it it has the scale it had before the polynomial
	// evaluation

	targetScale := ct.Scale
	for i := 0; i < evalModPoly.doubleAngle; i++ {
		qi := eval.params.Q()[evalModPoly.levelStart-evalModPoly.sinePoly.Depth()-evalModPoly.doubleAngle+i+1]
		targetScale = targetScale.Mul(rlwe.NewScale(qi))
		targetScale.Value.Sqrt(&targetScale.Value)
	}

	// Division by 1/2^r and change of variable for the Chebyshev evaluation
	if evalModPoly.sineType == CosDiscrete || evalModPoly.sineType == CosContinuous {
		offset := new(big.Float).Sub(evalModPoly.sinePoly.B, evalModPoly.sinePoly.A)
		offset.Mul(offset, new(big.Float).SetFloat64(evalModPoly.scFac))
		offset.Quo(new(big.Float).SetFloat64(-0.5), offset)
		eval.Add(ct, offset, ct)
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
		eval.Add(ct, -sqrt2pi, ct)
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
