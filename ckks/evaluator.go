package ckks

import (
	"errors"
	"math"
	"math/big"

	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/rlwe"
	"github.com/tuneinsight/lattigo/v4/rlwe/ringqp"
	"github.com/tuneinsight/lattigo/v4/utils"
)

// Evaluator is an interface implementing the methods to conduct homomorphic operations between ciphertext and/or plaintexts.
type Evaluator interface {
	// ========================
	// === Basic Arithmetic ===
	// ========================

	// Addition
	Add(ctIn *rlwe.Ciphertext, op1 rlwe.Operand, ctOut *rlwe.Ciphertext)
	AddNew(ctIn *rlwe.Ciphertext, op1 rlwe.Operand) (ctOut *rlwe.Ciphertext)

	// Subtraction
	Sub(ctIn *rlwe.Ciphertext, op1 rlwe.Operand, ctOut *rlwe.Ciphertext)
	SubNew(ctIn *rlwe.Ciphertext, op1 rlwe.Operand) (ctOut *rlwe.Ciphertext)

	// Negation
	Neg(ctIn *rlwe.Ciphertext, ctOut *rlwe.Ciphertext)
	NegNew(ctIn *rlwe.Ciphertext) (ctOut *rlwe.Ciphertext)

	// Constant Addition
	AddConstNew(ctIn *rlwe.Ciphertext, constant interface{}) (ctOut *rlwe.Ciphertext)
	AddConst(ctIn *rlwe.Ciphertext, constant interface{}, ctOut *rlwe.Ciphertext)

	// Constant Multiplication
	MultByConstNew(ctIn *rlwe.Ciphertext, constant interface{}) (ctOut *rlwe.Ciphertext)
	MultByConst(ctIn *rlwe.Ciphertext, constant interface{}, ctOut *rlwe.Ciphertext)
	MultByGaussianInteger(ctIn *rlwe.Ciphertext, cReal, cImag interface{}, ctOut *rlwe.Ciphertext)

	// Constant Multiplication with Addition
	MultByConstAndAdd(ctIn *rlwe.Ciphertext, constant interface{}, ctOut *rlwe.Ciphertext)
	MultByGaussianIntegerAndAdd(ctIn *rlwe.Ciphertext, cReal, cImag interface{}, ctOut *rlwe.Ciphertext)

	// Multiplication by the imaginary unit
	MultByiNew(ctIn *rlwe.Ciphertext) (ctOut *rlwe.Ciphertext)
	MultByi(ctIn *rlwe.Ciphertext, ctOut *rlwe.Ciphertext)
	DivByiNew(ctIn *rlwe.Ciphertext) (ctOut *rlwe.Ciphertext)
	DivByi(ctIn *rlwe.Ciphertext, ctOut *rlwe.Ciphertext)

	// Conjugation
	ConjugateNew(ctIn *rlwe.Ciphertext) (ctOut *rlwe.Ciphertext)
	Conjugate(ctIn *rlwe.Ciphertext, ctOut *rlwe.Ciphertext)

	// Multiplication
	Mul(ctIn *rlwe.Ciphertext, op1 rlwe.Operand, ctOut *rlwe.Ciphertext)
	MulNew(ctIn *rlwe.Ciphertext, op1 rlwe.Operand) (ctOut *rlwe.Ciphertext)
	MulRelin(ctIn *rlwe.Ciphertext, op1 rlwe.Operand, ctOut *rlwe.Ciphertext)
	MulRelinNew(ctIn *rlwe.Ciphertext, op1 rlwe.Operand) (ctOut *rlwe.Ciphertext)

	MulAndAdd(ctIn *rlwe.Ciphertext, op1 rlwe.Operand, ctOut *rlwe.Ciphertext)
	MulRelinAndAdd(ctIn *rlwe.Ciphertext, op1 rlwe.Operand, ctOut *rlwe.Ciphertext)

	// Slot Rotations
	RotateNew(ctIn *rlwe.Ciphertext, k int) (ctOut *rlwe.Ciphertext)
	Rotate(ctIn *rlwe.Ciphertext, k int, ctOut *rlwe.Ciphertext)
	RotateHoistedNew(ctIn *rlwe.Ciphertext, rotations []int) (ctOut map[int]*rlwe.Ciphertext)
	RotateHoisted(ctIn *rlwe.Ciphertext, rotations []int, ctOut map[int]*rlwe.Ciphertext)
	RotateHoistedNoModDownNew(level int, rotations []int, c0 *ring.Poly, c2DecompQP []ringqp.Poly) (cOut map[int]rlwe.CiphertextQP)

	// ===========================
	// === Advanced Arithmetic ===
	// ===========================

	// Polynomial evaluation
	EvaluatePoly(input interface{}, pol *Polynomial, targetScale rlwe.Scale) (ctOut *rlwe.Ciphertext, err error)
	EvaluatePolyVector(input interface{}, pols []*Polynomial, encoder Encoder, slotIndex map[int][]int, targetScale rlwe.Scale) (ctOut *rlwe.Ciphertext, err error)

	// Inversion
	InverseNew(ctIn *rlwe.Ciphertext, steps int) (ctOut *rlwe.Ciphertext, err error)

	// Linear Transformations
	LinearTransformNew(ctIn *rlwe.Ciphertext, linearTransform interface{}) (ctOut []*rlwe.Ciphertext)
	LinearTransform(ctIn *rlwe.Ciphertext, linearTransform interface{}, ctOut []*rlwe.Ciphertext)
	MultiplyByDiagMatrix(ctIn *rlwe.Ciphertext, matrix LinearTransform, c2DecompQP []ringqp.Poly, ctOut *rlwe.Ciphertext)
	MultiplyByDiagMatrixBSGS(ctIn *rlwe.Ciphertext, matrix LinearTransform, c2DecompQP []ringqp.Poly, ctOut *rlwe.Ciphertext)

	// Inner sum
	InnerSum(ctIn *rlwe.Ciphertext, batch, n int, ctOut *rlwe.Ciphertext)
	Average(ctIn *rlwe.Ciphertext, batch int, ctOut *rlwe.Ciphertext)

	// Replication (inverse of Inner sum)
	Replicate(ctIn *rlwe.Ciphertext, batch, n int, ctOut *rlwe.Ciphertext)

	// Trace
	Trace(ctIn *rlwe.Ciphertext, logSlots int, ctOut *rlwe.Ciphertext)
	TraceNew(ctIn *rlwe.Ciphertext, logSlots int) (ctOut *rlwe.Ciphertext)

	// =============================
	// === Ciphertext Management ===
	// =============================

	// Key-Switching
	SwitchKeysNew(ctIn *rlwe.Ciphertext, switchingKey *rlwe.SwitchingKey) (ctOut *rlwe.Ciphertext)
	SwitchKeys(ctIn *rlwe.Ciphertext, switchingKey *rlwe.SwitchingKey, ctOut *rlwe.Ciphertext)

	// Degree Management
	RelinearizeNew(ctIn *rlwe.Ciphertext) (ctOut *rlwe.Ciphertext)
	Relinearize(ctIn *rlwe.Ciphertext, ctOut *rlwe.Ciphertext)

	// Scale Management
	ScaleUpNew(ctIn *rlwe.Ciphertext, scale rlwe.Scale) (ctOut *rlwe.Ciphertext)
	ScaleUp(ctIn *rlwe.Ciphertext, scale rlwe.Scale, ctOut *rlwe.Ciphertext)
	SetScale(ctIn *rlwe.Ciphertext, scale rlwe.Scale)
	Rescale(ctIn *rlwe.Ciphertext, minScale rlwe.Scale, ctOut *rlwe.Ciphertext) (err error)

	// Level Management
	DropLevelNew(ctIn *rlwe.Ciphertext, levels int) (ctOut *rlwe.Ciphertext)
	DropLevel(ctIn *rlwe.Ciphertext, levels int)

	// ==============
	// === Others ===
	// ==============
	GetRLWEEvaluator() *rlwe.Evaluator
	BuffQ() [3]*ring.Poly
	BuffCt() *rlwe.Ciphertext
	ShallowCopy() Evaluator
	WithKey(rlwe.EvaluationKey) Evaluator
}

// evaluator is a struct that holds the necessary elements to execute the homomorphic operations between Ciphertexts and/or Plaintexts.
// It also holds a memory buffer used to store intermediate computations.
type evaluator struct {
	*evaluatorBase
	*evaluatorBuffers
	*rlwe.Evaluator
}

type evaluatorBase struct {
	params Parameters
}

type evaluatorBuffers struct {
	buffQ  [3]*ring.Poly    // Memory buffer in order: for MForm(c0), MForm(c1), c2
	buffCt *rlwe.Ciphertext // Memory buffer for ciphertexts that need to be scaled up (to be eventually removed)
}

// BuffQ returns a pointer to the internal memory buffer buffQ.
func (eval *evaluator) BuffQ() [3]*ring.Poly {
	return eval.buffQ
}

// BuffCt returns a pointer to the internal memory buffer buffCt.
func (eval *evaluator) BuffCt() *rlwe.Ciphertext {
	return eval.buffCt
}

func newEvaluatorBase(params Parameters) *evaluatorBase {
	ev := new(evaluatorBase)
	ev.params = params
	return ev
}

func newEvaluatorBuffers(evalBase *evaluatorBase) *evaluatorBuffers {
	buff := new(evaluatorBuffers)
	params := evalBase.params
	ringQ := params.RingQ()
	buff.buffQ = [3]*ring.Poly{ringQ.NewPoly(), ringQ.NewPoly(), ringQ.NewPoly()}
	buff.buffCt = NewCiphertext(params, 2, params.MaxLevel())
	return buff
}

// NewEvaluator creates a new Evaluator, that can be used to do homomorphic
// operations on the Ciphertexts and/or Plaintexts. It stores a memory buffer
// and Ciphertexts that will be used for intermediate values.
func NewEvaluator(params Parameters, evaluationKey rlwe.EvaluationKey) Evaluator {
	eval := new(evaluator)
	eval.evaluatorBase = newEvaluatorBase(params)
	eval.evaluatorBuffers = newEvaluatorBuffers(eval.evaluatorBase)
	eval.Evaluator = rlwe.NewEvaluator(params.Parameters, &evaluationKey)

	return eval
}

// GetRLWEEvaluator returns the underlying *rlwe.Evaluator.
func (eval *evaluator) GetRLWEEvaluator() *rlwe.Evaluator {
	return eval.Evaluator
}

func (eval *evaluator) PermuteNTTIndexesForKey(rtks *rlwe.RotationKeySet) *map[uint64][]uint64 {
	if rtks == nil {
		return &map[uint64][]uint64{}
	}
	PermuteNTTIndex := make(map[uint64][]uint64, len(rtks.Keys))
	for galEl := range rtks.Keys {
		PermuteNTTIndex[galEl] = eval.params.RingQ().PermuteNTTIndex(galEl)
	}
	return &PermuteNTTIndex
}

func (eval *evaluator) checkBinary(op0, op1, opOut rlwe.Operand, opOutMinDegree int) {
	if op0 == nil || op1 == nil || opOut == nil {
		panic("cannot checkBinary: rlwe.Operands cannot be nil")
	}

	if op0.Degree()+op1.Degree() == 0 {
		panic("cannot checkBinary: rlwe.Operands cannot be both plaintext")
	}

	if opOut.Degree() < opOutMinDegree {
		panic("cannot checkBinary: receiver rlwe.Operand degree is too small")
	}

	if !op0.El().IsNTT {
		panic("cannot checkBinary: op0 must be in NTT")
	}

	if !op1.El().IsNTT {
		panic("cannot checkBinary: op1 must be in NTT")
	}
}

func (eval *evaluator) newCiphertextBinary(op0, op1 rlwe.Operand) (ctOut *rlwe.Ciphertext) {

	maxDegree := utils.MaxInt(op0.Degree(), op1.Degree())
	minLevel := utils.MinInt(op0.Level(), op1.Level())

	return NewCiphertext(eval.params, maxDegree, minLevel)
}

// Add adds op1 to ctIn and returns the result in ctOut.
func (eval *evaluator) Add(ctIn *rlwe.Ciphertext, op1 rlwe.Operand, ctOut *rlwe.Ciphertext) {
	eval.checkBinary(ctIn, op1, ctOut, utils.MaxInt(ctIn.Degree(), op1.Degree()))
	eval.evaluateInPlace(ctIn, op1, ctOut, eval.params.RingQ().AddLvl)
}

// AddNew adds op1 to ctIn and returns the result in a newly created element.
func (eval *evaluator) AddNew(ctIn *rlwe.Ciphertext, op1 rlwe.Operand) (ctOut *rlwe.Ciphertext) {
	ctOut = eval.newCiphertextBinary(ctIn, op1)
	eval.Add(ctIn, op1, ctOut)
	return
}

// Sub subtracts op1 from ctIn and returns the result in ctOut.
func (eval *evaluator) Sub(ctIn *rlwe.Ciphertext, op1 rlwe.Operand, ctOut *rlwe.Ciphertext) {

	eval.checkBinary(ctIn, op1, ctOut, utils.MaxInt(ctIn.Degree(), op1.Degree()))

	eval.evaluateInPlace(ctIn, op1, ctOut, eval.params.RingQ().SubLvl)

	level := utils.MinInt(utils.MinInt(ctIn.Level(), op1.Level()), ctOut.Level())

	if ctIn.Degree() < op1.Degree() {
		for i := ctIn.Degree() + 1; i < op1.Degree()+1; i++ {
			eval.params.RingQ().NegLvl(level, ctOut.Value[i], ctOut.Value[i])
		}
	}

}

// SubNew subtracts op1 from ctIn and returns the result in a newly created element.
func (eval *evaluator) SubNew(ctIn *rlwe.Ciphertext, op1 rlwe.Operand) (ctOut *rlwe.Ciphertext) {
	ctOut = eval.newCiphertextBinary(ctIn, op1)
	eval.Sub(ctIn, op1, ctOut)
	return
}

func (eval *evaluator) evaluateInPlace(c0 *rlwe.Ciphertext, c1 rlwe.Operand, ctOut *rlwe.Ciphertext, evaluate func(int, *ring.Poly, *ring.Poly, *ring.Poly)) {

	var tmp0, tmp1 *rlwe.Ciphertext

	level := utils.MinInt(utils.MinInt(c0.Level(), c1.Level()), ctOut.Level())

	maxDegree := utils.MaxInt(c0.Degree(), c1.Degree())
	minDegree := utils.MinInt(c0.Degree(), c1.Degree())

	// Else resizes the receiver element
	ctOut.El().Resize(maxDegree, ctOut.Level())

	c0Scale := c0.GetScale().Float64()
	c1Scale := c1.GetScale().Float64()

	if ctOut.Level() > level {
		eval.DropLevel(ctOut, ctOut.Level()-utils.MinInt(c0.Level(), c1.Level()))
	}

	cmp := c0.GetScale().Cmp(c1.GetScale())

	// Checks whether or not the receiver element is the same as one of the input elements
	// and acts accordingly to avoid unnecessary element creation or element overwriting,
	// and scales properly the element before the evaluation.
	if ctOut == c0 {

		if cmp == 1 && math.Floor(c0Scale/c1Scale) > 1 {

			tmp1 = eval.buffCt.El()
			tmp1.Scale = ctOut.Scale

			eval.MultByConst(c1.El(), math.Floor(c0Scale/c1Scale), tmp1)

		} else if cmp == -1 && math.Floor(c1Scale/c0Scale) > 1 {

			eval.MultByConst(c0, math.Floor(c1Scale/c0Scale), c0)

			ctOut.Scale = c1.GetScale()

			tmp1 = c1.El()
		} else {
			tmp1 = c1.El()
		}

		tmp0 = c0.El()

	} else if ctOut == c1 {

		if cmp == 1 && math.Floor(c0Scale/c1Scale) > 1 {

			eval.MultByConst(c1.El(), math.Floor(c0Scale/c1Scale), ctOut)

			ctOut.Scale = c0.Scale

			tmp0 = c0.El()

		} else if cmp == -1 && math.Floor(c1Scale/c0Scale) > 1 {

			tmp0 = eval.buffCt.El()
			tmp0.Scale = ctOut.Scale

			eval.MultByConst(c0, math.Floor(c1Scale/c0Scale), tmp0)
		} else {
			tmp0 = c0.El()
		}

		tmp1 = c1.El()

	} else {

		if cmp == 1 && math.Floor(c0Scale/c1Scale) > 1 {

			tmp1 = eval.buffCt.El()

			tmp1.Scale = ctOut.Scale

			eval.MultByConst(c1.El(), math.Floor(c0Scale/c1Scale), tmp1)

			tmp0 = c0.El()

		} else if cmp == -1 && math.Floor(c1Scale/c0Scale) > 1 {

			tmp0 = eval.buffCt.El()

			tmp0.Scale = ctOut.Scale

			eval.MultByConst(c0, math.Floor(c1Scale/c0Scale), tmp0)

			tmp1 = c1.El()

		} else {
			tmp0 = c0.El()
			tmp1 = c1.El()
		}
	}

	for i := 0; i < minDegree+1; i++ {
		evaluate(level, tmp0.Value[i], tmp1.Value[i], ctOut.El().Value[i])
	}

	ctOut.MetaData = c0.MetaData
	ctOut.Scale = c0.Scale.Max(c1.GetScale())

	// If the inputs degrees differ, it copies the remaining degree on the receiver.
	// Also checks that the receiver is not one of the inputs to avoid unnecessary work.

	if c0.Degree() > c1.Degree() && tmp0 != ctOut.El() {
		for i := minDegree + 1; i < maxDegree+1; i++ {
			ring.CopyLvl(level, tmp0.Value[i], ctOut.El().Value[i])
		}
	} else if c1.Degree() > c0.Degree() && tmp1 != ctOut.El() {
		for i := minDegree + 1; i < maxDegree+1; i++ {
			ring.CopyLvl(level, tmp1.Value[i], ctOut.El().Value[i])
		}
	}
}

// Neg negates the value of ct0 and returns the result in ctOut.
func (eval *evaluator) Neg(ct0 *rlwe.Ciphertext, ctOut *rlwe.Ciphertext) {

	level := utils.MinInt(ct0.Level(), ctOut.Level())

	if ct0.Degree() != ctOut.Degree() {
		panic("cannot Negate: invalid receiver Ciphertext does not match input Ciphertext degree")
	}

	for i := range ct0.Value {
		eval.params.RingQ().NegLvl(level, ct0.Value[i], ctOut.Value[i])
	}

	ctOut.MetaData = ct0.MetaData
}

// NegNew negates ct0 and returns the result in a newly created element.
func (eval *evaluator) NegNew(ct0 *rlwe.Ciphertext) (ctOut *rlwe.Ciphertext) {
	ctOut = NewCiphertext(eval.params, ct0.Degree(), ct0.Level())
	eval.Neg(ct0, ctOut)
	return
}

// AddConstNew adds the input constant (which can be a uint64, int64, float64 or complex128) to ct0 and returns the result in a new element.
func (eval *evaluator) AddConstNew(ct0 *rlwe.Ciphertext, constant interface{}) (ctOut *rlwe.Ciphertext) {
	ctOut = ct0.CopyNew()
	eval.AddConst(ct0, constant, ctOut)
	return ctOut
}

func (eval *evaluator) getConstAndScale(level int, constant interface{}) (cReal, cImag, scale float64) {

	// Converts to float64 and determines if a scaling is required (which is the case if either real or imag have a rational part)
	scale = 1
	switch constant := constant.(type) {
	case complex128:
		cReal = real(constant)
		cImag = imag(constant)

		if cReal != 0 {
			valueInt := int64(cReal)
			valueFloat := cReal - float64(valueInt)

			if valueFloat != 0 {
				scale = float64(eval.params.RingQ().Tables[level].Modulus)
			}
		}

		if cImag != 0 {
			valueInt := int64(cImag)
			valueFloat := cImag - float64(valueInt)

			if valueFloat != 0 {
				scale = float64(eval.params.RingQ().Tables[level].Modulus)
			}
		}

	case float64:
		cReal = constant
		cImag = float64(0)

		if cReal != 0 {
			valueInt := int64(cReal)
			valueFloat := cReal - float64(valueInt)

			if valueFloat != 0 {
				scale = float64(eval.params.RingQ().Tables[level].Modulus)
			}
		}

	case *big.Float:
		cf64, _ := constant.Float64()
		return eval.getConstAndScale(level, cf64)

	case uint64:
		cReal = float64(constant)
		cImag = float64(0)

	case int64:
		cReal = float64(constant)
		cImag = float64(0)

	case int:
		cReal = float64(constant)
		cImag = float64(0)
	}

	if eval.params.RingType() == ring.ConjugateInvariant {
		cImag = float64(0)
	}

	return
}

// AddConst adds the input constant (which can be a uint64, int64, float64 or complex128) to ct0 and returns the result in ctOut.
func (eval *evaluator) AddConst(ct0 *rlwe.Ciphertext, constant interface{}, ctOut *rlwe.Ciphertext) {

	var level = utils.MinInt(ct0.Level(), ctOut.Level())
	var scaledConst, scaledConstReal, scaledConstImag, qi uint64

	cReal, cImag, _ := eval.getConstAndScale(level, constant)

	ringQ := eval.params.RingQ()

	cf64 := ctOut.Scale.Float64()

	ctOut.MetaData = ct0.MetaData

	// Component wise addition of the following vector to the ciphertext:
	// [a + b*psi_qi^2, ....., a + b*psi_qi^2, a - b*psi_qi^2, ...., a - b*psi_qi^2] mod Qi
	// [{                  N/2                }{                N/2               }]
	// Which is equivalent outside of the NTT domain to adding a to the first coefficient of ct0 and b to the N/2-th coefficient of ct0.
	for i := 0; i < level+1; i++ {
		scaledConstReal, scaledConstImag, scaledConst = 0, 0, 0
		qi = ringQ.Tables[i].Modulus

		if cReal != 0 {
			scaledConstReal = scaleUpExact(cReal, cf64, qi)
			scaledConst = scaledConstReal
		}

		if cImag != 0 {
			scaledConstImag = ring.MRed(scaleUpExact(cImag, cf64, qi), ringQ.Tables[i].RootsForward[1], qi, ringQ.Tables[i].MRedParams)
			scaledConst = ring.CRed(scaledConst+scaledConstImag, qi)
		}

		p0tmp := ct0.Value[0].Coeffs[i]
		p1tmp := ctOut.Value[0].Coeffs[i]

		ring.AddScalarVec(p0tmp[:ringQ.N()>>1], p1tmp[:ringQ.N()>>1], scaledConst, qi)

		if cImag != 0 {
			scaledConst = ring.CRed(scaledConstReal+(qi-scaledConstImag), qi)
		}

		ring.AddScalarVec(p0tmp[ringQ.N()>>1:], p1tmp[ringQ.N()>>1:], scaledConst, qi)
	}
}

// MultByConstAndAdd multiplies ct0 by the input constant, and adds it to the receiver element (it does not modify the input
// element), e.g., ctOut(x) = ctOut(x) + ct0(x) * (a+bi). This functions removes the need of storing the intermediate value c(x) * (a+bi).
// This function will modify the level and the scale of the receiver element depending on the level and the scale of the input
// element and the type of the constant. The level of the receiver element will be set to min(input.level, receiver.level).
// The scale of the receiver element will be set to the scale that the input element would have after the multiplication by the constant.
func (eval *evaluator) MultByConstAndAdd(ct0 *rlwe.Ciphertext, constant interface{}, ctOut *rlwe.Ciphertext) {

	var level = utils.MinInt(ct0.Level(), ctOut.Level())

	// Forces a drop of ctOut level to ct0 level
	if ctOut.Level() > level {
		eval.DropLevel(ctOut, ctOut.Level()-level)
	}

	cReal, cImag, scale := eval.getConstAndScale(level, constant)

	var scaledConst, scaledConstReal, scaledConstImag uint64

	c0f64 := ct0.Scale.Float64()
	c1f64 := ctOut.Scale.Float64()

	ringQ := eval.params.RingQ()

	// If a scaling would be required to multiply by the constant,
	// it equalizes scales such that the scales match in the end.
	if scale != 1 {

		// If ctOut scaling is smaller than ct0's scale + the default scaling,
		// then brings ctOut scale to ct0's scale.
		if c1f64 < c0f64*scale {

			if scale := math.Floor((scale * c0f64) / c1f64); scale > 1 {

				eval.MultByConst(ctOut, scale, ctOut)

			}

			ctOut.MetaData = ct0.MetaData
			ctOut.Scale = ct0.Scale.Mul(rlwe.NewScale(scale))

			// If ctOut.scale > ((a+bi)*scale)*ct0(x), then it sets the scale to
			// bring c(x)*scale to the level of ctOut(x) scale
		} else if c1f64 > c0f64*scale {
			scale = c1f64 / c0f64
		}

		// If no scaling is required, then it sets the appropriate scale such that
		// ct0(x)*scale matches ctOut(x) scale without modifying ct0(x) scale.
	} else {

		if c1f64 > c0f64 {

			scale = c1f64 / c0f64

		} else if c0f64 > c1f64 {

			if scale := math.Floor(c0f64 / c1f64); scale > 1 {
				eval.MultByConst(ctOut, scale, ctOut)
			}

			ctOut.MetaData = ct0.MetaData
			ctOut.Scale = ct0.Scale
		}
	}

	// Component-wise multiplication of the following vector to the ciphertext:
	// [a + b*psi_qi^2, ....., a + b*psi_qi^2, a - b*psi_qi^2, ...., a - b*psi_qi^2] mod Qi
	// [{                  N/2                }{                N/2               }]
	// Which is equivalent outside of the NTT domain to adding a to the first coefficient of ct0 and b to the N/2-th coefficient of ct0.
	for i := 0; i < level+1; i++ {

		Table := ringQ.Tables[i]

		qi := Table.Modulus
		mredParams := Table.MRedParams
		bredParams := Table.BRedParams

		scaledConstReal = 0
		scaledConstImag = 0
		scaledConst = 0

		if cReal != 0 {
			scaledConstReal = scaleUpExact(cReal, scale, qi)
			scaledConst = scaledConstReal
		}

		if cImag != 0 {
			scaledConstImag = scaleUpExact(cImag, scale, qi)
			scaledConstImag = ring.MRed(scaledConstImag, Table.RootsForward[1], qi, mredParams)
			scaledConst = ring.CRed(scaledConst+scaledConstImag, qi)
		}

		scaledConst = ring.MForm(scaledConst, qi, bredParams)

		for u := range ct0.Value {
			p0tmp := ct0.Value[u].Coeffs[i]
			p1tmp := ctOut.Value[u].Coeffs[i]
			ring.MulScalarMontgomeryAndAddVec(p0tmp[:ringQ.N()>>1], p1tmp[:ringQ.N()>>1], scaledConst, qi, mredParams)
		}

		if cImag != 0 {
			scaledConst = ring.CRed(scaledConstReal+(qi-scaledConstImag), qi)
			scaledConst = ring.MForm(scaledConst, qi, bredParams)
		}

		for u := range ct0.Value {
			p0tmp := ct0.Value[u].Coeffs[i]
			p1tmp := ctOut.Value[u].Coeffs[i]
			ring.MulScalarMontgomeryAndAddVec(p0tmp[ringQ.N()>>1:], p1tmp[ringQ.N()>>1:], scaledConst, qi, mredParams)
		}
	}
}

// MultByConstNew multiplies ct0 by the input constant and returns the result in a newly created element.
// The scale of the output element will depend on the scale of the input element and the constant (if the constant
// needs to be scaled (its rational part is not zero)). The constant can be a uint64, int64, float64 or complex128.
func (eval *evaluator) MultByConstNew(ct0 *rlwe.Ciphertext, constant interface{}) (ctOut *rlwe.Ciphertext) {
	ctOut = NewCiphertext(eval.params, ct0.Degree(), ct0.Level())
	eval.MultByConst(ct0, constant, ctOut)
	return
}

// MultByConst multiplies ct0 by the input constant and returns the result in ctOut.
// The scale of the output element will depend on the scale of the input element and the constant (if the constant
// needs to be scaled (its rational part is not zero)). The constant can be a uint64, int64, float64 or complex128.
func (eval *evaluator) MultByConst(ct0 *rlwe.Ciphertext, constant interface{}, ctOut *rlwe.Ciphertext) {

	var level = utils.MinInt(ct0.Level(), ctOut.Level())

	cReal, cImag, scale := eval.getConstAndScale(level, constant)

	// Component wise multiplication of the following vector with the ciphertext:
	// [a + b*psi_qi^2, ....., a + b*psi_qi^2, a - b*psi_qi^2, ...., a - b*psi_qi^2] mod Qi
	// [{                  N/2                }{                N/2               }]
	// Which is equivalent outside of the NTT domain to adding a to the first coefficient of ct0 and b to the N/2-th coefficient of ct0.
	ringQ := eval.params.RingQ()
	var scaledConst, scaledConstReal, scaledConstImag uint64
	for i := 0; i < level+1; i++ {

		Table := ringQ.Tables[i]

		qi := Table.Modulus
		bredParams := Table.BRedParams
		mredParams := Table.MRedParams

		scaledConstReal = 0
		scaledConstImag = 0
		scaledConst = 0

		if cReal != 0 {
			scaledConstReal = scaleUpExact(cReal, scale, qi)
			scaledConst = scaledConstReal
		}

		if cImag != 0 {
			scaledConstImag = scaleUpExact(cImag, scale, qi)
			scaledConstImag = ring.MRed(scaledConstImag, Table.RootsForward[1], qi, mredParams)
			scaledConst = ring.CRed(scaledConst+scaledConstImag, qi)
		}

		scaledConst = ring.MForm(scaledConst, qi, bredParams)

		for u := range ct0.Value {
			p0tmp := ct0.Value[u].Coeffs[i]
			p1tmp := ctOut.Value[u].Coeffs[i]
			ring.MulScalarMontgomeryVec(p0tmp[:ringQ.N()>>1], p1tmp[:ringQ.N()>>1], scaledConst, qi, mredParams)
		}

		if cImag != 0 {
			scaledConst = ring.CRed(scaledConstReal+(qi-scaledConstImag), qi)
			scaledConst = ring.MForm(scaledConst, qi, bredParams)
		}

		for u := range ct0.Value {
			p0tmp := ct0.Value[u].Coeffs[i]
			p1tmp := ctOut.Value[u].Coeffs[i]
			ring.MulScalarMontgomeryVec(p0tmp[ringQ.N()>>1:], p1tmp[ringQ.N()>>1:], scaledConst, qi, mredParams)
		}
	}

	ctOut.MetaData = ct0.MetaData
	ctOut.Scale = ct0.Scale.Mul(rlwe.NewScale(scale))
}

// MultByGaussianInteger multiples the ct0 by the gaussian integer cReal + i*cImag and returns the result on ctOut.
// Accepted types for cReal and cImag are uint64, int64 and big.Int.
func (eval *evaluator) MultByGaussianInteger(ct0 *rlwe.Ciphertext, cReal, cImag interface{}, ctOut *rlwe.Ciphertext) {

	ringQ := eval.params.RingQ()

	level := utils.MinInt(ct0.Level(), ctOut.Level())
	var scaledConst, scaledConstReal, scaledConstImag uint64

	ctOut.MetaData = ct0.MetaData

	for i := 0; i < level+1; i++ {

		Table := ringQ.Tables[i]

		qi := Table.Modulus
		bredParams := Table.BRedParams
		mredParams := Table.MRedParams

		scaledConstReal = interfaceMod(cReal, qi)

		if eval.params.RingType() != ring.ConjugateInvariant {
			scaledConstImag = interfaceMod(cImag, qi)
		}

		scaledConst = scaledConstReal

		if scaledConstImag != 0 {
			scaledConstImag = ring.MRed(scaledConstImag, Table.RootsForward[1], qi, mredParams)
			scaledConst = ring.CRed(scaledConst+scaledConstImag, qi)
		}

		scaledConst = ring.MForm(scaledConst, qi, bredParams)

		for u := range ct0.Value {
			p0tmp := ct0.Value[u].Coeffs[i]
			p1tmp := ctOut.Value[u].Coeffs[i]
			ring.MulScalarMontgomeryVec(p0tmp[:ringQ.N()>>1], p1tmp[:ringQ.N()>>1], scaledConst, qi, mredParams)
		}

		if cImag != 0 {
			scaledConst = ring.CRed(scaledConstReal+(qi-scaledConstImag), qi)
			scaledConst = ring.MForm(scaledConst, qi, bredParams)
		}

		for u := range ct0.Value {
			p0tmp := ct0.Value[u].Coeffs[i]
			p1tmp := ctOut.Value[u].Coeffs[i]
			ring.MulScalarMontgomeryVec(p0tmp[ringQ.N()>>1:], p1tmp[ringQ.N()>>1:], scaledConst, qi, mredParams)
		}
	}
}

// MultByGaussianIntegerAndAdd multiples the ct0 by the gaussian integer cReal + i*cImag and adds the result on ctOut.
// Accepted types for cReal and cImag are uint64, int64 and big.Int.
func (eval *evaluator) MultByGaussianIntegerAndAdd(ct0 *rlwe.Ciphertext, cReal, cImag interface{}, ctOut *rlwe.Ciphertext) {

	ringQ := eval.params.RingQ()

	level := utils.MinInt(ct0.Level(), ctOut.Level())
	var scaledConst, scaledConstReal, scaledConstImag uint64

	for i := 0; i < level+1; i++ {

		Table := ringQ.Tables[i]

		qi := Table.Modulus
		bredParams := Table.BRedParams
		mredParams := Table.MRedParams

		scaledConstReal = interfaceMod(cReal, qi)

		if eval.params.RingType() != ring.ConjugateInvariant {
			scaledConstImag = interfaceMod(cImag, qi)
		}

		scaledConst = scaledConstReal

		if scaledConstImag != 0 {
			scaledConstImag = ring.MRed(scaledConstImag, Table.RootsForward[1], qi, mredParams)
			scaledConst = ring.CRed(scaledConst+scaledConstImag, qi)
		}

		scaledConst = ring.MForm(scaledConst, qi, bredParams)

		for u := range ct0.Value {
			p0tmp := ct0.Value[u].Coeffs[i]
			p1tmp := ctOut.Value[u].Coeffs[i]
			ring.MulScalarMontgomeryAndAddVec(p0tmp[:ringQ.N()>>1], p1tmp[:ringQ.N()>>1], scaledConst, qi, mredParams)
		}

		if cImag != 0 {
			scaledConst = ring.CRed(scaledConstReal+(qi-scaledConstImag), qi)
			scaledConst = ring.MForm(scaledConst, qi, bredParams)
		}

		for u := range ct0.Value {
			p0tmp := ct0.Value[u].Coeffs[i]
			p1tmp := ctOut.Value[u].Coeffs[i]
			ring.MulScalarMontgomeryAndAddVec(p0tmp[ringQ.N()>>1:], p1tmp[ringQ.N()>>1:], scaledConst, qi, mredParams)
		}
	}
}

// MultByiNew multiplies ct0 by the imaginary number i, and returns the result in a newly created element.
// It does not change the scale.
func (eval *evaluator) MultByiNew(ct0 *rlwe.Ciphertext) (ctOut *rlwe.Ciphertext) {

	if eval.params.RingType() == ring.ConjugateInvariant {
		panic("cannot MultByiNew: method not supported when params.RingType() == ring.ConjugateInvariant")
	}

	ctOut = NewCiphertext(eval.params, 1, ct0.Level())
	eval.MultByi(ct0, ctOut)
	return ctOut
}

// MultByi multiplies ct0 by the imaginary number i, and returns the result in ctOut.
// It does not change the scale.
func (eval *evaluator) MultByi(ct0 *rlwe.Ciphertext, ctOut *rlwe.Ciphertext) {

	if eval.params.RingType() == ring.ConjugateInvariant {
		panic("cannot MultByi: method not supported when params.RingType() == ring.ConjugateInvariant")
	}

	var level = utils.MinInt(ct0.Level(), ctOut.Level())
	ctOut.MetaData = ct0.MetaData

	ringQ := eval.params.RingQ()

	// Equivalent to a product by the monomial x^(n/2) outside of the NTT domain
	for i := 0; i < level+1; i++ {

		Table := ringQ.Tables[i]

		qi := Table.Modulus
		mredParams := Table.MRedParams
		imag := Table.RootsForward[1] // Psi^2

		for u := range ctOut.Value {
			p0tmp := ct0.Value[u].Coeffs[i]
			p1tmp := ctOut.Value[u].Coeffs[i]
			ring.MulScalarMontgomeryVec(p0tmp[:ringQ.N()>>1], p1tmp[:ringQ.N()>>1], imag, qi, mredParams)
		}

		imag = qi - imag

		for u := range ctOut.Value {
			p0tmp := ct0.Value[u].Coeffs[i]
			p1tmp := ctOut.Value[u].Coeffs[i]
			ring.MulScalarMontgomeryVec(p0tmp[ringQ.N()>>1:], p1tmp[ringQ.N()>>1:], imag, qi, mredParams)
		}
	}
}

// DivByiNew multiplies ct0 by the imaginary number 1/i = -i, and returns the result in a newly created element.
// It does not change the scale.
func (eval *evaluator) DivByiNew(ct0 *rlwe.Ciphertext) (ctOut *rlwe.Ciphertext) {

	if eval.params.RingType() == ring.ConjugateInvariant {
		panic("cannot DivByiNew: method not supported when params.RingType() == ring.ConjugateInvariant")
	}

	ctOut = NewCiphertext(eval.params, 1, ct0.Level())
	eval.DivByi(ct0, ctOut)
	return
}

// DivByi multiplies ct0 by the imaginary number 1/i = -i, and returns the result in ctOut.
// It does not change the scale.
func (eval *evaluator) DivByi(ct0 *rlwe.Ciphertext, ctOut *rlwe.Ciphertext) {

	if eval.params.RingType() == ring.ConjugateInvariant {
		panic("cannot DivByi: method not supported when params.RingType() == ring.ConjugateInvariant")
	}

	var level = utils.MinInt(ct0.Level(), ctOut.Level())

	ringQ := eval.params.RingQ()

	ctOut.MetaData = ct0.MetaData

	// Equivalent to a product by the monomial x^(3*n/2) outside of the NTT domain
	for i := 0; i < level+1; i++ {

		Table := ringQ.Tables[i]

		qi := Table.Modulus
		mredParams := Table.MRedParams
		imag := qi - Table.RootsForward[1] // -Psi^2

		for u := range ctOut.Value {
			p0tmp := ct0.Value[u].Coeffs[i]
			p1tmp := ctOut.Value[u].Coeffs[i]
			ring.MulScalarMontgomeryVec(p0tmp[:ringQ.N()>>1], p1tmp[:ringQ.N()>>1], imag, qi, mredParams)
		}

		imag = qi - imag // Psi^2

		for u := range ctOut.Value {
			p0tmp := ct0.Value[u].Coeffs[i]
			p1tmp := ctOut.Value[u].Coeffs[i]
			ring.MulScalarMontgomeryVec(p0tmp[ringQ.N()>>1:], p1tmp[ringQ.N()>>1:], imag, qi, mredParams)
		}
	}
}

// ScaleUpNew multiplies ct0 by 2^scale and sets its scale to its previous scale
// plus 2^n. It returns the result in a newly created element.
func (eval *evaluator) ScaleUpNew(ct0 *rlwe.Ciphertext, scale rlwe.Scale) (ctOut *rlwe.Ciphertext) {
	ctOut = NewCiphertext(eval.params, ct0.Degree(), ct0.Level())
	eval.ScaleUp(ct0, scale, ctOut)
	return
}

// ScaleUp multiplies ct0 by 2^scale and sets its scale to its previous scale
// plus 2^n. It returns the result in ctOut.
func (eval *evaluator) ScaleUp(ct0 *rlwe.Ciphertext, scale rlwe.Scale, ctOut *rlwe.Ciphertext) {
	eval.MultByConst(ct0, scale.Uint64(), ctOut)
	ctOut.MetaData = ct0.MetaData
	ctOut.Scale = ct0.Scale.Mul(scale)
}

// SetScale sets the scale of the ciphertext to the input scale (consumes a level)
func (eval *evaluator) SetScale(ct *rlwe.Ciphertext, scale rlwe.Scale) {

	eval.MultByConst(ct, scale.Float64()/ct.Scale.Float64(), ct)
	if err := eval.Rescale(ct, scale, ct); err != nil {
		panic(err)
	}
	ct.Scale = scale
}

// DropLevelNew reduces the level of ct0 by levels and returns the result in a newly created element.
// No rescaling is applied during this procedure.
func (eval *evaluator) DropLevelNew(ct0 *rlwe.Ciphertext, levels int) (ctOut *rlwe.Ciphertext) {
	ctOut = ct0.CopyNew()
	eval.DropLevel(ctOut, levels)
	return
}

// DropLevel reduces the level of ct0 by levels and returns the result in ct0.
// No rescaling is applied during this procedure.
func (eval *evaluator) DropLevel(ct0 *rlwe.Ciphertext, levels int) {
	ct0.Resize(ct0.Degree(), ct0.Level()-levels)
}

// RescaleNew divides ct0 by the last modulus in the moduli chain, and repeats this
// procedure (consuming one level each time) until the scale reaches the original scale or before it goes below it, and returns the result
// in a newly created element. Since all the moduli in the moduli chain are generated to be close to the
// original scale, this procedure is equivalent to dividing the input element by the scale and adding
// some error.
// Returns an error if "threshold <= 0", ct.scale = 0, ct.Level() = 0, ct.IsNTT() != true
func (eval *evaluator) RescaleNew(ct0 *rlwe.Ciphertext, minScale rlwe.Scale) (ctOut *rlwe.Ciphertext, err error) {

	ctOut = NewCiphertext(eval.params, ct0.Degree(), ct0.Level())

	return ctOut, eval.Rescale(ct0, minScale, ctOut)
}

// Rescale divides ct0 by the last modulus in the moduli chain, and repeats this
// procedure (consuming one level each time) until the scale reaches the original scale or before it goes below it, and returns the result
// in ctOut. Since all the moduli in the moduli chain are generated to be close to the
// original scale, this procedure is equivalent to dividing the input element by the scale and adding
// some error.
// Returns an error if "minScale <= 0", ct.scale = 0, ct.Level() = 0, ct.IsNTT() != true or if ct.Leve() != ctOut.Level()
func (eval *evaluator) Rescale(ctIn *rlwe.Ciphertext, minScale rlwe.Scale, ctOut *rlwe.Ciphertext) (err error) {

	ringQ := eval.params.RingQ()

	if minScale.Cmp(rlwe.NewScale(0)) != 1 {
		return errors.New("cannot Rescale: minScale is <0")
	}

	minScale = minScale.Div(rlwe.NewScale(2))

	if ctIn.Scale.Cmp(rlwe.NewScale(0)) != 1 {
		return errors.New("cannot Rescale: ciphertext scale is <0")
	}

	if ctIn.Level() == 0 {
		return errors.New("cannot Rescale: input Ciphertext already at level 0")
	}

	if ctOut.Degree() != ctIn.Degree() {
		return errors.New("cannot Rescale: ctIn.Degree() != ctOut.Degree()")
	}

	ctOut.MetaData = ctIn.MetaData

	currentLevel := ctIn.Level()

	// Divides the scale by each moduli of the modulus chain as long as the scale isn't smaller than minScale/2
	// or until the output Level() would be zero
	var nbRescales int
	for currentLevel >= 0 {

		scale := ctOut.Scale.Div(rlwe.NewScale(ringQ.Tables[currentLevel].Modulus))

		if scale.Cmp(minScale) == -1 {
			break
		}

		ctOut.Scale = scale

		nbRescales++
		currentLevel--
	}

	if nbRescales > 0 {
		level := ctIn.Level()
		for i := range ctOut.Value {
			ringQ.DivRoundByLastModulusManyNTTLvl(level, nbRescales, ctIn.Value[i], eval.buffQ[0], ctOut.Value[i])
		}
		ctOut.Resize(ctOut.Degree(), level-nbRescales)
	} else {
		if ctIn != ctOut {
			ctOut.Copy(ctIn)
		}
	}

	return nil
}

// MulNew multiplies ctIn with op1 without relinearization and returns the result in a newly created element.
// The procedure will panic if either ctIn.Degree or op1.Degree > 1.
func (eval *evaluator) MulNew(ctIn *rlwe.Ciphertext, op1 rlwe.Operand) (ctOut *rlwe.Ciphertext) {
	ctOut = NewCiphertext(eval.params, ctIn.Degree()+op1.Degree(), utils.MinInt(ctIn.Level(), op1.Level()))
	eval.mulRelin(ctIn, op1, false, ctOut)
	return
}

// Mul multiplies ctIn with op1 without relinearization and returns the result in ctOut.
// The procedure will panic if either ctIn or op1 are have a degree higher than 1.
// The procedure will panic if ctOut.Degree != ctIn.Degree + op1.Degree.
func (eval *evaluator) Mul(ctIn *rlwe.Ciphertext, op1 rlwe.Operand, ctOut *rlwe.Ciphertext) {
	eval.mulRelin(ctIn, op1, false, ctOut)
}

// MulRelinNew multiplies ctIn with op1 with relinearization and returns the result in a newly created element.
// The procedure will panic if either ctIn.Degree or op1.Degree > 1.
// The procedure will panic if the evaluator was not created with an relinearization key.
func (eval *evaluator) MulRelinNew(ctIn *rlwe.Ciphertext, op1 rlwe.Operand) (ctOut *rlwe.Ciphertext) {
	ctOut = NewCiphertext(eval.params, 1, utils.MinInt(ctIn.Level(), op1.Level()))
	eval.mulRelin(ctIn, op1, true, ctOut)
	return
}

// MulRelin multiplies ctIn with op1 with relinearization and returns the result in ctOut.
// The procedure will panic if either ctIn.Degree or op1.Degree > 1.
// The procedure will panic if ctOut.Degree != ctIn.Degree + op1.Degree.
// The procedure will panic if the evaluator was not created with an relinearization key.
func (eval *evaluator) MulRelin(ctIn *rlwe.Ciphertext, op1 rlwe.Operand, ctOut *rlwe.Ciphertext) {
	eval.mulRelin(ctIn, op1, true, ctOut)
}

func (eval *evaluator) mulRelin(ctIn *rlwe.Ciphertext, op1 rlwe.Operand, relin bool, ctOut *rlwe.Ciphertext) {

	level := utils.MinInt(utils.MinInt(ctIn.Level(), op1.Level()), ctOut.Level())

	if ctIn.Degree()+op1.Degree() > 2 {
		panic("cannot MulRelin: the sum of the input elements' total degree cannot be larger than 2")
	}

	ctOut.MetaData = ctIn.MetaData
	ctOut.Scale = ctIn.Scale.Mul(op1.GetScale())

	ringQ := eval.params.RingQ()

	var c00, c01, c0, c1, c2 *ring.Poly

	// Case Ciphertext (x) Ciphertext
	if ctIn.Degree() == 1 && op1.Degree() == 1 {

		c00 = eval.buffQ[0]
		c01 = eval.buffQ[1]

		c0 = ctOut.Value[0]
		c1 = ctOut.Value[1]

		if !relin {
			ctOut.El().Resize(2, level)
			c2 = ctOut.Value[2]
		} else {
			ctOut.El().Resize(1, level)
			c2 = eval.buffQ[2]
		}

		eval.checkBinary(ctIn, op1, ctOut, ctOut.Degree())

		// Avoid overwriting if the second input is the output
		var tmp0, tmp1 *rlwe.Ciphertext
		if op1.El() == ctOut.El() {
			tmp0, tmp1 = op1.El(), ctIn.El()
		} else {
			tmp0, tmp1 = ctIn.El(), op1.El()
		}

		ringQ.MFormLvl(level, tmp0.Value[0], c00)
		ringQ.MFormLvl(level, tmp0.Value[1], c01)

		if ctIn.El() == op1.El() { // squaring case
			ringQ.MulCoeffsMontgomeryLvl(level, c00, tmp1.Value[0], c0) // c0 = c[0]*c[0]
			ringQ.MulCoeffsMontgomeryLvl(level, c01, tmp1.Value[1], c2) // c2 = c[1]*c[1]
			ringQ.MulCoeffsMontgomeryLvl(level, c00, tmp1.Value[1], c1) // c1 = 2*c[0]*c[1]
			ringQ.AddLvl(level, c1, c1, c1)

		} else { // regular case
			ringQ.MulCoeffsMontgomeryLvl(level, c00, tmp1.Value[0], c0) // c0 = c0[0]*c0[0]
			ringQ.MulCoeffsMontgomeryLvl(level, c01, tmp1.Value[1], c2) // c2 = c0[1]*c1[1]
			ringQ.MulCoeffsMontgomeryLvl(level, c00, tmp1.Value[1], c1)
			ringQ.MulCoeffsMontgomeryAndAddLvl(level, c01, tmp1.Value[0], c1) // c1 = c0[0]*c1[1] + c0[1]*c1[0]
		}

		if relin {

			if eval.Rlk == nil {
				panic("cannot MulRelin: relinearization key is missing")
			}

			tmpCt := &rlwe.Ciphertext{Value: []*ring.Poly{eval.BuffQP[1].Q, eval.BuffQP[2].Q}}
			tmpCt.IsNTT = true

			eval.GadgetProduct(level, c2, eval.Rlk.Keys[0].GadgetCiphertext, tmpCt)
			ringQ.AddLvl(level, c0, tmpCt.Value[0], ctOut.Value[0])
			ringQ.AddLvl(level, c1, tmpCt.Value[1], ctOut.Value[1])
		}

		// Case Plaintext (x) Ciphertext or Ciphertext (x) Plaintext
	} else {

		eval.checkBinary(ctIn, op1, ctOut, ctOut.Degree())

		var c0 *ring.Poly
		var c1 []*ring.Poly
		if ctIn.Degree() == 0 {
			c0 = eval.buffQ[0]
			ringQ.MFormLvl(level, ctIn.Value[0], c0)
			c1 = op1.El().Value

		} else {
			c0 = eval.buffQ[0]
			ringQ.MFormLvl(level, op1.El().Value[0], c0)
			c1 = ctIn.Value
		}

		ctOut.El().Resize(ctIn.Degree()+op1.Degree(), level)

		for i := range c1 {
			ringQ.MulCoeffsMontgomeryLvl(level, c0, c1[i], ctOut.Value[i])
		}
	}
}

// MulAndAdd multiplies ctIn with op1 without relinearization and adds the result on ctOut.
// User must ensure that ctOut.scale <= ctIn.scale * op1.scale.
// If ctOut.scale < ctIn.scale * op1.scale, then scales up ctOut before adding the result.
// The procedure will panic if either ctIn or op1 are have a degree higher than 1.
// The procedure will panic if ctOut.Degree != ctIn.Degree + op1.Degree.
// The procedure will panic if ctOut = ctIn or op1.
func (eval *evaluator) MulAndAdd(ctIn *rlwe.Ciphertext, op1 rlwe.Operand, ctOut *rlwe.Ciphertext) {
	eval.mulRelinAndAdd(ctIn, op1, false, ctOut)
}

// MulRelinAndAdd multiplies ctIn with op1 with relinearization and adds the result on ctOut.
// User must ensure that ctOut.scale <= ctIn.scale * op1.scale.
// If ctOut.scale < ctIn.scale * op1.scale, then scales up ctOut before adding the result.
// The procedure will panic if either ctIn.Degree or op1.Degree > 1.
// The procedure will panic if ctOut.Degree != ctIn.Degree + op1.Degree.
// The procedure will panic if the evaluator was not created with an relinearization key.
// The procedure will panic if ctOut = ctIn or op1.
func (eval *evaluator) MulRelinAndAdd(ctIn *rlwe.Ciphertext, op1 rlwe.Operand, ctOut *rlwe.Ciphertext) {
	eval.mulRelinAndAdd(ctIn, op1, true, ctOut)
}

func (eval *evaluator) mulRelinAndAdd(ctIn *rlwe.Ciphertext, op1 rlwe.Operand, relin bool, ctOut *rlwe.Ciphertext) {

	eval.checkBinary(ctIn, op1, ctOut, utils.MaxInt(ctIn.Degree(), op1.Degree()))

	level := utils.MinInt(utils.MinInt(ctIn.Level(), op1.Level()), ctOut.Level())

	if ctIn.Degree()+op1.Degree() > 2 {
		panic("cannot MulRelinAndAdd: the sum of the input elements' degree cannot be larger than 2")
	}

	if ctIn.El() == ctOut.El() || op1.El() == ctOut.El() {
		panic("cannot MulRelinAndAdd: ctOut must be different from op0 and op1")
	}

	c0f64 := ctIn.Scale.Float64()
	c1f64 := op1.GetScale().Float64()
	c2f64 := ctOut.Scale.Float64()

	resScale := c0f64 * c1f64

	if c2f64 < resScale {
		eval.MultByConst(ctOut, math.Round(resScale/c2f64), ctOut)
		ctOut.Scale = rlwe.NewScale(resScale)
	}

	ringQ := eval.params.RingQ()

	var c00, c01, c0, c1, c2 *ring.Poly

	// Case Ciphertext (x) Ciphertext
	if ctIn.Degree() == 1 && op1.Degree() == 1 {

		c00 = eval.buffQ[0]
		c01 = eval.buffQ[1]

		c0 = ctOut.Value[0]
		c1 = ctOut.Value[1]

		if !relin {
			ctOut.El().Resize(2, level)
			c2 = ctOut.Value[2]
		} else {
			// No resize here since we add on ctOut
			c2 = eval.buffQ[2]
		}

		tmp0, tmp1 := ctIn.El(), op1.El()

		ringQ.MFormLvl(level, tmp0.Value[0], c00)
		ringQ.MFormLvl(level, tmp0.Value[1], c01)

		ringQ.MulCoeffsMontgomeryAndAddLvl(level, c00, tmp1.Value[0], c0) // c0 += c[0]*c[0]
		ringQ.MulCoeffsMontgomeryAndAddLvl(level, c00, tmp1.Value[1], c1) // c1 += c[0]*c[1]
		ringQ.MulCoeffsMontgomeryAndAddLvl(level, c01, tmp1.Value[0], c1) // c1 += c[1]*c[0]

		if relin {

			if eval.Rlk == nil {
				panic("cannot MulRelinAndAdd: relinearization key is missing")
			}

			ringQ.MulCoeffsMontgomeryLvl(level, c01, tmp1.Value[1], c2) // c2 += c[1]*c[1]

			tmpCt := &rlwe.Ciphertext{Value: []*ring.Poly{eval.BuffQP[1].Q, eval.BuffQP[2].Q}}
			tmpCt.IsNTT = true

			eval.GadgetProduct(level, c2, eval.Rlk.Keys[0].GadgetCiphertext, tmpCt)
			ringQ.AddLvl(level, c0, tmpCt.Value[0], c0)
			ringQ.AddLvl(level, c1, tmpCt.Value[1], c1)
		} else {
			ringQ.MulCoeffsMontgomeryAndAddLvl(level, c01, tmp1.Value[1], c2) // c2 += c[1]*c[1]
		}

		// Case Plaintext (x) Ciphertext or Ciphertext (x) Plaintext
	} else {

		if ctOut.Degree() < ctIn.Degree() {
			ctOut.Resize(ctIn.Degree(), level)
		}

		c00 := eval.buffQ[0]

		ringQ.MFormLvl(level, op1.El().Value[0], c00)
		for i := range ctIn.Value {
			ringQ.MulCoeffsMontgomeryAndAddLvl(level, ctIn.Value[i], c00, ctOut.Value[i])
		}
	}
}

// RelinearizeNew applies the relinearization procedure on ct0 and returns the result in a newly
// created Ciphertext. The input Ciphertext must be of degree two.
func (eval *evaluator) RelinearizeNew(ct0 *rlwe.Ciphertext) (ctOut *rlwe.Ciphertext) {
	ctOut = NewCiphertext(eval.params, 1, ct0.Level())
	eval.Relinearize(ct0, ctOut)
	return
}

// SwitchKeysNew re-encrypts ct0 under a different key and returns the result in a newly created element.
// It requires a SwitchingKey, which is computed from the key under which the Ciphertext is currently encrypted,
// and the key under which the Ciphertext will be re-encrypted.
func (eval *evaluator) SwitchKeysNew(ct0 *rlwe.Ciphertext, switchingKey *rlwe.SwitchingKey) (ctOut *rlwe.Ciphertext) {
	ctOut = NewCiphertext(eval.params, ct0.Degree(), ct0.Level())
	eval.SwitchKeys(ct0, switchingKey, ctOut)
	return
}

// RotateNew rotates the columns of ct0 by k positions to the left, and returns the result in a newly created element.
// If the provided element is a Ciphertext, a key-switching operation is necessary and a rotation key for the specific rotation needs to be provided.
func (eval *evaluator) RotateNew(ct0 *rlwe.Ciphertext, k int) (ctOut *rlwe.Ciphertext) {
	ctOut = NewCiphertext(eval.params, ct0.Degree(), ct0.Level())
	eval.Rotate(ct0, k, ctOut)
	return
}

// Rotate rotates the columns of ct0 by k positions to the left and returns the result in ctOut.
// If the provided element is a Ciphertext, a key-switching operation is necessary and a rotation key for the specific rotation needs to be provided.
func (eval *evaluator) Rotate(ct0 *rlwe.Ciphertext, k int, ctOut *rlwe.Ciphertext) {
	eval.Automorphism(ct0, eval.params.GaloisElementForColumnRotationBy(k), ctOut)
}

// ConjugateNew conjugates ct0 (which is equivalent to a row rotation) and returns the result in a newly
// created element. If the provided element is a Ciphertext, a key-switching operation is necessary and a rotation key
// for the row rotation needs to be provided.
func (eval *evaluator) ConjugateNew(ct0 *rlwe.Ciphertext) (ctOut *rlwe.Ciphertext) {

	if eval.params.RingType() == ring.ConjugateInvariant {
		panic("cannot ConjugateNew: method is not supported when params.RingType() == ring.ConjugateInvariant")
	}

	ctOut = NewCiphertext(eval.params, ct0.Degree(), ct0.Level())
	eval.Conjugate(ct0, ctOut)
	return
}

// Conjugate conjugates ct0 (which is equivalent to a row rotation) and returns the result in ctOut.
// If the provided element is a Ciphertext, a key-switching operation is necessary and a rotation key for the row rotation needs to be provided.
func (eval *evaluator) Conjugate(ct0 *rlwe.Ciphertext, ctOut *rlwe.Ciphertext) {

	if eval.params.RingType() == ring.ConjugateInvariant {
		panic("cannot Conjugate: method is not supported when params.RingType() == ring.ConjugateInvariant")
	}

	eval.Automorphism(ct0, eval.params.GaloisElementForRowRotation(), ctOut)
}

func (eval *evaluator) RotateHoistedNoModDownNew(level int, rotations []int, c0 *ring.Poly, c2DecompQP []ringqp.Poly) (cOut map[int]rlwe.CiphertextQP) {
	ringQ := eval.params.RingQ()
	ringP := eval.params.RingP()
	cOut = make(map[int]rlwe.CiphertextQP)
	for _, i := range rotations {
		if i != 0 {
			cOut[i] = rlwe.CiphertextQP{Value: [2]ringqp.Poly{{Q: ringQ.NewPolyLvl(level), P: ringP.NewPoly()}, {Q: ringQ.NewPolyLvl(level), P: ringP.NewPoly()}}, MetaData: rlwe.MetaData{IsNTT: true}}
			eval.AutomorphismHoistedNoModDown(level, c0, c2DecompQP, eval.params.GaloisElementForColumnRotationBy(i), cOut[i])
		}
	}

	return
}

// ShallowCopy creates a shallow copy of this evaluator in which all the read-only data-structures are
// shared with the receiver and the temporary buffers are reallocated. The receiver and the returned
// Evaluators can be used concurrently.
func (eval *evaluator) ShallowCopy() Evaluator {
	return &evaluator{
		evaluatorBase:    eval.evaluatorBase,
		Evaluator:        eval.Evaluator.ShallowCopy(),
		evaluatorBuffers: newEvaluatorBuffers(eval.evaluatorBase),
	}
}

// WithKey creates a shallow copy of the receiver Evaluator for which the new EvaluationKey is evaluationKey
// and where the temporary buffers are shared. The receiver and the returned Evaluators cannot be used concurrently.
func (eval *evaluator) WithKey(evaluationKey rlwe.EvaluationKey) Evaluator {
	return &evaluator{
		Evaluator:        eval.Evaluator.WithKey(&evaluationKey),
		evaluatorBase:    eval.evaluatorBase,
		evaluatorBuffers: eval.evaluatorBuffers,
	}
}
