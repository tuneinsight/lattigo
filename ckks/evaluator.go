package ckks

import (
	"errors"
	"math"

	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/rlwe"
	"github.com/tuneinsight/lattigo/v4/rlwe/ringqp"
	"github.com/tuneinsight/lattigo/v4/utils"
)

// Operand is a common interface for Ciphertext and Plaintext types.
type Operand interface {
	El() *rlwe.Ciphertext
	Degree() int
	Level() int
	Scale() float64
	SetScale(float64)
}

// Evaluator is an interface implementing the methods to conduct homomorphic operations between ciphertext and/or plaintexts.
type Evaluator interface {
	// ========================
	// === Basic Arithmetic ===
	// ========================

	// Addition
	Add(ctIn *Ciphertext, op1 Operand, ctOut *Ciphertext)
	AddNoMod(ctIn *Ciphertext, op1 Operand, ctOut *Ciphertext)
	AddNew(ctIn *Ciphertext, op1 Operand) (ctOut *Ciphertext)
	AddNoModNew(ctIn *Ciphertext, op1 Operand) (ctOut *Ciphertext)

	// Subtraction
	Sub(ctIn *Ciphertext, op1 Operand, ctOut *Ciphertext)
	SubNoMod(ctIn *Ciphertext, op1 Operand, ctOut *Ciphertext)
	SubNew(ctIn *Ciphertext, op1 Operand) (ctOut *Ciphertext)
	SubNoModNew(ctIn *Ciphertext, op1 Operand) (ctOut *Ciphertext)

	// Negation
	Neg(ctIn *Ciphertext, ctOut *Ciphertext)
	NegNew(ctIn *Ciphertext) (ctOut *Ciphertext)

	// Constant Addition
	AddConstNew(ctIn *Ciphertext, constant interface{}) (ctOut *Ciphertext)
	AddConst(ctIn *Ciphertext, constant interface{}, ctOut *Ciphertext)

	// Constant Multiplication
	MultByConstNew(ctIn *Ciphertext, constant interface{}) (ctOut *Ciphertext)
	MultByConst(ctIn *Ciphertext, constant interface{}, ctOut *Ciphertext)
	MultByGaussianInteger(ctIn *Ciphertext, cReal, cImag interface{}, ctOut *Ciphertext)

	// Constant Multiplication with Addition
	MultByConstAndAdd(ctIn *Ciphertext, constant interface{}, ctOut *Ciphertext)
	MultByGaussianIntegerAndAdd(ctIn *Ciphertext, cReal, cImag interface{}, ctOut *Ciphertext)

	// Multiplication by the imaginary unit
	MultByiNew(ctIn *Ciphertext) (ctOut *Ciphertext)
	MultByi(ctIn *Ciphertext, ctOut *Ciphertext)
	DivByiNew(ctIn *Ciphertext) (ctOut *Ciphertext)
	DivByi(ctIn *Ciphertext, ctOut *Ciphertext)

	// Conjugation
	ConjugateNew(ctIn *Ciphertext) (ctOut *Ciphertext)
	Conjugate(ctIn *Ciphertext, ctOut *Ciphertext)

	// Multiplication
	Mul(ctIn *Ciphertext, op1 Operand, ctOut *Ciphertext)
	MulNew(ctIn *Ciphertext, op1 Operand) (ctOut *Ciphertext)
	MulRelin(ctIn *Ciphertext, op1 Operand, ctOut *Ciphertext)
	MulRelinNew(ctIn *Ciphertext, op1 Operand) (ctOut *Ciphertext)

	MulAndAdd(ctIn *Ciphertext, op1 Operand, ctOut *Ciphertext)
	MulRelinAndAdd(ctIn *Ciphertext, op1 Operand, ctOut *Ciphertext)

	// Slot Rotations
	RotateNew(ctIn *Ciphertext, k int) (ctOut *Ciphertext)
	Rotate(ctIn *Ciphertext, k int, ctOut *Ciphertext)
	RotateHoistedNew(ctIn *Ciphertext, rotations []int) (ctOut map[int]*Ciphertext)
	RotateHoisted(ctIn *Ciphertext, rotations []int, ctOut map[int]*Ciphertext)
	RotateHoistedNoModDownNew(level int, rotations []int, c0 *ring.Poly, c2DecompQP []ringqp.Poly) (cOut map[int][2]ringqp.Poly)

	// ===========================
	// === Advanced Arithmetic ===
	// ===========================

	// Multiplication by 2^{s}
	MulByPow2New(ctIn *Ciphertext, pow2 int) (ctOut *Ciphertext)
	MulByPow2(ctIn *Ciphertext, pow2 int, ctOut *Ciphertext)

	// Exponentiation
	PowerOf2(ctIn *Ciphertext, logPow2 int, ctOut *Ciphertext)
	Power(ctIn *Ciphertext, degree int, ctOut *Ciphertext)
	PowerNew(ctIn *Ciphertext, degree int) (ctOut *Ciphertext)

	// Polynomial evaluation
	EvaluatePoly(input interface{}, pol *Polynomial, targetScale float64) (ctOut *Ciphertext, err error)
	EvaluatePolyVector(input interface{}, pols []*Polynomial, encoder Encoder, slotIndex map[int][]int, targetScale float64) (ctOut *Ciphertext, err error)

	// Inversion
	InverseNew(ctIn *Ciphertext, steps int) (ctOut *Ciphertext)

	// Linear Transformations
	LinearTransformNew(ctIn *Ciphertext, linearTransform interface{}) (ctOut []*Ciphertext)
	LinearTransform(ctIn *Ciphertext, linearTransform interface{}, ctOut []*Ciphertext)
	MultiplyByDiagMatrix(ctIn *Ciphertext, matrix LinearTransform, c2DecompQP []ringqp.Poly, ctOut *Ciphertext)
	MultiplyByDiagMatrixBSGS(ctIn *Ciphertext, matrix LinearTransform, c2DecompQP []ringqp.Poly, ctOut *Ciphertext)

	// Inner sum
	InnerSumLog(ctIn *Ciphertext, batch, n int, ctOut *Ciphertext)
	InnerSum(ctIn *Ciphertext, batch, n int, ctOut *Ciphertext)
	Average(ctIn *Ciphertext, batch int, ctOut *Ciphertext)

	// Replication (inverse of Inner sum)
	ReplicateLog(ctIn *Ciphertext, batch, n int, ctOut *Ciphertext)
	Replicate(ctIn *Ciphertext, batch, n int, ctOut *Ciphertext)

	// Trace
	Trace(ctIn *Ciphertext, logSlots int, ctOut *Ciphertext)
	TraceNew(ctIn *Ciphertext, logSlots int) (ctOut *Ciphertext)

	// =============================
	// === Ciphertext Management ===
	// =============================

	// Key-Switching
	SwitchKeysNew(ctIn *Ciphertext, switchingKey *rlwe.SwitchingKey) (ctOut *Ciphertext)
	SwitchKeys(ctIn *Ciphertext, switchingKey *rlwe.SwitchingKey, ctOut *Ciphertext)

	// Degree Management
	RelinearizeNew(ctIn *Ciphertext) (ctOut *Ciphertext)
	Relinearize(ctIn *Ciphertext, ctOut *Ciphertext)

	// Scale Management
	ScaleUpNew(ctIn *Ciphertext, scale float64) (ctOut *Ciphertext)
	ScaleUp(ctIn *Ciphertext, scale float64, ctOut *Ciphertext)
	SetScale(ctIn *Ciphertext, scale float64)
	Rescale(ctIn *Ciphertext, minScale float64, ctOut *Ciphertext) (err error)

	// Level Management
	DropLevelNew(ctIn *Ciphertext, levels int) (ctOut *Ciphertext)
	DropLevel(ctIn *Ciphertext, levels int)

	// Modular Overflow Management
	ReduceNew(ctIn *Ciphertext) (ctOut *Ciphertext)
	Reduce(ctIn *Ciphertext, ctOut *Ciphertext) error

	// ==============
	// === Others ===
	// ==============
	GetRLWEEvaluator() *rlwe.Evaluator
	BuffQ() [3]*ring.Poly
	BuffCt() *Ciphertext
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
	buffQ  [3]*ring.Poly // Memory buffer in order: for MForm(c0), MForm(c1), c2
	buffCt *Ciphertext   // Memory buffer for ciphertexts that need to be scaled up (to be eventually removed)
}

// BuffQ returns a pointer to the internal memory buffer buffQ.
func (eval *evaluator) BuffQ() [3]*ring.Poly {
	return eval.buffQ
}

// BuffCt returns a pointer to the internal memory buffer buffCt.
func (eval *evaluator) BuffCt() *Ciphertext {
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
	buff.buffCt = NewCiphertext(params, 2, params.MaxLevel(), params.DefaultScale())
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

func (eval *evaluator) checkBinary(op0, op1, opOut Operand, opOutMinDegree int) {
	if op0 == nil || op1 == nil || opOut == nil {
		panic("cannot checkBinary: operands cannot be nil")
	}

	if op0.Degree()+op1.Degree() == 0 {
		panic("cannot checkBinary: operands cannot be both plaintext")
	}

	if opOut.Degree() < opOutMinDegree {
		panic("cannot checkBinary: receiver operand degree is too small")
	}

	for _, pol := range op0.El().Value {
		if !pol.IsNTT {
			panic("cannot checkBinary: op0 must be in NTT")
		}
	}

	for _, pol := range op1.El().Value {
		if !pol.IsNTT {
			panic("cannot checkBinary: op1 must be in NTT")
		}
	}
}

func (eval *evaluator) newCiphertextBinary(op0, op1 Operand) (ctOut *Ciphertext) {

	maxDegree := utils.MaxInt(op0.Degree(), op1.Degree())
	maxScale := utils.MaxFloat64(op0.Scale(), op1.Scale())
	minLevel := utils.MinInt(op0.Level(), op1.Level())

	return NewCiphertext(eval.params, maxDegree, minLevel, maxScale)
}

// Add adds op1 to ctIn and returns the result in ctOut.
func (eval *evaluator) Add(ctIn *Ciphertext, op1 Operand, ctOut *Ciphertext) {
	eval.checkBinary(ctIn, op1, ctOut, utils.MaxInt(ctIn.Degree(), op1.Degree()))
	eval.evaluateInPlace(ctIn, op1, ctOut, eval.params.RingQ().AddLvl)
}

// AddNoMod adds op1 to ctIn and returns the result in ctOut, without modular reduction.
func (eval *evaluator) AddNoMod(ctIn *Ciphertext, op1 Operand, ctOut *Ciphertext) {
	eval.checkBinary(ctIn, op1, ctOut, utils.MaxInt(ctIn.Degree(), op1.Degree()))
	eval.evaluateInPlace(ctIn, op1, ctOut, eval.params.RingQ().AddNoModLvl)
}

// AddNew adds op1 to ctIn and returns the result in a newly created element.
func (eval *evaluator) AddNew(ctIn *Ciphertext, op1 Operand) (ctOut *Ciphertext) {
	ctOut = eval.newCiphertextBinary(ctIn, op1)
	eval.Add(ctIn, op1, ctOut)
	return
}

// AddNoModNew adds op1 to ctIn without modular reduction, and returns the result in a newly created element.
func (eval *evaluator) AddNoModNew(ctIn *Ciphertext, op1 Operand) (ctOut *Ciphertext) {
	ctOut = eval.newCiphertextBinary(ctIn, op1)
	eval.AddNoMod(ctIn, op1, ctOut)
	return
}

// Sub subtracts op1 from ctIn and returns the result in ctOut.
func (eval *evaluator) Sub(ctIn *Ciphertext, op1 Operand, ctOut *Ciphertext) {

	eval.checkBinary(ctIn, op1, ctOut, utils.MaxInt(ctIn.Degree(), op1.Degree()))

	eval.evaluateInPlace(ctIn, op1, ctOut, eval.params.RingQ().SubLvl)

	level := utils.MinInt(utils.MinInt(ctIn.Level(), op1.Level()), ctOut.Level())

	if ctIn.Degree() < op1.Degree() {
		for i := ctIn.Degree() + 1; i < op1.Degree()+1; i++ {
			eval.params.RingQ().NegLvl(level, ctOut.Value[i], ctOut.Value[i])
		}
	}

}

// SubNoMod subtracts op1 from ctIn and returns the result in ctOut, without modular reduction.
func (eval *evaluator) SubNoMod(ctIn *Ciphertext, op1 Operand, ctOut *Ciphertext) {

	eval.checkBinary(ctIn, op1, ctOut, utils.MaxInt(ctIn.Degree(), op1.Degree()))

	eval.evaluateInPlace(ctIn, op1, ctOut, eval.params.RingQ().SubNoModLvl)

	level := utils.MinInt(utils.MinInt(ctIn.Level(), op1.Level()), ctOut.Level())

	if ctIn.Degree() < op1.Degree() {
		for i := ctIn.Degree() + 1; i < op1.Degree()+1; i++ {
			eval.params.RingQ().NegLvl(level, ctOut.Value[i], ctOut.Value[i])
		}
	}

}

// SubNew subtracts op1 from ctIn and returns the result in a newly created element.
func (eval *evaluator) SubNew(ctIn *Ciphertext, op1 Operand) (ctOut *Ciphertext) {
	ctOut = eval.newCiphertextBinary(ctIn, op1)
	eval.Sub(ctIn, op1, ctOut)
	return
}

// SubNoModNew subtracts op1 from ctIn without modular reduction, and returns the result in a newly created element.
func (eval *evaluator) SubNoModNew(ctIn *Ciphertext, op1 Operand) (ctOut *Ciphertext) {
	ctOut = eval.newCiphertextBinary(ctIn, op1)
	eval.SubNoMod(ctIn, op1, ctOut)
	return
}

func (eval *evaluator) evaluateInPlace(c0, c1, ctOut Operand, evaluate func(int, *ring.Poly, *ring.Poly, *ring.Poly)) {

	var tmp0, tmp1 *rlwe.Ciphertext

	level := utils.MinInt(utils.MinInt(c0.Level(), c1.Level()), ctOut.Level())

	maxDegree := utils.MaxInt(c0.Degree(), c1.Degree())
	minDegree := utils.MinInt(c0.Degree(), c1.Degree())

	// Else resizes the receiver element
	ctOut.El().Resize(maxDegree, ctOut.Level())

	c0Scale := c0.Scale()
	c1Scale := c1.Scale()
	ctOutScale := ctOut.Scale()

	if ctOut.Level() > level {
		eval.DropLevel(&Ciphertext{ctOut.El(), ctOutScale}, ctOut.Level()-utils.MinInt(c0.Level(), c1.Level()))
	}

	// Checks whether or not the receiver element is the same as one of the input elements
	// and acts accordingly to avoid unnecessary element creation or element overwriting,
	// and scales properly the element before the evaluation.
	if ctOut == c0 {

		if c0Scale > c1Scale && math.Floor(c0Scale/c1Scale) > 1 {

			tmp1 = eval.buffCt.El()

			eval.MultByConst(&Ciphertext{c1.El(), c1Scale}, math.Floor(c0Scale/c1Scale), &Ciphertext{tmp1, ctOutScale})

		} else if c1Scale > c0Scale && math.Floor(c1Scale/c0Scale) > 1 {

			eval.MultByConst(&Ciphertext{c0.El(), c0Scale}, math.Floor(c1Scale/c0Scale), &Ciphertext{c0.El(), c0Scale})

			ctOut.SetScale(c1Scale)

			tmp1 = c1.El()

		} else {

			tmp1 = c1.El()
		}

		tmp0 = c0.El()

	} else if ctOut == c1 {

		if c1Scale > c0Scale && math.Floor(c1Scale/c0Scale) > 1 {

			tmp0 = eval.buffCt.El()

			eval.MultByConst(&Ciphertext{c0.El(), c0Scale}, math.Floor(c1Scale/c0Scale), &Ciphertext{tmp0, ctOutScale})

		} else if c0Scale > c1Scale && math.Floor(c0Scale/c1Scale) > 1 {

			eval.MultByConst(&Ciphertext{c1.El(), c1Scale}, math.Floor(c0Scale/c1Scale), &Ciphertext{ctOut.El(), ctOutScale})

			ctOut.SetScale(c0Scale)

			tmp0 = c0.El()

		} else {

			tmp0 = c0.El()
		}

		tmp1 = c1.El()

	} else {

		if c1Scale > c0Scale && math.Floor(c1Scale/c0Scale) > 1 {

			tmp0 = eval.buffCt.El()

			eval.MultByConst(&Ciphertext{c0.El(), c0Scale}, math.Floor(c1Scale/c0Scale), &Ciphertext{tmp0, ctOutScale})

			tmp1 = c1.El()

		} else if c0Scale > c1Scale && math.Floor(c0Scale/c1Scale) > 1 {

			tmp1 = eval.buffCt.El()

			eval.MultByConst(&Ciphertext{c1.El(), c1Scale}, math.Floor(c0Scale/c1Scale), &Ciphertext{tmp1, ctOutScale})

			tmp0 = c0.El()

		} else {
			tmp0 = c0.El()
			tmp1 = c1.El()
		}
	}

	for i := 0; i < minDegree+1; i++ {
		evaluate(level, tmp0.Value[i], tmp1.Value[i], ctOut.El().Value[i])
	}

	ctOut.SetScale(utils.MaxFloat64(c0Scale, c1Scale))

	// If the inputs degrees differ, it copies the remaining degree on the receiver.
	// Also checks that the receiver is not one of the inputs to avoid unnecessary work.

	if c0.Degree() > c1.Degree() && tmp0 != ctOut.El() {
		for i := minDegree + 1; i < maxDegree+1; i++ {
			ring.CopyValuesLvl(level, tmp0.Value[i], ctOut.El().Value[i])
		}
	} else if c1.Degree() > c0.Degree() && tmp1 != ctOut.El() {
		for i := minDegree + 1; i < maxDegree+1; i++ {
			ring.CopyValuesLvl(level, tmp1.Value[i], ctOut.El().Value[i])
		}
	}
}

// Neg negates the value of ct0 and returns the result in ctOut.
func (eval *evaluator) Neg(ct0 *Ciphertext, ctOut *Ciphertext) {

	level := utils.MinInt(ct0.Level(), ctOut.Level())

	if ct0.Degree() != ctOut.Degree() {
		panic("cannot Negate: invalid receiver Ciphertext does not match input Ciphertext degree")
	}

	for i := range ct0.Value {
		eval.params.RingQ().NegLvl(level, ct0.Value[i], ctOut.Value[i])
	}

	ctOut.scale = ct0.scale
}

// NegNew negates ct0 and returns the result in a newly created element.
func (eval *evaluator) NegNew(ct0 *Ciphertext) (ctOut *Ciphertext) {
	ctOut = NewCiphertext(eval.params, ct0.Degree(), ct0.Level(), ct0.scale)
	eval.Neg(ct0, ctOut)
	return
}

// AddConstNew adds the input constant (which can be a uint64, int64, float64 or complex128) to ct0 and returns the result in a new element.
func (eval *evaluator) AddConstNew(ct0 *Ciphertext, constant interface{}) (ctOut *Ciphertext) {
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
				scale = float64(eval.params.RingQ().Modulus[level])
			}
		}

		if cImag != 0 {
			valueInt := int64(cImag)
			valueFloat := cImag - float64(valueInt)

			if valueFloat != 0 {
				scale = float64(eval.params.RingQ().Modulus[level])
			}
		}

	case float64:
		cReal = constant
		cImag = float64(0)

		if cReal != 0 {
			valueInt := int64(cReal)
			valueFloat := cReal - float64(valueInt)

			if valueFloat != 0 {
				scale = float64(eval.params.RingQ().Modulus[level])
			}
		}

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
func (eval *evaluator) AddConst(ct0 *Ciphertext, constant interface{}, ctOut *Ciphertext) {

	var level = utils.MinInt(ct0.Level(), ctOut.Level())
	var scaledConst, scaledConstReal, scaledConstImag, qi uint64

	cReal, cImag, _ := eval.getConstAndScale(level, constant)

	ringQ := eval.params.RingQ()

	ctOut.scale = ct0.scale

	// Component wise addition of the following vector to the ciphertext:
	// [a + b*psi_qi^2, ....., a + b*psi_qi^2, a - b*psi_qi^2, ...., a - b*psi_qi^2] mod Qi
	// [{                  N/2                }{                N/2               }]
	// Which is equivalent outside of the NTT domain to adding a to the first coefficient of ct0 and b to the N/2-th coefficient of ct0.
	for i := 0; i < level+1; i++ {
		scaledConstReal, scaledConstImag, scaledConst = 0, 0, 0
		qi = ringQ.Modulus[i]

		if cReal != 0 {
			scaledConstReal = scaleUpExact(cReal, ctOut.scale, qi)
			scaledConst = scaledConstReal
		}

		if cImag != 0 {
			scaledConstImag = ring.MRed(scaleUpExact(cImag, ctOut.scale, qi), ringQ.NttPsi[i][1], qi, ringQ.MredParams[i])
			scaledConst = ring.CRed(scaledConst+scaledConstImag, qi)
		}

		p0tmp := ct0.Value[0].Coeffs[i]
		p1tmp := ctOut.Value[0].Coeffs[i]

		ring.AddScalarVec(p0tmp[:ringQ.N>>1], p1tmp[:ringQ.N>>1], scaledConst, qi)

		if cImag != 0 {
			scaledConst = ring.CRed(scaledConstReal+(qi-scaledConstImag), qi)
		}

		ring.AddScalarVec(p0tmp[ringQ.N>>1:], p1tmp[ringQ.N>>1:], scaledConst, qi)
	}
}

// MultByConstAndAdd multiplies ct0 by the input constant, and adds it to the receiver element (it does not modify the input
// element), e.g., ctOut(x) = ctOut(x) + ct0(x) * (a+bi). This functions removes the need of storing the intermediate value c(x) * (a+bi).
// This function will modify the level and the scale of the receiver element depending on the level and the scale of the input
// element and the type of the constant. The level of the receiver element will be set to min(input.level, receiver.level).
// The scale of the receiver element will be set to the scale that the input element would have after the multiplication by the constant.
func (eval *evaluator) MultByConstAndAdd(ct0 *Ciphertext, constant interface{}, ctOut *Ciphertext) {

	var level = utils.MinInt(ct0.Level(), ctOut.Level())

	// Forces a drop of ctOut level to ct0 level
	if ctOut.Level() > level {
		eval.DropLevel(ctOut, ctOut.Level()-level)
	}

	cReal, cImag, scale := eval.getConstAndScale(level, constant)

	var scaledConst, scaledConstReal, scaledConstImag uint64

	ringQ := eval.params.RingQ()

	// If a scaling would be required to multiply by the constant,
	// it equalizes scales such that the scales match in the end.
	if scale != 1 {

		// If ctOut scaling is smaller than ct0's scale + the default scaling,
		// then brings ctOut scale to ct0's scale.
		if ctOut.scale < ct0.scale*scale {

			if scale := math.Floor((scale * ct0.scale) / ctOut.scale); scale > 1 {

				eval.MultByConst(ctOut, scale, ctOut)

			}

			ctOut.scale = scale * ct0.scale

			// If ctOut.scale > ((a+bi)*scale)*ct0(x), then it sets the scale to
			// bring c(x)*scale to the level of ctOut(x) scale
		} else if ctOut.scale > ct0.scale*scale {
			scale = ctOut.scale / ct0.scale
		}

		// If no scaling is required, then it sets the appropriate scale such that
		// ct0(x)*scale matches ctOut(x) scale without modifying ct0(x) scale.
	} else {

		if ctOut.scale > ct0.scale {

			scale = ctOut.scale / ct0.scale

		} else if ct0.scale > ctOut.scale {

			if scale := math.Floor(ct0.scale / ctOut.scale); scale > 1 {
				eval.MultByConst(ctOut, scale, ctOut)
			}

			ctOut.scale = ct0.scale
		}
	}

	// Component-wise multiplication of the following vector to the ciphertext:
	// [a + b*psi_qi^2, ....., a + b*psi_qi^2, a - b*psi_qi^2, ...., a - b*psi_qi^2] mod Qi
	// [{                  N/2                }{                N/2               }]
	// Which is equivalent outside of the NTT domain to adding a to the first coefficient of ct0 and b to the N/2-th coefficient of ct0.
	for i := 0; i < level+1; i++ {

		qi := ringQ.Modulus[i]
		mredParams := ringQ.MredParams[i]
		bredParams := ringQ.BredParams[i]

		scaledConstReal = 0
		scaledConstImag = 0
		scaledConst = 0

		if cReal != 0 {
			scaledConstReal = scaleUpExact(cReal, scale, qi)
			scaledConst = scaledConstReal
		}

		if cImag != 0 {
			scaledConstImag = scaleUpExact(cImag, scale, qi)
			scaledConstImag = ring.MRed(scaledConstImag, ringQ.NttPsi[i][1], qi, mredParams)
			scaledConst = ring.CRed(scaledConst+scaledConstImag, qi)
		}

		scaledConst = ring.MForm(scaledConst, qi, bredParams)

		for u := range ct0.Value {
			p0tmp := ct0.Value[u].Coeffs[i]
			p1tmp := ctOut.Value[u].Coeffs[i]
			ring.MulScalarMontgomeryAndAddVec(p0tmp[:ringQ.N>>1], p1tmp[:ringQ.N>>1], scaledConst, qi, mredParams)
		}

		if cImag != 0 {
			scaledConst = ring.CRed(scaledConstReal+(qi-scaledConstImag), qi)
			scaledConst = ring.MForm(scaledConst, qi, bredParams)
		}

		for u := range ct0.Value {
			p0tmp := ct0.Value[u].Coeffs[i]
			p1tmp := ctOut.Value[u].Coeffs[i]
			ring.MulScalarMontgomeryAndAddVec(p0tmp[ringQ.N>>1:], p1tmp[ringQ.N>>1:], scaledConst, qi, mredParams)
		}
	}
}

// MultByConstNew multiplies ct0 by the input constant and returns the result in a newly created element.
// The scale of the output element will depend on the scale of the input element and the constant (if the constant
// needs to be scaled (its rational part is not zero)). The constant can be a uint64, int64, float64 or complex128.
func (eval *evaluator) MultByConstNew(ct0 *Ciphertext, constant interface{}) (ctOut *Ciphertext) {
	ctOut = NewCiphertext(eval.params, ct0.Degree(), ct0.Level(), ct0.scale)
	eval.MultByConst(ct0, constant, ctOut)
	return
}

// MultByConst multiplies ct0 by the input constant and returns the result in ctOut.
// The scale of the output element will depend on the scale of the input element and the constant (if the constant
// needs to be scaled (its rational part is not zero)). The constant can be a uint64, int64, float64 or complex128.
func (eval *evaluator) MultByConst(ct0 *Ciphertext, constant interface{}, ctOut *Ciphertext) {

	var level = utils.MinInt(ct0.Level(), ctOut.Level())

	cReal, cImag, scale := eval.getConstAndScale(level, constant)

	// Component wise multiplication of the following vector with the ciphertext:
	// [a + b*psi_qi^2, ....., a + b*psi_qi^2, a - b*psi_qi^2, ...., a - b*psi_qi^2] mod Qi
	// [{                  N/2                }{                N/2               }]
	// Which is equivalent outside of the NTT domain to adding a to the first coefficient of ct0 and b to the N/2-th coefficient of ct0.
	ringQ := eval.params.RingQ()
	var scaledConst, scaledConstReal, scaledConstImag uint64
	for i := 0; i < level+1; i++ {

		qi := ringQ.Modulus[i]
		bredParams := ringQ.BredParams[i]
		mredParams := ringQ.MredParams[i]

		scaledConstReal = 0
		scaledConstImag = 0
		scaledConst = 0

		if cReal != 0 {
			scaledConstReal = scaleUpExact(cReal, scale, qi)
			scaledConst = scaledConstReal
		}

		if cImag != 0 {
			scaledConstImag = scaleUpExact(cImag, scale, qi)
			scaledConstImag = ring.MRed(scaledConstImag, ringQ.NttPsi[i][1], qi, mredParams)
			scaledConst = ring.CRed(scaledConst+scaledConstImag, qi)
		}

		scaledConst = ring.MForm(scaledConst, qi, bredParams)

		for u := range ct0.Value {
			p0tmp := ct0.Value[u].Coeffs[i]
			p1tmp := ctOut.Value[u].Coeffs[i]
			ring.MulScalarMontgomeryVec(p0tmp[:ringQ.N>>1], p1tmp[:ringQ.N>>1], scaledConst, qi, mredParams)
		}

		if cImag != 0 {
			scaledConst = ring.CRed(scaledConstReal+(qi-scaledConstImag), qi)
			scaledConst = ring.MForm(scaledConst, qi, bredParams)
		}

		for u := range ct0.Value {
			p0tmp := ct0.Value[u].Coeffs[i]
			p1tmp := ctOut.Value[u].Coeffs[i]
			ring.MulScalarMontgomeryVec(p0tmp[ringQ.N>>1:], p1tmp[ringQ.N>>1:], scaledConst, qi, mredParams)
		}
	}

	ctOut.scale = ct0.scale * scale
}

// MultByGaussianInteger multiples the ct0 by the gaussian integer cReal + i*cImag and returns the result on ctOut.
// Accepted types for cReal and cImag are uint64, int64 and big.Int.
func (eval *evaluator) MultByGaussianInteger(ct0 *Ciphertext, cReal, cImag interface{}, ctOut *Ciphertext) {

	ringQ := eval.params.RingQ()

	level := utils.MinInt(ct0.Level(), ctOut.Level())
	var scaledConst, scaledConstReal, scaledConstImag uint64

	ctOut.scale = ct0.scale

	for i := 0; i < level+1; i++ {

		qi := ringQ.Modulus[i]
		bredParams := ringQ.BredParams[i]
		mredParams := ringQ.MredParams[i]

		scaledConstReal = interfaceMod(cReal, qi)

		if eval.params.RingType() != ring.ConjugateInvariant {
			scaledConstImag = interfaceMod(cImag, qi)
		}

		scaledConst = scaledConstReal

		if scaledConstImag != 0 {
			scaledConstImag = ring.MRed(scaledConstImag, ringQ.NttPsi[i][1], qi, mredParams)
			scaledConst = ring.CRed(scaledConst+scaledConstImag, qi)
		}

		scaledConst = ring.MForm(scaledConst, qi, bredParams)

		for u := range ct0.Value {
			p0tmp := ct0.Value[u].Coeffs[i]
			p1tmp := ctOut.Value[u].Coeffs[i]
			ring.MulScalarMontgomeryVec(p0tmp[:ringQ.N>>1], p1tmp[:ringQ.N>>1], scaledConst, qi, mredParams)
		}

		if cImag != 0 {
			scaledConst = ring.CRed(scaledConstReal+(qi-scaledConstImag), qi)
			scaledConst = ring.MForm(scaledConst, qi, bredParams)
		}

		for u := range ct0.Value {
			p0tmp := ct0.Value[u].Coeffs[i]
			p1tmp := ctOut.Value[u].Coeffs[i]
			ring.MulScalarMontgomeryVec(p0tmp[ringQ.N>>1:], p1tmp[ringQ.N>>1:], scaledConst, qi, mredParams)
		}
	}
}

// MultByGaussianIntegerAndAdd multiples the ct0 by the gaussian integer cReal + i*cImag and adds the result on ctOut.
// Accepted types for cReal and cImag are uint64, int64 and big.Int.
func (eval *evaluator) MultByGaussianIntegerAndAdd(ct0 *Ciphertext, cReal, cImag interface{}, ctOut *Ciphertext) {

	ringQ := eval.params.RingQ()

	level := utils.MinInt(ct0.Level(), ctOut.Level())
	var scaledConst, scaledConstReal, scaledConstImag uint64

	for i := 0; i < level+1; i++ {

		qi := ringQ.Modulus[i]
		bredParams := ringQ.BredParams[i]
		mredParams := ringQ.MredParams[i]

		scaledConstReal = interfaceMod(cReal, qi)

		if eval.params.RingType() != ring.ConjugateInvariant {
			scaledConstImag = interfaceMod(cImag, qi)
		}

		scaledConst = scaledConstReal

		if scaledConstImag != 0 {
			scaledConstImag = ring.MRed(scaledConstImag, ringQ.NttPsi[i][1], qi, mredParams)
			scaledConst = ring.CRed(scaledConst+scaledConstImag, qi)
		}

		scaledConst = ring.MForm(scaledConst, qi, bredParams)

		for u := range ct0.Value {
			p0tmp := ct0.Value[u].Coeffs[i]
			p1tmp := ctOut.Value[u].Coeffs[i]
			ring.MulScalarMontgomeryAndAddVec(p0tmp[:ringQ.N>>1], p1tmp[:ringQ.N>>1], scaledConst, qi, mredParams)
		}

		if cImag != 0 {
			scaledConst = ring.CRed(scaledConstReal+(qi-scaledConstImag), qi)
			scaledConst = ring.MForm(scaledConst, qi, bredParams)
		}

		for u := range ct0.Value {
			p0tmp := ct0.Value[u].Coeffs[i]
			p1tmp := ctOut.Value[u].Coeffs[i]
			ring.MulScalarMontgomeryAndAddVec(p0tmp[ringQ.N>>1:], p1tmp[ringQ.N>>1:], scaledConst, qi, mredParams)
		}
	}
}

// MultByiNew multiplies ct0 by the imaginary number i, and returns the result in a newly created element.
// It does not change the scale.
func (eval *evaluator) MultByiNew(ct0 *Ciphertext) (ctOut *Ciphertext) {

	if eval.params.RingType() == ring.ConjugateInvariant {
		panic("cannot MultByiNew: method not supported when params.RingType() == ring.ConjugateInvariant")
	}

	ctOut = NewCiphertext(eval.params, 1, ct0.Level(), ct0.scale)
	eval.MultByi(ct0, ctOut)
	return ctOut
}

// MultByi multiplies ct0 by the imaginary number i, and returns the result in ctOut.
// It does not change the scale.
func (eval *evaluator) MultByi(ct0 *Ciphertext, ctOut *Ciphertext) {

	if eval.params.RingType() == ring.ConjugateInvariant {
		panic("cannot MultByi: method not supported when params.RingType() == ring.ConjugateInvariant")
	}

	var level = utils.MinInt(ct0.Level(), ctOut.Level())
	ctOut.scale = ct0.scale

	ringQ := eval.params.RingQ()

	var imag uint64

	// Equivalent to a product by the monomial x^(n/2) outside of the NTT domain
	for i := 0; i < level+1; i++ {

		qi := ringQ.Modulus[i]
		mredParams := ringQ.MredParams[i]

		imag = ringQ.NttPsi[i][1] // Psi^2

		for u := range ctOut.Value {
			p0tmp := ct0.Value[u].Coeffs[i]
			p1tmp := ctOut.Value[u].Coeffs[i]
			ring.MulScalarMontgomeryVec(p0tmp[:ringQ.N>>1], p1tmp[:ringQ.N>>1], imag, qi, mredParams)
		}

		imag = qi - imag

		for u := range ctOut.Value {
			p0tmp := ct0.Value[u].Coeffs[i]
			p1tmp := ctOut.Value[u].Coeffs[i]
			ring.MulScalarMontgomeryVec(p0tmp[ringQ.N>>1:], p1tmp[ringQ.N>>1:], imag, qi, mredParams)
		}
	}
}

// DivByiNew multiplies ct0 by the imaginary number 1/i = -i, and returns the result in a newly created element.
// It does not change the scale.
func (eval *evaluator) DivByiNew(ct0 *Ciphertext) (ctOut *Ciphertext) {

	if eval.params.RingType() == ring.ConjugateInvariant {
		panic("cannot DivByiNew: method not supported when params.RingType() == ring.ConjugateInvariant")
	}

	ctOut = NewCiphertext(eval.params, 1, ct0.Level(), ct0.scale)
	eval.DivByi(ct0, ctOut)
	return
}

// DivByi multiplies ct0 by the imaginary number 1/i = -i, and returns the result in ctOut.
// It does not change the scale.
func (eval *evaluator) DivByi(ct0 *Ciphertext, ctOut *Ciphertext) {

	if eval.params.RingType() == ring.ConjugateInvariant {
		panic("cannot DivByi: method not supported when params.RingType() == ring.ConjugateInvariant")
	}

	var level = utils.MinInt(ct0.Level(), ctOut.Level())

	ringQ := eval.params.RingQ()

	ctOut.scale = ct0.scale

	var imag uint64

	// Equivalent to a product by the monomial x^(3*n/2) outside of the NTT domain
	for i := 0; i < level+1; i++ {

		qi := ringQ.Modulus[i]
		mredParams := ringQ.MredParams[i]

		imag = qi - ringQ.NttPsi[i][1] // -Psi^2

		for u := range ctOut.Value {
			p0tmp := ct0.Value[u].Coeffs[i]
			p1tmp := ctOut.Value[u].Coeffs[i]
			ring.MulScalarMontgomeryVec(p0tmp[:ringQ.N>>1], p1tmp[:ringQ.N>>1], imag, qi, mredParams)
		}

		imag = ringQ.NttPsi[i][1] // Psi^2

		for u := range ctOut.Value {
			p0tmp := ct0.Value[u].Coeffs[i]
			p1tmp := ctOut.Value[u].Coeffs[i]
			ring.MulScalarMontgomeryVec(p0tmp[ringQ.N>>1:], p1tmp[ringQ.N>>1:], imag, qi, mredParams)
		}
	}
}

// ScaleUpNew multiplies ct0 by 2^scale and sets its scale to its previous scale
// plus 2^n. It returns the result in a newly created element.
func (eval *evaluator) ScaleUpNew(ct0 *Ciphertext, scale float64) (ctOut *Ciphertext) {
	ctOut = NewCiphertext(eval.params, ct0.Degree(), ct0.Level(), ct0.scale)
	eval.ScaleUp(ct0, scale, ctOut)
	return
}

// ScaleUp multiplies ct0 by 2^scale and sets its scale to its previous scale
// plus 2^n. It returns the result in ctOut.
func (eval *evaluator) ScaleUp(ct0 *Ciphertext, scale float64, ctOut *Ciphertext) {
	eval.MultByConst(ct0, uint64(scale), ctOut)
	ctOut.scale = ct0.scale * scale
}

// SetScale sets the scale of the ciphertext to the input scale (consumes a level)
func (eval *evaluator) SetScale(ct *Ciphertext, scale float64) {
	eval.MultByConst(ct, scale/ct.scale, ct)
	if err := eval.Rescale(ct, scale, ct); err != nil {
		panic(err)
	}
	ct.scale = scale
}

// MulByPow2New multiplies ct0 by 2^pow2 and returns the result in a newly created element.
func (eval *evaluator) MulByPow2New(ct0 *Ciphertext, pow2 int) (ctOut *Ciphertext) {
	ctOut = NewCiphertext(eval.params, ct0.Degree(), ct0.Level(), ct0.scale)
	eval.MulByPow2(ct0, pow2, ctOut)
	return
}

// MulByPow2 multiplies ct0 by 2^pow2 and returns the result in ctOut.
func (eval *evaluator) MulByPow2(ct0 *Ciphertext, pow2 int, ctOut *Ciphertext) {
	var level = utils.MinInt(ct0.Level(), ctOut.Level())
	ctOut.scale = ct0.scale
	for i := range ctOut.Value {
		eval.params.RingQ().MulByPow2Lvl(level, ct0.Value[i], pow2, ctOut.Value[i])
	}
}

// ReduceNew applies a modular reduction to ct0 and returns the result in a newly created element.
// To be used in conjunction with functions that do not apply modular reduction.
func (eval *evaluator) ReduceNew(ct0 *Ciphertext) (ctOut *Ciphertext) {

	ctOut = NewCiphertext(eval.params, ct0.Degree(), ct0.Level(), ct0.scale)

	_ = eval.Reduce(ct0, ctOut)

	return ctOut
}

// Reduce applies a modular reduction to ct0 and returns the result in ctOut.
// To be used in conjunction with functions that do not apply modular reduction.
func (eval *evaluator) Reduce(ct0 *Ciphertext, ctOut *Ciphertext) error {

	if ct0.Degree() != ctOut.Degree() {
		return errors.New("cannot Reduce: degrees of receiver Ciphertext and input Ciphertext do not match")
	}

	for i := range ct0.Value {
		eval.params.RingQ().ReduceLvl(utils.MinInt(ct0.Level(), ctOut.Level()), ct0.Value[i], ctOut.Value[i])
	}

	ctOut.scale = ct0.scale

	return nil
}

// DropLevelNew reduces the level of ct0 by levels and returns the result in a newly created element.
// No rescaling is applied during this procedure.
func (eval *evaluator) DropLevelNew(ct0 *Ciphertext, levels int) (ctOut *Ciphertext) {
	ctOut = ct0.CopyNew()
	eval.DropLevel(ctOut, levels)
	return
}

// DropLevel reduces the level of ct0 by levels and returns the result in ct0.
// No rescaling is applied during this procedure.
func (eval *evaluator) DropLevel(ct0 *Ciphertext, levels int) {
	ct0.Resize(ct0.Degree(), ct0.Level()-levels)
}

// RescaleNew divides ct0 by the last modulus in the moduli chain, and repeats this
// procedure (consuming one level each time) until the scale reaches the original scale or before it goes below it, and returns the result
// in a newly created element. Since all the moduli in the moduli chain are generated to be close to the
// original scale, this procedure is equivalent to dividing the input element by the scale and adding
// some error.
// Returns an error if "threshold <= 0", ct.scale = 0, ct.Level() = 0, ct.IsNTT() != true
func (eval *evaluator) RescaleNew(ct0 *Ciphertext, threshold float64) (ctOut *Ciphertext, err error) {

	ctOut = NewCiphertext(eval.params, ct0.Degree(), ct0.Level(), ct0.scale)

	return ctOut, eval.Rescale(ct0, threshold, ctOut)
}

// Rescale divides ct0 by the last modulus in the moduli chain, and repeats this
// procedure (consuming one level each time) until the scale reaches the original scale or before it goes below it, and returns the result
// in ctOut. Since all the moduli in the moduli chain are generated to be close to the
// original scale, this procedure is equivalent to dividing the input element by the scale and adding
// some error.
// Returns an error if "minScale <= 0", ct.scale = 0, ct.Level() = 0, ct.IsNTT() != true or if ct.Leve() != ctOut.Level()
func (eval *evaluator) Rescale(ctIn *Ciphertext, minScale float64, ctOut *Ciphertext) (err error) {

	ringQ := eval.params.RingQ()

	if minScale <= 0 {
		return errors.New("cannot Rescale: minScale is 0")
	}

	if ctIn.scale == 0 {
		return errors.New("cannot Rescale: ciphertext scale is 0")
	}

	if ctIn.Level() == 0 {
		return errors.New("cannot Rescale: input Ciphertext already at level 0")
	}

	if ctOut.Degree() != ctIn.Degree() {
		return errors.New("cannot Rescale: ctIn.Degree() != ctOut.Degree()")
	}

	ctOut.scale = ctIn.scale

	var nbRescales int
	// Divides the scale by each moduli of the modulus chain as long as the scale isn't smaller than minScale/2
	// or until the output Level() would be zero
	for ctIn.Level()-nbRescales >= 0 && ctOut.scale/float64(ringQ.Modulus[ctIn.Level()-nbRescales]) >= minScale/2 {
		ctOut.scale /= (float64(ringQ.Modulus[ctIn.Level()-nbRescales]))
		nbRescales++
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
func (eval *evaluator) MulNew(ctIn *Ciphertext, op1 Operand) (ctOut *Ciphertext) {
	ctOut = NewCiphertext(eval.params, ctIn.Degree()+op1.Degree(), utils.MinInt(ctIn.Level(), op1.Level()), 0)
	eval.mulRelin(ctIn, op1, false, ctOut)
	return
}

// Mul multiplies ctIn with op1 without relinearization and returns the result in ctOut.
// The procedure will panic if either ctIn or op1 are have a degree higher than 1.
// The procedure will panic if ctOut.Degree != ctIn.Degree + op1.Degree.
func (eval *evaluator) Mul(ctIn *Ciphertext, op1 Operand, ctOut *Ciphertext) {
	eval.mulRelin(ctIn, op1, false, ctOut)
}

// MulRelinNew multiplies ctIn with op1 with relinearization and returns the result in a newly created element.
// The procedure will panic if either ctIn.Degree or op1.Degree > 1.
// The procedure will panic if the evaluator was not created with an relinearization key.
func (eval *evaluator) MulRelinNew(ctIn *Ciphertext, op1 Operand) (ctOut *Ciphertext) {
	ctOut = NewCiphertext(eval.params, 1, utils.MinInt(ctIn.Level(), op1.Level()), 0)
	eval.mulRelin(ctIn, op1, true, ctOut)
	return
}

// MulRelin multiplies ctIn with op1 with relinearization and returns the result in ctOut.
// The procedure will panic if either ctIn.Degree or op1.Degree > 1.
// The procedure will panic if ctOut.Degree != ctIn.Degree + op1.Degree.
// The procedure will panic if the evaluator was not created with an relinearization key.
func (eval *evaluator) MulRelin(ctIn *Ciphertext, op1 Operand, ctOut *Ciphertext) {
	eval.mulRelin(ctIn, op1, true, ctOut)
}

func (eval *evaluator) mulRelin(ctIn *Ciphertext, op1 Operand, relin bool, ctOut *Ciphertext) {

	eval.checkBinary(ctIn, op1, ctOut, utils.MaxInt(ctIn.Degree(), op1.Degree()))

	level := utils.MinInt(utils.MinInt(ctIn.Level(), op1.Level()), ctOut.Level())

	if ctIn.Degree()+op1.Degree() > 2 {
		panic("cannot MulRelin: the sum of the input elements' total degree cannot be larger than 2")
	}

	ctOut.scale = ctIn.Scale() * op1.Scale()

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
			c2.IsNTT = true
			eval.GadgetProduct(level, c2, eval.Rlk.Keys[0].GadgetCiphertext, eval.BuffQP[1].Q, eval.BuffQP[2].Q)
			ringQ.AddLvl(level, c0, eval.BuffQP[1].Q, ctOut.Value[0])
			ringQ.AddLvl(level, c1, eval.BuffQP[2].Q, ctOut.Value[1])
		}

		// Case Plaintext (x) Ciphertext or Ciphertext (x) Plaintext
	} else {

		if ctOut.Degree() < ctIn.Degree() {
			ctOut.El().Resize(ctIn.Degree(), level)
		} else {
			ctOut.El().Resize(ctOut.Degree(), level)
		}

		c00 := eval.buffQ[0]

		ringQ.MFormLvl(level, op1.El().Value[0], c00)
		for i := range ctIn.Value {
			ringQ.MulCoeffsMontgomeryLvl(level, ctIn.Value[i], c00, ctOut.Value[i])
		}
	}
}

// MulAndAdd multiplies ctIn with op1 without relinearization and adds the result on ctOut.
// User must ensure that ctOut.scale <= ctIn.scale * op1.scale.
// If ctOut.scale < ctIn.scale * op1.scale, then scales up ctOut before adding the result.
// The procedure will panic if either ctIn or op1 are have a degree higher than 1.
// The procedure will panic if ctOut.Degree != ctIn.Degree + op1.Degree.
// The procedure will panic if ctOut = ctIn or op1.
func (eval *evaluator) MulAndAdd(ctIn *Ciphertext, op1 Operand, ctOut *Ciphertext) {
	eval.mulRelinAndAdd(ctIn, op1, false, ctOut)
}

// MulRelinAndAdd multiplies ctIn with op1 with relinearization and adds the result on ctOut.
// User must ensure that ctOut.scale <= ctIn.scale * op1.scale.
// If ctOut.scale < ctIn.scale * op1.scale, then scales up ctOut before adding the result.
// The procedure will panic if either ctIn.Degree or op1.Degree > 1.
// The procedure will panic if ctOut.Degree != ctIn.Degree + op1.Degree.
// The procedure will panic if the evaluator was not created with an relinearization key.
// The procedure will panic if ctOut = ctIn or op1.
func (eval *evaluator) MulRelinAndAdd(ctIn *Ciphertext, op1 Operand, ctOut *Ciphertext) {
	eval.mulRelinAndAdd(ctIn, op1, true, ctOut)
}

func (eval *evaluator) mulRelinAndAdd(ctIn *Ciphertext, op1 Operand, relin bool, ctOut *Ciphertext) {

	eval.checkBinary(ctIn, op1, ctOut, utils.MaxInt(ctIn.Degree(), op1.Degree()))

	level := utils.MinInt(utils.MinInt(ctIn.Level(), op1.Level()), ctOut.Level())

	if ctIn.Degree()+op1.Degree() > 2 {
		panic("cannot MulRelinAndAdd: the sum of the input elements' degree cannot be larger than 2")
	}

	if ctIn.El() == ctOut.El() || op1.El() == ctOut.El() {
		panic("cannot MulRelinAndAdd: ctOut must be different from op0 and op1")
	}

	resScale := ctIn.scale * op1.Scale()

	if ctOut.scale < resScale {
		eval.MultByConst(ctOut, math.Round(resScale/ctOut.scale), ctOut)
		ctOut.scale = resScale
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
			c2.IsNTT = true
			ringQ.MulCoeffsMontgomeryLvl(level, c01, tmp1.Value[1], c2) // c2 += c[1]*c[1]
			eval.GadgetProduct(level, c2, eval.Rlk.Keys[0].GadgetCiphertext, eval.BuffQP[1].Q, eval.BuffQP[2].Q)
			ringQ.AddLvl(level, c0, eval.BuffQP[1].Q, c0)
			ringQ.AddLvl(level, c1, eval.BuffQP[2].Q, c1)
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
func (eval *evaluator) RelinearizeNew(ct0 *Ciphertext) (ctOut *Ciphertext) {
	ctOut = NewCiphertext(eval.params, 1, ct0.Level(), ct0.scale)
	eval.Relinearize(ct0, ctOut)
	return
}

// Relinearize applies the relinearization procedure on ct0 and returns the result in ctOut. The input Ciphertext must be of degree two.
func (eval *evaluator) Relinearize(ct0 *Ciphertext, ctOut *Ciphertext) {
	eval.Evaluator.Relinearize(ct0.Ciphertext, ctOut.Ciphertext)
	ctOut.scale = ct0.scale
}

// SwitchKeysNew re-encrypts ct0 under a different key and returns the result in a newly created element.
// It requires a SwitchingKey, which is computed from the key under which the Ciphertext is currently encrypted,
// and the key under which the Ciphertext will be re-encrypted.
func (eval *evaluator) SwitchKeysNew(ct0 *Ciphertext, switchingKey *rlwe.SwitchingKey) (ctOut *Ciphertext) {
	ctOut = NewCiphertext(eval.params, ct0.Degree(), ct0.Level(), ct0.scale)
	eval.SwitchKeys(ct0, switchingKey, ctOut)
	return
}

// SwitchKeys re-encrypts ctIn under a different key and returns the result in ctOut.
// It requires a SwitchingKey, which is computed from the key under which the Ciphertext is currently encrypted,
// and the key under which the Ciphertext will be re-encrypted.
func (eval *evaluator) SwitchKeys(ctIn *Ciphertext, switchingKey *rlwe.SwitchingKey, ctOut *Ciphertext) {
	eval.Evaluator.SwitchKeys(ctIn.Ciphertext, switchingKey, ctOut.Ciphertext)
	ctOut.scale = ctIn.scale
}

// RotateNew rotates the columns of ct0 by k positions to the left, and returns the result in a newly created element.
// If the provided element is a Ciphertext, a key-switching operation is necessary and a rotation key for the specific rotation needs to be provided.
func (eval *evaluator) RotateNew(ct0 *Ciphertext, k int) (ctOut *Ciphertext) {
	ctOut = NewCiphertext(eval.params, ct0.Degree(), ct0.Level(), ct0.scale)
	eval.Rotate(ct0, k, ctOut)
	return
}

// Rotate rotates the columns of ct0 by k positions to the left and returns the result in ctOut.
// If the provided element is a Ciphertext, a key-switching operation is necessary and a rotation key for the specific rotation needs to be provided.
func (eval *evaluator) Rotate(ct0 *Ciphertext, k int, ctOut *Ciphertext) {
	eval.Automorphism(ct0.Ciphertext, eval.params.GaloisElementForColumnRotationBy(k), ctOut.Ciphertext)
	ctOut.scale = ct0.scale
}

// ConjugateNew conjugates ct0 (which is equivalent to a row rotation) and returns the result in a newly
// created element. If the provided element is a Ciphertext, a key-switching operation is necessary and a rotation key
// for the row rotation needs to be provided.
func (eval *evaluator) ConjugateNew(ct0 *Ciphertext) (ctOut *Ciphertext) {

	if eval.params.RingType() == ring.ConjugateInvariant {
		panic("cannot ConjugateNew: method is not supported when params.RingType() == ring.ConjugateInvariant")
	}

	ctOut = NewCiphertext(eval.params, ct0.Degree(), ct0.Level(), ct0.scale)
	eval.Conjugate(ct0, ctOut)
	return
}

// Conjugate conjugates ct0 (which is equivalent to a row rotation) and returns the result in ctOut.
// If the provided element is a Ciphertext, a key-switching operation is necessary and a rotation key for the row rotation needs to be provided.
func (eval *evaluator) Conjugate(ct0 *Ciphertext, ctOut *Ciphertext) {

	if eval.params.RingType() == ring.ConjugateInvariant {
		panic("cannot Conjugate: method is not supported when params.RingType() == ring.ConjugateInvariant")
	}

	eval.Automorphism(ct0.Ciphertext, eval.params.GaloisElementForRowRotation(), ctOut.Ciphertext)
	ctOut.scale = ct0.scale
}

func (eval *evaluator) RotateHoistedNoModDownNew(level int, rotations []int, c0 *ring.Poly, c2DecompQP []ringqp.Poly) (cOut map[int][2]ringqp.Poly) {
	ringQ := eval.params.RingQ()
	ringP := eval.params.RingP()
	cOut = make(map[int][2]ringqp.Poly)
	for _, i := range rotations {
		if i != 0 {
			cOut[i] = [2]ringqp.Poly{{Q: ringQ.NewPolyLvl(level), P: ringP.NewPoly()}, {Q: ringQ.NewPolyLvl(level), P: ringP.NewPoly()}}
			eval.AutomorphismHoistedNoModDown(level, c0, c2DecompQP, eval.params.GaloisElementForColumnRotationBy(i), cOut[i][0].Q, cOut[i][1].Q, cOut[i][0].P, cOut[i][1].P)
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
