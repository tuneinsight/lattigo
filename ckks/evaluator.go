package ckks

import (
	"errors"
	"fmt"
	"math/big"

	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/rlwe"
	"github.com/tuneinsight/lattigo/v4/rlwe/ringqp"
	"github.com/tuneinsight/lattigo/v4/utils"
	"github.com/tuneinsight/lattigo/v4/utils/bignum"
)

// Evaluator is an interface implementing the methods to conduct homomorphic operations between ciphertext and/or plaintexts.
type Evaluator interface {
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
	EvaluatePolyVector(input interface{}, pols []*bignum.Polynomial, encoder *Encoder, slotIndex map[int][]int, targetScale rlwe.Scale) (ctOut *rlwe.Ciphertext, err error)

	// GoldschmidtDivision
	GoldschmidtDivisionNew(ct *rlwe.Ciphertext, minValue, log2Targetprecision float64, btp rlwe.Bootstrapper) (ctInv *rlwe.Ciphertext, err error)

	// Linear Transformations
	LinearTransformNew(op0 *rlwe.Ciphertext, linearTransform interface{}) (ctOut []*rlwe.Ciphertext)
	LinearTransform(op0 *rlwe.Ciphertext, linearTransform interface{}, ctOut []*rlwe.Ciphertext)
	MultiplyByDiagMatrix(op0 *rlwe.Ciphertext, matrix LinearTransform, c2DecompQP []ringqp.Poly, ctOut *rlwe.Ciphertext)
	MultiplyByDiagMatrixBSGS(op0 *rlwe.Ciphertext, matrix LinearTransform, c2DecompQP []ringqp.Poly, ctOut *rlwe.Ciphertext)

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

	// ==============
	// === Others ===
	// ==============
	CheckBinary(op0, op1, opOut rlwe.Operand, opOutMinDegree int) (degree, level int)
	CheckUnary(op0, opOut rlwe.Operand) (degree, level int)
	GetRLWEEvaluator() *rlwe.Evaluator
	BuffQ() [3]*ring.Poly
	BuffCt() *rlwe.Ciphertext
	ShallowCopy() Evaluator
	WithKey(evk rlwe.EvaluationKeySetInterface) Evaluator
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
func NewEvaluator(params Parameters, evk rlwe.EvaluationKeySetInterface) Evaluator {
	eval := new(evaluator)
	eval.evaluatorBase = newEvaluatorBase(params)
	eval.evaluatorBuffers = newEvaluatorBuffers(eval.evaluatorBase)
	eval.Evaluator = rlwe.NewEvaluator(params.Parameters, evk)

	return eval
}

// GetRLWEEvaluator returns the underlying *rlwe.Evaluator.
func (eval *evaluator) GetRLWEEvaluator() *rlwe.Evaluator {
	return eval.Evaluator
}

// Add adds op1 to op0 and returns the result in op2.
func (eval *evaluator) Add(op0 *rlwe.Ciphertext, op1 interface{}, op2 *rlwe.Ciphertext) {

	switch op1 := op1.(type) {
	case rlwe.Operand:
		_, level := eval.CheckBinary(op0, op1, op2, utils.MaxInt(op0.Degree(), op1.Degree()))
		eval.evaluateInPlace(level, op0, op1, op2, eval.params.RingQ().AtLevel(level).Add)
	default:
		level := utils.MinInt(op0.Level(), op2.Level())
		RNSReal, RNSImag := bigComplexToRNSScalar(eval.params.RingQ().AtLevel(level), &op0.Scale.Value, bignum.ToComplex(op1, eval.params.DefaultPrecision()))
		op2.Resize(op0.Degree(), level)
		eval.evaluateWithScalar(level, op0.Value[:1], RNSReal, RNSImag, op2.Value[:1], eval.params.RingQ().AtLevel(level).AddDoubleRNSScalar)
	}
}

// AddNew adds op1 to op0 and returns the result in a newly created element op2.
func (eval *evaluator) AddNew(op0 *rlwe.Ciphertext, op1 interface{}) (op2 *rlwe.Ciphertext) {
	op2 = op0.CopyNew()
	eval.Add(op2, op1, op2)
	return
}

// Sub subtracts op1 from op0 and returns the result in op2.
func (eval *evaluator) Sub(op0 *rlwe.Ciphertext, op1 interface{}, op2 *rlwe.Ciphertext) {

	switch op1 := op1.(type) {
	case rlwe.Operand:
		_, level := eval.CheckBinary(op0, op1, op2, utils.MaxInt(op0.Degree(), op1.Degree()))

		eval.evaluateInPlace(level, op0, op1, op2, eval.params.RingQ().AtLevel(level).Sub)

		if op0.Degree() < op1.Degree() {
			for i := op0.Degree() + 1; i < op1.Degree()+1; i++ {
				eval.params.RingQ().AtLevel(level).Neg(op2.Value[i], op2.Value[i])
			}
		}
	default:
		level := utils.MinInt(op0.Level(), op2.Level())
		RNSReal, RNSImag := bigComplexToRNSScalar(eval.params.RingQ().AtLevel(level), &op0.Scale.Value, bignum.ToComplex(op1, eval.params.DefaultPrecision()))
		op2.Resize(op0.Degree(), level)
		eval.evaluateWithScalar(level, op0.Value[:1], RNSReal, RNSImag, op2.Value[:1], eval.params.RingQ().AtLevel(level).SubDoubleRNSScalar)
	}
}

// SubNew subtracts op1 from op0 and returns the result in a newly created element op2.
func (eval *evaluator) SubNew(op0 *rlwe.Ciphertext, op1 interface{}) (op2 *rlwe.Ciphertext) {
	op2 = op0.CopyNew()
	eval.Sub(op2, op1, op2)
	return
}

func (eval *evaluator) evaluateInPlace(level int, c0 *rlwe.Ciphertext, c1 rlwe.Operand, ctOut *rlwe.Ciphertext, evaluate func(*ring.Poly, *ring.Poly, *ring.Poly)) {

	var tmp0, tmp1 *rlwe.Ciphertext

	maxDegree := utils.Max(c0.Degree(), c1.Degree())
	minDegree := utils.Min(c0.Degree(), c1.Degree())

	// Else resizes the receiver element
	ctOut.El().Resize(maxDegree, ctOut.Level())

	c0Scale := c0.GetMetaData().Scale
	c1Scale := c1.GetMetaData().Scale

	if ctOut.Level() > level {
		eval.DropLevel(ctOut, ctOut.Level()-utils.Min(c0.Level(), c1.Level()))
	}

	cmp := c0.GetMetaData().Scale.Cmp(c1.GetMetaData().Scale)

	// Checks whether or not the receiver element is the same as one of the input elements
	// and acts accordingly to avoid unnecessary element creation or element overwriting,
	// and scales properly the element before the evaluation.
	if ctOut == c0 {

		if cmp == 1 {

			ratioFlo := c0Scale.Div(c1Scale).Value

			ratioInt, _ := ratioFlo.Int(nil)

			if ratioInt.Cmp(new(big.Int).SetUint64(0)) == 1 {

				tmp1 = rlwe.NewCiphertextAtLevelFromPoly(level, eval.buffCt.Value[:c1.Degree()+1])
				tmp1.MetaData = ctOut.MetaData

				eval.Mul(c1.El(), ratioInt, tmp1)
			}

		} else if cmp == -1 {

			ratioFlo := c1Scale.Div(c0Scale).Value

			ratioInt, _ := ratioFlo.Int(nil)

			if ratioInt.Cmp(new(big.Int).SetUint64(0)) == 1 {

				eval.Mul(c0, ratioInt, c0)

				ctOut.Scale = c1.GetMetaData().Scale

				tmp1 = c1.El()
			}

		} else {
			tmp1 = &rlwe.Ciphertext{OperandQ: *c1.El()}
		}

		tmp0 = &rlwe.Ciphertext{OperandQ: *c0.El()}

	} else if ctOut == c1 {

		if cmp == 1 {

			ratioFlo := c0Scale.Div(c1Scale).Value

			ratioInt, _ := ratioFlo.Int(nil)

			if ratioInt.Cmp(new(big.Int).SetUint64(0)) == 1 {
				eval.Mul(c1.El(), ratioInt, ctOut)

				ctOut.Scale = c0.Scale

				tmp0 = c0.El()
			}

		} else if cmp == -1 {

			ratioFlo := c1Scale.Div(c0Scale).Value

			ratioInt, _ := ratioFlo.Int(nil)

			if ratioInt.Cmp(new(big.Int).SetUint64(0)) == 1 {
				// Will avoid resizing on the output
				tmp0 = rlwe.NewCiphertextAtLevelFromPoly(level, eval.buffCt.Value[:c0.Degree()+1])
				tmp0.MetaData = ctOut.MetaData

				eval.Mul(c0, ratioInt, tmp0)
			}

		} else {
			tmp0 = &rlwe.Ciphertext{OperandQ: *c0.El()}
		}

		tmp1 = &rlwe.Ciphertext{OperandQ: *c1.El()}

	} else {

		if cmp == 1 {

			ratioFlo := c0Scale.Div(c1Scale).Value

			ratioInt, _ := ratioFlo.Int(nil)

			if ratioInt.Cmp(new(big.Int).SetUint64(0)) == 1 {
				// Will avoid resizing on the output
				tmp1 = rlwe.NewCiphertextAtLevelFromPoly(level, eval.buffCt.Value[:c1.Degree()+1])
				tmp1.MetaData = ctOut.MetaData

				eval.Mul(c1.El(), ratioInt, tmp1)

				tmp0 = c0.El()
			}

		} else if cmp == -1 {

			ratioFlo := c1Scale.Div(c0Scale).Value

			ratioInt, _ := ratioFlo.Int(nil)

			if ratioInt.Cmp(new(big.Int).SetUint64(0)) == 1 {

				tmp0 = rlwe.NewCiphertextAtLevelFromPoly(level, eval.buffCt.Value[:c0.Degree()+1])
				tmp0.MetaData = ctOut.MetaData

				eval.Mul(c0, ratioInt, tmp0)

				tmp1 = c1.El()

			}

		} else {
			tmp0 = &rlwe.Ciphertext{OperandQ: *c0.El()}
			tmp1 = &rlwe.Ciphertext{OperandQ: *c1.El()}
		}
	}

	for i := 0; i < minDegree+1; i++ {
		evaluate(tmp0.Value[i], tmp1.Value[i], ctOut.El().Value[i])
	}

	scale := c0.Scale.Max(c1.GetMetaData().Scale)

	ctOut.MetaData = c0.MetaData
	ctOut.Scale = scale

	// If the inputs degrees differ, it copies the remaining degree on the receiver.
	// Also checks that the receiver is not one of the inputs to avoid unnecessary work.

	if c0.Degree() > c1.Degree() && &tmp0.OperandQ != ctOut.El() {
		for i := minDegree + 1; i < maxDegree+1; i++ {
			ring.Copy(tmp0.Value[i], ctOut.El().Value[i])
		}
	} else if c1.Degree() > c0.Degree() && &tmp1.OperandQ != ctOut.El() {
		for i := minDegree + 1; i < maxDegree+1; i++ {
			ring.Copy(tmp1.Value[i], ctOut.El().Value[i])
		}
	}
}

func (eval *evaluator) evaluateWithScalar(level int, p0 []*ring.Poly, RNSReal, RNSImag ring.RNSScalar, p1 []*ring.Poly, evaluate func(*ring.Poly, ring.RNSScalar, ring.RNSScalar, *ring.Poly)) {

	// Component wise operation with the following vector:
	// [a + b*psi_qi^2, ....., a + b*psi_qi^2, a - b*psi_qi^2, ...., a - b*psi_qi^2] mod Qi
	// [{                  N/2                }{                N/2               }]
	// Which is equivalent outside of the NTT domain to evaluating a to the first coefficient of ct0 and b to the N/2-th coefficient of ct0.
	for i, s := range eval.params.RingQ().SubRings[:level+1] {
		RNSImag[i] = ring.MRed(RNSImag[i], s.RootsForward[1], s.Modulus, s.MRedConstant)
		RNSReal[i], RNSImag[i] = ring.CRed(RNSReal[i]+RNSImag[i], s.Modulus), ring.CRed(RNSReal[i]+s.Modulus-RNSImag[i], s.Modulus)
	}

	for i := range p0 {
		evaluate(p0[i], RNSReal, RNSImag, p1[i])
	}
}

// ScaleUpNew multiplies ct0 by scale and sets its scale to its previous scale times scale returns the result in ctOut.
func (eval *evaluator) ScaleUpNew(ct0 *rlwe.Ciphertext, scale rlwe.Scale) (ctOut *rlwe.Ciphertext) {
	ctOut = NewCiphertext(eval.params, ct0.Degree(), ct0.Level())
	eval.ScaleUp(ct0, scale, ctOut)
	return
}

// ScaleUp multiplies ct0 by scale and sets its scale to its previous scale times scale returns the result in ctOut.
func (eval *evaluator) ScaleUp(ct0 *rlwe.Ciphertext, scale rlwe.Scale, ctOut *rlwe.Ciphertext) {
	eval.Mul(ct0, scale.Uint64(), ctOut)
	ctOut.MetaData = ct0.MetaData
	ctOut.Scale = ct0.Scale.Mul(scale)
}

// SetScale sets the scale of the ciphertext to the input scale (consumes a level).
func (eval *evaluator) SetScale(ct *rlwe.Ciphertext, scale rlwe.Scale) {
	ratioFlo := scale.Div(ct.Scale).Value
	eval.Mul(ct, &ratioFlo, ct)
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
func (eval *evaluator) Rescale(op0 *rlwe.Ciphertext, minScale rlwe.Scale, ctOut *rlwe.Ciphertext) (err error) {

	if minScale.Cmp(rlwe.NewScale(0)) != 1 {
		return errors.New("cannot Rescale: minScale is <0")
	}

	minScale = minScale.Div(rlwe.NewScale(2))

	if op0.Scale.Cmp(rlwe.NewScale(0)) != 1 {
		return errors.New("cannot Rescale: ciphertext scale is <0")
	}

	if op0.Level() == 0 {
		return errors.New("cannot Rescale: input Ciphertext already at level 0")
	}

	if ctOut.Degree() != op0.Degree() {
		return errors.New("cannot Rescale: op0.Degree() != ctOut.Degree()")
	}

	ctOut.MetaData = op0.MetaData

	newLevel := op0.Level()

	ringQ := eval.params.RingQ().AtLevel(op0.Level())

	// Divides the scale by each moduli of the modulus chain as long as the scale isn't smaller than minScale/2
	// or until the output Level() would be zero
	var nbRescales int
	for newLevel >= 0 {

		scale := ctOut.Scale.Div(rlwe.NewScale(ringQ.SubRings[newLevel].Modulus))

		if scale.Cmp(minScale) == -1 {
			break
		}

		ctOut.Scale = scale

		nbRescales++
		newLevel--
	}

	if nbRescales > 0 {
		for i := range ctOut.Value {
			ringQ.DivRoundByLastModulusManyNTT(nbRescales, op0.Value[i], eval.buffQ[0], ctOut.Value[i])
		}
		ctOut.Resize(ctOut.Degree(), newLevel)
	} else {
		if op0 != ctOut {
			ctOut.Copy(op0)
		}
	}

	return nil
}

// MulNew multiplies op0 with op1 without relinearization and returns the result in a newly created element op2.
//
// op1.(type) can be rlwe.Operand, complex128, float64, int, int64, uint64. *big.Float, *big.Int or *ring.Complex.
//
// If op1.(type) == rlwe.Operand:
// - The procedure will panic if either op0.Degree or op1.Degree > 1.
func (eval *evaluator) MulNew(op0 *rlwe.Ciphertext, op1 interface{}) (op2 *rlwe.Ciphertext) {
	op2 = op0.CopyNew()
	eval.Mul(op2, op1, op2)
	return
}

// Mul multiplies op0 with op1 without relinearization and returns the result in ctOut.
//
// op1.(type) can be rlwe.Operand, complex128, float64, int, int64, uint64. *big.Float, *big.Int or *ring.Complex.
//
// If op1.(type) == rlwe.Operand:
// - The procedure will panic if either op0 or op1 are have a degree higher than 1.
// - The procedure will panic if op2.Degree != op0.Degree + op1.Degree.
func (eval *evaluator) Mul(op0 *rlwe.Ciphertext, op1 interface{}, op2 *rlwe.Ciphertext) {
	switch op1 := op1.(type) {
	case rlwe.Operand:
		eval.mulRelin(op0, op1, false, op2)
	default:
		level := utils.MinInt(op0.Level(), op2.Level())
		op2.Resize(op0.Degree(), level)

		ringQ := eval.params.RingQ().AtLevel(level)

		cmplxBig := bignum.ToComplex(op1, eval.params.DefaultPrecision())

		var scale rlwe.Scale

		if cmplxBig.IsInt() {
			scale = rlwe.NewScale(1)
		} else {
			scale = rlwe.NewScale(ringQ.SubRings[level].Modulus)

			for i := 1; i < eval.params.DefaultScaleModuliRatio(); i++ {
				scale = scale.Mul(rlwe.NewScale(ringQ.SubRings[level-i].Modulus))
			}
		}

		RNSReal, RNSImag := bigComplexToRNSScalar(ringQ, &scale.Value, cmplxBig)

		eval.evaluateWithScalar(level, op0.Value, RNSReal, RNSImag, op2.Value, ringQ.MulDoubleRNSScalar)
		op2.MetaData = op0.MetaData
		op2.Scale = op0.Scale.Mul(scale)
	}
}

// MulRelinNew multiplies op0 with op1 with relinearization and returns the result in a newly created element.
// The procedure will panic if either op0.Degree or op1.Degree > 1.
// The procedure will panic if the evaluator was not created with an relinearization key.
func (eval *evaluator) MulRelinNew(op0 *rlwe.Ciphertext, op1 rlwe.Operand) (ctOut *rlwe.Ciphertext) {
	ctOut = NewCiphertext(eval.params, 1, utils.MinInt(op0.Level(), op1.Level()))
	eval.mulRelin(op0, op1, true, ctOut)
	return
}

// MulRelin multiplies op0 with op1 with relinearization and returns the result in ctOut.
// The procedure will panic if either op0.Degree or op1.Degree > 1.
// The procedure will panic if ctOut.Degree != op0.Degree + op1.Degree.
// The procedure will panic if the evaluator was not created with an relinearization key.
func (eval *evaluator) MulRelin(op0 *rlwe.Ciphertext, op1 rlwe.Operand, ctOut *rlwe.Ciphertext) {
	eval.mulRelin(op0, op1, true, ctOut)
}

func (eval *evaluator) mulRelin(op0 *rlwe.Ciphertext, op1 rlwe.Operand, relin bool, ctOut *rlwe.Ciphertext) {

	if op0.Degree()+op1.Degree() > 2 {
		panic("cannot MulRelin: the sum of the input elements' total degree cannot be larger than 2")
	}

	ctOut.MetaData = op0.MetaData
	ctOut.Scale = op0.Scale.Mul(op1.GetMetaData().Scale)

	var c00, c01, c0, c1, c2 *ring.Poly

	// Case Ciphertext (x) Ciphertext
	if op0.Degree() == 1 && op1.Degree() == 1 {

		_, level := eval.CheckBinary(op0, op1, ctOut, ctOut.Degree())

		ringQ := eval.params.RingQ().AtLevel(level)

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
		var tmp0, tmp1 *rlwe.OperandQ
		if op1.El() == ctOut.El() {
			tmp0, tmp1 = op1.El(), op0.El()
		} else {
			tmp0, tmp1 = op0.El(), op1.El()
		}

		ringQ.MForm(tmp0.Value[0], c00)
		ringQ.MForm(tmp0.Value[1], c01)

		if op0.El() == op1.El() { // squaring case
			ringQ.MulCoeffsMontgomery(c00, tmp1.Value[0], c0) // c0 = c[0]*c[0]
			ringQ.MulCoeffsMontgomery(c01, tmp1.Value[1], c2) // c2 = c[1]*c[1]
			ringQ.MulCoeffsMontgomery(c00, tmp1.Value[1], c1) // c1 = 2*c[0]*c[1]
			ringQ.Add(c1, c1, c1)

		} else { // regular case
			ringQ.MulCoeffsMontgomery(c00, tmp1.Value[0], c0) // c0 = c0[0]*c0[0]
			ringQ.MulCoeffsMontgomery(c01, tmp1.Value[1], c2) // c2 = c0[1]*c1[1]
			ringQ.MulCoeffsMontgomery(c00, tmp1.Value[1], c1)
			ringQ.MulCoeffsMontgomeryThenAdd(c01, tmp1.Value[0], c1) // c1 = c0[0]*c1[1] + c0[1]*c1[0]
		}

		if relin {

			var rlk *rlwe.RelinearizationKey
			var err error
			if rlk, err = eval.CheckAndGetRelinearizationKey(); err != nil {
				panic(fmt.Errorf("cannot relinearize: %w", err))
			}

			tmpCt := &rlwe.Ciphertext{}
			tmpCt.Value = []*ring.Poly{eval.BuffQP[1].Q, eval.BuffQP[2].Q}
			tmpCt.IsNTT = true

			eval.GadgetProduct(level, c2, &rlk.GadgetCiphertext, tmpCt)
			ringQ.Add(c0, tmpCt.Value[0], ctOut.Value[0])
			ringQ.Add(c1, tmpCt.Value[1], ctOut.Value[1])
		}

		// Case Plaintext (x) Ciphertext or Ciphertext (x) Plaintext
	} else {

		_, level := eval.CheckBinary(op0, op1, ctOut, ctOut.Degree())

		ringQ := eval.params.RingQ().AtLevel(level)

		var c0 *ring.Poly
		var c1 []*ring.Poly
		if op0.Degree() == 0 {
			c0 = eval.buffQ[0]
			ringQ.MForm(op0.Value[0], c0)
			c1 = op1.El().Value

		} else {
			c0 = eval.buffQ[0]
			ringQ.MForm(op1.El().Value[0], c0)
			c1 = op0.Value
		}

		ctOut.El().Resize(op0.Degree()+op1.Degree(), level)

		for i := range c1 {
			ringQ.MulCoeffsMontgomery(c0, c1[i], ctOut.Value[i])
		}
	}
}

// MulThenAdd evaluate op2 = op2 + op0 * op1.
//
// op1.(type) can be rlwe.Operand, complex128, float64, int, int64, uint64. *big.Float, *big.Int or *ring.Complex.
//
// If op1.(type) is complex128, float64, int, int64, uint64. *big.Float, *big.Int or *ring.Complex:
//
// This function will not modify op0 but will multiply op2 by Q[min(op0.Level(), op2.Level())] if:
// - op0.Scale == op2.Scale
// - constant is not a Gaussian integer.
//
// If op0.Scale == op2.Scale, and constant is not a Gaussian integer, then the constant will be scaled by
// Q[min(op0.Level(), op2.Level())] else if op2.Scale > op0.Scale, the constant will be scaled by op2.Scale/op0.Scale.
//
// To correctly use this function, make sure that either op0.Scale == op2.Scale or
// op2.Scale = op0.Scale * Q[min(op0.Level(), op2.Level())].
//
// If op1.(type) is rlwe.Operand, the multiplication is carried outwithout relinearization and:
//
// This function will panic if op0.Scale > op2.Scale.
// User must ensure that op2.scale <= op0.scale * op1.scale.
// If op2.scale < op0.scale * op1.scale, then scales up op2 before adding the result.
// Additionally, the procedure will panic if:
// - either op0 or op1 are have a degree higher than 1.
// - op2.Degree != op0.Degree + op1.Degree.
// - op2 = op0 or op1.
func (eval *evaluator) MulThenAdd(op0 *rlwe.Ciphertext, op1 interface{}, op2 *rlwe.Ciphertext) {
	switch op1 := op1.(type) {
	case rlwe.Operand:
		eval.mulRelinThenAdd(op0, op1, false, op2)
	default:
		var level = utils.MinInt(op0.Level(), op2.Level())

		ringQ := eval.params.RingQ().AtLevel(level)

		op2.Resize(op2.Degree(), level)

		cmplxBig := bignum.ToComplex(op1, eval.params.DefaultPrecision())

		var scaleRLWE rlwe.Scale

		// If op0 and op2 scales are identical, but the op1 is not a Gaussian integer then multiplies op2 by scaleRLWE.
		// This ensures noiseless addition with op2 = scaleRLWE * op2 + op0 * round(scalar * scaleRLWE).
		if cmp := op0.Scale.Cmp(op2.Scale); cmp == 0 {

			if cmplxBig.IsInt() {
				scaleRLWE = rlwe.NewScale(1)
			} else {
				scaleRLWE = rlwe.NewScale(ringQ.SubRings[level].Modulus)

				for i := 1; i < eval.params.DefaultScaleModuliRatio(); i++ {
					scaleRLWE = scaleRLWE.Mul(rlwe.NewScale(ringQ.SubRings[level-i].Modulus))
				}

				scaleInt := new(big.Int)
				scaleRLWE.Value.Int(scaleInt)
				eval.Mul(op2, scaleInt, op2)
				op2.Scale = op2.Scale.Mul(scaleRLWE)
			}

		} else if cmp == -1 { // op2.Scale > op0.Scale then the scaling factor for op1 becomes the quotient between the two scales
			scaleRLWE = op2.Scale.Div(op0.Scale)
		} else {
			panic("MulThenAdd: op0.Scale > op2.Scale is not supported")
		}

		RNSReal, RNSImag := bigComplexToRNSScalar(ringQ, &scaleRLWE.Value, cmplxBig)

		eval.evaluateWithScalar(level, op0.Value, RNSReal, RNSImag, op2.Value, ringQ.MulDoubleRNSScalarThenAdd)
	}
}

// MulRelinThenAdd multiplies op0 with op1 with relinearization and adds the result on op2.
// User must ensure that op2.scale <= op0.scale * op1.scale.
// If op2.scale < op0.scale * op1.scale, then scales up op2 before adding the result.
// The procedure will panic if either op0.Degree or op1.Degree > 1.
// The procedure will panic if op2.Degree != op0.Degree + op1.Degree.
// The procedure will panic if the evaluator was not created with an relinearization key.
// The procedure will panic if op2 = op0 or op1.
func (eval *evaluator) MulRelinThenAdd(op0 *rlwe.Ciphertext, op1 rlwe.Operand, op2 *rlwe.Ciphertext) {
	eval.mulRelinThenAdd(op0, op1, true, op2)
}

func (eval *evaluator) mulRelinThenAdd(op0 *rlwe.Ciphertext, op1 rlwe.Operand, relin bool, op2 *rlwe.Ciphertext) {

	_, level := eval.CheckBinary(op0, op1, op2, utils.MaxInt(op0.Degree(), op1.Degree()))

	if op0.Degree()+op1.Degree() > 2 {
		panic("cannot MulRelinThenAdd: the sum of the input elements' degree cannot be larger than 2")
	}

	if op0.El() == op2.El() || op1.El() == op2.El() {
		panic("cannot MulRelinThenAdd: op2 must be different from op0 and op1")
	}

	resScale := op0.Scale.Mul(op1.GetMetaData().Scale)

	if op2.Scale.Cmp(resScale) == -1 {
		ratio := resScale.Div(op2.Scale)
		// Only scales up if int(ratio) >= 2
		if ratio.Float64() >= 2.0 {
			eval.Mul(op2, &ratio.Value, op2)
			op2.Scale = resScale
		}
	}

	ringQ := eval.params.RingQ().AtLevel(level)

	var c00, c01, c0, c1, c2 *ring.Poly

	// Case Ciphertext (x) Ciphertext
	if op0.Degree() == 1 && op1.Degree() == 1 {

		c00 = eval.buffQ[0]
		c01 = eval.buffQ[1]

		c0 = op2.Value[0]
		c1 = op2.Value[1]

		if !relin {
			op2.El().Resize(2, level)
			c2 = op2.Value[2]
		} else {
			// No resize here since we add on op2
			c2 = eval.buffQ[2]
		}

		tmp0, tmp1 := op0.El(), op1.El()

		ringQ.MForm(tmp0.Value[0], c00)
		ringQ.MForm(tmp0.Value[1], c01)

		ringQ.MulCoeffsMontgomeryThenAdd(c00, tmp1.Value[0], c0) // c0 += c[0]*c[0]
		ringQ.MulCoeffsMontgomeryThenAdd(c00, tmp1.Value[1], c1) // c1 += c[0]*c[1]
		ringQ.MulCoeffsMontgomeryThenAdd(c01, tmp1.Value[0], c1) // c1 += c[1]*c[0]

		if relin {

			var rlk *rlwe.RelinearizationKey
			var err error
			if rlk, err = eval.CheckAndGetRelinearizationKey(); err != nil {
				panic(fmt.Errorf("cannot relinearize: %w", err))
			}

			ringQ.MulCoeffsMontgomery(c01, tmp1.Value[1], c2) // c2 += c[1]*c[1]

			tmpCt := &rlwe.Ciphertext{}
			tmpCt.Value = []*ring.Poly{eval.BuffQP[1].Q, eval.BuffQP[2].Q}
			tmpCt.IsNTT = true

			eval.GadgetProduct(level, c2, &rlk.GadgetCiphertext, tmpCt)
			ringQ.Add(c0, tmpCt.Value[0], c0)
			ringQ.Add(c1, tmpCt.Value[1], c1)
		} else {
			ringQ.MulCoeffsMontgomeryThenAdd(c01, tmp1.Value[1], c2) // c2 += c[1]*c[1]
		}

		// Case Plaintext (x) Ciphertext or Ciphertext (x) Plaintext
	} else {

		if op2.Degree() < op0.Degree() {
			op2.Resize(op0.Degree(), level)
		}

		c00 := eval.buffQ[0]

		ringQ.MForm(op1.El().Value[0], c00)
		for i := range op0.Value {
			ringQ.MulCoeffsMontgomeryThenAdd(op0.Value[i], c00, op2.Value[i])
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

// ApplyEvaluationKeyNew applies the rlwe.EvaluationKey on ct0 and returns the result on a new ciphertext ctOut.
func (eval *evaluator) ApplyEvaluationKeyNew(ct0 *rlwe.Ciphertext, evk *rlwe.EvaluationKey) (ctOut *rlwe.Ciphertext) {
	ctOut = NewCiphertext(eval.params, ct0.Degree(), ct0.Level())
	eval.ApplyEvaluationKey(ct0, evk, ctOut)
	return
}

// RotateNew rotates the columns of ct0 by k positions to the left, and returns the result in a newly created element.
// The method will panic if the evaluator hasn't been given an evaluation key set with the appropriate GaloisKey.
func (eval *evaluator) RotateNew(ct0 *rlwe.Ciphertext, k int) (ctOut *rlwe.Ciphertext) {
	ctOut = NewCiphertext(eval.params, ct0.Degree(), ct0.Level())
	eval.Rotate(ct0, k, ctOut)
	return
}

// Rotate rotates the columns of ct0 by k positions to the left and returns the result in ctOut.
// The method will panic if the evaluator hasn't been given an evaluation key set with the appropriate GaloisKey.
func (eval *evaluator) Rotate(ct0 *rlwe.Ciphertext, k int, ctOut *rlwe.Ciphertext) {
	eval.Automorphism(ct0, eval.params.GaloisElementForColumnRotationBy(k), ctOut)
}

// ConjugateNew conjugates ct0 (which is equivalent to a row rotation) and returns the result in a newly created element.
// The method will panic if the evaluator hasn't been given an evaluation key set with the appropriate GaloisKey.
func (eval *evaluator) ConjugateNew(ct0 *rlwe.Ciphertext) (ctOut *rlwe.Ciphertext) {

	if eval.params.RingType() == ring.ConjugateInvariant {
		panic("cannot ConjugateNew: method is not supported when params.RingType() == ring.ConjugateInvariant")
	}

	ctOut = NewCiphertext(eval.params, ct0.Degree(), ct0.Level())
	eval.Conjugate(ct0, ctOut)
	return
}

// Conjugate conjugates ct0 (which is equivalent to a row rotation) and returns the result in ctOut.
// The method will panic if the evaluator hasn't been given an evaluation key set with the appropriate GaloisKey.
func (eval *evaluator) Conjugate(ct0 *rlwe.Ciphertext, ctOut *rlwe.Ciphertext) {

	if eval.params.RingType() == ring.ConjugateInvariant {
		panic("cannot Conjugate: method is not supported when params.RingType() == ring.ConjugateInvariant")
	}

	eval.Automorphism(ct0, eval.params.GaloisElementForRowRotation(), ctOut)
}

func (eval *evaluator) RotateHoistedLazyNew(level int, rotations []int, ct *rlwe.Ciphertext, c2DecompQP []ringqp.Poly) (cOut map[int]*rlwe.OperandQP) {
	cOut = make(map[int]*rlwe.OperandQP)
	for _, i := range rotations {
		if i != 0 {
			cOut[i] = rlwe.NewOperandQP(eval.params.Parameters, 1, level, eval.params.MaxLevelP())
			eval.AutomorphismHoistedLazy(level, ct, c2DecompQP, eval.params.GaloisElementForColumnRotationBy(i), cOut[i])
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
func (eval *evaluator) WithKey(evk rlwe.EvaluationKeySetInterface) Evaluator {
	return &evaluator{
		Evaluator:        eval.Evaluator.WithKey(evk),
		evaluatorBase:    eval.evaluatorBase,
		evaluatorBuffers: eval.evaluatorBuffers,
	}
}
