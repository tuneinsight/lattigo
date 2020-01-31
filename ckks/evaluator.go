package ckks

import (
	"errors"
	"github.com/ldsec/lattigo/ring"
	"github.com/ldsec/lattigo/utils"
	"math"
)

// Evaluator is an interface implementing the methodes to conduct homomorphic operations between ciphertext and/or plaintexts.
type Evaluator interface {
	Add(op0, op1 Operand, ctOut *Ciphertext)
	AddNoMod(op0, op1 Operand, ctOut *Ciphertext)
	AddNew(op0, op1 Operand) (ctOut *Ciphertext)
	AddNoModNew(op0, op1 Operand) (ctOut *Ciphertext)
	Sub(op0, op1 Operand, ctOut *Ciphertext)
	SubNoMod(op0, op1 Operand, ctOut *Ciphertext)
	SubNew(op0, op1 Operand) (ctOut *Ciphertext)
	SubNoModNew(op0, op1 Operand) (ctOut *Ciphertext)
	Neg(ct0 *Ciphertext, ctOut *Ciphertext)
	NegNew(ct0 *Ciphertext) (ctOut *Ciphertext)
	AddConstNew(ct0 *Ciphertext, constant interface{}) (ctOut *Ciphertext)
	AddConst(ct0 *Ciphertext, constant interface{}, ctOut *Ciphertext)
	MultByConstAndAdd(ct0 *Ciphertext, constant interface{}, ctOut *Ciphertext)
	MultByConstNew(ct0 *Ciphertext, constant interface{}) (ctOut *Ciphertext)
	MultByConst(ct0 *Ciphertext, constant interface{}, ctOut *Ciphertext)
	MultByiNew(ct0 *Ciphertext) (ctOut *Ciphertext)
	MultByi(ct0 *Ciphertext, ct1 *Ciphertext)
	DivByiNew(ct0 *Ciphertext) (ctOut *Ciphertext)
	DivByi(ct0 *Ciphertext, ct1 *Ciphertext)
	ScaleUpNew(ct0 *Ciphertext, scale float64) (ctOut *Ciphertext)
	ScaleUp(ct0 *Ciphertext, scale float64, ctOut *Ciphertext)
	SetScale(ct *Ciphertext, scale float64)
	MulByPow2New(ct0 *Ciphertext, pow2 uint64) (ctOut *Ciphertext)
	MulByPow2(ct0 *ckksElement, pow2 uint64, ctOut *ckksElement)
	ReduceNew(ct0 *Ciphertext) (ctOut *Ciphertext)
	Reduce(ct0 *Ciphertext, ctOut *Ciphertext) error
	DropLevelNew(ct0 *Ciphertext, levels uint64) (ctOut *Ciphertext)
	DropLevel(ct0 *Ciphertext, levels uint64) (err error)
	Rescale(ct0 *Ciphertext, threshold float64, c1 *Ciphertext) (err error)
	RescaleMany(ct0 *Ciphertext, nbRescales uint64, c1 *Ciphertext) (err error)
	MulRelinNew(op0, op1 Operand, evakey *EvaluationKey) (ctOut *Ciphertext)
	MulRelin(op0, op1 Operand, evakey *EvaluationKey, ctOut *Ciphertext)
	RelinearizeNew(ct0 *Ciphertext, evakey *EvaluationKey) (ctOut *Ciphertext)
	Relinearize(ct0 *Ciphertext, evakey *EvaluationKey, ctOut *Ciphertext)
	SwitchKeysNew(ct0 *Ciphertext, switchingKey *SwitchingKey) (ctOut *Ciphertext)
	SwitchKeys(ct0 *Ciphertext, switchingKey *SwitchingKey, ctOut *Ciphertext)
	RotateColumnsNew(ct0 *Ciphertext, k uint64, evakey *RotationKeys) (ctOut *Ciphertext)
	RotateColumns(ct0 *Ciphertext, k uint64, evakey *RotationKeys, ctOut *Ciphertext)
	RotateHoisted(ctIn *Ciphertext, rotations []uint64, rotkeys *RotationKeys) (cOut map[uint64]*Ciphertext)
	ConjugateNew(ct0 *Ciphertext, evakey *RotationKeys) (ctOut *Ciphertext)
	Conjugate(ct0 *Ciphertext, evakey *RotationKeys, ctOut *Ciphertext)
	PowerOf2(el0 *Ciphertext, logPow2 uint64, evakey *EvaluationKey, elOut *Ciphertext)
	PowerNew(op *Ciphertext, degree uint64, evakey *EvaluationKey) (opOut *Ciphertext)
	Power(ct0 *Ciphertext, degree uint64, evakey *EvaluationKey, res *Ciphertext)
	InverseNew(ct0 *Ciphertext, steps uint64, evakey *EvaluationKey) (res *Ciphertext)
	EvaluatePolyFast(ct *Ciphertext, coeffs interface{}, evakey *EvaluationKey) (res *Ciphertext)
	EvaluatePolyEco(ct *Ciphertext, coeffs interface{}, evakey *EvaluationKey) (res *Ciphertext)
	EvaluateChebyFast(ct *Ciphertext, cheby *ChebyshevInterpolation, evakey *EvaluationKey) (res *Ciphertext)
	EvaluateChebyEco(ct *Ciphertext, cheby *ChebyshevInterpolation, evakey *EvaluationKey) (res *Ciphertext)
	EvaluateChebyEcoSpecial(ct *Ciphertext, n complex128, cheby *ChebyshevInterpolation, evakey *EvaluationKey) (res *Ciphertext)
	EvaluateChebyFastSpecial(ct *Ciphertext, n complex128, cheby *ChebyshevInterpolation, evakey *EvaluationKey) (res *Ciphertext)
	Bootstrapp(ct *Ciphertext, bootcontext *BootContext) (res *Ciphertext)
	BootstrappBetterSine(ct *Ciphertext, bootcontext *BootContext) (res *Ciphertext)
}

// evaluator is a struct that holds the necessary elements to execute the homomorphic operations between Ciphertexts and/or Plaintexts.
// It also holds a small memory pool used to store intermediate computations.
type evaluator struct {
	params      *Parameters
	ckksContext *Context
	ringpool    [6]*ring.Poly

	poolQ [4]*ring.Poly
	poolP [3]*ring.Poly

	ctxpool *Ciphertext

	baseconverter *ring.FastBasisExtender
	decomposer    *ring.Decomposer
}

// NewEvaluator creates a new Evaluator, that can be used to do homomorphic
// operations on the Ciphertexts and/or Plaintexts. It stores a small pool of polynomials
// and Ciphertexts that will be used for intermediate values.
func NewEvaluator(params *Parameters) Evaluator {

	if !params.isValid {
		panic("cannot newEvaluator: parameters are invalid (check if the generation was done properly)")
	}

	ckksContext := newContext(params)
	q := ckksContext.contextQ
	p := ckksContext.contextP

	return &evaluator{
		params:        params.Copy(),
		ckksContext:   ckksContext,
		ringpool:      [6]*ring.Poly{q.NewPoly(), q.NewPoly(), q.NewPoly(), q.NewPoly(), q.NewPoly(), q.NewPoly()},
		poolQ:         [4]*ring.Poly{q.NewPoly(), q.NewPoly(), q.NewPoly(), q.NewPoly()},
		poolP:         [3]*ring.Poly{p.NewPoly(), p.NewPoly(), p.NewPoly()},
		ctxpool:       NewCiphertext(params, 1, params.MaxLevel(), params.Scale),
		baseconverter: ring.NewFastBasisExtender(q, p),
		decomposer:    ring.NewDecomposer(q.Modulus, p.Modulus),
	}
}

func (eval *evaluator) getElemAndCheckBinary(op0, op1, opOut Operand, opOutMinDegree uint64) (el0, el1, elOut *ckksElement) {
	if op0 == nil || op1 == nil || opOut == nil {
		panic("operands cannot be nil")
	}

	if op0.Degree()+op1.Degree() == 0 {
		panic("operands cannot be both plaintext")
	}

	if opOut.Degree() < opOutMinDegree {
		panic("receiver operand degree is too small")
	}
	el0, el1, elOut = op0.Element(), op1.Element(), opOut.Element()
	return // TODO: more checks on elements
}

func (eval *evaluator) getElemAndCheckUnary(op0, opOut Operand, opOutMinDegree uint64) (el0, elOut *ckksElement) {
	if op0 == nil || opOut == nil {
		panic("operand cannot be nil")
	}

	if op0.Degree() == 0 {
		panic("operand cannot be plaintext")
	}

	if opOut.Degree() < opOutMinDegree {
		panic("receiver operand degree is too small")
	}
	el0, elOut = op0.Element(), opOut.Element()
	return // TODO: more checks on elements
}

func (eval *evaluator) newCiphertextBinary(op0, op1 Operand) (ctOut *Ciphertext) {

	maxDegree := utils.MaxUint64(op0.Degree(), op1.Degree())
	maxScale := utils.MaxFloat64(op0.Scale(), op1.Scale())
	minLevel := utils.MinUint64(op0.Level(), op1.Level())

	return NewCiphertext(eval.params, maxDegree, minLevel, maxScale)
}

// Add adds op0 to op1 and returns the result in ctOut.
func (eval *evaluator) Add(op0, op1 Operand, ctOut *Ciphertext) {
	el0, el1, elOut := eval.getElemAndCheckBinary(op0, op1, ctOut, utils.MaxUint64(op0.Degree(), op1.Degree()))
	eval.evaluateInPlace(el0, el1, elOut, eval.ckksContext.contextQ.AddLvl)
}

// AddNoMod adds op0 to op1 and returns the result in ctOut, without modular reduction.
func (eval *evaluator) AddNoMod(op0, op1 Operand, ctOut *Ciphertext) {
	el0, el1, elOut := eval.getElemAndCheckBinary(op0, op1, ctOut, utils.MaxUint64(op0.Degree(), op1.Degree()))
	eval.evaluateInPlace(el0, el1, elOut, eval.ckksContext.contextQ.AddNoModLvl)
}

// AddNew adds op0 to op1 and returns the result in a newly created element.
func (eval *evaluator) AddNew(op0, op1 Operand) (ctOut *Ciphertext) {
	ctOut = eval.newCiphertextBinary(op0, op1)
	eval.Add(op0, op1, ctOut)
	return
}

// AddNoModNew adds op0 to op1 without modular reduction, and returns the result in a newly created element.
func (eval *evaluator) AddNoModNew(op0, op1 Operand) (ctOut *Ciphertext) {
	ctOut = eval.newCiphertextBinary(op0, op1)
	eval.AddNoMod(op0, op1, ctOut)
	return
}

// Sub subtracts op1 from op0 and returns the result in ctOut.
func (eval *evaluator) Sub(op0, op1 Operand, ctOut *Ciphertext) {

	el0, el1, elOut := eval.getElemAndCheckBinary(op0, op1, ctOut, utils.MaxUint64(op0.Degree(), op1.Degree()))

	eval.evaluateInPlace(el0, el1, elOut, eval.ckksContext.contextQ.SubLvl)

	level := utils.MinUint64(utils.MinUint64(el0.Level(), el1.Level()), elOut.Level())

	if el0.Degree() < el1.Degree() {
		for i := el0.Degree() + 1; i < el1.Degree()+1; i++ {
			eval.ckksContext.contextQ.NegLvl(level, elOut.Value()[i], elOut.Value()[i])
		}
	}

}

// SubNoMod subtracts op1 from op0 and returns the result in ctOut, without modular reduction.
func (eval *evaluator) SubNoMod(op0, op1 Operand, ctOut *Ciphertext) {

	el0, el1, elOut := eval.getElemAndCheckBinary(op0, op1, ctOut, utils.MaxUint64(op0.Degree(), op1.Degree()))

	eval.evaluateInPlace(el0, el1, elOut, eval.ckksContext.contextQ.SubNoModLvl)

	level := utils.MinUint64(utils.MinUint64(el0.Level(), el1.Level()), elOut.Level())

	if el0.Degree() < el1.Degree() {
		for i := el0.Degree() + 1; i < el1.Degree()+1; i++ {
			eval.ckksContext.contextQ.NegLvl(level, elOut.Value()[i], elOut.Value()[i])
		}
	}

}

// SubNew subtracts op1 from op0 and returns the result in a newly created element.
func (eval *evaluator) SubNew(op0, op1 Operand) (ctOut *Ciphertext) {
	ctOut = eval.newCiphertextBinary(op0, op1)
	eval.Sub(op0, op1, ctOut)
	return
}

// SubNoModNew subtracts op1 from op0 without modular reduction, and returns the result in a newly created element.
func (eval *evaluator) SubNoModNew(op0, op1 Operand) (ctOut *Ciphertext) {
	ctOut = eval.newCiphertextBinary(op0, op1)
	eval.SubNoMod(op0, op1, ctOut)
	return
}

func (eval *evaluator) evaluateInPlace(c0, c1, ctOut *ckksElement, evaluate func(uint64, *ring.Poly, *ring.Poly, *ring.Poly)) {

	var tmp0, tmp1 *ckksElement // TODO : use eval mem pool

	level := utils.MinUint64(utils.MinUint64(c0.Level(), c1.Level()), ctOut.Level())

	maxDegree := utils.MaxUint64(c0.Degree(), c1.Degree())
	minDegree := utils.MinUint64(c0.Degree(), c1.Degree())

	// Else resizes the receiver element
	ctOut.Resize(eval.params, maxDegree)
	eval.DropLevel(ctOut.Ciphertext(), ctOut.Level()-utils.MinUint64(c0.Level(), c1.Level()))

	// Checks whether or not the receiver element is the same as one of the input elements
	// and acts accordingly to avoid unnecessary element creation or element overwriting,
	// and scales properly the element before the evaluation.
	if ctOut == c0 {

		if c0.Scale() > c1.Scale() {

			tmp1 = eval.ctxpool.Element()

			if uint64(c0.Scale()/c1.Scale()) != 0 {
				eval.MultByConst(c1.Ciphertext(), uint64(c0.Scale()/c1.Scale()), tmp1.Ciphertext())
			}

		} else if c1.Scale() > c0.Scale() {

			if uint64(c1.Scale()/c0.Scale()) != 0 {
				eval.MultByConst(c0.Ciphertext(), uint64(c1.Scale()/c0.Scale()), c0.Ciphertext())
			}

			c0.SetScale(c1.Scale())

			tmp1 = c1

		} else {

			tmp1 = c1
		}

		tmp0 = c0

	} else if ctOut == c1 {

		if c1.Scale() > c0.Scale() {

			tmp0 = eval.ctxpool.Element()
			if uint64(c1.Scale()/c0.Scale()) != 0 {
				eval.MultByConst(c0.Ciphertext(), uint64(c1.Scale()/c0.Scale()), tmp0.Ciphertext())
			}

		} else if c0.Scale() > c1.Scale() {

			if uint64(c0.Scale()/c1.Scale()) != 0 {
				eval.MultByConst(c1.Ciphertext(), uint64(c0.Scale()/c1.Scale()), ctOut.Ciphertext())
			}

			ctOut.SetScale(c0.Scale())

			tmp0 = c0

		} else {

			tmp0 = c0
		}

		tmp1 = c1

	} else {

		if c1.Scale() > c0.Scale() {

			tmp0 = eval.ctxpool.Element()

			if uint64(c1.Scale()/c0.Scale()) != 0 {
				eval.MultByConst(c0.Ciphertext(), uint64(c1.Scale()/c0.Scale()), tmp0.Ciphertext())
			}

			tmp1 = c1

		} else if c0.Scale() > c1.Scale() {

			tmp1 = eval.ctxpool.Element()

			if uint64(c0.Scale()/c1.Scale()) != 0 {
				eval.MultByConst(c1.Ciphertext(), uint64(c0.Scale()/c1.Scale()), tmp1.Ciphertext())
			}

			tmp0 = c0

		} else {
			tmp0 = c0
			tmp1 = c1
		}
	}

	for i := uint64(0); i < minDegree+1; i++ {
		evaluate(level, tmp0.Value()[i], tmp1.Value()[i], ctOut.Value()[i])
	}

	ctOut.SetScale(utils.MaxFloat64(c0.Scale(), c1.Scale()))

	// If the inputs degrees differ, it copies the remaining degree on the receiver.
	// Also checks that the receiver is not one of the inputs to avoid unnecessary work.

	if c0.Degree() > c1.Degree() && tmp0 != ctOut {
		for i := minDegree + 1; i < maxDegree+1; i++ {
			eval.ckksContext.contextQ.CopyLvl(level, tmp0.Value()[i], ctOut.Value()[i])
		}
	} else if c1.Degree() > c0.Degree() && tmp1 != ctOut {
		for i := minDegree + 1; i < maxDegree+1; i++ {
			eval.ckksContext.contextQ.CopyLvl(level, tmp1.Value()[i], ctOut.Value()[i])
		}
	}
}

// Neg negates the value of ct0 and returns the result in ctOut.
func (eval *evaluator) Neg(ct0 *Ciphertext, ctOut *Ciphertext) {

	level := utils.MinUint64(ct0.Level(), ctOut.Level())

	if ct0.Degree() != ctOut.Degree() {
		panic("cannot Negate: invalid receiver Ciphertext does not match input Ciphertext degree")
	}

	for i := range ct0.value {
		eval.ckksContext.contextQ.NegLvl(level, ct0.value[i], ctOut.Value()[i])
	}
}

// NegNew negates ct0 and returns the result in a newly created element.
func (eval *evaluator) NegNew(ct0 *Ciphertext) (ctOut *Ciphertext) {
	ctOut = NewCiphertext(eval.params, ct0.Degree(), ct0.Level(), ct0.Scale())
	eval.Neg(ct0, ctOut)
	return
}

// AddConstNew adds the input constant (which can be a uint64, int64, float64 or complex128) to ct0 and returns the result in a new element.
func (eval *evaluator) AddConstNew(ct0 *Ciphertext, constant interface{}) (ctOut *Ciphertext) {
	ctOut = ct0.CopyNew().Ciphertext()
	eval.AddConst(ct0, constant, ctOut)
	return ctOut
}

// AddConst adds the input constant (which can be a uint64, int64, float64 or complex128) to ct0 and returns the result in ctOut.
func (eval *evaluator) AddConst(ct0 *Ciphertext, constant interface{}, ctOut *Ciphertext) {

	var level uint64

	level = utils.MinUint64(ct0.Level(), ctOut.Level())

	var cReal, cImag float64

	switch constant.(type) {
	case complex128:
		cReal = real(constant.(complex128))
		cImag = imag(constant.(complex128))

	case float64:
		cReal = constant.(float64)
		cImag = float64(0)

	case uint64:
		cReal = float64(constant.(uint64))
		cImag = float64(0)

	case int64:
		cReal = float64(constant.(int64))
		cImag = float64(0)

	case int:
		cReal = float64(constant.(int))
		cImag = float64(0)
	}

	var scaledConst, scaledConstReal, scaledConstImag uint64

	context := eval.ckksContext.contextQ

	// Component wise addition of the following vector to the ciphertext:
	// [a + b*psi_qi^2, ....., a + b*psi_qi^2, a - b*psi_qi^2, ...., a - b*psi_qi^2] mod Qi
	// [{                  N/2                }{                N/2               }]
	// Which is equivalent outside of the NTT domain to adding a to the first coefficient of ct0 and b to the N/2-th coefficient of ct0.
	var qi uint64
	for i := uint64(0); i < level+1; i++ {
		scaledConstReal = 0
		scaledConstImag = 0
		scaledConst = 0

		qi = context.Modulus[i]

		if cReal != 0 {
			scaledConstReal = scaleUpExact(cReal, ctOut.Scale(), qi)
			scaledConst = scaledConstReal
		}

		if cImag != 0 {
			scaledConstImag = ring.MRed(scaleUpExact(cImag, ctOut.Scale(), qi), context.GetNttPsi()[i][1], qi, context.GetMredParams()[i])
			scaledConst = ring.CRed(scaledConst+scaledConstImag, qi)
		}

		p1tmp := ctOut.Value()[0].Coeffs[i]
		p0tmp := ct0.value[0].Coeffs[i]

		for j := uint64(0); j < eval.ckksContext.n>>1; j++ {
			p1tmp[j] = ring.CRed(p0tmp[j]+scaledConst, qi)
		}

		if cImag != 0 {
			scaledConst = ring.CRed(scaledConstReal+(qi-scaledConstImag), qi)
		}

		for j := eval.ckksContext.n >> 1; j < eval.ckksContext.n; j++ {
			p1tmp[j] = ring.CRed(p0tmp[j]+scaledConst, qi)
		}
	}
}

// MultByConstAndAdd multiplies ct0 by the input constant, and adds it to the receiver element (it does not modify the input
// element), e.g., ctOut(x) = ctOut(x) + ct0(x) * (a+bi). This functions removes the need of storing the intermediate value c(x) * (a+bi).
// This function will modify the level and the scale of the receiver element depending on the level and the scale of the input
// element and the type of the constant. The level of the receiver element will be set to min(input.level, receiver.level).
// The scale of the receiver element will be set to the scale that the input element would have after the multiplication by the constant.
func (eval *evaluator) MultByConstAndAdd(ct0 *Ciphertext, constant interface{}, ctOut *Ciphertext) {

	var level uint64

	level = utils.MinUint64(ct0.Level(), ctOut.Level())

	// Forces a drop of ctOut level to ct0 level
	if ctOut.Level() > level {
		eval.DropLevel(ctOut, ctOut.Level()-level)
	}

	var cReal, cImag float64
	var scale float64

	// Converts to float64 and determines if a scaling is required (which is the case if either real or imag have a rational part)
	scale = 1
	switch constant.(type) {
	case complex128:
		cReal = real(constant.(complex128))
		cImag = imag(constant.(complex128))

		if cReal != 0 {
			valueInt := int64(cReal)
			valueFloat := cReal - float64(valueInt)

			if valueFloat != 0 {
				scale = float64(eval.ckksContext.contextQ.Modulus[level])
			}
		}

		if cImag != 0 {
			valueInt := int64(cImag)
			valueFloat := cImag - float64(valueInt)

			if valueFloat != 0 {
				scale = float64(eval.ckksContext.contextQ.Modulus[level])
			}
		}

	case float64:
		cReal = constant.(float64)
		cImag = float64(0)

		if cReal != 0 {
			valueInt := int64(cReal)
			valueFloat := cReal - float64(valueInt)

			if valueFloat != 0 {
				scale = float64(eval.ckksContext.contextQ.Modulus[level])
			}
		}

	case uint64:
		cReal = float64(constant.(uint64))
		cImag = float64(0)

	case int64:
		cReal = float64(constant.(int64))
		cImag = float64(0)

	case int:
		cReal = float64(constant.(int))
		cImag = float64(0)
	}

	var scaledConst, scaledConstReal, scaledConstImag uint64

	context := eval.ckksContext.contextQ

	// If a scaling would be required to multiply by the constant,
	// it equalizes scales such that the scales match in the end.
	if scale != 1 {

		// If ctOut scaling is smaller than ct0's scale + the default scaling,
		// then brings ctOut scale to ct0's scale.
		if ctOut.Scale() < ct0.Scale()*scale {

			if uint64((scale*ct0.Scale())/ctOut.Scale()) != 0 {

				eval.MultByConst(ctOut, uint64((scale*ct0.Scale())/ctOut.Scale()), ctOut)

			}

			ctOut.SetScale(scale * ct0.Scale())

			// If ctOut.Scale() > ((a+bi)*scale)*ct0(x), then it sets the scale to
			// bring c(x)*scale to the level of ctOut(x) scale
		} else if ctOut.Scale() > ct0.Scale()*scale {
			scale = ctOut.Scale() / ct0.Scale()
		}

		// If no scaling is required, then it sets the appropriate scale such that
		// ct0(x)*scale matches ctOut(x) scale without modifying ct0(x) scale.
	} else {

		if ctOut.Scale() > ct0.Scale() {

			scale = ctOut.Scale() / ct0.Scale()

		} else if ct0.Scale() > ctOut.Scale() {

			if uint64(ct0.Scale()/ctOut.Scale()) != 0 {
				eval.MultByConst(ctOut, ct0.Scale()/ctOut.Scale(), ctOut)
			}

			ctOut.SetScale(ct0.Scale())
		}
	}

	// Component-wise multiplication of the following vector to the ciphertext:
	// [a + b*psi_qi^2, ....., a + b*psi_qi^2, a - b*psi_qi^2, ...., a - b*psi_qi^2] mod Qi
	// [{                  N/2                }{                N/2               }]
	// Which is equivalent outside of the NTT domain to adding a to the first coefficient of ct0 and b to the N/2-th coefficient of ct0.
	for i := uint64(0); i < level+1; i++ {

		qi := context.Modulus[i]
		mredParams := context.GetMredParams()[i]
		bredParams := context.GetBredParams()[i]

		scaledConstReal = 0
		scaledConstImag = 0
		scaledConst = 0

		if cReal != 0 {
			scaledConstReal = scaleUpExact(cReal, scale, qi)
			scaledConst = scaledConstReal
		}

		if cImag != 0 {
			scaledConstImag = scaleUpExact(cImag, scale, qi)
			scaledConstImag = ring.MRed(scaledConstImag, context.GetNttPsi()[i][1], qi, mredParams)
			scaledConst = ring.CRed(scaledConst+scaledConstImag, qi)
		}

		scaledConst = ring.MForm(scaledConst, qi, bredParams)

		for u := range ct0.Value() {
			p0tmp := ct0.Value()[u].Coeffs[i]
			p1tmp := ctOut.Value()[u].Coeffs[i]
			for j := uint64(0); j < eval.ckksContext.n>>1; j++ {
				p1tmp[j] = ring.CRed(p1tmp[j]+ring.MRed(p0tmp[j], scaledConst, qi, mredParams), qi)
			}
		}

		if cImag != 0 {
			scaledConst = ring.CRed(scaledConstReal+(qi-scaledConstImag), qi)
			scaledConst = ring.MForm(scaledConst, qi, bredParams)
		}

		for u := range ct0.Value() {
			p0tmp := ct0.Value()[u].Coeffs[i]
			p1tmp := ctOut.Value()[u].Coeffs[i]
			for j := eval.ckksContext.n >> 1; j < eval.ckksContext.n; j++ {
				p1tmp[j] = ring.CRed(p1tmp[j]+ring.MRed(p0tmp[j], scaledConst, qi, mredParams), qi)
			}
		}
	}
}

// MultByConstNew multiplies ct0 by the input constant and returns the result in a newly created element.
// The scale of the output element will depend on the scale of the input element and the constant (if the constant
// needs to be scaled (its rational part is not zero)). The constant can be a uint64, int64, float64 or complex128.
func (eval *evaluator) MultByConstNew(ct0 *Ciphertext, constant interface{}) (ctOut *Ciphertext) {
	ctOut = NewCiphertext(eval.params, ct0.Degree(), ct0.Level(), ct0.Scale())
	eval.MultByConst(ct0, constant, ctOut)
	return
}

// MultByConst multiplies ct0 by the input constant and returns the result in ctOut.
// The scale of the output element will depend on the scale of the input element and the constant (if the constant
// needs to be scaled (its rational part is not zero)). The constant can be a uint64, int64, float64 or complex128.
func (eval *evaluator) MultByConst(ct0 *Ciphertext, constant interface{}, ctOut *Ciphertext) {

	var level uint64

	level = utils.MinUint64(ct0.Level(), ctOut.Level())

	var cReal, cImag float64
	var scale float64

	// Converts to float64 and determines if a scaling is required (which is the case if either real or imag have a rational part)
	scale = 1
	switch constant.(type) {
	case complex128:
		cReal = real(constant.(complex128))
		cImag = imag(constant.(complex128))

		if cReal != 0 {
			valueInt := int64(cReal)
			valueFloat := cReal - float64(valueInt)

			if valueFloat != 0 {
				scale = float64(eval.ckksContext.contextQ.Modulus[level])
			}
		}

		if cImag != 0 {
			valueInt := int64(cImag)
			valueFloat := cImag - float64(valueInt)

			if valueFloat != 0 {
				scale = float64(eval.ckksContext.contextQ.Modulus[level])
			}
		}

	case float64:
		cReal = constant.(float64)
		cImag = float64(0)

		if cReal != 0 {
			valueInt := int64(cReal)
			valueFloat := cReal - float64(valueInt)

			if valueFloat != 0 {
				scale = float64(eval.ckksContext.contextQ.Modulus[level])
			}
		}

	case uint64:
		cReal = float64(constant.(uint64))
		cImag = float64(0)

	case int64:
		cReal = float64(constant.(int64))
		cImag = float64(0)

	case int:
		cReal = float64(constant.(int))
		cImag = float64(0)
	}

	// Component wise multiplication of the following vector with the ciphertext:
	// [a + b*psi_qi^2, ....., a + b*psi_qi^2, a - b*psi_qi^2, ...., a - b*psi_qi^2] mod Qi
	// [{                  N/2                }{                N/2               }]
	// Which is equivalent outside of the NTT domain to adding a to the first coefficient of ct0 and b to the N/2-th coefficient of ct0.
	context := eval.ckksContext.contextQ
	var scaledConst, scaledConstReal, scaledConstImag uint64
	for i := uint64(0); i < level+1; i++ {

		qi := context.Modulus[i]
		bredParams := context.GetBredParams()[i]
		mredParams := context.GetMredParams()[i]

		scaledConstReal = 0
		scaledConstImag = 0
		scaledConst = 0

		if cReal != 0 {
			scaledConstReal = scaleUpExact(cReal, scale, qi)
			scaledConst = scaledConstReal
		}

		if cImag != 0 {
			scaledConstImag = scaleUpExact(cImag, scale, qi)
			scaledConstImag = ring.MRed(scaledConstImag, context.GetNttPsi()[i][1], qi, mredParams)
			scaledConst = ring.CRed(scaledConst+scaledConstImag, qi)
		}

		scaledConst = ring.MForm(scaledConst, qi, bredParams)

		for u := range ct0.Value() {
			p0tmp := ct0.Value()[u].Coeffs[i]
			p1tmp := ctOut.Value()[u].Coeffs[i]
			for j := uint64(0); j < eval.ckksContext.n>>1; j++ {
				p1tmp[j] = ring.MRed(p0tmp[j], scaledConst, qi, mredParams)
			}
		}

		if cImag != 0 {
			scaledConst = ring.CRed(scaledConstReal+(qi-scaledConstImag), qi)
			scaledConst = ring.MForm(scaledConst, qi, bredParams)
		}

		for u := range ct0.Value() {
			p0tmp := ct0.Value()[u].Coeffs[i]
			p1tmp := ctOut.Value()[u].Coeffs[i]
			for j := eval.ckksContext.n >> 1; j < eval.ckksContext.n; j++ {
				p1tmp[j] = ring.MRed(p0tmp[j], scaledConst, qi, mredParams)
			}
		}
	}

	ctOut.SetScale(ct0.Scale() * scale)
}

// MultByiNew multiplies ct0 by the imaginary number i, and returns the result in a newly created element.
// It does not change the scale.
func (eval *evaluator) MultByiNew(ct0 *Ciphertext) (ctOut *Ciphertext) {
	ctOut = NewCiphertext(eval.params, 1, ct0.Level(), ct0.Scale())
	eval.MultByi(ct0, ctOut)
	return ctOut
}

// MultByi multiplies ct0 by the imaginary number i, and returns the result in ctOut.
// It does not change the scale.
func (eval *evaluator) MultByi(ct0 *Ciphertext, ctOut *Ciphertext) {

	var level uint64

	level = utils.MinUint64(ct0.Level(), ctOut.Level())

	context := eval.ckksContext.contextQ

	var imag uint64

	// Equivalent to a product by the monomial x^(n/2) outside of the NTT domain
	for i := uint64(0); i < level+1; i++ {

		qi := context.Modulus[i]
		mredParams := context.GetMredParams()[i]

		imag = context.GetNttPsi()[i][1] // Psi^2

		for u := range ctOut.value {
			p0tmp := ct0.value[u].Coeffs[i]
			p1tmp := ctOut.value[u].Coeffs[i]
			for j := uint64(0); j < context.N>>1; j++ {
				p1tmp[j] = ring.MRed(p0tmp[j], imag, qi, mredParams)
			}
		}

		imag = qi - imag

		for u := range ctOut.value {
			p0tmp := ct0.value[u].Coeffs[i]
			p1tmp := ctOut.value[u].Coeffs[i]
			for j := context.N >> 1; j < context.N; j++ {
				p1tmp[j] = ring.MRed(p0tmp[j], imag, qi, mredParams)

			}
		}
	}
}

// DivByiNew multiplies ct0 by the imaginary number 1/i = -i, and returns the result in a newly created element.
// It does not change the scale.
func (eval *evaluator) DivByiNew(ct0 *Ciphertext) (ctOut *Ciphertext) {
	ctOut = NewCiphertext(eval.params, 1, ct0.Level(), ct0.Scale())
	eval.DivByi(ct0, ctOut)
	return
}

// DivByi multiplies ct0 by the imaginary number 1/i = -i, and returns the result in ctOut.
// It does not change the scale.
func (eval *evaluator) DivByi(ct0 *Ciphertext, ctOut *Ciphertext) {

	var level uint64

	level = utils.MinUint64(ct0.Level(), ctOut.Level())

	context := eval.ckksContext.contextQ

	var imag uint64

	// Equivalent to a product by the monomial x^(3*n/2) outside of the NTT domain
	for i := uint64(0); i < level+1; i++ {

		qi := context.Modulus[i]
		mredParams := context.GetMredParams()[i]

		imag = qi - context.GetNttPsi()[i][1] // -Psi^2

		for u := range ctOut.value {
			p0tmp := ct0.value[u].Coeffs[i]
			p1tmp := ctOut.value[u].Coeffs[i]
			for j := uint64(0); j < context.N>>1; j++ {
				p1tmp[j] = ring.MRed(p0tmp[j], imag, qi, mredParams)
			}
		}

		imag = context.GetNttPsi()[i][1] // Psi^2

		for u := range ctOut.value {
			p0tmp := ct0.value[u].Coeffs[i]
			p1tmp := ctOut.value[u].Coeffs[i]
			for j := context.N >> 1; j < context.N; j++ {
				p1tmp[j] = ring.MRed(p0tmp[j], imag, qi, mredParams)
			}
		}
	}
}

// ScaleUpNew multiplies ct0 by 2^scale and sets its scale to its previous scale
// plus 2^n. It returns the result in a newly created element.
func (eval *evaluator) ScaleUpNew(ct0 *Ciphertext, scale float64) (ctOut *Ciphertext) {
	ctOut = NewCiphertext(eval.params, ct0.Degree(), ct0.Level(), ct0.Scale())
	eval.ScaleUp(ct0, scale, ctOut)
	return
}

// ScaleUp multiplies ct0 by 2^scale and sets its scale to its previous scale
// plus 2^n. It returns the result in ctOut.
func (eval *evaluator) ScaleUp(ct0 *Ciphertext, scale float64, ctOut *Ciphertext) {
	eval.MultByConst(ct0, uint64(scale), ctOut)
	ctOut.SetScale(ct0.Scale() * scale)
	return
}

func (evaluator *evaluator) SetScale(ct *Ciphertext, scale float64) {

	var tmp float64

	tmp = evaluator.params.Scale

	evaluator.params.Scale = scale

	evaluator.MultByConst(ct, scale/ct.Scale(), ct)

	if err = evaluator.Rescale(ct, scale, ct); err != nil {
		panic(err)
	}

	ct.SetScale(scale)

	evaluator.params.Scale = tmp
}

// MulByPow2New multiplies ct0 by 2^pow2 and returns the result in a newly created element.
func (eval *evaluator) MulByPow2New(ct0 *Ciphertext, pow2 uint64) (ctOut *Ciphertext) {
	ctOut = NewCiphertext(eval.params, ct0.Degree(), ct0.Level(), ct0.Scale())
	eval.MulByPow2(ct0.Element(), pow2, ctOut.Element())
	return
}

// MulByPow2 multiplies ct0 by 2^pow2 and returns the result in ctOut.
func (eval *evaluator) MulByPow2(ct0 *ckksElement, pow2 uint64, ctOut *ckksElement) {
	var level uint64
	level = utils.MinUint64(ct0.Level(), ctOut.Level())
	for i := range ctOut.Value() {
		eval.ckksContext.contextQ.MulByPow2Lvl(level, ct0.value[i], pow2, ctOut.Value()[i])
	}
}

// ReduceNew applies a modular reduction to ct0 and returns the result in a newly created element.
// To be used in conjunction with functions that do not apply modular reduction.
func (eval *evaluator) ReduceNew(ct0 *Ciphertext) (ctOut *Ciphertext) {

	ctOut = NewCiphertext(eval.params, ct0.Degree(), ct0.Level(), ct0.Scale())

	_ = eval.Reduce(ct0, ctOut)

	return ctOut
}

// Reduce applies a modular reduction to ct0 and returns the result in ctOut.
// To be used in conjunction with functions that do not apply modular reduction.
func (eval *evaluator) Reduce(ct0 *Ciphertext, ctOut *Ciphertext) error {

	if ct0.Degree() != ctOut.Degree() {
		return errors.New("cannot Reduce: degrees of receiver Ciphertext and input Ciphertext do not match")
	}

	for i := range ct0.value {
		eval.ckksContext.contextQ.ReduceLvl(utils.MinUint64(ct0.Level(), ctOut.Level()), ct0.value[i], ctOut.value[i])
	}

	return nil
}

// DropLevelNew reduces the level of ct0 by levels and returns the result in a newly created element.
// No rescaling is applied during this procedure.
func (eval *evaluator) DropLevelNew(ct0 *Ciphertext, levels uint64) (ctOut *Ciphertext) {
	ctOut = ct0.CopyNew().Ciphertext()
	eval.DropLevel(ctOut, levels)
	return
}

// DropLevel reduces the level of ct0 by levels and returns the result in ct0.
// No rescaling is applied during this procedure.
func (eval *evaluator) DropLevel(ct0 *Ciphertext, levels uint64) (err error) {

	if ct0.Level() == 0 {
		return errors.New("cannot DropLevel: Ciphertext already at level 0")
	}

	level := ct0.Level()

	for i := range ct0.value {
		ct0.value[i].Coeffs = ct0.value[i].Coeffs[:level+1-levels]
	}

	return nil
}

// RescaleNew divides ct0 by the last modulus in the moduli chain, and repeats this
// procedure (consuming one level each time) until the scale reaches the original scale or before it goes below it, and returns the result
// in a newly created element. Since all the moduli in the moduli chain are generated to be close to the
// original scale, this procedure is equivalent to dividing the input element by the scale and adding
// some error.
func (eval *evaluator) RescaleNew(ct0 *Ciphertext, threshold float64) (ctOut *Ciphertext, err error) {

	ctOut = NewCiphertext(eval.params, ct0.Degree(), ct0.Level(), ct0.Scale())

	return ctOut, eval.Rescale(ct0, threshold, ctOut)
}

// Rescale divides ct0 by the last modulus in the moduli chain, and repeats this
// procedure (consuming one level each time) until the scale reaches the original scale or before it goes below it, and returns the result
// in ctOut. Since all the moduli in the moduli chain are generated to be close to the
// original scale, this procedure is equivalent to dividing the input element by the scale and adding
// some error.
func (eval *evaluator) Rescale(ct0 *Ciphertext, threshold float64, ctOut *Ciphertext) (err error) {

	ringContext := eval.ckksContext.contextQ

	if ct0.Level() == 0 {
		return errors.New("cannot Rescale: input Ciphertext already at level 0")
	}

	if ct0.Level() != ctOut.Level() {
		panic("cannot Rescale: degrees of receiver Ciphertext and input Ciphertext do not match")
	}

	if ct0.Scale() >= (threshold*float64(ringContext.Modulus[ctOut.Level()]))/2 {

		if !ct0.IsNTT() {
			panic("cannot Rescale: input Ciphertext not in NTT")
		}

		ctOut.Copy(ct0.Element())

		for ctOut.Scale() >= (threshold*float64(ringContext.Modulus[ctOut.Level()]))/2 && ctOut.Level() != 0 {

			ctOut.DivScale(float64(ringContext.Modulus[ctOut.Level()]))

			for i := range ctOut.Value() {
				eval.ckksContext.contextQ.DivRoundByLastModulusNTT(ctOut.Value()[i])
			}

		}

	} else {
		ctOut.Copy(ct0.Element())
	}

	return nil
}

// RescaleMany applies Rescale several times in a row on the input Ciphertext.
func (eval *evaluator) RescaleMany(ct0 *Ciphertext, nbRescales uint64, ctOut *Ciphertext) (err error) {

	if ct0.Level() < nbRescales {
		return errors.New("cannot RescaleMany: input Ciphertext level too low")
	}

	if ct0.Level() != ctOut.Level() {
		panic("cannot RescaleMany: degrees of receiver Ciphertext and input Ciphertext do not match")
	}

	if !ct0.IsNTT() {
		panic("cannot RescaleMany: input Ciphertext not in NTT")
	}

	ctOut.Copy(ct0.Element())

	for i := uint64(0); i < nbRescales; i++ {
		ctOut.DivScale(float64(eval.ckksContext.contextQ.Modulus[ctOut.Level()-i]))
	}

	for i := range ctOut.Value() {
		eval.ckksContext.contextQ.DivRoundByLastModulusManyNTT(ctOut.Value()[i], nbRescales)
	}

	return nil
}

// MulRelinNew multiplies ct0 by ct1 and returns the result in a newly created element. The new scale is
// the multiplication between the scales of the input elements (addition when the scale is represented in log2). An evaluation
// key can be provided to apply a relinearization step to reduce the degree of the output element. This evaluation key is only
// required when the two input elements are Ciphertexts. If no evaluation key is provided and the input elements are two Ciphertexts,
// the resulting Ciphertext will be of degree two. This function only accepts Plaintexts (degree zero) and/or Ciphertexts of degree one.
func (eval *evaluator) MulRelinNew(op0, op1 Operand, evakey *EvaluationKey) (ctOut *Ciphertext) {

	ctOut = NewCiphertext(eval.params, 1, utils.MinUint64(op0.Level(), op1.Level()), op0.Scale()+op1.Scale())
	eval.MulRelin(op0, op1, evakey, ctOut)

	return ctOut
}

// MulRelin multiplies ct0 by ct1 and returns the result in ctOut. The new scale is
// the multiplication between the scales of the input elements (addition when the scale is represented in log2). An evaluation
// key can be provided to apply a relinearization step to reduce the degree of the output element. This evaluation key is only
// required when the two input elements are Ciphertexts. If no evaluation key is provided and the input elements are two Ciphertexts,
// the resulting Ciphertext will be of degree two. This function only accepts Plaintexts (degree zero) and/or Ciphertexts of degree one.
func (eval *evaluator) MulRelin(op0, op1 Operand, evakey *EvaluationKey, ctOut *Ciphertext) {

	el0, el1, elOut := eval.getElemAndCheckBinary(op0, op1, ctOut, utils.MaxUint64(op0.Degree(), op1.Degree()))

	level := utils.MinUint64(utils.MinUint64(el0.Level(), el1.Level()), elOut.Level())

	if ctOut.Level() > level {
		eval.DropLevel(elOut.Ciphertext(), elOut.Level()-level)
	}

	if el0.Degree() > 1 || el1.Degree() > 1 {
		panic("cannot MulRelin: input elements must be of degree 0 or 1")
	}

	if !el0.IsNTT() {
		panic("cannot MulRelin: op0 must be in NTT")
	}

	if !el1.IsNTT() {
		panic("cannot MulRelin: op1 must be in NTT")
	}

	elOut.SetScale(el0.Scale() * el1.Scale())

	context := eval.ckksContext.contextQ

	var c00, c01, c0, c1, c2 *ring.Poly

	// Case Ciphertext (x) Ciphertext
	if el0.Degree()+el1.Degree() == 2 {

		c00 = eval.ringpool[0]
		c01 = eval.ringpool[1]

		// If the receiver Ciphertext is neither of the inputs,
		// we can write directly on it.
		if elOut != el0 && elOut != el1 {

			c0 = elOut.value[0]
			c1 = elOut.value[1]

			// If the evaluation key is nil and we can write directly on the receiver, then
			// it resizes the Ciphertext to a degree 2 Ciphertext
			if evakey == nil {

				elOut.Resize(eval.params, 2)
				c2 = elOut.value[2]

				// If there is however an evaluation key, then
				// we still use the mempool for the third element
			} else {

				c2 = eval.ringpool[4]
			}

			// If the receiver Ciphertext either one of the inputs,
			// then it makes use of the mempool for the three elements
		} else {

			c0 = eval.ringpool[2]
			c1 = eval.ringpool[3]
			c2 = eval.ringpool[4]
		}

		context.MFormLvl(level, el0.value[0], c00)
		context.MFormLvl(level, el0.value[1], c01)

		if el0 == el1 { // squaring case

			context.MulCoeffsMontgomeryLvl(level, c00, el1.value[0], c0) // c0 = c[0]*c[0]
			context.MulCoeffsMontgomeryLvl(level, c00, el1.value[1], c1) // c1 = 2*c[0]*c[1]
			context.AddLvl(level, c1, c1, c1)
			context.MulCoeffsMontgomeryLvl(level, c01, el1.value[1], c2) // c2 = c[1]*c[1]

		} else { // regular case

			context.MulCoeffsMontgomeryLvl(level, c00, el1.value[0], c0) // c0 = c0[0]*c0[0]
			context.MulCoeffsMontgomeryLvl(level, c00, el1.value[1], c1)
			context.MulCoeffsMontgomeryAndAddNoModLvl(level, c01, el1.value[0], c1) // c1 = c0[0]*c1[1] + c0[1]*c1[0]
			context.MulCoeffsMontgomeryLvl(level, c01, el1.value[1], c2)            // c2 = c0[1]*c1[1]
		}

		// Relinearize if a key was provided
		if evakey != nil {

			context.CopyLvl(level, c0, elOut.value[0])
			context.CopyLvl(level, c1, elOut.value[1])

			eval.switchKeysInPlace(c2, evakey.evakey, elOut.Ciphertext())

		} else { // Or copy the result on the output Ciphertext if it was one of the inputs
			if elOut == el0 || elOut == el1 {
				elOut.Resize(eval.params, 2)
				context.CopyLvl(level, c0, elOut.value[0])
				context.CopyLvl(level, c1, elOut.value[1])
				context.CopyLvl(level, c2, elOut.value[2])
			}
		}

		// Case Plaintext (x) Ciphertext or Ciphertext (x) Plaintext
	} else {

		var tmp0, tmp1 *ckksElement

		if el0.Degree() == 1 {
			tmp0, tmp1 = el1, el0
		} else {
			tmp0, tmp1 = el0, el1
		}

		c00 := eval.ringpool[0]
		c00.Zero()

		context.MFormLvl(level, tmp0.value[0], c00)
		context.MulCoeffsMontgomeryLvl(level, c00, tmp1.value[0], elOut.value[0])
		context.MulCoeffsMontgomeryLvl(level, c00, tmp1.value[1], elOut.value[1])
	}
}

// RelinearizeNew applies the relinearization procedure on ct0 and returns the result in a newly
// created Ciphertext. The input Ciphertext must be of degree two.
func (eval *evaluator) RelinearizeNew(ct0 *Ciphertext, evakey *EvaluationKey) (ctOut *Ciphertext) {
	ctOut = NewCiphertext(eval.params, 1, ct0.Level(), ct0.Scale())
	eval.Relinearize(ct0, evakey, ctOut)
	return
}

// Relinearize applies the relinearization procedure on ct0 and returns the result in ctOut. The input Ciphertext must be of degree two.
func (eval *evaluator) Relinearize(ct0 *Ciphertext, evakey *EvaluationKey, ctOut *Ciphertext) {
	if ct0.Degree() != 2 {
		panic("cannot Relinearize: input Ciphertext is not of degree 2")
	}

	if ctOut != ct0 {
		ctOut.SetScale(ct0.Scale())
	}

	level := utils.MinUint64(ct0.Level(), ctOut.Level())
	context := eval.ckksContext.contextQ

	context.CopyLvl(level, ct0.value[0], ctOut.value[0])
	context.CopyLvl(level, ct0.value[1], ctOut.value[1])

	eval.switchKeysInPlace(ct0.value[2], evakey.evakey, ctOut)

	ctOut.Resize(eval.params, 1)
}

// SwitchKeysNew re-encrypts ct0 under a different key and returns the result in a newly created element.
// It requires a SwitchingKey, which is computed from the key under which the Ciphertext is currently encrypted,
// and the key under which the Ciphertext will be re-encrypted.
func (eval *evaluator) SwitchKeysNew(ct0 *Ciphertext, switchingKey *SwitchingKey) (ctOut *Ciphertext) {
	ctOut = NewCiphertext(eval.params, ct0.Degree(), ct0.Level(), ct0.Scale())
	eval.SwitchKeys(ct0, switchingKey, ctOut)
	return
}

// SwitchKeys re-encrypts ct0 under a different key and returns the result in ctOut.
// It requires a SwitchingKey, which is computed from the key under which the Ciphertext is currently encrypted,
// and the key under which the Ciphertext will be re-encrypted.
func (eval *evaluator) SwitchKeys(ct0 *Ciphertext, switchingKey *SwitchingKey, ctOut *Ciphertext) {

	if ct0.Degree() != 1 || ctOut.Degree() != 1 {
		panic("cannot SwitchKeys: input and output Ciphertext must be of degree 1")
	}

	level := utils.MinUint64(ct0.Level(), ctOut.Level())
	context := eval.ckksContext.contextQ

	context.CopyLvl(level, ct0.value[0], ctOut.value[0])
	context.CopyLvl(level, ct0.value[1], ctOut.value[1])

	eval.switchKeysInPlace(ct0.value[1], switchingKey, ctOut)
}

// RotateColumnsNew rotates the columns of ct0 by k positions to the left, and returns the result in a newly created element.
// If the provided element is a Ciphertext, a key-switching operation is necessary and a rotation key for the specific rotation needs to be provided.
func (eval *evaluator) RotateColumnsNew(ct0 *Ciphertext, k uint64, evakey *RotationKeys) (ctOut *Ciphertext) {
	ctOut = NewCiphertext(eval.params, ct0.Degree(), ct0.Level(), ct0.Scale())
	eval.RotateColumns(ct0, k, evakey, ctOut)
	return
}

// RotateColumns rotates the columns of ct0 by k positions to the left and returns the result in ctOut.
// If the provided element is a Ciphertext, a key-switching operation is necessary and a rotation key for the specific rotation needs to be provided.
func (eval *evaluator) RotateColumns(ct0 *Ciphertext, k uint64, evakey *RotationKeys, ctOut *Ciphertext) {

	if ct0.Degree() != 1 || ctOut.Degree() != 1 {
		panic("cannot RotateColumns: input and output Ciphertext must be of degree 1")
	}

	k &= ((eval.ckksContext.n >> 1) - 1)

	if k == 0 {

		ctOut.Copy(ct0.Element())

	} else {

		ctOut.SetScale(ct0.Scale())

		// It checks in the RotationKeys if the corresponding rotation has been generated
		if evakey.evakeyRotColLeft[k] != nil {

			eval.permuteNTT(ct0, evakey.permuteNTTLeftIndex[k], evakey.evakeyRotColLeft[k], ctOut)

		} else {

			// If not, it checks if the left and right pow2 rotations have been generated
			hasPow2Rotations := true
			for i := uint64(1); i < eval.ckksContext.n>>1; i <<= 1 {
				if evakey.evakeyRotColLeft[i] == nil || evakey.evakeyRotColRight[i] == nil {
					hasPow2Rotations = false
					break
				}
			}

			// If yes, it computes the least amount of rotation between left and right required to apply the demanded rotation
			if hasPow2Rotations {

				if utils.HammingWeight64(k) <= utils.HammingWeight64((eval.ckksContext.n>>1)-k) {
					eval.rotateColumnsLPow2(ct0, k, evakey, ctOut)
				} else {
					eval.rotateColumnsRPow2(ct0, (eval.ckksContext.n>>1)-k, evakey, ctOut)
				}

				// Otherwise, it returns an error indicating that the keys have not been generated
			} else {
				panic("cannot RotateColumns: specific rotation and pow2 rotations have not been generated")
			}
		}
	}
}

// RotateHoisted takes an input Ciphertext and a list of rotations and returns a map of Ciphertext, where each element of the map is the input Ciphertext
// rotation by one element of the list. It is much faster than sequential calls to RotateColumns.
func (eval *evaluator) RotateHoisted(ct0 *Ciphertext, rotations []uint64, rotkeys *RotationKeys) (cOut map[uint64]*Ciphertext) {

	// Pre-computation for rotations using hoisting
	contextQ := eval.ckksContext.contextQ
	contextP := eval.ckksContext.contextP

	c2NTT := ct0.value[1]
	c2InvNTT := contextQ.NewPoly()
	contextQ.InvNTTLvl(ct0.Level(), c2NTT, c2InvNTT)

	alpha := eval.params.Alpha()
	beta := uint64(math.Ceil(float64(ct0.Level()+1) / float64(alpha)))

	c2QiQDecomp := make([]*ring.Poly, beta)
	c2QiPDecomp := make([]*ring.Poly, beta)

	for i := uint64(0); i < beta; i++ {
		c2QiQDecomp[i] = contextQ.NewPoly()
		c2QiPDecomp[i] = contextP.NewPoly()
		eval.decomposeAndSplitNTT(ct0.Level(), i, c2NTT, c2InvNTT, c2QiQDecomp[i], c2QiPDecomp[i])
	}

	cOut = make(map[uint64]*Ciphertext)

	for _, i := range rotations {

		i &= ((eval.ckksContext.n >> 1) - 1)

		if i == 0 {
			cOut[i] = ct0.CopyNew().Ciphertext()
		} else {
			cOut[i] = NewCiphertext(eval.params, 1, ct0.Level(), ct0.Scale())
			eval.switchKeyHoisted(ct0, c2QiQDecomp, c2QiPDecomp, i, rotkeys, cOut[i])
		}
	}

	return
}

func (eval *evaluator) switchKeyHoisted(ct0 *Ciphertext, c2QiQDecomp, c2QiPDecomp []*ring.Poly, k uint64, evakey *RotationKeys, ctOut *Ciphertext) {

	if ct0.Degree() != 1 || ctOut.Degree() != 1 {
		panic("cannot switchKeyHoisted: input and output Ciphertext must be of degree 1")
	}

	k &= (1 << (eval.ckksContext.logN - 1)) - 1

	if evakey.permuteNTTLeftIndex[k] == nil {
		panic("cannot switchKeyHoisted: specific rotation has not been generated")
	}

	ctOut.SetScale(ct0.Scale())

	var level, reduce uint64

	level = ctOut.Level()

	contextQ := eval.ckksContext.contextQ
	contextP := eval.ckksContext.contextP

	if ct0 != ctOut {
		ring.PermuteNTTWithIndex(ct0.value[0], evakey.permuteNTTLeftIndex[k], eval.ringpool[0])
		contextQ.CopyLvl(level, eval.ringpool[0], ctOut.value[0])
		ring.PermuteNTTWithIndex(ct0.value[1], evakey.permuteNTTLeftIndex[k], eval.ringpool[0])
		contextQ.CopyLvl(level, eval.ringpool[0], ctOut.value[1])
	} else {
		ring.PermuteNTTWithIndex(ct0.value[0], evakey.permuteNTTLeftIndex[k], ctOut.value[0])
		ring.PermuteNTTWithIndex(ct0.value[1], evakey.permuteNTTLeftIndex[k], ctOut.value[1])
	}

	for i := range eval.poolQ {
		eval.poolQ[i].Zero()
	}

	for i := range eval.poolP {
		eval.poolP[i].Zero()
	}

	c2QiQPermute := eval.poolQ[0]
	c2QiPPermute := eval.poolP[0]

	pool2Q := eval.poolQ[1]
	pool2P := eval.poolP[1]

	pool3Q := eval.poolQ[2]
	pool3P := eval.poolP[2]

	reduce = 0

	alpha := eval.params.Alpha()
	beta := uint64(math.Ceil(float64(level+1) / float64(alpha)))

	// Key switching with CRT decomposition for the Qi
	for i := uint64(0); i < beta; i++ {

		ring.PermuteNTTWithIndex(c2QiQDecomp[i], evakey.permuteNTTLeftIndex[k], c2QiQPermute)
		ring.PermuteNTTWithIndex(c2QiPDecomp[i], evakey.permuteNTTLeftIndex[k], c2QiPPermute)

		contextQ.MulCoeffsMontgomeryAndAddNoModLvl(level, evakey.evakeyRotColLeft[k].evakey[i][0], c2QiQPermute, pool2Q)
		contextQ.MulCoeffsMontgomeryAndAddNoModLvl(level, evakey.evakeyRotColLeft[k].evakey[i][1], c2QiQPermute, pool3Q)

		// We continue with the key-switch primes.
		for j, keysindex := uint64(0), eval.ckksContext.levels; j < uint64(len(contextP.Modulus)); j, keysindex = j+1, keysindex+1 {

			pj := contextP.Modulus[j]
			mredParams := contextP.GetMredParams()[j]

			key0 := evakey.evakeyRotColLeft[k].evakey[i][0].Coeffs[keysindex]
			key1 := evakey.evakeyRotColLeft[k].evakey[i][1].Coeffs[keysindex]
			p2tmp := pool2P.Coeffs[j]
			p3tmp := pool3P.Coeffs[j]
			c2tmp := c2QiPPermute.Coeffs[j]

			for y := uint64(0); y < contextP.N; y++ {
				p2tmp[y] += ring.MRed(key0[y], c2tmp[y], pj, mredParams)
				p3tmp[y] += ring.MRed(key1[y], c2tmp[y], pj, mredParams)
			}
		}

		if reduce&7 == 1 {
			contextQ.ReduceLvl(level, pool2Q, pool2Q)
			contextQ.ReduceLvl(level, pool3Q, pool3Q)
			contextP.Reduce(pool2P, pool2P)
			contextP.Reduce(pool3P, pool3P)
		}

		reduce++
	}

	if (reduce-1)&7 != 1 {
		contextQ.ReduceLvl(level, pool2Q, pool2Q)
		contextQ.ReduceLvl(level, pool3Q, pool3Q)
		contextP.Reduce(pool2P, pool2P)
		contextP.Reduce(pool3P, pool3P)
	}

	// Independent of context (parameter: level)
	// Computes pool2Q = pool2Q/pool2P and pool3Q = pool3Q/pool3P
	eval.baseconverter.ModDownSplitedNTTPQ(level, pool2Q, pool2P, pool2Q)
	eval.baseconverter.ModDownSplitedNTTPQ(level, pool3Q, pool3P, pool3Q)

	// Independent of context (parameter: level)
	contextQ.AddLvl(level, ctOut.value[0], pool2Q, ctOut.value[0])
	contextQ.AddLvl(level, ctOut.value[1], pool3Q, ctOut.value[1])
}

func (eval *evaluator) rotateColumnsLPow2(ct0 *Ciphertext, k uint64, evakey *RotationKeys, ctOut *Ciphertext) {
	eval.rotateColumnsPow2(ct0, k, evakey.permuteNTTLeftIndex, evakey.evakeyRotColLeft, ctOut)
}

func (eval *evaluator) rotateColumnsRPow2(ct0 *Ciphertext, k uint64, evakey *RotationKeys, ctOut *Ciphertext) {
	eval.rotateColumnsPow2(ct0, k, evakey.permuteNTTRightIndex, evakey.evakeyRotColRight, ctOut)
}

func (eval *evaluator) rotateColumnsPow2(ct0 *Ciphertext, k uint64, permuteNTTIndex map[uint64][]uint64, evakeyRotCol map[uint64]*SwitchingKey, ctOut *Ciphertext) {

	var evakeyIndex uint64

	evakeyIndex = 1

	level := utils.MinUint64(ct0.Level(), ctOut.Level())
	context := eval.ckksContext.contextQ

	context.CopyLvl(level, ct0.value[0], ctOut.value[0])
	context.CopyLvl(level, ct0.value[1], ctOut.value[1])

	for k > 0 {

		if k&1 == 1 {

			eval.permuteNTT(ctOut, permuteNTTIndex[evakeyIndex], evakeyRotCol[evakeyIndex], ctOut)
		}

		evakeyIndex <<= 1
		k >>= 1
	}
}

// ConjugateNew conjugates ct0 (which is equivalent to a row rotation) and returns the result in a newly
// created element. If the provided element is a Ciphertext, a key-switching operation is necessary and a rotation key
// for the row rotation needs to be provided.
func (eval *evaluator) ConjugateNew(ct0 *Ciphertext, evakey *RotationKeys) (ctOut *Ciphertext) {
	ctOut = NewCiphertext(eval.params, ct0.Degree(), ct0.Level(), ct0.Scale())
	eval.Conjugate(ct0, evakey, ctOut)
	return
}

// Conjugate conjugates ct0 (which is equivalent to a row rotation) and returns the result in ctOut.
// If the provided element is a Ciphertext, a key-switching operation is necessary and a rotation key for the row rotation needs to be provided.
func (eval *evaluator) Conjugate(ct0 *Ciphertext, evakey *RotationKeys, ctOut *Ciphertext) {

	if ct0.Degree() != 1 || ctOut.Degree() != 1 {
		panic("cannot Conjugate: input and output Ciphertext must be of degree 1")
	}

	if evakey.evakeyConjugate == nil {
		panic("cannot Conjugate: rows rotation key not generated")
	}

	ctOut.SetScale(ct0.Scale())

	eval.permuteNTT(ct0, evakey.permuteNTTConjugateIndex, evakey.evakeyConjugate, ctOut)
}

func (eval *evaluator) permuteNTT(ct0 *Ciphertext, index []uint64, evakey *SwitchingKey, ctOut *Ciphertext) {

	var el0, el1 *ring.Poly

	if ct0 != ctOut {
		el0, el1 = ctOut.value[0], ctOut.value[1]
	} else {
		el0, el1 = eval.ringpool[0], eval.ringpool[1]
	}

	ring.PermuteNTTWithIndex(ct0.value[0], index, el0)
	ring.PermuteNTTWithIndex(ct0.value[1], index, el1)

	level := utils.MinUint64(ct0.Level(), ctOut.Level())
	context := eval.ckksContext.contextQ

	context.CopyLvl(level, el0, ctOut.value[0])
	context.CopyLvl(level, el1, ctOut.value[1])

	eval.switchKeysInPlace(ctOut.value[1], evakey, ctOut)
}

// switchKeysInPlace applies the general key-switching procedure of the form [c0 + cx*evakey[0], c1 + cx*evakey[1]]
func (eval *evaluator) switchKeysInPlace(cx *ring.Poly, evakey *SwitchingKey, ctOut *Ciphertext) {
	var level, reduce uint64

	level = ctOut.Level()

	contextQ := eval.ckksContext.contextQ
	contextP := eval.ckksContext.contextP

	for i := range eval.poolQ {
		eval.poolQ[i].Zero()
	}

	for i := range eval.poolP {
		eval.poolP[i].Zero()
	}

	c2QiQ := eval.poolQ[0]
	c2QiP := eval.poolP[0]

	pool2Q := eval.poolQ[1]
	pool2P := eval.poolP[1]

	pool3Q := eval.poolQ[2]
	pool3P := eval.poolP[2]

	c2 := eval.poolQ[3]

	// We switch the element on which the switching key operation will be conducted out of the NTT domain

	contextQ.InvNTTLvl(level, cx, c2)

	reduce = 0

	alpha := eval.params.Alpha()
	beta := uint64(math.Ceil(float64(level+1) / float64(alpha)))

	// Key switching with CRT decomposition for the Qi
	for i := uint64(0); i < beta; i++ {

		eval.decomposeAndSplitNTT(level, i, cx, c2, c2QiQ, c2QiP)

		contextQ.MulCoeffsMontgomeryAndAddNoModLvl(level, evakey.evakey[i][0], c2QiQ, pool2Q)
		contextQ.MulCoeffsMontgomeryAndAddNoModLvl(level, evakey.evakey[i][1], c2QiQ, pool3Q)

		// We continue with the key-switch primes.
		for j, keysindex := uint64(0), eval.ckksContext.levels; j < uint64(len(contextP.Modulus)); j, keysindex = j+1, keysindex+1 {

			pj := contextP.Modulus[j]
			mredParams := contextP.GetMredParams()[j]

			key0 := evakey.evakey[i][0].Coeffs[keysindex]
			key1 := evakey.evakey[i][1].Coeffs[keysindex]
			c2tmp := c2QiP.Coeffs[j]
			p2tmp := pool2P.Coeffs[j]
			p3tmp := pool3P.Coeffs[j]

			for y := uint64(0); y < contextP.N; y++ {
				p2tmp[y] += ring.MRed(key0[y], c2tmp[y], pj, mredParams)
				p3tmp[y] += ring.MRed(key1[y], c2tmp[y], pj, mredParams)
			}
		}

		if reduce&7 == 1 {
			contextQ.ReduceLvl(level, pool2Q, pool2Q)
			contextQ.ReduceLvl(level, pool3Q, pool3Q)
			contextP.Reduce(pool2P, pool2P)
			contextP.Reduce(pool3P, pool3P)
		}

		reduce++
	}

	//Independent of context (parameter: level)
	if (reduce-1)&7 != 1 {
		contextQ.ReduceLvl(level, pool2Q, pool2Q)
		contextQ.ReduceLvl(level, pool3Q, pool3Q)
		contextP.Reduce(pool2P, pool2P)
		contextP.Reduce(pool3P, pool3P)
	}

	//Independent of context (parameter: level)
	// Computes pool2Q = pool2Q/pool2P and pool3Q = pool3Q/pool3P
	eval.baseconverter.ModDownSplitedNTTPQ(level, pool2Q, pool2P, pool2Q)
	eval.baseconverter.ModDownSplitedNTTPQ(level, pool3Q, pool3P, pool3Q)

	//Independent of context (parameter: level)
	contextQ.AddLvl(level, ctOut.value[0], pool2Q, ctOut.value[0])
	contextQ.AddLvl(level, ctOut.value[1], pool3Q, ctOut.value[1])
}

// decomposeAndSplitNTT decomposes the input polynomial into the target CRT basis.
func (eval *evaluator) decomposeAndSplitNTT(level, beta uint64, c2NTT, c2InvNTT, c2QiQ, c2QiP *ring.Poly) {

	contextQ := eval.ckksContext.contextQ
	contextP := eval.ckksContext.contextP

	eval.decomposer.DecomposeAndSplit(level, beta, c2InvNTT, c2QiQ, c2QiP)

	p0idxst := beta * eval.params.Alpha()
	p0idxed := p0idxst + eval.decomposer.Xalpha()[beta]

	// c2_qi = cx mod qi mod qi
	for x := uint64(0); x < level+1; x++ {

		qi := contextQ.Modulus[x]
		nttPsi := contextQ.GetNttPsi()[x]
		bredParams := contextQ.GetBredParams()[x]
		mredParams := contextQ.GetMredParams()[x]

		if p0idxst <= x && x < p0idxed {
			p0tmp := c2NTT.Coeffs[x]
			p1tmp := c2QiQ.Coeffs[x]
			for j := uint64(0); j < contextQ.N; j++ {
				p1tmp[j] = p0tmp[j]
			}
		} else {
			ring.NTT(c2QiQ.Coeffs[x], c2QiQ.Coeffs[x], contextQ.N, nttPsi, qi, mredParams, bredParams)
		}
	}
	// c2QiP = c2 mod qi mod pj
	contextP.NTT(c2QiP, c2QiP)
}
