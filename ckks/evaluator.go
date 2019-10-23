package ckks

import (
	"errors"
	"github.com/ldsec/lattigo/ring"
	"math"
)

// Evaluator is a struct holding the necessary elements to operates the homomorphic operations between ciphertext and/or plaintexts.
// It also holds a small memory pool used to store intermediate computations.
type Evaluator struct {
	ckkscontext   *CkksContext
	ringpool      [6]*ring.Poly
	keyswitchpool [4]*ring.Poly

	poolQ [3]*ring.Poly
	poolP [3]*ring.Poly

	ctxpool     *Ciphertext
	rescalepool []uint64

	baseconverter *ring.FastBasisExtender
	decomposer    *ring.ArbitraryDecomposer
}

// NewEvaluator creates a new Evaluator, that can be used to do homomorphic
// operations on the ciphertexts and/or plaintexts. It stores a small pool of polynomials
// and ciphertexts that will be used for intermediate values.
func (ckkscontext *CkksContext) NewEvaluator() (evaluator *Evaluator) {

	evaluator = new(Evaluator)
	evaluator.ckkscontext = ckkscontext

	for i := 0; i < 5; i++ {
		evaluator.ringpool[i] = evaluator.ckkscontext.contextQ.NewPoly()
	}

	for i := 0; i < 4; i++ {
		evaluator.keyswitchpool[i] = evaluator.ckkscontext.contextKeys.NewPoly()
	}

	for i := 0 ; i < 3 ; i++{
		evaluator.poolQ[i] = evaluator.ckkscontext.contextQ.NewPoly()
		evaluator.poolP[i] = evaluator.ckkscontext.contextP.NewPoly()
	}

	evaluator.rescalepool = make([]uint64, ckkscontext.n)

	evaluator.ctxpool = ckkscontext.NewCiphertext(1, ckkscontext.levels-1, ckkscontext.scale)

	evaluator.baseconverter = ring.NewFastBasisExtender(evaluator.ckkscontext.contextLevel[ckkscontext.levels-1].Modulus, evaluator.ckkscontext.specialprimes)

	evaluator.decomposer = ring.NewArbitraryDecomposer(ckkscontext.moduli, ckkscontext.specialprimes)

	return evaluator
}

func (evaluator *Evaluator) getElemAndCheckBinary(op0, op1, opOut Operand, opOutMinDegree uint64) (el0, el1, elOut *ckksElement, err error) {
	if op0 == nil || op1 == nil || opOut == nil {
		return nil, nil, nil, errors.New("operands cannot be nil")
	}

	if op0.Degree()+op1.Degree() == 0 {
		return nil, nil, nil, errors.New("operands cannot be both plaintext")
	}

	if opOut.Degree() < opOutMinDegree {
		return nil, nil, nil, errors.New("receiver operand degree is too small")
	}
	el0, el1, elOut = op0.Element(), op1.Element(), opOut.Element()
	return // TODO: more checks on elements
}

func (evaluator *Evaluator) getElemAndCheckUnary(op0, opOut Operand, opOutMinDegree uint64) (el0, elOut *ckksElement, err error) {
	if op0 == nil || opOut == nil {
		return nil, nil, errors.New("operand cannot be nil")
	}

	if op0.Degree() == 0 {
		return nil, nil, errors.New("operand cannot be plaintext")
	}

	if opOut.Degree() < opOutMinDegree {
		return nil, nil, errors.New("receiver operand degree is too small")
	}
	el0, elOut = op0.Element(), opOut.Element()
	return // TODO: more checks on elements
}

func (evaluator *Evaluator) newCiphertextBinary(op0, op1 Operand) (ctOut *Ciphertext) {

	maxDegree := max([]uint64{op0.Degree(), op1.Degree()})
	maxScale := maxflo([]float64{op0.Scale(), op1.Scale()})
	minLevel := min([]uint64{op0.Level(), op1.Level()})

	return evaluator.ckkscontext.NewCiphertext(maxDegree, minLevel, maxScale)
}

// Add adds op0 to op1 and returns the result on ctOut.
func (evaluator *Evaluator) Add(op0, op1 Operand, ctOut *Ciphertext) (err error) {

	el0, el1, elOut, err := evaluator.getElemAndCheckBinary(op0, op1, ctOut, max([]uint64{op0.Degree(), op1.Degree()}))
	if err != nil {
		return err
	}

	return evaluator.evaluateInPlace(el0, el1, elOut, evaluator.ckkscontext.contextLevel[min([]uint64{el0.Level(), el1.Level(), elOut.Level()})].Add)
}

// AddNoMod adds op0 to op1 and returns the result on ctOut, without modular reduction.
func (evaluator *Evaluator) AddNoMod(op0, op1 Operand, ctOut *Ciphertext) (err error) {

	el0, el1, elOut, err := evaluator.getElemAndCheckBinary(op0, op1, ctOut, max([]uint64{op0.Degree(), op1.Degree()}))
	if err != nil {
		return err
	}

	return evaluator.evaluateInPlace(el0, el1, elOut, evaluator.ckkscontext.contextLevel[min([]uint64{el0.Level(), el1.Level(), elOut.Level()})].AddNoMod)
}

// Add adds op0 to op1 and returns the result on a newly created element.
func (evaluator *Evaluator) AddNew(op0, op1 Operand) (ctOut *Ciphertext, err error) {

	ctOut = evaluator.newCiphertextBinary(op0, op1)

	return ctOut, evaluator.Add(op0, op1, ctOut)
}

// Add adds op0 to op1 without modular reduction, and returns the result on a newly created element.
func (evaluator *Evaluator) AddNoModNew(op0, op1 Operand) (ctOut *Ciphertext, err error) {

	ctOut = evaluator.newCiphertextBinary(op0, op1)

	return ctOut, evaluator.AddNoMod(op0, op1, ctOut)
}

// Sub subtracts op0 to op1 and returns the result on ctOut.
func (evaluator *Evaluator) Sub(op0, op1 Operand, ctOut *Ciphertext) (err error) {

	el0, el1, elOut, err := evaluator.getElemAndCheckBinary(op0, op1, ctOut, max([]uint64{op0.Degree(), op1.Degree()}))
	if err != nil {
		return err
	}

	minLevel := min([]uint64{el0.Level(), el1.Level(), elOut.Level()})

	if err = evaluator.evaluateInPlace(el0, el1, elOut, evaluator.ckkscontext.contextLevel[minLevel].Sub); err != nil {
		return err
	}

	if el0.Degree() < el1.Degree() {
		for i := el0.Degree() + 1; i < el1.Degree()+1; i++ {
			evaluator.ckkscontext.contextLevel[minLevel].Neg(elOut.Value()[i], elOut.Value()[i])
		}
	}

	return nil
}

// SubNoMod subtracts op0 to op1 and returns the result on ctOut, without modular reduction.
func (evaluator *Evaluator) SubNoMod(op0, op1 Operand, ctOut *Ciphertext) (err error) {

	el0, el1, elOut, err := evaluator.getElemAndCheckBinary(op0, op1, ctOut, max([]uint64{op0.Degree(), op1.Degree()}))
	if err != nil {
		return err
	}

	minLevel := min([]uint64{el0.Level(), el1.Level(), elOut.Level()})

	if err = evaluator.evaluateInPlace(el0, el1, elOut, evaluator.ckkscontext.contextLevel[minLevel].SubNoMod); err != nil {
		return err
	}

	if el0.Degree() < el1.Degree() {
		for i := el0.Degree() + 1; i < el1.Degree()+1; i++ {
			evaluator.ckkscontext.contextLevel[minLevel].Neg(elOut.Value()[i], elOut.Value()[i])
		}
	}

	return nil
}

// SubNew subtracts op0 to op1 and returns the result on a newly created element.
func (evaluator *Evaluator) SubNew(op0, op1 Operand) (ctOut *Ciphertext, err error) {

	ctOut = evaluator.newCiphertextBinary(op0, op1)

	return ctOut, evaluator.Sub(op0, op1, ctOut)
}

// SubNoModNew subtracts op0 to op1 without modular reduction, and returns the result on a newly created element.
func (evaluator *Evaluator) SubNoModNew(op0, op1 Operand) (ctOut *Ciphertext, err error) {

	ctOut = evaluator.newCiphertextBinary(op0, op1)

	return ctOut, evaluator.SubNoMod(op0, op1, ctOut)
}

func (evaluator *Evaluator) evaluateInPlace(c0, c1, ctOut *ckksElement, evaluate func(*ring.Poly, *ring.Poly, *ring.Poly)) (err error) {

	var tmp0, tmp1 *ckksElement // TODO : use evaluator mem pool

	maxDegree := max([]uint64{c0.Degree(), c1.Degree()})
	minDegree := min([]uint64{c0.Degree(), c1.Degree()})

	// Else resizes the receiver element
	ctOut.Resize(evaluator.ckkscontext, maxDegree)
	evaluator.DropLevel(ctOut, ctOut.Level()-min([]uint64{c0.Level(), c1.Level()}))

	// Checks wether or not the receiver element is the same as one of the input elements
	// and acts accordingly to avoid unnecessary element creation or element overwriting,
	// and scales properly the element before the evaluation.
	if ctOut == c0 {

		if c0.Scale() > c1.Scale() {

			tmp1 = evaluator.ctxpool.Element()

			if uint64(c0.Scale()/c1.Scale()) != 0 {
				if err = evaluator.MultConst(c1.Ciphertext(), uint64(c0.Scale()/c1.Scale()), tmp1.Ciphertext()); err != nil {
					return err
				}
			}

		} else if c1.Scale() > c0.Scale() {

			if uint64(c1.Scale()/c0.Scale()) != 0 {
				evaluator.MultConst(c0.Ciphertext(), uint64(c1.Scale()/c0.Scale()), c0.Ciphertext())
			}

			c0.SetScale(c1.Scale())

			tmp1 = c1

		} else {

			tmp1 = c1
		}

		tmp0 = c0

	} else if ctOut == c1 {

		if c1.Scale() > c0.Scale() {

			tmp0 = evaluator.ctxpool.Element()
			if uint64(c1.Scale()/c0.Scale()) != 0 {
				if err = evaluator.MultConst(c0.Ciphertext(), uint64(c1.Scale()/c0.Scale()), tmp0.Ciphertext()); err != nil {
					return err
				}
			}

		} else if c0.Scale() > c1.Scale() {

			if uint64(c0.Scale()/c1.Scale()) != 0 {
				evaluator.MultConst(c1.Ciphertext(), uint64(c0.Scale()/c1.Scale()), ctOut.Ciphertext())
			}

			ctOut.SetScale(c0.Scale())

			tmp0 = c0

		} else {

			tmp0 = c0
		}

		tmp1 = c1

	} else {

		if c1.Scale() > c0.Scale() {

			tmp0 = evaluator.ctxpool.Element()

			if uint64(c1.Scale()/c0.Scale()) != 0 {
				if err = evaluator.MultConst(c0.Ciphertext(), uint64(c1.Scale()/c0.Scale()), tmp0.Ciphertext()); err != nil {
					return err
				}
			}

			tmp1 = c1

		} else if c0.Scale() > c1.Scale() {

			tmp1 = evaluator.ctxpool.Element()

			if uint64(c0.Scale()/c1.Scale()) != 0 {
				if err = evaluator.MultConst(c1.Ciphertext(), uint64(c0.Scale()/c1.Scale()), tmp1.Ciphertext()); err != nil {
					return err
				}
			}

			tmp0 = c0

		} else {
			tmp0 = c0
			tmp1 = c1
		}
	}

	for i := uint64(0); i < minDegree+1; i++ {
		evaluate(tmp0.Value()[i], tmp1.Value()[i], ctOut.Value()[i])
	}

	ctOut.SetScale(maxflo([]float64{c0.Scale(), c1.Scale()}))

	// If the inputs degree differ, copies the remaining degree on the receiver
	// Also checks that the receiver is ont one of the inputs to avoid unnecessary work.

	context := evaluator.ckkscontext.contextLevel[ctOut.Level()]

	if c0.Degree() > c1.Degree() && tmp0 != ctOut {
		for i := minDegree + 1; i < maxDegree+1; i++ {
			context.Copy(tmp0.Value()[i], ctOut.Value()[i])
		}
	} else if c1.Degree() > c0.Degree() && tmp1 != ctOut {
		for i := minDegree + 1; i < maxDegree+1; i++ {
			context.Copy(tmp1.Value()[i], ctOut.Value()[i])
		}
	}

	return nil
}

// Neg negates the ct0 and returns the result on ctOut.
func (evaluator *Evaluator) Neg(ct0 *Ciphertext, ctOut *Ciphertext) (err error) {

	minLevel := min([]uint64{ct0.Level(), ctOut.Level()})

	if ct0.Degree() != ctOut.Degree() {
		return errors.New("cannot negate -> invalid receiver ciphertext does not match input ciphertext degree")
	}

	for i := range ct0.value {
		evaluator.ckkscontext.contextLevel[minLevel].Neg(ct0.value[i], ctOut.Value()[i])
	}

	return nil
}

// Neg negates ct0 and returns the result on a newly created element.
func (evaluator *Evaluator) NegNew(ct0 *Ciphertext) (ctOut *Ciphertext) {

	ctOut = evaluator.ckkscontext.NewCiphertext(ct0.Degree(), ct0.Level(), ct0.Scale())

	_ = evaluator.Neg(ct0, ctOut)

	return
}

// ExtractImagNew sets the real part of ct0 to the imaginary part of ct0 and sets the imaginary part of ct0 to zero, and returns the result on a new element.
// ex. f(a + b*i) = b. Requires a rotationkey for which the conjugate key has been generated. Scale is increased by one.
func (evaluator *Evaluator) ExtractImagNew(ct0 *Ciphertext, evakey *RotationKey) (ctOut *Ciphertext, err error) {

	ctOut = evaluator.ckkscontext.NewCiphertext(ct0.Degree(), ct0.Level(), ct0.Scale())

	return ctOut, evaluator.ExtractImag(ct0, evakey, ctOut)
}

// ExtractImag sets the real part of ct0 to the imaginary part of ct0 and sets the imaginary part of ct0 to zero, and returns the result on ctOut.
// ex. f(a + b*i) = b. Requires a rotationkey for which the conjugate key has been generated. Scale is increased by one.
func (evaluator *Evaluator) ExtractImag(ct0 *Ciphertext, evakey *RotationKey, ctOut *Ciphertext) (err error) {

	if err = evaluator.Conjugate(ct0, evakey, evaluator.ctxpool); err != nil {
		return err
	}

	if err = evaluator.MultByi(evaluator.ctxpool, evaluator.ctxpool); err != nil {
		return err
	}

	if err = evaluator.DivByi(ct0, ctOut); err != nil {
		return err
	}

	if err = evaluator.Add(ctOut, evaluator.ctxpool, ctOut); err != nil {
		return err
	}

	ctOut.SetScale(ctOut.Scale() * 2)

	return nil
}

// SwapRealImagNew swaps the real and imaginary parts of ct0 and returns the result on a newly created element, ex.
// f(a + b*i) = b + a * i. Requires a rotationkey for which the conjugate key has been generated.
func (evaluator *Evaluator) SwapRealImagNew(ct0 *Ciphertext, evakey *RotationKey) (ctOut *Ciphertext, err error) {

	ctOut = evaluator.ckkscontext.NewCiphertext(ct0.Degree(), ct0.Level(), ct0.Scale())

	return ctOut, evaluator.SwapRealImag(ct0, evakey, ctOut)
}

// SwapRealImagNew swaps the real and imaginary parts of ct0 and returns the result on ctOut, ex.
// f(a + b*i) = b + a * i. Requires a rotationkey for which the conjugate key has been generated.
func (evaluator *Evaluator) SwapRealImag(ct0 *Ciphertext, evakey *RotationKey, ctOut *Ciphertext) (err error) {

	if err = evaluator.DivByi(ct0, ctOut); err != nil {
		return err
	}

	if err = evaluator.Conjugate(ctOut, evakey, ctOut); err != nil {
		return err
	}

	return nil
}

// RemoveRealNew sets the real part of ct0 to zero and returns the result on a newly created element, ex. f(a + b*i) = b*i.
// Requires a rotationkey for which the conjugate key has been generated. Scale is increased by one.
func (evaluator *Evaluator) RemoveRealNew(ct0 *Ciphertext, evakey *RotationKey) (ctOut *Ciphertext, err error) {

	ctOut = evaluator.ckkscontext.NewCiphertext(ct0.Degree(), ct0.Level(), ct0.Scale())

	return ctOut, evaluator.RemoveReal(ct0, evakey, ctOut)
}

// RemoveReal sets the real part of ct0 to zero and returns the result on ctOut, ex. f(a + b*i) = b*i.
// Requires a rotationkey for which the conjugate key has been generated. Scale is increased by one.
func (evaluator *Evaluator) RemoveReal(ct0 *Ciphertext, evakey *RotationKey, ctOut *Ciphertext) (err error) {

	if ct0 != ctOut {

		if err = evaluator.Conjugate(ct0, evakey, ctOut); err != nil {
			return err
		}
		if err = evaluator.Sub(ct0, ctOut, ctOut); err != nil {
			return err
		}

	} else {

		if err = evaluator.Conjugate(ct0, evakey, evaluator.ctxpool); err != nil {
			return err
		}

		if err = evaluator.Sub(ctOut, evaluator.ctxpool, ctOut); err != nil {
			return err
		}
	}

	ctOut.SetScale(ctOut.Scale() * 2)

	return nil
}

// RemoveImagNew sets the imaginary part of ct0 to zero and returns the result on a newly created element, ex. f(a + b*i) = a.
// Requires a rotationkey for which the conjugate key has been generated. Scale is increased by one.
func (evaluator *Evaluator) RemoveImagNew(ct0 *Ciphertext, evakey *RotationKey) (ctOut *Ciphertext, err error) {

	ctOut = evaluator.ckkscontext.NewCiphertext(ct0.Degree(), ct0.Level(), ct0.Scale())

	return ctOut, evaluator.RemoveImag(ct0, evakey, ctOut)
}

// RemoveImag sets the imaginary part of ct0 to zero and returns the result on ctOut, ex. f(a + b*i) = a.
// Requires a rotationkey for which the conjugate key has been generated. Scale is increased by one.
func (evaluator *Evaluator) RemoveImag(ct0 *Ciphertext, evakey *RotationKey, ctOut *Ciphertext) (err error) {

	if ct0 != ctOut {

		if err = evaluator.Conjugate(ct0, evakey, ctOut); err != nil {
			return err
		}
		if err = evaluator.Add(ct0, ctOut, ctOut); err != nil {
			return err
		}

	} else {

		if err = evaluator.Conjugate(ct0, evakey, evaluator.ctxpool); err != nil {
			return err
		}

		if err = evaluator.Add(evaluator.ctxpool, ctOut, ctOut); err != nil {
			return err
		}

	}

	ctOut.SetScale(ctOut.Scale() * 2)

	return nil
}

// AddConstNew adds the input constant (which can be an uint64, int64, float64 or complex128) to ct0 and returns the result on a new element.
func (evaluator *Evaluator) AddConstNew(ct0 *Ciphertext, constant interface{}) (ctOut *Ciphertext) {
	ctOut = ct0.CopyNew().Ciphertext()
	_ = evaluator.AddConst(ct0, constant, ctOut)
	return ctOut
}

// AddConstNew adds the input constant (which can be an uint64, int64, float64 or complex128) to ct0 and returns the result on ctOut.
func (evaluator *Evaluator) AddConst(ct0 *Ciphertext, constant interface{}, ctOut *Ciphertext) (err error) {

	var level uint64

	level = min([]uint64{ct0.Level(), ctOut.Level()})

	var c_real, c_imag float64

	switch constant.(type) {
	case complex128:
		c_real = real(constant.(complex128))
		c_imag = imag(constant.(complex128))

	case float64:
		c_real = constant.(float64)
		c_imag = float64(0)

	case uint64:
		c_real = float64(constant.(uint64))
		c_imag = float64(0)

	case int64:
		c_real = float64(constant.(int64))
		c_imag = float64(0)

	case int:
		c_real = float64(constant.(int))
		c_imag = float64(0)
	}

	var scaled_const, scaled_const_real, scaled_const_imag uint64

	context := evaluator.ckkscontext.contextLevel[level]

	// Component wise addition of the following vector to the ciphertext :
	// [a + b*psi_qi^2, ....., a + b*psi_qi^2, a - b*psi_qi^2, ...., a - b*psi_qi^2] mod Qi
	// [{                  N/2                }{                N/2               }]
	// Which is equivalent outside of the NTT domain of adding a to the first coefficient of ct0 and b to the N/2th coefficient of ct0.
	for i := uint64(0); i < level+1; i++ {
		scaled_const_real = 0
		scaled_const_imag = 0
		scaled_const = 0

		if c_real != 0 {
			scaled_const_real = scaleUpExact(c_real, ct0.Scale(), evaluator.ckkscontext.moduli[i])
			scaled_const = scaled_const_real
		}

		if c_imag != 0 {
			scaled_const_imag = ring.MRed(scaleUpExact(c_imag, ct0.Scale(), evaluator.ckkscontext.moduli[i]), context.GetNttPsi()[i][1], context.Modulus[i], context.GetMredParams()[i])
			scaled_const = ring.CRed(scaled_const+scaled_const_imag, context.Modulus[i])
		}

		for j := uint64(0); j < evaluator.ckkscontext.n>>1; j++ {
			ctOut.Value()[0].Coeffs[i][j] = ring.CRed(ct0.value[0].Coeffs[i][j]+scaled_const, evaluator.ckkscontext.moduli[i])
		}

		if c_imag != 0 {
			scaled_const = ring.CRed(scaled_const_real+(context.Modulus[i]-scaled_const_imag), context.Modulus[i])
		}

		for j := evaluator.ckkscontext.n >> 1; j < evaluator.ckkscontext.n; j++ {
			ctOut.Value()[0].Coeffs[i][j] = ring.CRed(ct0.value[0].Coeffs[i][j]+scaled_const, evaluator.ckkscontext.moduli[i])
		}
	}

	return nil
}

// MultByConstAndAdd multiplies ct0 by the input constant, and adds it to the receiver element (does not modify the input
// element), ex. ctOut(x) = ctOut(x) + ct0(x) * (a+bi). This functions  removes the need of storing the intermediate value c(x) * (a+bi).
// This function will modifie the level and the scale of the receiver element depending on the level and the scale of the input
// element and the type of the constant. The level of the receiver element will be set to min(input.level, receiver.level).
// The scale of the receiver element will be set to the scale that the input element would have after the multiplication by the constant.
func (evaluator *Evaluator) MultByConstAndAdd(ct0 *Ciphertext, constant interface{}, ctOut *Ciphertext) (err error) {

	var level uint64

	level = min([]uint64{ct0.Level(), ctOut.Level()})

	// Forces a drop of cOut level to ct0 level
	if ctOut.Level() > level {
		evaluator.DropLevel(ctOut.Element(), ctOut.Level()-level)
	}

	var c_real, c_imag float64
	var scale float64

	// Converts to float64 and determines if a scale is required (which is the case if either real or imag has a rational part)
	scale = 1
	switch constant.(type) {
	case complex128:
		c_real = real(constant.(complex128))
		c_imag = imag(constant.(complex128))

		if c_real != 0 {
			value_int := int64(c_real)
			value_float := c_real - float64(value_int)

			if value_float != 0 {
				scale = evaluator.ckkscontext.scale
			}
		}

		if c_imag != 0 {
			value_int := int64(c_imag)
			value_float := c_imag - float64(value_int)

			if value_float != 0 {
				scale = evaluator.ckkscontext.scale
			}
		}

	case float64:
		c_real = constant.(float64)
		c_imag = float64(0)

		if c_real != 0 {
			value_int := int64(c_real)
			value_float := c_real - float64(value_int)

			if value_float != 0 {
				scale = evaluator.ckkscontext.scale
			}
		}

	case uint64:
		c_real = float64(constant.(uint64))
		c_imag = float64(0)

	case int64:
		c_real = float64(constant.(int64))
		c_imag = float64(0)

	case int:
		c_real = float64(constant.(int))
		c_imag = float64(0)
	}

	var scaled_const, scaled_const_real, scaled_const_imag uint64

	context := evaluator.ckkscontext.contextLevel[level]

	// If a scaling will be required to multiply by the constant,
	// equalizes scales such that the scales match in the end.
	if scale != 1 {

		// If ctOut scaling is smaller than ct0's scale + the default scaling,
		// then brings ctOut scale to ct0's scale.
		if ctOut.Scale() < ct0.Scale()*evaluator.ckkscontext.scale {

			if uint64((evaluator.ckkscontext.scale*ct0.Scale())/ctOut.Scale()) != 0 {

				evaluator.MultConst(ctOut, uint64((evaluator.ckkscontext.scale*ct0.Scale())/ctOut.Scale()), ctOut)

			}

			ctOut.SetScale(evaluator.ckkscontext.scale * ct0.Scale())

			// If ctOut.Scale() > ((a+bi)*scale)*ct0(x), then sets the scale to
			// bring c(x)*scale to the level of ctOut(x) scale
		} else if ctOut.Scale() > ct0.Scale()*evaluator.ckkscontext.scale {
			scale = ctOut.Scale() / ct0.Scale()
		}

		// If no scaling is required, the sets the appropriate scale such that
		// ct0(x)*scale matches ctOut(x) scale without modifiying ct0(x) scale.
	} else {

		if ctOut.Scale() > ct0.Scale() {

			scale = ctOut.Scale() / ct0.Scale()

		} else if ct0.Scale() > ctOut.Scale() {

			if uint64(ct0.Scale()/ctOut.Scale()) != 0 {
				evaluator.MultConst(ctOut, ct0.Scale()/ctOut.Scale(), ctOut)
			}

			ctOut.SetScale(ct0.Scale())
		}
	}

	// Component wise multiplication of the following vector to the ciphertext :
	// [a + b*psi_qi^2, ....., a + b*psi_qi^2, a - b*psi_qi^2, ...., a - b*psi_qi^2] mod Qi
	// [{                  N/2                }{                N/2               }]
	// Which is equivalent outside of the NTT domain of adding a to the first coefficient of ct0 and b to the N/2th coefficient of ct0.
	for i := uint64(0); i < level+1; i++ {

		scaled_const_real = 0
		scaled_const_imag = 0
		scaled_const = 0

		if c_real != 0 {
			scaled_const_real = scaleUpExact(c_real, scale, evaluator.ckkscontext.moduli[i])
			scaled_const = scaled_const_real
		}

		if c_imag != 0 {
			scaled_const_imag = scaleUpExact(c_imag, scale, evaluator.ckkscontext.moduli[i])
			scaled_const_imag = ring.MRed(scaled_const_imag, context.GetNttPsi()[i][1], context.Modulus[i], context.GetMredParams()[i])
			scaled_const = ring.CRed(scaled_const+scaled_const_imag, context.Modulus[i])
		}

		scaled_const = ring.MForm(scaled_const, context.Modulus[i], context.GetBredParams()[i])

		for u := range ct0.Value() {
			for j := uint64(0); j < evaluator.ckkscontext.n>>1; j++ {
				ctOut.Value()[u].Coeffs[i][j] = ring.CRed(ctOut.Value()[u].Coeffs[i][j]+ring.MRed(ct0.Value()[u].Coeffs[i][j], scaled_const, context.Modulus[i], context.GetMredParams()[i]), context.Modulus[i])
			}
		}

		if c_imag != 0 {
			scaled_const = ring.CRed(scaled_const_real+(context.Modulus[i]-scaled_const_imag), context.Modulus[i])
			scaled_const = ring.MForm(scaled_const, context.Modulus[i], context.GetBredParams()[i])
		}

		for u := range ct0.Value() {
			for j := evaluator.ckkscontext.n >> 1; j < evaluator.ckkscontext.n; j++ {
				ctOut.Value()[u].Coeffs[i][j] = ring.CRed(ctOut.Value()[u].Coeffs[i][j]+ring.MRed(ct0.Value()[u].Coeffs[i][j], scaled_const, context.Modulus[i], context.GetMredParams()[i]), context.Modulus[i])
			}
		}
	}

	return nil

}

// MultConstNew multiplies ct0 by the input constant and returns the result on a newly created element.
// The scale of the output element will depend on the scale of the input element and the constant (if the constant
// needs to be scaled (its rational part is not zero)). The constant can be an uint64, int64, float64 or complex128.
func (evaluator *Evaluator) MultConstNew(ct0 *Ciphertext, constant interface{}) (ctOut *Ciphertext, err error) {

	ctOut = evaluator.ckkscontext.NewCiphertext(ct0.Degree(), ct0.Level(), ct0.Scale())

	return ctOut, evaluator.MultConst(ct0, constant, ctOut)
}

// MultConstNew multiplies ct0 by the input constant and returns the result on ctOut.
// The scale of the output element will depend on the scale of the input element and the constant (if the constant
// needs to be scaled (its rational part is not zero)). The constant can be an uint64, int64, float64 or complex128.
func (evaluator *Evaluator) MultConst(ct0 *Ciphertext, constant interface{}, ctOut *Ciphertext) (err error) {

	var level uint64

	level = min([]uint64{ct0.Level(), ctOut.Level()})

	var c_real, c_imag float64
	var scale float64

	// Converts to float64 and determines if a scale is required (which is the case if either real or imag has a rational part)
	scale = 1
	switch constant.(type) {
	case complex128:
		c_real = real(constant.(complex128))
		c_imag = imag(constant.(complex128))

		if c_real != 0 {
			value_int := int64(c_real)
			value_float := c_real - float64(value_int)

			if value_float != 0 {
				scale = evaluator.ckkscontext.scale
			}
		}

		if c_imag != 0 {
			value_int := int64(c_imag)
			value_float := c_imag - float64(value_int)

			if value_float != 0 {
				scale = evaluator.ckkscontext.scale
			}
		}

	case float64:
		c_real = constant.(float64)
		c_imag = float64(0)

		if c_real != 0 {
			value_int := int64(c_real)
			value_float := c_real - float64(value_int)

			if value_float != 0 {
				scale = evaluator.ckkscontext.scale
			}
		}

	case uint64:
		c_real = float64(constant.(uint64))
		c_imag = float64(0)

	case int64:
		c_real = float64(constant.(int64))
		c_imag = float64(0)

	case int:
		c_real = float64(constant.(int))
		c_imag = float64(0)
	}

	// Component wise multiplication of the following vector to the ciphertext :
	// [a + b*psi_qi^2, ....., a + b*psi_qi^2, a - b*psi_qi^2, ...., a - b*psi_qi^2] mod Qi
	// [{                  N/2                }{                N/2               }]
	// Which is equivalent outside of the NTT domain of adding a to the first coefficient of ct0 and b to the N/2th coefficient of ct0.
	context := evaluator.ckkscontext.contextLevel[level]
	var scaled_const, scaled_const_real, scaled_const_imag uint64
	for i := uint64(0); i < level+1; i++ {

		scaled_const_real = 0
		scaled_const_imag = 0
		scaled_const = 0

		if c_real != 0 {
			scaled_const_real = scaleUpExact(c_real, scale, evaluator.ckkscontext.moduli[i])
			scaled_const = scaled_const_real
		}

		if c_imag != 0 {
			scaled_const_imag = scaleUpExact(c_imag, scale, evaluator.ckkscontext.moduli[i])
			scaled_const_imag = ring.MRed(scaled_const_imag, context.GetNttPsi()[i][1], context.Modulus[i], context.GetMredParams()[i])
			scaled_const = ring.CRed(scaled_const+scaled_const_imag, context.Modulus[i])
		}

		scaled_const = ring.MForm(scaled_const, context.Modulus[i], context.GetBredParams()[i])

		for u := range ct0.Value() {
			for j := uint64(0); j < evaluator.ckkscontext.n>>1; j++ {
				ctOut.Value()[u].Coeffs[i][j] = ring.MRed(ct0.Value()[u].Coeffs[i][j], scaled_const, context.Modulus[i], context.GetMredParams()[i])
			}
		}

		if c_imag != 0 {
			scaled_const = ring.CRed(scaled_const_real+(context.Modulus[i]-scaled_const_imag), context.Modulus[i])
			scaled_const = ring.MForm(scaled_const, context.Modulus[i], context.GetBredParams()[i])
		}

		for u := range ct0.Value() {
			for j := evaluator.ckkscontext.n >> 1; j < evaluator.ckkscontext.n; j++ {
				ctOut.Value()[u].Coeffs[i][j] = ring.MRed(ct0.Value()[u].Coeffs[i][j], scaled_const, context.Modulus[i], context.GetMredParams()[i])
			}
		}
	}

	ctOut.SetScale(ct0.Scale() * scale)

	return nil
}

// MultByiNew multiplies ct0 by the imaginary number i, and returns the result on a newly created element.
// Does not change the scale.
func (evaluator *Evaluator) MultByiNew(ct0 *Ciphertext) (ctOut *Ciphertext, err error) {

	ctOut = evaluator.ckkscontext.NewCiphertext(1, ct0.Level(), ct0.Scale())

	return ctOut, evaluator.MultByi(ct0, ctOut)
}

// MultByiNew multiplies ct0 by the imaginary number i, and returns the result on c1.
// Does not change the scale.
func (evaluator *Evaluator) MultByi(ct0 *Ciphertext, ctOut *Ciphertext) (err error) {

	var level uint64

	level = min([]uint64{ct0.Level(), ctOut.Level()})

	context := evaluator.ckkscontext.contextLevel[level]

	var imag uint64

	// Equivalent to a mult by monomial x^(n/2) outside of the NTT domain
	for i := uint64(0); i < level+1; i++ {

		imag = context.GetNttPsi()[i][1] // Psi^2

		for u := range ctOut.value {
			for j := uint64(0); j < evaluator.ckkscontext.n>>1; j++ {
				ctOut.value[u].Coeffs[i][j] = ring.MRed(ct0.value[u].Coeffs[i][j], imag, context.Modulus[i], context.GetMredParams()[i])
			}
		}

		imag = context.Modulus[i] - imag

		for u := range ctOut.value {
			for j := evaluator.ckkscontext.n >> 1; j < evaluator.ckkscontext.n; j++ {
				ctOut.value[u].Coeffs[i][j] = ring.MRed(ct0.value[u].Coeffs[i][j], imag, context.Modulus[i], context.GetMredParams()[i])

			}
		}
	}

	return nil
}

// DivByiNew multiplies ct0 by the imaginary number 1/i = -i, and returns the result on a newly created element.
// Does not change the scale.
func (evaluator *Evaluator) DivByiNew(ct0 *Ciphertext) (ctOut *Ciphertext, err error) {

	ctOut = evaluator.ckkscontext.NewCiphertext(1, ct0.Level(), ct0.Scale())

	return ctOut, evaluator.DivByi(ct0, ctOut)
}

// DivByi multiplies ct0 by the imaginary number 1/i = -i, and returns the result on c1.
// Does not change the scale.
func (evaluator *Evaluator) DivByi(ct0 *Ciphertext, c1 *Ciphertext) (err error) {

	var level uint64

	level = min([]uint64{ct0.Level(), c1.Level()})

	context := evaluator.ckkscontext.contextLevel[level]

	var imag uint64

	// Equivalent to a mult by monomial x^(3*n/2) outside of the NTT domain
	for i := uint64(0); i < level+1; i++ {

		imag = context.Modulus[i] - context.GetNttPsi()[i][1] // -Psi^2

		for u := range c1.value {
			for j := uint64(0); j < evaluator.ckkscontext.n>>1; j++ {
				c1.value[u].Coeffs[i][j] = ring.MRed(ct0.value[u].Coeffs[i][j], imag, context.Modulus[i], context.GetMredParams()[i])
			}
		}

		imag = context.GetNttPsi()[i][1] // Psi^2

		for u := range c1.value {
			for j := evaluator.ckkscontext.n >> 1; j < evaluator.ckkscontext.n; j++ {
				c1.value[u].Coeffs[i][j] = ring.MRed(ct0.value[u].Coeffs[i][j], imag, context.Modulus[i], context.GetMredParams()[i])
			}
		}
	}

	return nil

}

// ScaleUpNew multiplies ct0 by 2^scale and sets its scale to its previous scale
// plus 2^n. Returns the result on a newly created element.
func (evaluator *Evaluator) ScaleUpNew(ct0 *Ciphertext, scale float64) (ctOut *Ciphertext, err error) {

	ctOut = evaluator.ckkscontext.NewCiphertext(ct0.Degree(), ct0.Level(), ct0.Scale())

	return ctOut, evaluator.ScaleUp(ct0, scale, ctOut)
}

// ScaleUpNew multiplies ct0 by 2^scale and sets its scale to its previous scale
// plus 2^n. Returns the result on ctOut.
func (evaluator *Evaluator) ScaleUp(ct0 *Ciphertext, scale float64, ctOut *Ciphertext) (err error) {

	if err = evaluator.MultConst(ct0, uint64(scale), ctOut); err != nil {
		return err
	}

	ctOut.SetScale(ct0.Scale() * scale)
	return nil
}

// MutByPow2New multiplies the ct0 by 2^pow2 and returns the result on a newly created element.
func (evaluator *Evaluator) MulByPow2New(ct0 *Ciphertext, pow2 uint64) (ctOut *Ciphertext, err error) {

	ctOut = evaluator.ckkscontext.NewCiphertext(ct0.Degree(), ct0.Level(), ct0.Scale())

	return ctOut, evaluator.MulByPow2(ct0.Element(), pow2, ctOut.Element())
}

// MutByPow2New multiplies ct0 by 2^pow2 and returns the result on ctOut.
func (evaluator *Evaluator) MulByPow2(ct0 *ckksElement, pow2 uint64, ctOut *ckksElement) (err error) {

	var level uint64

	level = min([]uint64{ct0.Level(), ctOut.Level()})

	for i := range ctOut.Value() {
		evaluator.ckkscontext.contextLevel[level].MulByPow2(ct0.value[i], pow2, ctOut.Value()[i])
	}

	return nil
}

// Reduce applies a modular reduction ct0 and returns the result on a newly created element.
// To be used in conjonction with function not applying modular reduction.
func (evaluator *Evaluator) ReduceNew(ct0 *Ciphertext) (ctOut *Ciphertext) {

	ctOut = evaluator.ckkscontext.NewCiphertext(ct0.Degree(), ct0.Level(), ct0.Scale())

	_ = evaluator.Reduce(ct0, ctOut)

	return ctOut
}

// Reduce applies a modular reduction ct0 and returns the result on ctOut.
// To be used in conjonction with function not applying modular reduction.
func (evaluator *Evaluator) Reduce(ct0 *Ciphertext, ctOut *Ciphertext) error {

	if ct0.Degree() != ctOut.Degree() {
		return errors.New("cannot reduce -> receiver ciphertext does not match input ciphertext degree")
	}

	for i := range ct0.value {
		evaluator.ckkscontext.contextLevel[ct0.Level()].Reduce(ct0.value[i], ctOut.Value()[i])
	}

	return nil
}

// DropLevel reduces the level of ct0 by levels and returns the result on a newly created element.
// No rescaling is applied during this procedure.
func (evaluator *Evaluator) DropLevelNew(ct0 *ckksElement, levels uint64) (ctOut *ckksElement, err error) {

	ctOut = ct0.CopyNew()

	return ctOut, evaluator.DropLevel(ctOut, levels)
}

// DropLevel reduces the level of ct0 by levels and returns the result on ct0.
// No rescaling is applied during this procedure.
func (evaluator *Evaluator) DropLevel(ct0 *ckksElement, levels uint64) error {

	if ct0.Level() == 0 {
		return errors.New("cannot drop level -> ciphertext already at level 0")
	}

	level := ct0.Level()

	for i := range ct0.value {
		ct0.value[i].Coeffs = ct0.value[i].Coeffs[:level+1-levels]
	}

	for i := uint64(0); i < levels; i++ {
		ct0.CurrentModulus().DivRound(ct0.CurrentModulus(), ring.NewUint(evaluator.ckkscontext.moduli[level-i]))
	}

	return nil
}

// RescaleNew divides ct0 by the last modulus in the modulus chain, repeats this
// procedure (each time consuming a level) until the scale reaches the original scale or would go below it, and returns the result
// on a newly created element. Since all the moduli in the modulus chain are generated to be close to the
// original scale, this procedure is equivalement to dividing the input element by the scale and adding
// some error.
func (evaluator *Evaluator) RescaleNew(ct0 *Ciphertext, threshold float64) (ctOut *Ciphertext, err error) {

	ctOut = evaluator.ckkscontext.NewCiphertext(ct0.Degree(), ct0.Level(), ct0.Scale())

	return ctOut, evaluator.Rescale(ct0, threshold, ctOut)
}

// RescaleNew divides ct0 by the last modulus in the modulus chain, repeats this
// procedure (each time consuming a level) until the scale reaches the original scale or would go below it, and returns the result
// on c1. Since all the moduli in the modulus chain are generated to be close to the
// original scale, this procedure is equivalement to dividing the input element by the scale and adding
// some error.
func (evaluator *Evaluator) Rescale(ct0 *Ciphertext, threshold float64, c1 *Ciphertext) (err error) {

	if ct0.Level() == 0 {
		return errors.New("cannot rescale -> input ciphertext already at level 0")
	}

	if ct0.Level() != c1.Level() {
		return errors.New("cannot rescale -> reciever ciphertext does not match input ciphertext level")
	}

	if ct0.Scale() >= (threshold*evaluator.ckkscontext.scalechain[c1.Level()])/2 {

		if !ct0.IsNTT() {
			return errors.New("cannot rescale -> input ciphertext not in NTT")
		}

		c1.Copy(ct0.Element())

		for c1.Scale() >= (threshold*evaluator.ckkscontext.scalechain[c1.Level()])/2 && c1.Level() != 0 {

			c1.DivScale(evaluator.ckkscontext.scalechain[c1.Level()])

			c1.CurrentModulus().DivRound(c1.CurrentModulus(), ring.NewUint(evaluator.ckkscontext.moduli[c1.Level()]))

			for i := range c1.Value() {
				rescale(evaluator, c1.Value()[i], c1.Value()[i])
			}

		}

	} else {
		c1.Copy(ct0.Element())
	}

	return nil
}

// Performes a modulus switching : starts with a base Q = {q0, q1, ..., qL}, ends up with a base Q = {q0, q1, ..., qL-1}.
// The output is also divided by the removed modulus, i.e. by qL.
func rescale(evaluator *Evaluator, p0, p1 *ring.Poly) {

	// To save NTT transforms, we keep all the polynomials except the one to remove in the NTT domain.
	// We then reduce this polynomial by the respective modulus of each other polynomials, and transform
	// it back to their respective NTT domain. We then finaly apply the modulus switching. This brings the
	// total number of NTT transforms from 2n-1 to n.

	level := len(p0.Coeffs) - 1

	var Qi, InvQl uint64

	context := evaluator.ckkscontext.contextLevel[level]
	mredParams := context.GetMredParams()

	p_tmp := evaluator.ckkscontext.contextLevel[level].NewPoly()

	ring.InvNTT(p0.Coeffs[level], p0.Coeffs[level], context.N, context.GetNttPsiInv()[level], context.GetNttNInv()[level], context.Modulus[level], context.GetMredParams()[level])

	for i := 0; i < level; i++ {

		ring.NTT(p0.Coeffs[level], p_tmp.Coeffs[0], context.N, context.GetNttPsi()[i], context.Modulus[i], context.GetMredParams()[i], context.GetBredParams()[i])

		Qi = evaluator.ckkscontext.moduli[i]
		InvQl = evaluator.ckkscontext.rescaleParams[level-1][i]

		for j := uint64(0); j < evaluator.ckkscontext.n; j++ {
			p1.Coeffs[i][j] += Qi - ring.BRedAdd(p_tmp.Coeffs[0][j], Qi, context.GetBredParams()[i]) // x[i] - x[-1]
			p1.Coeffs[i][j] = ring.MRed(p1.Coeffs[i][j], InvQl, Qi, mredParams[i])                   // (x[i] - x[-1]) * InvQl
		}
	}

	p1.Coeffs = p1.Coeffs[:level]
}

// MulRelinNew multiplies ct0 by ct1 and returns the result on a newly created element. The new scale is
// the multiplication between scales of the input elements (addition when the scale is represented in log2). An evaluation
// key can be provided to apply a relinearization step and reduce the degree of the output element. This evaluation key is only
// required when the two inputs elements are ciphertexts. If not evaluationkey is provided and the input elements are two ciphertexts,
// the resulting ciphertext will be of degree two. This function only accepts plaintexts (degree zero) and/or ciphertexts of degree one.
func (evaluator *Evaluator) MulRelinNew(op0, op1 Operand, evakey *EvaluationKey) (ctOut *Ciphertext, err error) {

	ctOut = evaluator.ckkscontext.NewCiphertext(1, min([]uint64{op0.Level(), op1.Level()}), op0.Scale()+op1.Scale())

	return ctOut, evaluator.MulRelin(op0, op1, evakey, ctOut)
}

// MulRelinNew multiplies ct0 by ct1 and returns the result on ctOut. The new scale is
// the multiplication between scales of the input elements (addition when the scale is represented in log2). An evaluation
// key can be provided to apply a relinearization step and reduce the degree of the output element. This evaluation key is only
// required when the two inputs elements are ciphertexts. If not evaluationkey is provided and the input elements are two ciphertexts,
// the resulting ciphertext will be of degree two. This function only accepts plaintexts (degree zero) and/or ciphertexts of degree one.
func (evaluator *Evaluator) MulRelin(op0, op1 Operand, evakey *EvaluationKey, ctOut *Ciphertext) error {

	el0, el1, elOut, err := evaluator.getElemAndCheckBinary(op0, op1, ctOut, max([]uint64{op0.Degree(), op1.Degree()}))
	if err != nil {
		return err
	}

	minLevel := min([]uint64{el0.Level(), el1.Level(), elOut.Level()})

	if ctOut.Level() > minLevel {
		evaluator.DropLevel(elOut, elOut.Level()-minLevel)
	}

	if el0.Degree() > 1 || el1.Degree() > 1 {
		return errors.New("cannont mul -> input element and output element must be of degree 0 or 1")
	}

	if !el0.IsNTT() {
		return errors.New("cannot mul -> op0 must be in NTT to multiply")
	}

	if !el1.IsNTT() {
		return errors.New("cannot mul -> op1 must be in NTT to multiply")
	}

	elOut.SetScale(el0.Scale() * el1.Scale())

	context := evaluator.ckkscontext.contextLevel[minLevel]

	var c_00, c_01, c0, c1, c2 *ring.Poly

	// Case Ciphertext (x) Ciphertext
	if el0.Degree()+el1.Degree() == 2 {

		c_00 = evaluator.ringpool[0]
		c_01 = evaluator.ringpool[1]

		// If the receiver ciphertext is neither of the inptus,
		// we can write directly on it.
		if elOut != el0 && elOut != el1 {

			c0 = elOut.value[0]
			c1 = elOut.value[1]

			// If the evaluation key is nil and we can write directly on the receiver, then
			// resizes the cipher text to a degree 2 ciphertext
			if evakey == nil {

				elOut.Resize(evaluator.ckkscontext, 2)
				c2 = elOut.value[2]

				// If there is however an evaluation key, then
				// we still use the mempool for the third element
			} else {

				c2 = evaluator.ringpool[4]
			}

			// If the receiver ciphertext either one of the inputs,
			// then makes use of the mempool for the three elements
		} else {

			c0 = evaluator.ringpool[2]
			c1 = evaluator.ringpool[3]
			c2 = evaluator.ringpool[4]
		}

		context.MForm(el0.value[0], c_00)
		context.MForm(el0.value[1], c_01)

		if el0 == el1 { // squaring case

			context.MulCoeffsMontgomery(c_00, el1.value[0], c0) // c0 = c[0]*c[0]
			context.MulCoeffsMontgomery(c_00, el1.value[1], c1) // c1 = 2*c[0]*c[1]
			context.Add(c1, c1, c1)
			context.MulCoeffsMontgomery(c_01, el1.value[1], c2) // c2 = c[1]*c[1]

		} else { // regular case

			context.MulCoeffsMontgomery(c_00, el1.value[0], c0) // c0 = c0[0]*c0[0]
			context.MulCoeffsMontgomery(c_00, el1.value[1], c1)
			context.MulCoeffsMontgomeryAndAddNoMod(c_01, el1.value[0], c1) // c1 = c0[0]*c1[1] + c0[1]*c1[0]
			context.MulCoeffsMontgomery(c_01, el1.value[1], c2)            // c2 = c0[1]*c1[1]
		}

		// Relinearize if a key was provided
		if evakey != nil {

			context.Copy(c0, elOut.value[0])
			context.Copy(c1, elOut.value[1])

			evaluator.switchKeysInPlace(c2, evakey.evakey, elOut.Ciphertext())

		} else { // Or copies the result on the output ciphertext if it was one of the inputs
			if elOut == el0 || elOut == el1 {
				elOut.Resize(evaluator.ckkscontext, 2)
				context.Copy(c0, elOut.value[0])
				context.Copy(c1, elOut.value[1])
				context.Copy(c2, elOut.value[2])
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

		c_00 := evaluator.ringpool[0]
		c_00.Zero()

		context.MForm(tmp0.value[0], c_00)
		context.MulCoeffsMontgomery(c_00, tmp1.value[0], elOut.value[0])
		context.MulCoeffsMontgomery(c_00, tmp1.value[1], elOut.value[1])
	}

	return nil
}

// RelinearizeNew applies the relinearization procedure on ct0 and returns the result on a newly
// created ciphertext. Requires the input ciphertext to be of degree two.
func (evaluator *Evaluator) RelinearizeNew(ct0 *Ciphertext, evakey *EvaluationKey) (ctOut *Ciphertext, err error) {

	ctOut = evaluator.ckkscontext.NewCiphertext(1, ct0.Level(), ct0.Scale())

	return ctOut, evaluator.Relinearize(ct0, evakey, ctOut)
}

// RelinearizeNew applies the relinearization procedure on ct0 and returns the result on ctOut. Requires the input ciphertext to be of degree two.
func (evaluator *Evaluator) Relinearize(ct0 *Ciphertext, evakey *EvaluationKey, ctOut *Ciphertext) (err error) {
	if ct0.Degree() != 2 {
		return errors.New("cannot relinearize -> input is not of degree 2")
	}

	if ctOut != ct0 {
		ctOut.SetScale(ct0.Scale())
	}

	context := evaluator.ckkscontext.contextLevel[min([]uint64{ct0.Level(), ctOut.Level()})]
	context.Copy(ct0.value[0], ctOut.value[0])
	context.Copy(ct0.value[1], ctOut.value[1])

	evaluator.switchKeysInPlace(ct0.value[2], evakey.evakey, ctOut)

	ctOut.Resize(evaluator.ckkscontext, 1)

	return nil
}

// Switchkeys re-encrypts ct0 under a different key and returns the result on a newly created element.
// Requires a switchinkey, which is computed from the key under which the ciphertext is currently encrypted,
// and the key under which the ciphertext will be re-encrypted.
func (evaluator *Evaluator) SwitchKeysNew(ct0 *Ciphertext, switchingKey *SwitchingKey) (ctOut *Ciphertext, err error) {

	ctOut = evaluator.ckkscontext.NewCiphertext(ct0.Degree(), ct0.Level(), ct0.Scale())

	return ctOut, evaluator.SwitchKeys(ct0, switchingKey, ctOut)
}

// Switchkeys re-encrypts ct0 under a different key and returns the result on ctOut.
// Requires a switchinkey, which is computed from the key under which the ciphertext is currently encrypted,
// and the key under which the ciphertext will be re-encrypted.
func (evaluator *Evaluator) SwitchKeys(ct0 *Ciphertext, switchingKey *SwitchingKey, ctOut *Ciphertext) error {

	if ct0.Degree() != 1 || ctOut.Degree() != 1 {
		return errors.New("cannot switchkeys -> input and output ciphertext must be of degree 1")
	}

	context := evaluator.ckkscontext.contextLevel[min([]uint64{ct0.Level(), ctOut.Level()})]

	context.Copy(ct0.value[0], ctOut.value[0])
	context.Copy(ct0.value[1], ctOut.value[1])

	evaluator.switchKeysInPlace(ct0.value[1], switchingKey, ctOut)

	return nil
}

// RotateColumnsNew rotates the columns of ct0 by k position to the left, and returns the result on a newly created element.
// If the provided element is a ciphertext, a keyswitching operation is necessary and a rotation key for the specific rotation needs to be provided.
func (evaluator *Evaluator) RotateColumnsNew(ct0 *Ciphertext, k uint64, evakey *RotationKey) (ctOut *Ciphertext, err error) {

	ctOut = evaluator.ckkscontext.NewCiphertext(ct0.Degree(), ct0.Level(), ct0.Scale())

	return ctOut, evaluator.RotateColumns(ct0, k, evakey, ctOut)
}

// RotateColumns rotates the columns of ct0 by k position to the left and returns the result on the provided receiver.
// If the provided element is a ciphertext, a keyswitching operation is necessary and a rotation key for the specific rotation needs to be provided.
func (evaluator *Evaluator) RotateColumns(ct0 *Ciphertext, k uint64, evakey *RotationKey, ctOut *Ciphertext) (err error) {

	if ct0.Degree() != 1 || ctOut.Degree() != 1 {
		return errors.New("cannot rotate -> input and output ciphertext must be of degree 1")
	}

	k &= ((evaluator.ckkscontext.n >> 1) - 1)

	if k == 0 {
		ctOut.Copy(ct0.Element())
		return nil
	}

	// Looks in the rotationkey if the corresponding rotation has been generated
	if evakey.evakey_rot_col_L[k] != nil {

		evaluator.permuteNTT(ct0, evaluator.ckkscontext.galElRotColLeft[k], evakey.evakey_rot_col_L[k], ctOut)
		return nil

	} else {

		// If not looks if the left and right pow2 rotations have been generated
		has_pow2_rotations := true
		for i := uint64(1); i < evaluator.ckkscontext.n>>1; i <<= 1 {
			if evakey.evakey_rot_col_L[i] == nil || evakey.evakey_rot_col_R[i] == nil {
				has_pow2_rotations = false
				break
			}
		}

		// If yes, computes the least amount of rotation between left and right required to apply the demanded rotation
		if has_pow2_rotations {

			if hammingWeight64(k) <= hammingWeight64((evaluator.ckkscontext.n>>1)-k) {
				evaluator.rotateColumnsLPow2(ct0, k, evakey, ctOut)
			} else {
				evaluator.rotateColumnsRPow2(ct0, (evaluator.ckkscontext.n>>1)-k, evakey, ctOut)
			}

			return nil

			// Else returns an error indicating that the keys have not been generated
		} else {
			return errors.New("cannot rotate -> specific rotation and pow2 rotations have not been generated")
		}
	}
}

func (evaluator *Evaluator) rotateColumnsLPow2(ct0 *Ciphertext, k uint64, evakey *RotationKey, ctOut *Ciphertext) {
	evaluator.rotateColumnsPow2(ct0, evaluator.ckkscontext.gen, k, evakey.evakey_rot_col_L, ctOut)
}

func (evaluator *Evaluator) rotateColumnsRPow2(ct0 *Ciphertext, k uint64, evakey *RotationKey, ctOut *Ciphertext) {
	evaluator.rotateColumnsPow2(ct0, evaluator.ckkscontext.genInv, k, evakey.evakey_rot_col_R, ctOut)
}

func (evaluator *Evaluator) rotateColumnsPow2(ct0 *Ciphertext, generator, k uint64, evakey_rot_col map[uint64]*SwitchingKey, ctOut *Ciphertext) {

	var mask, evakey_index uint64

	mask = (evaluator.ckkscontext.n << 1) - 1

	evakey_index = 1

	evaluator.ckkscontext.contextLevel[ctOut.Level()].Copy(ct0.value[0], ctOut.value[0])
	evaluator.ckkscontext.contextLevel[ctOut.Level()].Copy(ct0.value[1], ctOut.value[1])

	for k > 0 {

		if k&1 == 1 {

			evaluator.permuteNTT(ctOut, generator, evakey_rot_col[evakey_index], ctOut)
		}

		generator *= generator
		generator &= mask

		evakey_index <<= 1
		k >>= 1
	}
}

// ConjugateNew conjugates ct0 (which is equivalement to a row rotation) and returns the result on a newly
// created element. If the provided element is a ciphertext, a keyswitching operation is necessary and a rotation key
// for the row rotation needs to be provided.
func (evaluator *Evaluator) ConjugateNew(ct0 *Ciphertext, evakey *RotationKey) (ctOut *Ciphertext, err error) {

	ctOut = evaluator.ckkscontext.NewCiphertext(ct0.Degree(), ct0.Level(), ct0.Scale())

	return ctOut, evaluator.Conjugate(ct0, evakey, ctOut)
}

// ConjugateNew conjugates c0 (which is equivalement to a row rotation) and returns the result on c1.
// If the provided element is a ciphertext, a keyswitching operation is necessary and a rotation key for the row rotation needs to be provided.
func (evaluator *Evaluator) Conjugate(ct0 *Ciphertext, evakey *RotationKey, ctOut *Ciphertext) (err error) {

	if ct0.Degree() != 1 || ctOut.Degree() != 1 {
		return errors.New("cannot rotate -> input and output ciphertext must be of degree 1")
	}

	if evakey.evakey_rot_row == nil {
		return errors.New("cannot rotate -> : rows rotation key not generated")
	}

	evaluator.permuteNTT(ct0, evaluator.ckkscontext.galElRotRow, evakey.evakey_rot_row, ctOut)

	return
}

func (evaluator *Evaluator) permuteNTT(ct0 *Ciphertext, generator uint64, evakey *SwitchingKey, ctOut *Ciphertext) {

	var el0, el1 *ring.Poly

	if ct0 != ctOut {
		el0, el1 = ctOut.value[0], ctOut.value[1]
	} else {
		el0, el1 = evaluator.ringpool[0], evaluator.ringpool[1]
	}

	ring.PermuteNTT(ct0.value[0], generator, el0)
	ring.PermuteNTT(ct0.value[1], generator, el1)

	evaluator.ckkscontext.contextLevel[ctOut.Level()].Copy(el0, ctOut.value[0])
	evaluator.ckkscontext.contextLevel[ctOut.Level()].Copy(el1, ctOut.value[1])

	evaluator.switchKeysInPlace(ctOut.value[1], evakey, ctOut)
}

// Applies the general keyswitching procedure of the form [ctOut[0] + cx*evakey[0], ctOut[0] + cx*evakey[1]]

// Applies the general keyswitching procedure of the form [c0 + cx*evakey[0], c1 + cx*evakey[1]]
func (evaluator *Evaluator) switchKeysInPlace(cx *ring.Poly, evakey *SwitchingKey, ctOut *Ciphertext) {

	var level, reduce uint64

	level = ctOut.Level()

	contextQ := evaluator.ckkscontext.contextQ
	contextP := evaluator.ckkscontext.contextP

	for i := range evaluator.poolQ {
		evaluator.poolQ[i].Zero()
	}

	for i := range evaluator.poolP {
		evaluator.poolP[i].Zero()
	}

	c2_qiQ := evaluator.poolQ[0]
	c2_qiP := evaluator.poolP[0]

	//c2_tmp := contextkeys.NewPoly()
	c2 := evaluator.keyswitchpool[0]

	pool2Q := evaluator.poolQ[1]
	pool2P := evaluator.poolP[1]

	pool3Q := evaluator.poolQ[2]
	pool3P := evaluator.poolP[2]

	// We switch the element on which the switching key operation will be conducted out of the NTT domain

	//Independant of context (parameter : level)
	contextQ.InvNTTLvl(level, cx, c2)

	reduce = 0

	alpha := evaluator.ckkscontext.alpha
	beta := uint64(math.Ceil(float64(level+1) / float64(alpha)))

	c2_qi_ntt := make([]uint64, contextQ.N)

	// Key switching with crt decomposition for the Qi
	for i := uint64(0); i < beta; i++ {

		evaluator.decomposer.DecomposeAndSplit(level, i, c2, c2_qiQ, c2_qiP)

		p0idxst := i * evaluator.ckkscontext.alpha
		p0idxed := p0idxst + evaluator.decomposer.Xalpha()[i]

		//Independant of context (parameter : level)
		// We start by applying the keyswitch to the moduli of the ciphertext's context
		for x := uint64(0); x < level+1; x++ {

			qi := contextQ.Modulus[x]
			nttPsi := contextQ.GetNttPsi()[x]
			bredParams := contextQ.GetBredParams()[x]
			mredParams := contextQ.GetMredParams()[x]

			// c2_qi = cx mod qi mod qi

			// if the modulus index is equal to the one that had be to be decomposed,
			// we simply copy the coefficients of cx, which are already in NTT
			if p0idxst <= x && x < p0idxed {
				for j := uint64(0); j < contextQ.N; j++ {
					c2_qi_ntt[j] = cx.Coeffs[x][j]
				}
			} else {
				ring.NTT(c2_qiQ.Coeffs[x], c2_qi_ntt, contextQ.N, nttPsi, qi, mredParams, bredParams)
			}

			for y := uint64(0); y < contextQ.N; y++ {
				pool2Q.Coeffs[x][y] += ring.MRed(evakey.evakey[i][0].Coeffs[x][y], c2_qi_ntt[y], qi, mredParams)
				pool3Q.Coeffs[x][y] += ring.MRed(evakey.evakey[i][1].Coeffs[x][y], c2_qi_ntt[y], qi, mredParams)
			}
		}

		// And continue with the special primes, however since the keys are stored with all the moduli, we need to omit
		// the moduli between the current level and the special primes
		for j, keysindex := uint64(0), evaluator.ckkscontext.levels; j < uint64(len(evaluator.ckkscontext.specialprimes)); j, keysindex = j+1, keysindex+1 {

			pj := contextP.Modulus[j]
			nttPsi := contextP.GetNttPsi()[j]
			bredParams := contextP.GetBredParams()[j]
			mredParams := contextP.GetMredParams()[j]

			// c2_qi = cx mod qi mod pj
			ring.NTT(c2_qiP.Coeffs[j], c2_qi_ntt, contextP.N, nttPsi, pj, mredParams, bredParams)

			for y := uint64(0); y < contextP.N; y++ {
				pool2P.Coeffs[j][y] += ring.MRed(evakey.evakey[i][0].Coeffs[keysindex][y], c2_qi_ntt[y], pj, mredParams)
				pool3P.Coeffs[j][y] += ring.MRed(evakey.evakey[i][1].Coeffs[keysindex][y], c2_qi_ntt[y], pj, mredParams)
			}
		}

		//Independant of context (parameter : level)
		if reduce&1 == 1 {
			contextQ.ReduceLvl(level, pool2Q, pool2Q)
			contextQ.ReduceLvl(level, pool3Q, pool3Q)

			contextP.Reduce(pool2P, pool2P)
			contextP.Reduce(pool3P, pool3P)
		}

		reduce++
	}

	//Independant of context (parameter : level)
	if (reduce-1)&1 != 1 {
		contextQ.ReduceLvl(level, pool2Q, pool2Q)
		contextQ.ReduceLvl(level, pool3Q, pool3Q)

		contextP.Reduce(pool2P, pool2P)
		contextP.Reduce(pool3P, pool3P)
	}

	//Independant of context (parameter : level)
	evaluator.baseconverter.ModDownSplitedNTT(contextQ, contextP, evaluator.ckkscontext.rescaleParamsKeys, level, pool2Q, pool2P, pool2Q, evaluator.keyswitchpool[0])
	evaluator.baseconverter.ModDownSplitedNTT(contextQ, contextP, evaluator.ckkscontext.rescaleParamsKeys, level, pool3Q, pool3P, pool3Q, evaluator.keyswitchpool[0])

	//Independant of context (parameter : level)
	contextQ.AddLvl(level, ctOut.value[0], pool2Q, ctOut.value[0])
	contextQ.AddLvl(level, ctOut.value[1], pool3Q, ctOut.value[1])
}
