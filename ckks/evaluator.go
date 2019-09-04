package ckks

import (
	"errors"
	"github.com/lca1/lattigo/ring"
)

// Evaluator is a struct holding the necessary elements to operates the homomorphic operations between ciphertext and/or plaintexts.
// It also holds a small memory pool used to store intermediate computations.
type Evaluator struct {
	ckkscontext *CkksContext
	ringpool    [7]*ring.Poly
	ctxpool     *Ciphertext
}

// NewEvaluator creates a new Evaluator, that can be used to do homomorphic
// operations on the ciphertexts and/or plaintexts. It stores a small pool of polynomials
// and ciphertexts that will be used for intermediate values.
func (ckkscontext *CkksContext) NewEvaluator() (evaluator *Evaluator) {

	evaluator = new(Evaluator)
	evaluator.ckkscontext = ckkscontext

	context := ckkscontext.contextLevel[ckkscontext.levels-1]

	for i := 0; i < 7; i++ {
		evaluator.ringpool[i] = context.NewPoly()
	}

	evaluator.ctxpool = ckkscontext.NewCiphertext(1, ckkscontext.levels-1, ckkscontext.logScale)

	return evaluator
}

// Add adds c1 to c0 and returns the result on cOut.
func (evaluator *Evaluator) Add(c0 CkksElement, c1 CkksElement, cOut CkksElement) (err error) {

	minLevel := min([]uint64{c0.Level(), c1.Level(), cOut.Level()})

	if err = evaluateInPlace(evaluator, c0, c1, cOut, evaluator.ckkscontext.contextLevel[minLevel].Add); err != nil {
		return err
	}

	return nil
}

// AddNoMod adds c1 to c0 and returns the result on cOut, without modular reduction.
func (evaluator *Evaluator) AddNoMod(c0 *Ciphertext, c1 CkksElement, cOut *Ciphertext) (err error) {

	minlevel := min([]uint64{c0.Level(), c1.Level(), cOut.Level()})

	if err = evaluateInPlace(evaluator, c0, c1, cOut, evaluator.ckkscontext.contextLevel[minlevel].AddNoMod); err != nil {
		return err
	}
	return nil
}

// Add adds c1 to c0 and returns the result on a newly created element.
func (evaluator *Evaluator) AddNew(c0 CkksElement, c1 CkksElement) (cOut CkksElement, err error) {

	minlevel := min([]uint64{c0.Level(), c1.Level()})

	if cOut, err = evaluateNew(c0, c1, evaluator.ckkscontext.contextLevel[minlevel].Add, evaluator.ckkscontext, minlevel); err != nil {
		return nil, err
	}

	return
}

// Add adds c1 to c0 without modular reduction, and returns the result on a newly created element.
func (evaluator *Evaluator) AddNoModNew(c0 CkksElement, c1 CkksElement) (cOut CkksElement, err error) {

	minlevel := min([]uint64{c0.Level(), c1.Level()})

	if cOut, err = evaluateNew(c0, c1, evaluator.ckkscontext.contextLevel[minlevel].AddNoMod, evaluator.ckkscontext, minlevel); err != nil {
		return nil, err
	}

	return
}

// Sub subtracts c1 to c0 and returns the result on cOut.
func (evaluator *Evaluator) Sub(c0 CkksElement, c1 CkksElement, cOut CkksElement) (err error) {

	minlevel := min([]uint64{c0.Level(), c1.Level(), cOut.Level()})

	if err = evaluateInPlace(evaluator, c0, c1, cOut, evaluator.ckkscontext.contextLevel[minlevel].Sub); err != nil {
		return err
	}

	if c0.Degree() < c1.Degree() {
		for i := c0.Degree() + 1; i < c1.Degree()+1; i++ {
			evaluator.ckkscontext.contextLevel[minlevel].Neg(cOut.Value()[i], cOut.Value()[i])
		}
	}

	return nil
}

// SubNoMod subtracts c1 to c0 and returns the result on cOut, without modular reduction.
func (evaluator *Evaluator) SubNoMod(c0 CkksElement, c1 CkksElement, cOut CkksElement) (err error) {

	minlevel := min([]uint64{c0.Level(), c1.Level(), cOut.Level()})

	if err = evaluateInPlace(evaluator, c0, c1, cOut, evaluator.ckkscontext.contextLevel[minlevel].SubNoMod); err != nil {
		return err
	}

	if c0.Degree() < c1.Degree() {
		for i := c0.Degree() + 1; i < c1.Degree()+1; i++ {
			evaluator.ckkscontext.contextLevel[minlevel].Neg(cOut.Value()[i], cOut.Value()[i])
		}
	}

	return nil
}

// SubNew subtracts c1 to c0 and returns the result on a newly created element.
func (evaluator *Evaluator) SubNew(c0 CkksElement, c1 CkksElement) (cOut CkksElement, err error) {

	minlevel := min([]uint64{c0.Level(), c1.Level()})

	if cOut, err = evaluateNew(c0, c1, evaluator.ckkscontext.contextLevel[minlevel].Sub, evaluator.ckkscontext, minlevel); err != nil {
		return nil, err
	}

	if c0.Degree() < c1.Degree() {
		for i := c0.Degree() + 1; i < c1.Degree()+1; i++ {
			evaluator.ckkscontext.contextLevel[minlevel].Neg(cOut.Value()[i], cOut.Value()[i])
		}
	}

	return
}

// SubNoModNew subtracts c1 to c0 without modular reduction, and returns the result on a newly created element.
func (evaluator *Evaluator) SubNoModNew(c0 CkksElement, c1 CkksElement) (cOut CkksElement, err error) {

	minlevel := min([]uint64{c0.Level(), c1.Level()})

	if cOut, err = evaluateNew(c0, c1, evaluator.ckkscontext.contextLevel[minlevel].SubNoMod, evaluator.ckkscontext, minlevel); err != nil {
		return nil, err
	}

	if c0.Degree() < c1.Degree() {
		for i := c0.Degree() + 1; i < c1.Degree()+1; i++ {
			evaluator.ckkscontext.contextLevel[minlevel].Neg(cOut.Value()[i], cOut.Value()[i])
		}
	}

	return
}

func evaluateInPlace(evaluator *Evaluator, c0, c1, cOut CkksElement, evaluate func(*ring.Poly, *ring.Poly, *ring.Poly)) (err error) {

	var tmp0, tmp1 CkksElement // TODO : use evaluator mem pool

	maxDegree := max([]uint64{c0.Degree(), c1.Degree()})
	minDegree := min([]uint64{c0.Degree(), c1.Degree()})

	// Checks the validity of the receiver element
	if cOut.Degree() == 0 && cOut.Degree() < maxDegree {
		return errors.New("cannot evaluate(c0, c1 cOut) -> cOut is a plaintext (or an invalid ciphertext of degree 0) while c1 and/or c2 are ciphertexts of degree >= 1")
	} else {
		// Else resizes the receiver element
		cOut.Resize(evaluator.ckkscontext, maxDegree)
		evaluator.DropLevel(cOut, cOut.Level()-min([]uint64{c0.Level(), c1.Level()}))
	}

	// Checks wether or not the receiver element is the same as one of the input elements
	// and acts accordingly to avoid unnecessary element creation or element overwriting,
	// and scales properly the element before the evaluation.
	if cOut == c0 {

		if c0.Scale() > c1.Scale() {

			tmp1 = evaluator.ctxpool
			if err = evaluator.MulByPow2(c1, c0.Scale()-c1.Scale(), tmp1); err != nil {
				return err
			}

		} else if c1.Scale() > c0.Scale() {

			evaluator.MulByPow2(c0, c1.Scale()-c0.Scale(), c0)
			c0.SetScale(c1.Scale())

			tmp1 = c1

		} else {

			tmp1 = c1
		}

		tmp0 = c0

	} else if cOut == c1 {

		if c1.Scale() > c0.Scale() {
			tmp0 = evaluator.ctxpool
			if err = evaluator.MulByPow2(c0, c1.Scale()-c0.Scale(), tmp0); err != nil {
				return err
			}

		} else if c0.Scale() > c1.Scale() {

			evaluator.MulByPow2(c1, c0.Scale()-c1.Scale(), cOut)
			cOut.SetScale(c0.Scale())

			tmp0 = c0

		} else {

			tmp0 = c0
		}

		tmp1 = c1

	} else {

		if c1.Scale() > c0.Scale() {
			tmp0 = evaluator.ctxpool
			if err = evaluator.MulByPow2(c0, c1.Scale()-c0.Scale(), tmp0); err != nil {
				return err
			}
			tmp1 = c1

		} else if c0.Scale() > c1.Scale() {

			tmp0 = c0
			if err = evaluator.MulByPow2(c1, c0.Scale()-c1.Scale(), tmp1); err != nil {
				return err
			}
		} else {
			tmp0 = c0
			tmp1 = c1
		}
	}

	for i := uint64(0); i < minDegree+1; i++ {
		evaluate(tmp0.Value()[i], tmp1.Value()[i], cOut.Value()[i])
	}

	cOut.SetScale(max([]uint64{c0.Scale(), c1.Scale()}))

	// If the inputs degree differ, copies the remaining degree on the receiver
	// Also checks that the receiver is ont one of the inputs to avoid unnecessary work.

	if c0.Degree() > c1.Degree() && tmp0 != cOut {
		for i := minDegree + 1; i < maxDegree+1; i++ {
			tmp0.Value()[i].Copy(cOut.Value()[i])
		}
	} else if c1.Degree() > c0.Degree() && tmp1 != cOut {
		for i := minDegree + 1; i < maxDegree+1; i++ {
			tmp1.Value()[i].Copy(cOut.Value()[i])
		}
	}

	return nil
}

func evaluateNew(c0 CkksElement, c1 CkksElement, evaluate func(*ring.Poly, *ring.Poly, *ring.Poly), ckkscontext *CkksContext, level uint64) (cOut CkksElement, err error) {

	if c0.Degree() >= c1.Degree() {

		cOut = c0.CopyNew()

		for i := range c1.Value() {
			evaluate(cOut.Value()[i], c1.Value()[i], cOut.Value()[i])
		}

	} else {

		cOut = c1.CopyNew()

		for i := range c0.Value() {
			evaluate(cOut.Value()[i], c0.Value()[i], cOut.Value()[i])
		}
	}

	return cOut, nil
}

// Neg negates the c0 and returns the result on cOut.
func (evaluator *Evaluator) Neg(c0 CkksElement, cOut CkksElement) (err error) {

	minlevel := min([]uint64{c0.Level(), cOut.Level()})

	if c0.Degree() != cOut.Degree() {
		return errors.New("error : invalid receiver ciphertext (degree not equal to input ciphertext")
	}

	for i := range c0.Value() {
		evaluator.ckkscontext.contextLevel[minlevel].Neg(c0.Value()[i], cOut.Value()[i])
	}

	return nil
}

// Neg negates c0 and returns the result on a newly created element.
func (evaluator *Evaluator) NegNew(c0 CkksElement) (cOut CkksElement) {

	if c0.Degree() == 0 {
		cOut = evaluator.ckkscontext.NewPlaintext(c0.Level(), c0.Scale())
	} else {
		cOut = evaluator.ckkscontext.NewCiphertext(c0.Degree(), c0.Level(), c0.Scale())
	}

	for i := range c0.Value() {
		evaluator.ckkscontext.contextLevel[c0.Level()].Neg(c0.Value()[i], cOut.Value()[i])
	}

	return cOut
}

// ExtractImagNew sets the real part of c0 to the imaginary part of c0 and sets the imaginary part of c0 to zero, and returns the result on a new element.
// ex. f(a + b*i) = b. Requires a rotationkey for which the conjugate key has been generated. Scale is increased by one.
func (evaluator *Evaluator) ExtractImagNew(c0 *Ciphertext, evakey *RotationKey) (cOut *Ciphertext, err error) {

	cOut = evaluator.ckkscontext.NewCiphertext(c0.Degree(), c0.Level(), c0.Scale())

	if err = evaluator.ExtractImag(c0, evakey, cOut); err != nil {
		return nil, err
	}

	return cOut, nil
}

// ExtractImag sets the real part of c0 to the imaginary part of c0 and sets the imaginary part of c0 to zero, and returns the result on cOut.
// ex. f(a + b*i) = b. Requires a rotationkey for which the conjugate key has been generated. Scale is increased by one.
func (evaluator *Evaluator) ExtractImag(c0 *Ciphertext, evakey *RotationKey, cOut *Ciphertext) (err error) {

	if err = evaluator.Conjugate(c0, evakey, evaluator.ctxpool); err != nil {
		return err
	}

	if err = evaluator.MultByi(evaluator.ctxpool, evaluator.ctxpool); err != nil {
		return err
	}

	if err = evaluator.DivByi(c0, cOut); err != nil {
		return err
	}

	if err = evaluator.Add(cOut, evaluator.ctxpool, cOut); err != nil {
		return err
	}

	cOut.SetScale(cOut.Scale() + 1)

	return nil
}

// SwapRealImagNew swaps the real and imaginary parts of c0and returns the result on a newly created element, ex.
// f(a + b*i) = b + a * i. Requires a rotationkey for which the conjugate key has been generated.
func (evaluator *Evaluator) SwapRealImagNew(c0 *Ciphertext, evakey *RotationKey) (cOut *Ciphertext, err error) {

	cOut = evaluator.ckkscontext.NewCiphertext(c0.Degree(), c0.Level(), c0.Scale())

	if err = evaluator.SwapRealImag(c0, evakey, cOut); err != nil {
		return nil, err
	}

	return cOut, nil
}

// SwapRealImagNew swaps the real and imaginary parts of c0 and returns the result on cOut, ex.
// f(a + b*i) = b + a * i. Requires a rotationkey for which the conjugate key has been generated.
func (evaluator *Evaluator) SwapRealImag(c0 *Ciphertext, evakey *RotationKey, cOut *Ciphertext) (err error) {

	if err = evaluator.DivByi(c0, cOut); err != nil {
		return err
	}

	if err = evaluator.Conjugate(cOut, evakey, cOut); err != nil {
		return err
	}

	return nil
}

// RemoveRealNew sets the real part of c0 to zero and returns the result on a newly created element, ex. f(a + b*i) = b*i.
// Requires a rotationkey for which the conjugate key has been generated. Scale is increased by one.
func (evaluator *Evaluator) RemoveRealNew(c0 *Ciphertext, evakey *RotationKey) (cOut *Ciphertext, err error) {

	cOut = evaluator.ckkscontext.NewCiphertext(c0.Degree(), c0.Level(), c0.Scale())

	if err = evaluator.RemoveReal(c0, evakey, cOut); err != nil {
		return nil, err
	}

	return cOut, nil
}

// RemoveReal sets the real part of c0 to zero and returns the result on cOut, ex. f(a + b*i) = b*i.
// Requires a rotationkey for which the conjugate key has been generated. Scale is increased by one.
func (evaluator *Evaluator) RemoveReal(c0 *Ciphertext, evakey *RotationKey, cOut *Ciphertext) (err error) {

	if c0 != cOut {

		if err = evaluator.Conjugate(c0, evakey, cOut); err != nil {
			return err
		}
		if err = evaluator.Sub(c0, cOut, cOut); err != nil {
			return err
		}

	} else {

		if err = evaluator.Conjugate(c0, evakey, evaluator.ctxpool); err != nil {
			return err
		}

		if err = evaluator.Sub(cOut, evaluator.ctxpool, cOut); err != nil {
			return err
		}
	}

	cOut.SetScale(cOut.Scale() + 1)

	return nil
}

// RemoveImagNew sets the imaginary part of c0 to zero and returns the result on a newly created element, ex. f(a + b*i) = a.
// Requires a rotationkey for which the conjugate key has been generated. Scale is increased by one.
func (evaluator *Evaluator) RemoveImagNew(c0 *Ciphertext, evakey *RotationKey) (cOut *Ciphertext, err error) {

	cOut = evaluator.ckkscontext.NewCiphertext(c0.Degree(), c0.Level(), c0.Scale())

	if err = evaluator.RemoveImag(c0, evakey, cOut); err != nil {
		return nil, err
	}

	return cOut, nil
}

// RemoveImag sets the imaginary part of c0 to zero and returns the result on cOut, ex. f(a + b*i) = a.
// Requires a rotationkey for which the conjugate key has been generated. Scale is increased by one.
func (evaluator *Evaluator) RemoveImag(c0 *Ciphertext, evakey *RotationKey, cOut *Ciphertext) (err error) {

	if c0 != cOut {

		if err = evaluator.Conjugate(c0, evakey, cOut); err != nil {
			return err
		}
		if err = evaluator.Add(c0, cOut, cOut); err != nil {
			return err
		}

	} else {

		if err = evaluator.Conjugate(c0, evakey, evaluator.ctxpool); err != nil {
			return err
		}

		if err = evaluator.Add(evaluator.ctxpool, cOut, cOut); err != nil {
			return err
		}

	}

	cOut.SetScale(cOut.Scale() + 1)

	return nil
}

// AddConstNew adds the input constant (which can be an uint64, int64, float64 or complex128) to c0 and returns the result on a new element.
func (evaluator *Evaluator) AddConstNew(c0 CkksElement, constant interface{}) (cOut CkksElement) {
	cOut = c0.CopyNew()
	evaluator.AddConst(c0, constant, cOut)
	return cOut
}

// AddConstNew adds the input constant (which can be an uint64, int64, float64 or complex128) to c0 and returns the result on cOut.
func (evaluator *Evaluator) AddConst(c0 CkksElement, constant interface{}, cOut CkksElement) (err error) {

	var level uint64

	if level, err = checkLevels([]CkksElement{c0, cOut}); err != nil {
		return err
	}

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
	// Which is equivalent outside of the NTT domain of adding a to the first coefficient of c0 and b to the N/2th coefficient of c0.
	for i := uint64(0); i < level+1; i++ {
		scaled_const_real = 0
		scaled_const_imag = 0
		scaled_const = 0

		if c_real != 0 {
			scaled_const_real = scaleUp(c_real, c0.Scale(), evaluator.ckkscontext.moduli[i])
			scaled_const = scaled_const_real
		}

		if c_imag != 0 {
			scaled_const_imag = ring.MRed(scaleUp(c_imag, c0.Scale(), evaluator.ckkscontext.moduli[i]), context.GetNttPsi()[i][1], context.Modulus[i], context.GetMredParams()[i])
			scaled_const = ring.CRed(scaled_const+scaled_const_imag, context.Modulus[i])
		}

		for j := uint64(0); j < evaluator.ckkscontext.n>>1; j++ {
			cOut.Value()[0].Coeffs[i][j] = ring.CRed(c0.Value()[0].Coeffs[i][j]+scaled_const, evaluator.ckkscontext.moduli[i])
		}

		if c_imag != 0 {
			scaled_const = ring.CRed(scaled_const_real+(context.Modulus[i]-scaled_const_imag), context.Modulus[i])
		}

		for j := evaluator.ckkscontext.n >> 1; j < evaluator.ckkscontext.n; j++ {
			cOut.Value()[0].Coeffs[i][j] = ring.CRed(c0.Value()[0].Coeffs[i][j]+scaled_const, evaluator.ckkscontext.moduli[i])
		}
	}

	return nil
}

// MultByConstAndAdd multiplies c0 by the input constant, and adds it to the receiver element (does not modify the input
// element), ex. cOut(x) = cOut(x) + c0(x) * (a+bi). This functions  removes the need of storing the intermediate value c(x) * (a+bi).
// This function will modifie the level and the scale of the receiver element depending on the level and the scale of the input
// element and the type of the constant. The level of the receiver element will be set to min(input.level, receiver.level).
// The scale of the receiver element will be set to the scale that the input element would have after the multiplication by the constant.
func (evaluator *Evaluator) MultByConstAndAdd(c0 CkksElement, constant interface{}, cOut CkksElement) (err error) {

	var level uint64

	level = min([]uint64{c0.Level(), cOut.Level()})

	// Forces a drop of cOut level to c0 level
	if cOut.Level() > level {
		evaluator.DropLevel(cOut, cOut.Level()-level)
	}

	var c_real, c_imag float64
	var scale uint64

	// Converts to float64 and determines if a scale is required (which is the case if either real or imag has a rational part)
	scale = 0
	switch constant.(type) {
	case complex128:
		c_real = real(constant.(complex128))
		c_imag = imag(constant.(complex128))

		if c_real != 0 {
			value_int := int64(c_real)
			value_float := c_real - float64(value_int)

			if value_float != 0 {
				scale = evaluator.ckkscontext.logScale
			}
		}

		if c_imag != 0 {
			value_int := int64(c_imag)
			value_float := c_imag - float64(value_int)

			if value_float != 0 {
				scale = evaluator.ckkscontext.logScale
			}
		}

	case float64:
		c_real = constant.(float64)
		c_imag = float64(0)

		if c_real != 0 {
			value_int := int64(c_real)
			value_float := c_real - float64(value_int)

			if value_float != 0 {
				scale = evaluator.ckkscontext.logScale
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
	if scale != 0 {

		// If cOut scaling is smaller than c0's scale + the default scaling,
		// then brings cOut scale to c0's scale.
		if cOut.Scale() < c0.Scale()+evaluator.ckkscontext.logScale {

			evaluator.MulByPow2(cOut, (evaluator.ckkscontext.logScale+c0.Scale())-cOut.Scale(), cOut)
			cOut.SetScale(c0.Scale() + evaluator.ckkscontext.logScale)

			// If cOut.Scale() > ((a+bi)*scale)*c0(x), then sets the scale to
			// bring c(x)*scale to the level of cOut(x) scale
		} else if cOut.Scale() > c0.Scale()+evaluator.ckkscontext.logScale {
			scale = cOut.Scale() - c0.Scale()
		}

		// If no scaling is required, the sets the appropriate scale such that
		// c0(x)*scale matches cOut(x) scale without modifiying c0(x) scale.
	} else {

		if cOut.Scale() > c0.Scale() {

			scale = cOut.Scale() - c0.Scale()

		} else if c0.Scale() > cOut.Scale() {

			evaluator.MulByPow2(cOut, c0.Scale()-cOut.Scale(), cOut)
			cOut.SetScale(c0.Scale())
		}
	}

	// Component wise multiplication of the following vector to the ciphertext :
	// [a + b*psi_qi^2, ....., a + b*psi_qi^2, a - b*psi_qi^2, ...., a - b*psi_qi^2] mod Qi
	// [{                  N/2                }{                N/2               }]
	// Which is equivalent outside of the NTT domain of adding a to the first coefficient of c0 and b to the N/2th coefficient of c0.
	for i := uint64(0); i < level+1; i++ {

		scaled_const_real = 0
		scaled_const_imag = 0
		scaled_const = 0

		if c_real != 0 {
			scaled_const_real = scaleUp(c_real, scale, evaluator.ckkscontext.moduli[i])
			scaled_const = scaled_const_real
		}

		if c_imag != 0 {
			scaled_const_imag = scaleUp(c_imag, scale, evaluator.ckkscontext.moduli[i])
			scaled_const_imag = ring.MRed(scaled_const_imag, context.GetNttPsi()[i][1], context.Modulus[i], context.GetMredParams()[i])
			scaled_const = ring.CRed(scaled_const+scaled_const_imag, context.Modulus[i])
		}

		scaled_const = ring.MForm(scaled_const, context.Modulus[i], context.GetBredParams()[i])

		for u := range c0.Value() {
			for j := uint64(0); j < evaluator.ckkscontext.n>>1; j++ {
				cOut.Value()[u].Coeffs[i][j] = ring.CRed(cOut.Value()[u].Coeffs[i][j]+ring.MRed(c0.Value()[u].Coeffs[i][j], scaled_const, context.Modulus[i], context.GetMredParams()[i]), context.Modulus[i])
			}
		}

		if c_imag != 0 {
			scaled_const = ring.CRed(scaled_const_real+(context.Modulus[i]-scaled_const_imag), context.Modulus[i])
			scaled_const = ring.MForm(scaled_const, context.Modulus[i], context.GetBredParams()[i])
		}

		for u := range c0.Value() {
			for j := evaluator.ckkscontext.n >> 1; j < evaluator.ckkscontext.n; j++ {
				cOut.Value()[u].Coeffs[i][j] = ring.CRed(cOut.Value()[u].Coeffs[i][j]+ring.MRed(c0.Value()[u].Coeffs[i][j], scaled_const, context.Modulus[i], context.GetMredParams()[i]), context.Modulus[i])
			}
		}
	}

	return nil

}

// MultConstNew multiplies c0 by the input constant and returns the result on a newly created element.
// The scale of the output element will depend on the scale of the input element and the constant (if the constant
// needs to be scaled (its rational part is not zero)). The constant can be an uint64, int64, float64 or complex128.
func (evaluator *Evaluator) MultConstNew(c0 CkksElement, constant interface{}) (cOut CkksElement, err error) {

	cOut = evaluator.ckkscontext.NewCiphertext(1, c0.Level(), c0.Scale())

	if err = evaluator.MultConst(c0, constant, cOut); err != nil {
		return nil, err
	}

	return cOut, nil
}

// MultConstNew multiplies c0 by the input constant and returns the result on cOut.
// The scale of the output element will depend on the scale of the input element and the constant (if the constant
// needs to be scaled (its rational part is not zero)). The constant can be an uint64, int64, float64 or complex128.
func (evaluator *Evaluator) MultConst(c0 CkksElement, constant interface{}, cOut CkksElement) (err error) {

	var level uint64

	level = min([]uint64{c0.Level(), cOut.Level()})

	var c_real, c_imag float64
	var scale uint64

	// Converts to float64 and determines if a scale is required (which is the case if either real or imag has a rational part)
	scale = 0
	switch constant.(type) {
	case complex128:
		c_real = real(constant.(complex128))
		c_imag = imag(constant.(complex128))

		if c_real != 0 {
			value_int := int64(c_real)
			value_float := c_real - float64(value_int)

			if value_float != 0 {
				scale = evaluator.ckkscontext.logScale
			}
		}

		if c_imag != 0 {
			value_int := int64(c_imag)
			value_float := c_imag - float64(value_int)

			if value_float != 0 {
				scale = evaluator.ckkscontext.logScale
			}
		}

	case float64:
		c_real = constant.(float64)
		c_imag = float64(0)

		if c_real != 0 {
			value_int := int64(c_real)
			value_float := c_real - float64(value_int)

			if value_float != 0 {
				scale = evaluator.ckkscontext.logScale
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
	// Which is equivalent outside of the NTT domain of adding a to the first coefficient of c0 and b to the N/2th coefficient of c0.
	context := evaluator.ckkscontext.contextLevel[level]
	var scaled_const, scaled_const_real, scaled_const_imag uint64
	for i := uint64(0); i < level+1; i++ {

		scaled_const_real = 0
		scaled_const_imag = 0
		scaled_const = 0

		if c_real != 0 {
			scaled_const_real = scaleUp(c_real, scale, evaluator.ckkscontext.moduli[i])
			scaled_const = scaled_const_real
		}

		if c_imag != 0 {
			scaled_const_imag = scaleUp(c_imag, scale, evaluator.ckkscontext.moduli[i])
			scaled_const_imag = ring.MRed(scaled_const_imag, context.GetNttPsi()[i][1], context.Modulus[i], context.GetMredParams()[i])
			scaled_const = ring.CRed(scaled_const+scaled_const_imag, context.Modulus[i])
		}

		scaled_const = ring.MForm(scaled_const, context.Modulus[i], context.GetBredParams()[i])

		for u := range c0.Value() {
			for j := uint64(0); j < evaluator.ckkscontext.n>>1; j++ {
				cOut.Value()[u].Coeffs[i][j] = ring.MRed(c0.Value()[u].Coeffs[i][j], scaled_const, context.Modulus[i], context.GetMredParams()[i])
			}
		}

		if c_imag != 0 {
			scaled_const = ring.CRed(scaled_const_real+(context.Modulus[i]-scaled_const_imag), context.Modulus[i])
			scaled_const = ring.MForm(scaled_const, context.Modulus[i], context.GetBredParams()[i])
		}

		for u := range c0.Value() {
			for j := evaluator.ckkscontext.n >> 1; j < evaluator.ckkscontext.n; j++ {
				cOut.Value()[u].Coeffs[i][j] = ring.MRed(c0.Value()[u].Coeffs[i][j], scaled_const, context.Modulus[i], context.GetMredParams()[i])
			}
		}
	}

	cOut.SetScale(c0.Scale() + scale)

	return nil
}

// MultByiNew multiplies c0 by the imaginary number i, and returns the result on a newly created element.
// Does not change the scale.
func (evaluator *Evaluator) MultByiNew(c0 CkksElement) (c1 CkksElement, err error) {
	c1 = evaluator.ckkscontext.NewCiphertext(1, c0.Level(), c0.Scale())

	if err = evaluator.MultByi(c0, c1); err != nil {
		return nil, err
	}

	return c1, nil
}

// MultByiNew multiplies c0 by the imaginary number i, and returns the result on c1.
// Does not change the scale.
func (evaluator *Evaluator) MultByi(c0 CkksElement, c1 CkksElement) (err error) {

	var level uint64

	if level, err = checkLevels([]CkksElement{c0, c1}); err != nil {
		return err
	}

	context := evaluator.ckkscontext.contextLevel[level]

	var imag uint64

	// Equivalent to a mult by monomial x^(n/2) outside of the NTT domain
	for i := uint64(0); i < level+1; i++ {

		imag = context.GetNttPsi()[i][1] // Psi^2

		for u := range c1.Value() {
			for j := uint64(0); j < evaluator.ckkscontext.n>>1; j++ {
				c1.Value()[u].Coeffs[i][j] = ring.MRed(c0.Value()[u].Coeffs[i][j], imag, context.Modulus[i], context.GetMredParams()[i])
			}
		}

		imag = context.Modulus[i] - imag

		for u := range c1.Value() {
			for j := evaluator.ckkscontext.n >> 1; j < evaluator.ckkscontext.n; j++ {
				c1.Value()[u].Coeffs[i][j] = ring.MRed(c0.Value()[u].Coeffs[i][j], imag, context.Modulus[i], context.GetMredParams()[i])

			}
		}
	}

	return nil
}

// DivByiNew multiplies c0 by the imaginary number 1/i = -i, and returns the result on a newly created element.
// Does not change the scale.
func (evaluator *Evaluator) DivByiNew(c0 CkksElement) (c1 CkksElement, err error) {
	c1 = evaluator.ckkscontext.NewCiphertext(1, c0.Level(), c0.Scale())

	if err = evaluator.DivByi(c0, c1); err != nil {
		return nil, err
	}

	return c1, nil
}

// DivByi multiplies c0 by the imaginary number 1/i = -i, and returns the result on c1.
// Does not change the scale.
func (evaluator *Evaluator) DivByi(c0 CkksElement, c1 CkksElement) (err error) {

	var level uint64

	if level, err = checkLevels([]CkksElement{c0, c1}); err != nil {
		return err
	}

	context := evaluator.ckkscontext.contextLevel[level]

	var imag uint64

	// Equivalent to a mult by monomial x^(3*n/2) outside of the NTT domain
	for i := uint64(0); i < level+1; i++ {

		imag = context.Modulus[i] - context.GetNttPsi()[i][1] // -Psi^2

		for u := range c1.Value() {
			for j := uint64(0); j < evaluator.ckkscontext.n>>1; j++ {
				c1.Value()[u].Coeffs[i][j] = ring.MRed(c0.Value()[u].Coeffs[i][j], imag, context.Modulus[i], context.GetMredParams()[i])
			}
		}

		imag = context.GetNttPsi()[i][1] // Psi^2

		for u := range c1.Value() {
			for j := evaluator.ckkscontext.n >> 1; j < evaluator.ckkscontext.n; j++ {
				c1.Value()[u].Coeffs[i][j] = ring.MRed(c0.Value()[u].Coeffs[i][j], imag, context.Modulus[i], context.GetMredParams()[i])
			}
		}
	}

	return nil

}

// ScaleUpNew multiplies c0 by 2^scale and sets its scale to its previous scale
// plus 2^n. Returns the result on a newly created element.
func (evaluator *Evaluator) ScaleUpNew(c0 CkksElement, scale uint64) (cOut CkksElement, err error) {

	if cOut, err = evaluator.MulByPow2New(c0, scale); err != nil {
		return nil, err
	}

	cOut.SetScale(cOut.Scale() + scale)
	return cOut, nil
}

// ScaleUpNew multiplies c0 by 2^scale and sets its scale to its previous scale
// plus 2^n. Returns the result on cOut.
func (evaluator *Evaluator) ScaleUp(c0 CkksElement, scale uint64, cOut CkksElement) (err error) {
	if err = evaluator.MulByPow2(c0, scale, cOut); err != nil {
		return err
	}
	cOut.SetScale(c0.Scale() + scale)
	return nil
}

// MutByPow2New multiplies the c0 by 2^pow2 and returns the result on a newly created element.
func (evaluator *Evaluator) MulByPow2New(c0 CkksElement, pow2 uint64) (cOut CkksElement, err error) {

	cOut = evaluator.ckkscontext.NewCiphertext(1, c0.Level(), c0.Scale())

	if err = evaluator.MulByPow2(c0, pow2, cOut); err != nil {
		return nil, err
	}

	return cOut, nil
}

// MutByPow2New multiplies c0 by 2^pow2 and returns the result on cOut.
func (evaluator *Evaluator) MulByPow2(c0 CkksElement, pow2 uint64, cOut CkksElement) (err error) {

	var level uint64

	if level, err = checkLevels([]CkksElement{c0, cOut}); err != nil {
		return err
	}

	for i := range cOut.Value() {
		evaluator.ckkscontext.contextLevel[level].MulByPow2(c0.Value()[i], pow2, cOut.Value()[i])
	}

	return nil
}

// Reduce applies a modular reduction c0 and returns the result on a newly created element.
// To be used in conjonction with function not applying modular reduction.
func (evaluator *Evaluator) ReduceNew(c0 CkksElement) (cOut CkksElement) {

	if c0.Degree() == 0 {
		cOut = evaluator.ckkscontext.NewPlaintext(c0.Level(), c0.Scale())
	} else {
		cOut = evaluator.ckkscontext.NewCiphertext(c0.Degree(), c0.Level(), c0.Scale())
	}

	for i := range c0.Value() {
		evaluator.ckkscontext.contextLevel[c0.Level()].Reduce(c0.Value()[i], cOut.Value()[i])
	}

	return nil
}

// Reduce applies a modular reduction c0 and returns the result on cOut.
// To be used in conjonction with function not applying modular reduction.
func (evaluator *Evaluator) Reduce(c0 CkksElement, cOut CkksElement) error {

	if c0.Degree() != cOut.Degree() {
		return errors.New("error : invalide ciphertext receiver (degree doesn't match c0.Degree")
	}

	for i := range c0.Value() {
		evaluator.ckkscontext.contextLevel[c0.Level()].Reduce(c0.Value()[i], cOut.Value()[i])
	}

	return nil
}

// DropLevel reduces the level of c0 by levels and returns the result on a newly created element.
// No rescaling is applied during this procedure.
func (evaluator *Evaluator) DropLevelNew(c0 CkksElement, levels uint64) (cOut CkksElement, err error) {

	if c0.Level() == 0 {
		return nil, errors.New("cannot drop level -> element already at level 0")
	}

	cOut = c0.CopyNew()
	evaluator.DropLevel(cOut, levels)

	return cOut, nil
}

// DropLevel reduces the level of c0 by levels and returns the result on c0.
// No rescaling is applied during this procedure.
func (evaluator *Evaluator) DropLevel(c0 CkksElement, levels uint64) error {

	if c0.Level() == 0 {
		return errors.New("error : cannot drop level, ciphertext already at level 0")
	}

	level := c0.Level()

	for i := range c0.Value() {
		c0.Value()[i].Coeffs = c0.Value()[i].Coeffs[:level+1-levels]
	}

	for i := uint64(0); i < levels; i++ {
		c0.CurrentModulus().DivRound(c0.CurrentModulus(), ring.NewUint(evaluator.ckkscontext.moduli[level-i]))
	}

	return nil
}

// RescaleNew divides c0 by the last modulus in the modulus chain, repeats this
// procedure (each time consuming a level) until the scale reaches the original scale or would go below it, and returns the result
// on a newly created element. Since all the moduli in the modulus chain are generated to be close to the
// original scale, this procedure is equivalement to dividing the input element by the scale and adding
// some error.
func (evaluator *Evaluator) RescaleNew(c0 CkksElement) (cOut CkksElement, err error) {

	if c0.Level() == 0 {
		return nil, errors.New("can't rescale, ciphertext already at level 0")
	}

	cOut = evaluator.ckkscontext.NewCiphertext(c0.Degree(), c0.Level(), c0.Scale())
	evaluator.Rescale(c0, cOut)

	return cOut, nil
}

// RescaleNew divides c0 by the last modulus in the modulus chain, repeats this
// procedure (each time consuming a level) until the scale reaches the original scale or would go below it, and returns the result
// on c1. Since all the moduli in the modulus chain are generated to be close to the
// original scale, this procedure is equivalement to dividing the input element by the scale and adding
// some error.
func (evaluator *Evaluator) Rescale(c0, c1 CkksElement) (err error) {

	if c0.Level() != c1.Level() {
		return errors.New("invalid receiver : ciphertexts not on the same level")
	}

	if c0.Scale() >= evaluator.ckkscontext.logQ+evaluator.ckkscontext.logScale {

		if c0.Level() == 0 {
			return errors.New("can't rescale, ciphertext already at level 0")
		}

		if !c0.IsNTT() {
			return errors.New("ciphertext not in NTT")
		}

		if c0 != c1 {
			c0.Copy(c1.(*Ciphertext)) // TODO : make copy work for both plaintext and ciphertext
		}

		for c1.Scale() >= evaluator.ckkscontext.logScale+evaluator.ckkscontext.logQ && c1.Level() > 0 {

			c1.CurrentModulus().DivRound(c1.CurrentModulus(), ring.NewUint(evaluator.ckkscontext.moduli[c1.Level()]))

			for i := range c1.Value() {
				rescale(evaluator, c1.Value()[i], c1.Value()[i])
			}

			c1.SetScale(c1.Scale() - evaluator.ckkscontext.logQ)
		}

	} else {
		if c0 != c1 {
			c0.Copy(c1.(*Ciphertext))
		}
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

	p_tmp := evaluator.ringpool[0]

	ring.InvNTT(p0.Coeffs[level], p0.Coeffs[level], context.N, context.GetNttPsiInv()[level], context.GetNttNInv()[level], context.Modulus[level], context.GetMredParams()[level])

	for i := 0; i < level; i++ {

		ring.NTT(p0.Coeffs[level], p_tmp.Coeffs[0], context.N, context.GetNttPsi()[i], context.Modulus[i], context.GetMredParams()[i], context.GetBredParams()[i])

		Qi = evaluator.ckkscontext.moduli[i]
		InvQl = evaluator.ckkscontext.rescalParams[level-1][i]

		for j := uint64(0); j < evaluator.ckkscontext.n; j++ {
			p1.Coeffs[i][j] = ring.CRed(p0.Coeffs[i][j]+(Qi-p_tmp.Coeffs[0][j]), Qi) // x[i] - x[-1]
			p1.Coeffs[i][j] = ring.MRed(p1.Coeffs[i][j], InvQl, Qi, mredParams[i])   // (x[i] - x[-1]) * InvQl
		}
	}

	p1.Coeffs = p1.Coeffs[:level]
}

// MulRelinNew multiplies ct0 by ct1 and returns the result on a newly created element. The new scale is
// the multiplication between scales of the input elements (addition when the scale is represented in log2). An evaluation
// key can be provided to apply a relinearization step and reduce the degree of the output element. This evaluation key is only
// required when the two inputs elements are ciphertexts. If not evaluationkey is provided and the input elements are two ciphertexts,
// the resulting ciphertext will be of degree two. This function only accepts plaintexts (degree zero) and/or ciphertexts of degree one.
func (evaluator *Evaluator) MulRelinNew(ct0, ct1 CkksElement, evakey *EvaluationKey) (cOut CkksElement, err error) {

	minlevel := min([]uint64{ct0.Level(), ct1.Level()})

	if ct0.Degree()+ct1.Degree() == 0 {
		cOut = evaluator.ckkscontext.NewPlaintext(minlevel, ct0.Scale()+ct1.Scale())
	} else {
		cOut = evaluator.ckkscontext.NewCiphertext(1, minlevel, ct0.Scale()+ct1.Scale())
	}

	if err = evaluator.MulRelin(ct0, ct1, evakey, cOut); err != nil {
		return nil, err
	}

	return cOut, nil
}

// MulRelinNew multiplies ct0 by ct1 and returns the result on cOut. The new scale is
// the multiplication between scales of the input elements (addition when the scale is represented in log2). An evaluation
// key can be provided to apply a relinearization step and reduce the degree of the output element. This evaluation key is only
// required when the two inputs elements are ciphertexts. If not evaluationkey is provided and the input elements are two ciphertexts,
// the resulting ciphertext will be of degree two. This function only accepts plaintexts (degree zero) and/or ciphertexts of degree one.
func (evaluator *Evaluator) MulRelin(ct0, ct1 CkksElement, evakey *EvaluationKey, cOut CkksElement) error {

	minlevel := min([]uint64{ct0.Level(), ct1.Level()})

	if cOut.Level() > minlevel {
		evaluator.DropLevel(cOut, cOut.Level()-minlevel)
	}

	if ct0.Degree() > 1 || ct1.Degree() > 1 {
		return errors.New("cannont mul -> input ciphertexts and output ciphertext must be of degree 0 or 1")
	}

	if ct0.Degree()+ct1.Degree() > 0 && cOut.Degree() == 0 {
		return errors.New("cannot mul -> receiver is a plaintext while one of the inputs is a ciphertext")
	}

	if !ct0.IsNTT() {
		return errors.New("cannot mul -> ct0 must be in NTT to multiply")
	}

	if !ct1.IsNTT() {
		return errors.New("cannot mul -> ct1 must be in NTT to multiply")
	}

	if ct0 == cOut {
		cOut.SetScale(cOut.Scale() + ct1.Scale())
	} else if ct1 == cOut {
		cOut.SetScale(cOut.Scale() + ct0.Scale())
	} else {
		cOut.SetScale(ct0.Scale() + ct1.Scale())
	}

	context := evaluator.ckkscontext.contextLevel[minlevel]

	var c_00, c_01, c0, c1, c2 *ring.Poly

	// Case Ciphertext (x) Ciphertext
	if ct0.Degree()+ct1.Degree() == 2 {

		c_00 = evaluator.ringpool[0]
		c_01 = evaluator.ringpool[1]

		// If the receiver ciphertext is neither of the inptus,
		// we can write directly on it.
		if cOut != ct0 && cOut != ct1 {

			c0 = cOut.Value()[0]
			c1 = cOut.Value()[1]

			// If the evaluation key is nil and we can write directly on the receiver, then
			// resizes the cipher text to a degree 2 ciphertext
			if evakey == nil {

				cOut.Resize(evaluator.ckkscontext, 2)
				c2 = cOut.Value()[2]

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

		context.MForm(ct0.Value()[0], c_00)
		context.MForm(ct0.Value()[1], c_01)

		if ct0 == ct1 { // squaring case

			context.MulCoeffsMontgomery(c_00, ct1.Value()[0], c0) // c0 = c[0]*c[0]
			context.MulCoeffsMontgomery(c_00, ct1.Value()[1], c1) // c1 = 2*c[0]*c[1]
			context.Add(c1, c1, c1)
			context.MulCoeffsMontgomery(c_01, ct1.Value()[1], c2) // c2 = c[1]*c[1]

		} else { // regular case

			c1a := evaluator.ringpool[5]
			c1b := evaluator.ringpool[6]

			context.Add(c_00, c_01, c1a)
			context.Add(ct1.Value()[0], ct1.Value()[1], c1b)

			context.MulCoeffsMontgomery(c_00, ct1.Value()[0], c0) // c0 = c0[0]*c1[0]
			context.MulCoeffsMontgomery(c1a, c1b, c1)             // c1 = c0[0]*c1[0] + c0[0]*c1[1] + c0[1]*c1[0] + c0[1]*c1[1]
			context.MulCoeffsMontgomery(c_01, ct1.Value()[1], c2) // c2 = c0[1]*c1[1]

			context.Sub(c1, c0, c1) // c2 = c0[0]*c1[1] + c0[1]*c1[0] + c0[1]*c1[1]
			context.Sub(c1, c2, c1) // c2 = c0[0]*c1[1] + c0[1]*c1[0]

		}

		// Relinearize if a key was provided
		if evakey != nil {

			switchKeys(evaluator, c0, c1, c2, evakey.evakey, cOut.(*Ciphertext).value[0], cOut.(*Ciphertext).value[1])

		} else { // Or copies the result on the output ciphertext if it was one of the inputs
			if cOut == ct0 || cOut == ct1 {
				cOut.Resize(evaluator.ckkscontext, 2)
				c0.Copy(cOut.Value()[0])
				c1.Copy(cOut.Value()[1])
				c2.Copy(cOut.Value()[2])
			}
		}

		// Case Plaintext (x) Ciphertext or Ciphertext (x) Plaintext
	} else if ct0.Degree()+ct1.Degree() == 1 {

		var tmp0, tmp1 CkksElement

		if ct0.Degree() == 1 {
			tmp0, tmp1 = ct1, ct0
		} else {
			tmp0, tmp1 = ct0, ct1
		}

		c_00 := evaluator.ringpool[0]
		c_00.Zero()

		context.MForm(tmp0.Value()[0], c_00)
		context.MulCoeffsMontgomery(c_00, tmp1.Value()[0], cOut.Value()[0])
		context.MulCoeffsMontgomery(c_00, tmp1.Value()[1], cOut.Value()[1])

		// Case Plaintext (x) Plaintext
	} else {
		context.MulCoeffs(ct0.Value()[0], ct1.Value()[0], cOut.Value()[0])
	}

	return nil
}

// RelinearizeNew applies the relinearization procedure on cIn and returns the result on a newly
// created ciphertext. Requires the input ciphertext to be of degree two.
func (evaluator *Evaluator) RelinearizeNew(cIn *Ciphertext, evakey *EvaluationKey) (cOut *Ciphertext, err error) {

	if cIn.Degree() < 2 || cIn.Degree() > 2 {
		return nil, errors.New("cannot relinearize -> input is not of degree 2")
	}

	cOut = evaluator.ckkscontext.NewCiphertext(1, cIn.Level(), cIn.Scale())

	switchKeys(evaluator, cIn.value[0], cIn.value[1], cIn.value[2], evakey.evakey, cOut.value[0], cOut.value[1])

	return
}

// RelinearizeNew applies the relinearization procedure on cIn and returns the result on cOut. Requires the input ciphertext to be of degree two.
func (evaluator *Evaluator) Relinearize(cIn *Ciphertext, evakey *EvaluationKey, cOut *Ciphertext) (err error) {
	if cIn.Degree() != 2 {
		return errors.New("cannot relinearize -> input is not of degree 2")
	}

	if cOut != cIn {
		cOut.SetScale(cIn.Scale())
	}
	switchKeys(evaluator, cIn.value[0], cIn.value[1], cIn.value[2], evakey.evakey, cOut.value[0], cOut.value[1])
	cOut.Resize(evaluator.ckkscontext, 1)

	return nil
}

// Switchkeys re-encrypts cIn under a different key and returns the result on a newly created element.
// Requires a switchinkey, which is computed from the key under which the ciphertext is currently encrypted,
// and the key under which the ciphertext will be re-encrypted.
func (evaluator *Evaluator) SwitchKeysNew(cIn *Ciphertext, switchingKey *SwitchingKey) (cOut *Ciphertext, err error) {

	if cIn.Degree() != 1 {
		return nil, errors.New("error : ciphertext must be of degree 1 to allow key switching")
	}

	cOut = evaluator.ckkscontext.NewCiphertext(cIn.Degree(), cIn.Level(), cIn.Scale())

	switchKeys(evaluator, cIn.value[0], cIn.value[1], cIn.value[1], switchingKey, cOut.value[0], cOut.value[1])

	return cOut, nil
}

// Switchkeys re-encrypts cIn under a different key and returns the result on cOut.
// Requires a switchinkey, which is computed from the key under which the ciphertext is currently encrypted,
// and the key under which the ciphertext will be re-encrypted.
func (evaluator *Evaluator) SwitchKeys(cIn *Ciphertext, switchingKey *SwitchingKey, cOut *Ciphertext) error {

	if cIn.Degree() != 1 {
		return errors.New("error : ciphertext must be of degree 1 to allow key switching")
	}

	if cOut.Degree() != 1 {
		return errors.New("error : receiver ciphertext must be of degree 1 to allow key switching")
	}

	switchKeys(evaluator, cIn.value[0], cIn.value[1], cIn.value[1], switchingKey, cOut.value[0], cOut.value[1])

	return nil
}

// RotateColumnsNew rotates the columns of c0 by k position to the left, and returns the result on a newly created element.
// If the provided element is a ciphertext, a keyswitching operation is necessary and a rotation key for the specific rotation needs to be provided.
func (evaluator *Evaluator) RotateColumnsNew(c0 CkksElement, k uint64, evakey *RotationKey) (cOut CkksElement, err error) {

	k &= ((evaluator.ckkscontext.n >> 1) - 1)

	if k == 0 {

		cOut = c0.CopyNew()

		return cOut, nil
	}

	if c0.Degree() > 1 {
		return nil, errors.New("cannot rotate -> input and or output degree not 0 or 1")
	}

	if c0.Degree() > 0 {
		cOut = evaluator.ckkscontext.NewCiphertext(c0.Degree(), c0.Level(), c0.Scale())
	} else {
		cOut = evaluator.ckkscontext.NewPlaintext(c0.Level(), c0.Scale())
	}

	if err = evaluator.RotateColumns(c0, k, evakey, cOut); err != nil {
		return nil, err
	}

	return cOut, nil
}

// RotateColumns rotates the columns of c0 by k position to the left and returns the result on the provided receiver.
// If the provided element is a ciphertext, a keyswitching operation is necessary and a rotation key for the specific rotation needs to be provided.
func (evaluator *Evaluator) RotateColumns(c0 CkksElement, k uint64, evakey *RotationKey, c1 CkksElement) (err error) {

	k &= ((evaluator.ckkscontext.n >> 1) - 1)

	if k == 0 {
		if c0 != c1 {
			if err = c0.Copy(c1); err != nil {
				return err
			}
		}

		return nil
	}

	if c1.Degree() != c0.Degree() {
		return errors.New("cannot rotate -> receiver degree doesn't match input degree ")
	}

	if c0.Degree() > 1 {
		return errors.New("cannot rotate -> input and or output degree not 0 or 1")
	}

	// Looks in the rotationkey if the corresponding rotation has been generated
	if evakey.evakey_rot_col_L[k] != nil {

		if c0.Degree() == 0 {

			if c1 != c0 {

				ring.PermuteNTT(c0.Value()[0], evaluator.ckkscontext.galElRotColLeft[k], c1.Value()[0])

			} else {

				ring.PermuteNTT(c0.Value()[0], evaluator.ckkscontext.galElRotColLeft[k], evaluator.ringpool[0])

				evaluator.ringpool[0].Copy(c1.Value()[0])
			}

			return nil

		} else {

			if c1 != c0 {

				ring.PermuteNTT(c0.Value()[0], evaluator.ckkscontext.galElRotColLeft[k], c1.Value()[0])
				ring.PermuteNTT(c0.Value()[1], evaluator.ckkscontext.galElRotColLeft[k], c1.Value()[1])

			} else {

				ring.PermuteNTT(c0.Value()[0], evaluator.ckkscontext.galElRotColLeft[k], evaluator.ringpool[0])
				ring.PermuteNTT(c0.Value()[1], evaluator.ckkscontext.galElRotColLeft[k], evaluator.ringpool[1])

				evaluator.ringpool[0].Copy(c1.Value()[0])
				evaluator.ringpool[1].Copy(c1.Value()[1])
			}

			switchKeys(evaluator, c1.Value()[0], c1.Value()[1], c1.Value()[1], evakey.evakey_rot_col_L[k], c1.(*Ciphertext).value[0], c1.(*Ciphertext).value[1])

			return nil
		}
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
				rotateColumnsLPow2(evaluator, c0, k, evakey, c1)
			} else {
				rotateColumnsRPow2(evaluator, c0, (evaluator.ckkscontext.n>>1)-k, evakey, c1)
			}

			return nil

			// Else returns an error indicating that the keys have not been generated
		} else {
			return errors.New("cannot rotate -> specific rotation and pow2 rotations have not been generated")
		}
	}
}

func rotateColumnsLPow2(evaluator *Evaluator, c0 CkksElement, k uint64, evakey *RotationKey, c1 CkksElement) {
	rotateColumnsPow2(evaluator, c0, evaluator.ckkscontext.gen, k, evakey.evakey_rot_col_L, c1)
}

func rotateColumnsRPow2(evaluator *Evaluator, c0 CkksElement, k uint64, evakey *RotationKey, c1 CkksElement) {
	rotateColumnsPow2(evaluator, c0, evaluator.ckkscontext.genInv, k, evakey.evakey_rot_col_R, c1)
}

func rotateColumnsPow2(evaluator *Evaluator, c0 CkksElement, generator, k uint64, evakey_rot_col map[uint64]*SwitchingKey, c1 CkksElement) {

	var mask, evakey_index uint64

	mask = (evaluator.ckkscontext.n << 1) - 1

	evakey_index = 1

	c0.Copy(c1)

	for k > 0 {

		if k&1 == 1 {

			if c1.Degree() == 0 {

				ring.PermuteNTT(c1.Value()[0], generator, evaluator.ringpool[0])

				evaluator.ringpool[0].Copy(c1.Value()[0])

			} else {

				ring.PermuteNTT(c1.Value()[0], generator, evaluator.ringpool[0])
				ring.PermuteNTT(c1.Value()[1], generator, evaluator.ringpool[1])

				evaluator.ringpool[0].Copy(c1.Value()[0])
				evaluator.ringpool[1].Copy(c1.Value()[1])

				switchKeys(evaluator, c1.Value()[0], c1.Value()[1], c1.Value()[1], evakey_rot_col[evakey_index], c1.(*Ciphertext).value[0], c1.(*Ciphertext).value[1])
			}
		}

		generator *= generator
		generator &= mask

		evakey_index <<= 1
		k >>= 1
	}
}

// ConjugateNew conjugates c0 (which is equivalement to a row rotation) and returns the result on a newly
// created element. If the provided element is a ciphertext, a keyswitching operation is necessary and a rotation key
// for the row rotation needs to be provided.
func (evaluator *Evaluator) ConjugateNew(c0 CkksElement, evakey *RotationKey) (cOut CkksElement, err error) {

	if c0.Degree() > 1 {
		return nil, errors.New("cannot rotate -> input and or output degree not 0 or 1")
	}

	if c0.Degree() == 1 {
		cOut = evaluator.ckkscontext.NewCiphertext(c0.Degree(), c0.Level(), c0.Scale())
	} else {
		cOut = evaluator.ckkscontext.NewPlaintext(c0.Level(), c0.Scale())
	}

	if err = evaluator.Conjugate(c0, evakey, cOut); err != nil {
		return nil, err
	}

	return cOut, nil

}

// ConjugateNew conjugates c0 (which is equivalement to a row rotation) and returns the result on c1.
// If the provided element is a ciphertext, a keyswitching operation is necessary and a rotation key for the row rotation needs to be provided.
func (evaluator *Evaluator) Conjugate(c0 CkksElement, evakey *RotationKey, c1 CkksElement) error {

	if c0.Degree() > 1 {
		return errors.New("cannot rotate -> input degree not 0 or 1")
	}

	if c1.Degree() != c0.Degree() {
		return errors.New("cannot rotate -> receiver degree doesn't match input degree ")
	}

	if c1.Degree() == 0 {

		cTmp0 := evaluator.ringpool[0]
		ring.PermuteNTT(c0.Value()[0], evaluator.ckkscontext.galElRotRow, cTmp0)
		cTmp0.Copy(c1.Value()[0])

	} else {

		if evakey.evakey_rot_row == nil {
			return errors.New("error : rows rotation key not generated")
		}

		cTmp0 := evaluator.ringpool[0]
		cTmp1 := evaluator.ringpool[1]

		ring.PermuteNTT(c0.Value()[0], evaluator.ckkscontext.galElRotRow, cTmp0)
		ring.PermuteNTT(c0.Value()[1], evaluator.ckkscontext.galElRotRow, cTmp1)

		switchKeys(evaluator, cTmp0, cTmp1, cTmp1, evakey.evakey_rot_row, c1.(*Ciphertext).value[0], c1.(*Ciphertext).value[1])
	}

	return nil
}

// Applies the general keyswitching procedure of the form [c0 + cx*evakey[0], c1 + cx*evakey[1]]
func switchKeys(evaluator *Evaluator, c0, c1, cx *ring.Poly, evakey *SwitchingKey, c0Out, c1Out *ring.Poly) {

	var level, mask, reduce, bitLog uint64

	level = uint64(len(c0Out.Coeffs)) - 1
	context := evaluator.ckkscontext.contextLevel[level]

	c2_qi_w := evaluator.ringpool[5]
	c2 := evaluator.ringpool[4]
	context.InvNTT(cx, c2)

	if c0 != c0Out {
		context.Copy(c0, c0Out)
	}

	if c1 != c1Out {
		context.Copy(c1, c1Out)
	}

	mask = uint64(((1 << evakey.bitDecomp) - 1))

	reduce = 0

	for i := range context.Modulus {

		bitLog = uint64(len(evakey.evakey[i]))

		for j := uint64(0); j < bitLog; j++ {
			//c2_qi_w = (c2_qi_w >> (w*z)) & (w-1)
			for u := uint64(0); u < evaluator.ckkscontext.n; u++ {
				for v := range context.Modulus {
					c2_qi_w.Coeffs[v][u] = (c2.Coeffs[i][u] >> (j * evakey.bitDecomp)) & mask
				}
			}

			context.NTT(c2_qi_w, c2_qi_w)

			context.MulCoeffsMontgomeryAndAddNoMod(evakey.evakey[i][j][0], c2_qi_w, c0Out)
			context.MulCoeffsMontgomeryAndAddNoMod(evakey.evakey[i][j][1], c2_qi_w, c1Out)

			if reduce&7 == 7 {
				context.Reduce(c0Out, c0Out)
				context.Reduce(c1Out, c1Out)
			}

			reduce += 1
		}
	}

	if (reduce-1)&7 != 7 {
		context.Reduce(c0Out, c0Out)
		context.Reduce(c1Out, c1Out)
	}
}
