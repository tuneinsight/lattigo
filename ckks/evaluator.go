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

// Evaluator is a struct that holds the necessary elements to execute the homomorphic operations between Ciphertexts and/or Plaintexts.
// It also holds a memory buffer used to store intermediate computations.
type Evaluator struct {
	parameters Parameters
	*Encoder
	*evaluatorBuffers
	*rlwe.Evaluator
}

// NewEvaluator creates a new Evaluator, that can be used to do homomorphic
// operations on the Ciphertexts and/or Plaintexts. It stores a memory buffer
// and Ciphertexts that will be used for intermediate values.
func NewEvaluator(parameters Parameters, evk rlwe.EvaluationKeySetInterface) *Evaluator {
	return &Evaluator{
		parameters:       parameters,
		Encoder:          NewEncoder(parameters),
		evaluatorBuffers: newEvaluatorBuffers(parameters),
		Evaluator:        rlwe.NewEvaluator(parameters.Parameters, evk),
	}
}

type evaluatorBuffers struct {
	buffQ [3]*ring.Poly // Memory buffer in order: for MForm(c0), MForm(c1), c2
}

// BuffQ returns a pointer to the internal memory buffer buffQ.
func (eval *Evaluator) BuffQ() [3]*ring.Poly {
	return eval.buffQ
}

func newEvaluatorBuffers(parameters Parameters) *evaluatorBuffers {
	buff := new(evaluatorBuffers)
	ringQ := parameters.RingQ()
	buff.buffQ = [3]*ring.Poly{ringQ.NewPoly(), ringQ.NewPoly(), ringQ.NewPoly()}
	return buff
}

// Add adds op1 to op0 and returns the result in op2.
func (eval *Evaluator) Add(op0 *rlwe.Ciphertext, op1 interface{}, op2 *rlwe.Ciphertext) {

	switch op1 := op1.(type) {
	case rlwe.Operand:

		// Checks operand validity and retrieves minimum level
		_, level := eval.CheckBinary(op0.El(), op1.El(), op2.El(), utils.Max(op0.Degree(), op1.Degree()))

		// Generic inplace evaluation
		eval.evaluateInPlace(level, op0, op1.El(), op2, eval.parameters.RingQ().AtLevel(level).Add)

	case complex128, float64, int, int64, uint, uint64, *big.Int, *big.Float, *bignum.Complex:

		// Retrieves minimum level
		level := utils.Min(op0.Level(), op2.Level())

		// Resizes output to minimum level
		op2.Resize(op0.Degree(), level)

		// Convertes the scalar to a complex RNS scalar
		RNSReal, RNSImag := bigComplexToRNSScalar(eval.parameters.RingQ().AtLevel(level), &op0.Scale.Value, bignum.ToComplex(op1, eval.parameters.DefaultPrecision()))

		// Generic inplace evaluation
		eval.evaluateWithScalar(level, op0.Value[:1], RNSReal, RNSImag, op2.Value[:1], eval.parameters.RingQ().AtLevel(level).AddDoubleRNSScalar)

		// Copies the metadata on the output
		op2.MetaData = op0.MetaData

	case []complex128, []float64, []*big.Float, []*bignum.Complex:

		// Retrieves minimum level
		level := utils.Min(op0.Level(), op2.Level())

		// Resizes output to minimum level
		op2.Resize(op0.Degree(), level)

		// Instantiates new plaintext from buffer
		pt := rlwe.NewPlaintextAtLevelFromPoly(level, eval.buffQ[0])
		pt.MetaData = op0.MetaData // Sets the metadata, notably matches scalses

		// Encodes the vector on the plaintext
		if err := eval.Encoder.Encode(op1, pt); err != nil {
			panic(err)
		}

		// Generic in place evaluation
		eval.evaluateInPlace(level, op0, pt.El(), op2, eval.parameters.RingQ().AtLevel(level).Add)
	default:
		panic(fmt.Errorf("op1.(type) must be rlwe.Operand, complex128, float64, int, int64, uint, uint64, *big.Int, *big.Float, *bignum.Complex, []complex128, []float64, []*big.Float or []*bignum.Complex, but is %T", op1))
	}
}

// AddNew adds op1 to op0 and returns the result in a newly created element op2.
func (eval *Evaluator) AddNew(op0 *rlwe.Ciphertext, op1 interface{}) (op2 *rlwe.Ciphertext) {
	op2 = op0.CopyNew()
	eval.Add(op2, op1, op2)
	return
}

// Sub subtracts op1 from op0 and returns the result in op2.
func (eval *Evaluator) Sub(op0 *rlwe.Ciphertext, op1 interface{}, op2 *rlwe.Ciphertext) {

	switch op1 := op1.(type) {
	case rlwe.Operand:

		// Checks operand validity and retrieves minimum level
		_, level := eval.CheckBinary(op0.El(), op1.El(), op2.El(), utils.Max(op0.Degree(), op1.Degree()))

		// Generic inplace evaluation
		eval.evaluateInPlace(level, op0, op1.El(), op2, eval.parameters.RingQ().AtLevel(level).Sub)

		// Negates high degree ciphertext coefficients if the degree of the second operand is larger than the first operand
		if op0.Degree() < op1.Degree() {
			for i := op0.Degree() + 1; i < op1.Degree()+1; i++ {
				eval.parameters.RingQ().AtLevel(level).Neg(op2.Value[i], op2.Value[i])
			}
		}
	case complex128, float64, int, int64, uint, uint64, *big.Int, *big.Float, *bignum.Complex:

		// Retrieves minimum level
		level := utils.Min(op0.Level(), op2.Level())

		// Resizes output to minimum level
		op2.Resize(op0.Degree(), level)

		// Convertes the scalar to a complex RNS scalar
		RNSReal, RNSImag := bigComplexToRNSScalar(eval.parameters.RingQ().AtLevel(level), &op0.Scale.Value, bignum.ToComplex(op1, eval.parameters.DefaultPrecision()))

		// Generic inplace evaluation
		eval.evaluateWithScalar(level, op0.Value[:1], RNSReal, RNSImag, op2.Value[:1], eval.parameters.RingQ().AtLevel(level).SubDoubleRNSScalar)

		// Copies the metadata on the output
		op2.MetaData = op0.MetaData

	case []complex128, []float64, []*big.Float, []*bignum.Complex:

		// Retrieves minimum level
		level := utils.Min(op0.Level(), op2.Level())

		// Resizes output to minimum level
		op2.Resize(op0.Degree(), level)

		// Instantiates new plaintext from buffer
		pt := rlwe.NewPlaintextAtLevelFromPoly(level, eval.buffQ[0])
		pt.MetaData = op0.MetaData

		// Encodes the vector on the plaintext
		if err := eval.Encoder.Encode(op1, pt); err != nil {
			panic(err)
		}

		// Generic inplace evaluation
		eval.evaluateInPlace(level, op0, pt.El(), op2, eval.parameters.RingQ().AtLevel(level).Sub)

	default:
		panic(fmt.Errorf("op1.(type) must be rlwe.Operand, complex128, float64, int, int64, uint, uint64, *big.Int, *big.Float, *bignum.Complex, []complex128, []float64, []*big.Float or []*bignum.Complex, but is %T", op1))
	}
}

// SubNew subtracts op1 from op0 and returns the result in a newly created element op2.
func (eval *Evaluator) SubNew(op0 *rlwe.Ciphertext, op1 interface{}) (op2 *rlwe.Ciphertext) {
	op2 = op0.CopyNew()
	eval.Sub(op2, op1, op2)
	return
}

func (eval *Evaluator) evaluateInPlace(level int, c0 *rlwe.Ciphertext, c1 *rlwe.OperandQ, ctOut *rlwe.Ciphertext, evaluate func(*ring.Poly, *ring.Poly, *ring.Poly)) {

	var tmp0, tmp1 *rlwe.Ciphertext

	maxDegree := utils.Max(c0.Degree(), c1.Degree())
	minDegree := utils.Min(c0.Degree(), c1.Degree())

	// Else resizes the receiver element
	ctOut.El().Resize(maxDegree, ctOut.Level())

	c0Scale := c0.Scale
	c1Scale := c1.Scale

	if ctOut.Level() > level {
		eval.DropLevel(ctOut, ctOut.Level()-utils.Min(c0.Level(), c1.Level()))
	}

	cmp := c0.Scale.Cmp(c1.Scale)

	// Checks whether or not the receiver element is the same as one of the input elements
	// and acts accordingly to avoid unnecessary element creation or element overwriting,
	// and scales properly the element before the evaluation.
	if ctOut == c0 {

		if cmp == 1 {

			ratioFlo := c0Scale.Div(c1Scale).Value

			ratioInt, _ := ratioFlo.Int(nil)

			if ratioInt.Cmp(new(big.Int).SetUint64(0)) == 1 {

				tmp1 = rlwe.NewCiphertextAtLevelFromPoly(level, eval.BuffCt.Value[:c1.Degree()+1])
				tmp1.MetaData = ctOut.MetaData

				eval.Mul(&rlwe.Ciphertext{OperandQ: *c1}, ratioInt, tmp1)
			}

		} else if cmp == -1 {

			ratioFlo := c1Scale.Div(c0Scale).Value

			ratioInt, _ := ratioFlo.Int(nil)

			if ratioInt.Cmp(new(big.Int).SetUint64(0)) == 1 {

				eval.Mul(c0, ratioInt, c0)

				ctOut.Scale = c1.Scale

				tmp1 = &rlwe.Ciphertext{OperandQ: *c1}
			}

		} else {
			tmp1 = &rlwe.Ciphertext{OperandQ: *c1}
		}

		tmp0 = c0

	} else if &ctOut.OperandQ == c1 {

		if cmp == 1 {

			ratioFlo := c0Scale.Div(c1Scale).Value

			ratioInt, _ := ratioFlo.Int(nil)

			if ratioInt.Cmp(new(big.Int).SetUint64(0)) == 1 {
				eval.Mul(&rlwe.Ciphertext{OperandQ: *c1}, ratioInt, ctOut)

				ctOut.Scale = c0.Scale

				tmp0 = c0
			}

		} else if cmp == -1 {

			ratioFlo := c1Scale.Div(c0Scale).Value

			ratioInt, _ := ratioFlo.Int(nil)

			if ratioInt.Cmp(new(big.Int).SetUint64(0)) == 1 {
				// Will avoid resizing on the output
				tmp0 = rlwe.NewCiphertextAtLevelFromPoly(level, eval.BuffCt.Value[:c0.Degree()+1])
				tmp0.MetaData = ctOut.MetaData

				eval.Mul(c0, ratioInt, tmp0)
			}

		} else {
			tmp0 = c0
		}

		tmp1 = &rlwe.Ciphertext{OperandQ: *c1}

	} else {

		if cmp == 1 {

			ratioFlo := c0Scale.Div(c1Scale).Value

			ratioInt, _ := ratioFlo.Int(nil)

			if ratioInt.Cmp(new(big.Int).SetUint64(0)) == 1 {
				// Will avoid resizing on the output
				tmp1 = rlwe.NewCiphertextAtLevelFromPoly(level, eval.BuffCt.Value[:c1.Degree()+1])
				tmp1.MetaData = ctOut.MetaData

				eval.Mul(&rlwe.Ciphertext{OperandQ: *c1}, ratioInt, tmp1)

				tmp0 = c0
			}

		} else if cmp == -1 {

			ratioFlo := c1Scale.Div(c0Scale).Value

			ratioInt, _ := ratioFlo.Int(nil)

			if ratioInt.Cmp(new(big.Int).SetUint64(0)) == 1 {

				tmp0 = rlwe.NewCiphertextAtLevelFromPoly(level, eval.BuffCt.Value[:c0.Degree()+1])
				tmp0.MetaData = ctOut.MetaData

				eval.Mul(c0, ratioInt, tmp0)

				tmp1 = &rlwe.Ciphertext{OperandQ: *c1}

			}

		} else {
			tmp0 = c0
			tmp1 = &rlwe.Ciphertext{OperandQ: *c1}
		}
	}

	for i := 0; i < minDegree+1; i++ {
		evaluate(tmp0.Value[i], tmp1.Value[i], ctOut.El().Value[i])
	}

	scale := c0.Scale.Max(c1.Scale)

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

func (eval *Evaluator) evaluateWithScalar(level int, p0 []*ring.Poly, RNSReal, RNSImag ring.RNSScalar, p1 []*ring.Poly, evaluate func(*ring.Poly, ring.RNSScalar, ring.RNSScalar, *ring.Poly)) {

	// Component wise operation with the following vector:
	// [a + b*psi_qi^2, ....., a + b*psi_qi^2, a - b*psi_qi^2, ...., a - b*psi_qi^2] mod Qi
	// [{                  N/2                }{                N/2               }]
	// Which is equivalent outside of the NTT domain to evaluating a to the first coefficient of ct0 and b to the N/2-th coefficient of ct0.
	for i, s := range eval.parameters.RingQ().SubRings[:level+1] {
		RNSImag[i] = ring.MRed(RNSImag[i], s.RootsForward[1], s.Modulus, s.MRedConstant)
		RNSReal[i], RNSImag[i] = ring.CRed(RNSReal[i]+RNSImag[i], s.Modulus), ring.CRed(RNSReal[i]+s.Modulus-RNSImag[i], s.Modulus)
	}

	for i := range p0 {
		evaluate(p0[i], RNSReal, RNSImag, p1[i])
	}
}

// ScaleUpNew multiplies ct0 by scale and sets its scale to its previous scale times scale returns the result in ctOut.
func (eval *Evaluator) ScaleUpNew(ct0 *rlwe.Ciphertext, scale rlwe.Scale) (ctOut *rlwe.Ciphertext) {
	ctOut = NewCiphertext(eval.parameters, ct0.Degree(), ct0.Level())
	eval.ScaleUp(ct0, scale, ctOut)
	return
}

// ScaleUp multiplies ct0 by scale and sets its scale to its previous scale times scale returns the result in ctOut.
func (eval *Evaluator) ScaleUp(ct0 *rlwe.Ciphertext, scale rlwe.Scale, ctOut *rlwe.Ciphertext) {
	eval.Mul(ct0, scale.Uint64(), ctOut)
	ctOut.MetaData = ct0.MetaData
	ctOut.Scale = ct0.Scale.Mul(scale)
}

// SetScale sets the scale of the ciphertext to the input scale (consumes a level).
func (eval *Evaluator) SetScale(ct *rlwe.Ciphertext, scale rlwe.Scale) {
	ratioFlo := scale.Div(ct.Scale).Value
	eval.Mul(ct, &ratioFlo, ct)
	if err := eval.Rescale(ct, scale, ct); err != nil {
		panic(err)
	}
	ct.Scale = scale
}

// DropLevelNew reduces the level of ct0 by levels and returns the result in a newly created element.
// No rescaling is applied during this procedure.
func (eval *Evaluator) DropLevelNew(ct0 *rlwe.Ciphertext, levels int) (ctOut *rlwe.Ciphertext) {
	ctOut = ct0.CopyNew()
	eval.DropLevel(ctOut, levels)
	return
}

// DropLevel reduces the level of ct0 by levels and returns the result in ct0.
// No rescaling is applied during this procedure.
func (eval *Evaluator) DropLevel(ct0 *rlwe.Ciphertext, levels int) {
	ct0.Resize(ct0.Degree(), ct0.Level()-levels)
}

// RescaleNew divides ct0 by the last modulus in the moduli chain, and repeats this
// procedure (consuming one level each time) until the scale reaches the original scale or before it goes below it, and returns the result
// in a newly created element. Since all the moduli in the moduli chain are generated to be close to the
// original scale, this procedure is equivalent to dividing the input element by the scale and adding
// some error.
// Returns an error if "threshold <= 0", ct.scale = 0, ct.Level() = 0, ct.IsNTT() != true
func (eval *Evaluator) RescaleNew(ct0 *rlwe.Ciphertext, minScale rlwe.Scale) (ctOut *rlwe.Ciphertext, err error) {

	ctOut = NewCiphertext(eval.parameters, ct0.Degree(), ct0.Level())

	return ctOut, eval.Rescale(ct0, minScale, ctOut)
}

// Rescale divides ct0 by the last modulus in the moduli chain, and repeats this
// procedure (consuming one level each time) until the scale reaches the original scale or before it goes below it, and returns the result
// in ctOut. Since all the moduli in the moduli chain are generated to be close to the
// original scale, this procedure is equivalent to dividing the input element by the scale and adding
// some error.
// Returns an error if "minScale <= 0", ct.scale = 0, ct.Level() = 0, ct.IsNTT() != true or if ct.Leve() != ctOut.Level()
func (eval *Evaluator) Rescale(op0 *rlwe.Ciphertext, minScale rlwe.Scale, ctOut *rlwe.Ciphertext) (err error) {

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

	ringQ := eval.parameters.RingQ().AtLevel(op0.Level())

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
func (eval *Evaluator) MulNew(op0 *rlwe.Ciphertext, op1 interface{}) (op2 *rlwe.Ciphertext) {
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
func (eval *Evaluator) Mul(op0 *rlwe.Ciphertext, op1 interface{}, op2 *rlwe.Ciphertext) {
	switch op1 := op1.(type) {
	case rlwe.Operand:

		// Generic in place evaluation
		eval.mulRelin(op0, op1.El(), false, op2)

	case complex128, float64, int, int64, uint, uint64, *big.Int, *big.Float, *bignum.Complex:

		// Retrieves the minimum level
		level := utils.Min(op0.Level(), op2.Level())

		// Resizes output to minimum level
		op2.Resize(op0.Degree(), level)

		// Convertes the scalar to a *bignum.Complex
		cmplxBig := bignum.ToComplex(op1, eval.parameters.DefaultPrecision())

		// Gets the ring at the target level
		ringQ := eval.parameters.RingQ().AtLevel(level)

		var scale rlwe.Scale
		if cmplxBig.IsInt() {
			scale = rlwe.NewScale(1) // Scalar is a GaussianInteger, thus no scaling required
		} else {
			scale = rlwe.NewScale(ringQ.SubRings[level].Modulus) // Current modulus scaling factor

			// If DefaultScalingFactor > 2^60, then multiple moduli are used per single rescale
			// thus continues multiplying the scale with the appropriate number of moduli
			for i := 1; i < eval.parameters.DefaultScaleModuliRatio(); i++ {
				scale = scale.Mul(rlwe.NewScale(ringQ.SubRings[level-i].Modulus))
			}
		}

		// Convertes the *bignum.Complex to a complex RNS scalar
		RNSReal, RNSImag := bigComplexToRNSScalar(ringQ, &scale.Value, cmplxBig)

		// Generic in place evaluation
		eval.evaluateWithScalar(level, op0.Value, RNSReal, RNSImag, op2.Value, ringQ.MulDoubleRNSScalar)

		// Copies the metadata on the output
		op2.MetaData = op0.MetaData
		op2.Scale = op0.Scale.Mul(scale) // updates the scaling factor

	case []complex128, []float64, []*big.Float, []*bignum.Complex:

		// Retrieves minimum level
		level := utils.Min(op0.Level(), op2.Level())

		// Resizes output to minimum level
		op2.Resize(op0.Degree(), level)

		// Gets the ring at the target level
		ringQ := eval.parameters.RingQ().AtLevel(level)

		// Instantiates new plaintext from buffer
		pt := rlwe.NewPlaintextAtLevelFromPoly(level, eval.buffQ[0])
		pt.MetaData = op0.MetaData
		pt.Scale = rlwe.NewScale(ringQ.SubRings[level].Modulus)

		// If DefaultScalingFactor > 2^60, then multiple moduli are used per single rescale
		// thus continues multiplying the scale with the appropriate number of moduli
		for i := 1; i < eval.parameters.DefaultScaleModuliRatio(); i++ {
			pt.Scale = pt.Scale.Mul(rlwe.NewScale(ringQ.SubRings[level-i].Modulus))
		}

		// Encodes the vector on the plaintext
		if err := eval.Encoder.Encode(op1, pt); err != nil {
			panic(err)
		}

		// Generic in place evaluation
		eval.mulRelin(op0, pt.El(), false, op2)
	default:
		panic(fmt.Errorf("op1.(type) must be rlwe.Operand, complex128, float64, int, int64, uint, uint64, *big.Int, *big.Float, *bignum.Complex, []complex128, []float64, []*big.Float or []*bignum.Complex, but is %T", op1))
	}
}

// MulRelinNew multiplies op0 with op1 with relinearization and returns the result in a newly created element.
// The procedure will panic if either op0.Degree or op1.Degree > 1.
// The procedure will panic if the evaluator was not created with an relinearization key.
func (eval *Evaluator) MulRelinNew(op0 *rlwe.Ciphertext, op1 interface{}) (ctOut *rlwe.Ciphertext) {

	switch op1 := op1.(type) {
	case rlwe.Operand:
		ctOut = NewCiphertext(eval.parameters, 1, utils.Min(op0.Level(), op1.Level()))
		eval.mulRelin(op0, op1.El(), true, ctOut)
	default:
		ctOut = NewCiphertext(eval.parameters, 1, op0.Level())
		eval.Mul(op0, op1, ctOut)
	}
	return
}

// MulRelin multiplies op0 with op1 with relinearization and returns the result in ctOut.
// The procedure will panic if either op0.Degree or op1.Degree > 1.
// The procedure will panic if ctOut.Degree != op0.Degree + op1.Degree.
// The procedure will panic if the evaluator was not created with an relinearization key.
func (eval *Evaluator) MulRelin(op0 *rlwe.Ciphertext, op1 interface{}, ctOut *rlwe.Ciphertext) {
	switch op1 := op1.(type) {
	case rlwe.Operand:
		eval.mulRelin(op0, op1.El(), true, ctOut)
	default:
		eval.Mul(op0, op1, ctOut)
	}
}

func (eval *Evaluator) mulRelin(op0 *rlwe.Ciphertext, op1 *rlwe.OperandQ, relin bool, ctOut *rlwe.Ciphertext) {

	if op0.Degree()+op1.Degree() > 2 {
		panic("cannot MulRelin: the sum of the input elements' total degree cannot be larger than 2")
	}

	ctOut.MetaData = op0.MetaData
	ctOut.Scale = op0.Scale.Mul(op1.Scale)

	var c00, c01, c0, c1, c2 *ring.Poly

	// Case Ciphertext (x) Ciphertext
	if op0.Degree() == 1 && op1.Degree() == 1 {

		_, level := eval.CheckBinary(op0.El(), op1.El(), ctOut.El(), ctOut.Degree())

		ringQ := eval.parameters.RingQ().AtLevel(level)

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

		_, level := eval.CheckBinary(op0.El(), op1.El(), ctOut.El(), ctOut.Degree())

		ringQ := eval.parameters.RingQ().AtLevel(level)

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
// If op1.(type) is []complex128, []float64, []*big.Float or []*bignum.Complex:
// - If op2.Scale == op0.Scale, op1 will be encoded and scaled by Q[min(op0.Level(), op2.Level())]
// - If op2.Scale > op0.Scale, op1 will be encoded ans scaled by op2.Scale/op1.Scale.
// Then the method will recurse with op1 given as rlwe.Operand.
//
// If op1.(type) is rlwe.Operand, the multiplication is carried outwithout relinearization and:
//
// This function will panic if op0.Scale > op2.Scale and user must ensure that op2.scale <= op0.scale * op1.scale.
// If op2.scale < op0.scale * op1.scale, then scales up op2 before adding the result.
// Additionally, the procedure will panic if:
// - either op0 or op1 are have a degree higher than 1.
// - op2.Degree != op0.Degree + op1.Degree.
// - op2 = op0 or op1.
func (eval *Evaluator) MulThenAdd(op0 *rlwe.Ciphertext, op1 interface{}, op2 *rlwe.Ciphertext) {
	switch op1 := op1.(type) {
	case rlwe.Operand:

		// Generic in place evaluation
		eval.mulRelinThenAdd(op0, op1.El(), false, op2)
	case complex128, float64, int, int64, uint, uint64, *big.Int, *big.Float, *bignum.Complex:

		// Retrieves the minimum level
		level := utils.Min(op0.Level(), op2.Level())

		// Resizes the output to the minimum level
		op2.Resize(op2.Degree(), level)

		// Gets the ring at the minimum level
		ringQ := eval.parameters.RingQ().AtLevel(level)

		// Convertes the scalar to a *bignum.Complex
		cmplxBig := bignum.ToComplex(op1, eval.parameters.DefaultPrecision())

		var scaleRLWE rlwe.Scale

		// If op0 and op2 scales are identical, but the op1 is not a Gaussian integer then multiplies op2 by scaleRLWE.
		// This ensures noiseless addition with op2 = scaleRLWE * op2 + op0 * round(scalar * scaleRLWE).
		if cmp := op0.Scale.Cmp(op2.Scale); cmp == 0 {

			if cmplxBig.IsInt() {
				scaleRLWE = rlwe.NewScale(1)
			} else {
				scaleRLWE = rlwe.NewScale(ringQ.SubRings[level].Modulus)

				for i := 1; i < eval.parameters.DefaultScaleModuliRatio(); i++ {
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
	case []complex128, []float64, []*big.Float, []*bignum.Complex:

		// Retrieves minimum level
		level := utils.Min(op0.Level(), op2.Level())

		// Resizes output to minimum level
		op2.Resize(op0.Degree(), level)

		// Gets the ring at the target level
		ringQ := eval.parameters.RingQ().AtLevel(level)

		var scaleRLWE rlwe.Scale
		if cmp := op0.Scale.Cmp(op2.Scale); cmp == 0 { // If op0 and op2 scales are identical then multiplies op2 by scaleRLWE.

			scaleRLWE = rlwe.NewScale(ringQ.SubRings[level].Modulus)

			for i := 1; i < eval.parameters.DefaultScaleModuliRatio(); i++ {
				scaleRLWE = scaleRLWE.Mul(rlwe.NewScale(ringQ.SubRings[level-i].Modulus))
			}

			scaleInt := new(big.Int)
			scaleRLWE.Value.Int(scaleInt)
			eval.Mul(op2, scaleInt, op2)
			op2.Scale = op2.Scale.Mul(scaleRLWE)

		} else if cmp == -1 { // op2.Scale > op0.Scale then the scaling factor for op1 becomes the quotient between the two scales
			scaleRLWE = op2.Scale.Div(op0.Scale)
		} else {
			panic("MulThenAdd: op0.Scale > op2.Scale is not supported")
		}

		// Instantiates new plaintext from buffer
		pt := rlwe.NewPlaintextAtLevelFromPoly(level, eval.buffQ[0])
		pt.MetaData = op0.MetaData
		pt.Scale = scaleRLWE

		// Encodes the vector on the plaintext
		if err := eval.Encoder.Encode(op1, pt); err != nil {
			panic(err)
		}

		// Generic in place evaluation
		eval.mulRelinThenAdd(op0, pt.El(), false, op2)

	default:
		panic(fmt.Errorf("op1.(type) must be rlwe.Operand, complex128, float64, int, int64, uint, uint64, *big.Int, *big.Float, *bignum.Complex, []complex128, []float64, []*big.Float or []*bignum.Complex, but is %T", op1))
	}
}

// MulRelinThenAdd multiplies op0 with op1 with relinearization and adds the result on op2.
// User must ensure that op2.scale <= op0.scale * op1.scale.
// If op2.scale < op0.scale * op1.scale, then scales up op2 before adding the result.
// The procedure will panic if either op0.Degree or op1.Degree > 1.
// The procedure will panic if op2.Degree != op0.Degree + op1.Degree.
// The procedure will panic if the evaluator was not created with an relinearization key.
// The procedure will panic if op2 = op0 or op1.
func (eval *Evaluator) MulRelinThenAdd(op0 *rlwe.Ciphertext, op1 rlwe.Operand, op2 *rlwe.Ciphertext) {
	eval.mulRelinThenAdd(op0, op1.El(), true, op2)
}

func (eval *Evaluator) mulRelinThenAdd(op0 *rlwe.Ciphertext, op1 *rlwe.OperandQ, relin bool, op2 *rlwe.Ciphertext) {

	_, level := eval.CheckBinary(op0.El(), op1.El(), op2.El(), utils.Max(op0.Degree(), op1.Degree()))

	if op0.Degree()+op1.Degree() > 2 {
		panic("cannot MulRelinThenAdd: the sum of the input elements' degree cannot be larger than 2")
	}

	if op0.El() == op2.El() || op1.El() == op2.El() {
		panic("cannot MulRelinThenAdd: op2 must be different from op0 and op1")
	}

	resScale := op0.Scale.Mul(op1.Scale)

	if op2.Scale.Cmp(resScale) == -1 {
		ratio := resScale.Div(op2.Scale)
		// Only scales up if int(ratio) >= 2
		if ratio.Float64() >= 2.0 {
			eval.Mul(op2, &ratio.Value, op2)
			op2.Scale = resScale
		}
	}

	ringQ := eval.parameters.RingQ().AtLevel(level)

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
func (eval *Evaluator) RelinearizeNew(ct0 *rlwe.Ciphertext) (ctOut *rlwe.Ciphertext) {
	ctOut = NewCiphertext(eval.parameters, 1, ct0.Level())
	eval.Relinearize(ct0, ctOut)
	return
}

// ApplyEvaluationKeyNew applies the rlwe.EvaluationKey on ct0 and returns the result on a new ciphertext ctOut.
func (eval *Evaluator) ApplyEvaluationKeyNew(ct0 *rlwe.Ciphertext, evk *rlwe.EvaluationKey) (ctOut *rlwe.Ciphertext) {
	ctOut = NewCiphertext(eval.parameters, ct0.Degree(), ct0.Level())
	eval.ApplyEvaluationKey(ct0, evk, ctOut)
	return
}

// RotateNew rotates the columns of ct0 by k positions to the left, and returns the result in a newly created element.
// The method will panic if the evaluator hasn't been given an evaluation key set with the appropriate GaloisKey.
func (eval *Evaluator) RotateNew(ct0 *rlwe.Ciphertext, k int) (ctOut *rlwe.Ciphertext) {
	ctOut = NewCiphertext(eval.parameters, ct0.Degree(), ct0.Level())
	eval.Rotate(ct0, k, ctOut)
	return
}

// Rotate rotates the columns of ct0 by k positions to the left and returns the result in ctOut.
// The method will panic if the evaluator hasn't been given an evaluation key set with the appropriate GaloisKey.
func (eval *Evaluator) Rotate(ct0 *rlwe.Ciphertext, k int, ctOut *rlwe.Ciphertext) {
	eval.Automorphism(ct0, eval.parameters.GaloisElement(k), ctOut)
}

// ConjugateNew conjugates ct0 (which is equivalent to a row rotation) and returns the result in a newly created element.
// The method will panic if the evaluator hasn't been given an evaluation key set with the appropriate GaloisKey.
func (eval *Evaluator) ConjugateNew(ct0 *rlwe.Ciphertext) (ctOut *rlwe.Ciphertext) {

	if eval.parameters.RingType() == ring.ConjugateInvariant {
		panic("cannot ConjugateNew: method is not supported when parameters.RingType() == ring.ConjugateInvariant")
	}

	ctOut = NewCiphertext(eval.parameters, ct0.Degree(), ct0.Level())
	eval.Conjugate(ct0, ctOut)
	return
}

// Conjugate conjugates ct0 (which is equivalent to a row rotation) and returns the result in ctOut.
// The method will panic if the evaluator hasn't been given an evaluation key set with the appropriate GaloisKey.
func (eval *Evaluator) Conjugate(ct0 *rlwe.Ciphertext, ctOut *rlwe.Ciphertext) {

	if eval.parameters.RingType() == ring.ConjugateInvariant {
		panic("cannot Conjugate: method is not supported when parameters.RingType() == ring.ConjugateInvariant")
	}

	eval.Automorphism(ct0, eval.parameters.GaloisElementInverse(), ctOut)
}

// RotateHoistedNew takes an input Ciphertext and a list of rotations and returns a map of Ciphertext, where each element of the map is the input Ciphertext
// rotation by one element of the list. It is much faster than sequential calls to Rotate.
func (eval *Evaluator) RotateHoistedNew(ctIn *rlwe.Ciphertext, rotations []int) (ctOut map[int]*rlwe.Ciphertext) {
	ctOut = make(map[int]*rlwe.Ciphertext)
	for _, i := range rotations {
		ctOut[i] = NewCiphertext(eval.parameters, 1, ctIn.Level())
	}
	eval.RotateHoisted(ctIn, rotations, ctOut)
	return
}

// RotateHoisted takes an input Ciphertext and a list of rotations and populates a map of pre-allocated Ciphertexts,
// where each element of the map is the input Ciphertext rotation by one element of the list.
// It is much faster than sequential calls to Rotate.
func (eval *Evaluator) RotateHoisted(ctIn *rlwe.Ciphertext, rotations []int, ctOut map[int]*rlwe.Ciphertext) {
	levelQ := ctIn.Level()
	eval.DecomposeNTT(levelQ, eval.parameters.MaxLevelP(), eval.parameters.PCount(), ctIn.Value[1], ctIn.IsNTT, eval.BuffDecompQP)
	for _, i := range rotations {
		eval.AutomorphismHoisted(levelQ, ctIn, eval.BuffDecompQP, eval.parameters.GaloisElement(i), ctOut[i])
	}
}

func (eval *Evaluator) RotateHoistedLazyNew(level int, rotations []int, ct *rlwe.Ciphertext, c2DecompQP []ringqp.Poly) (cOut map[int]*rlwe.OperandQP) {
	cOut = make(map[int]*rlwe.OperandQP)
	for _, i := range rotations {
		if i != 0 {
			cOut[i] = rlwe.NewOperandQP(eval.parameters.Parameters, 1, level, eval.parameters.MaxLevelP())
			eval.AutomorphismHoistedLazy(level, ct, c2DecompQP, eval.parameters.GaloisElement(i), cOut[i])
		}
	}

	return
}

// Parameters returns the Parametrs of the underlying struct as an rlwe.ParametersInterface.
func (eval *Evaluator) Parameters() rlwe.ParametersInterface {
	return eval.parameters
}

// ShallowCopy creates a shallow copy of this evaluator in which all the read-only data-structures are
// shared with the receiver and the temporary buffers are reallocated. The receiver and the returned
// Evaluators can be used concurrently.
func (eval *Evaluator) ShallowCopy() *Evaluator {
	return &Evaluator{
		parameters:       eval.parameters,
		Encoder:          NewEncoder(eval.parameters),
		Evaluator:        eval.Evaluator.ShallowCopy(),
		evaluatorBuffers: newEvaluatorBuffers(eval.parameters),
	}
}

// WithKey creates a shallow copy of the receiver Evaluator for which the new EvaluationKey is evaluationKey
// and where the temporary buffers are shared. The receiver and the returned Evaluators cannot be used concurrently.
func (eval *Evaluator) WithKey(evk rlwe.EvaluationKeySetInterface) *Evaluator {
	return &Evaluator{
		Evaluator:        eval.Evaluator.WithKey(evk),
		parameters:       eval.parameters,
		Encoder:          eval.Encoder,
		evaluatorBuffers: eval.evaluatorBuffers,
	}
}
