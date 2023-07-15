package ckks

import (
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
func NewEvaluator(parameters Parameters, evk rlwe.EvaluationKeySet) *Evaluator {
	return &Evaluator{
		parameters:       parameters,
		Encoder:          NewEncoder(parameters),
		evaluatorBuffers: newEvaluatorBuffers(parameters),
		Evaluator:        rlwe.NewEvaluator(parameters.Parameters, evk),
	}
}

type evaluatorBuffers struct {
	buffQ [3]ring.Poly // Memory buffer in order: for MForm(c0), MForm(c1), c2
}

// BuffQ returns a pointer to the internal memory buffer buffQ.
func (eval Evaluator) BuffQ() [3]ring.Poly {
	return eval.buffQ
}

func newEvaluatorBuffers(parameters Parameters) *evaluatorBuffers {
	buff := new(evaluatorBuffers)
	ringQ := parameters.RingQ()
	buff.buffQ = [3]ring.Poly{ringQ.NewPoly(), ringQ.NewPoly(), ringQ.NewPoly()}
	return buff
}

// Add adds op1 to op0 and returns the result in opOut.
// The following types are accepted for op1:
// - rlwe.OperandInterface[ring.Poly]
// - complex128, float64, int, int64, uint, uint64, *big.Int, *big.Float, *bignum.Complex
// - []complex128, []float64, []*big.Float or []*bignum.Complex
// Passing an invalid type will return an error.
func (eval Evaluator) Add(op0 *rlwe.Ciphertext, op1 interface{}, opOut *rlwe.Ciphertext) (err error) {

	switch op1 := op1.(type) {
	case rlwe.OperandInterface[ring.Poly]:

		// Checks operand validity and retrieves minimum level
		_, level, err := eval.InitOutputBinaryOp(op0.El(), op1.El(), utils.Max(op0.Degree(), op1.Degree()), opOut.El())

		if err != nil {
			return err
		}

		// Generic inplace evaluation
		eval.evaluateInPlace(level, op0, op1.El(), opOut, eval.parameters.RingQ().AtLevel(level).Add)

	case complex128, float64, int, int64, uint, uint64, *big.Int, *big.Float, *bignum.Complex:

		// Retrieves minimum level
		level := utils.Min(op0.Level(), opOut.Level())

		// Resizes output to minimum level
		opOut.Resize(op0.Degree(), level)

		// Convertes the scalar to a complex RNS scalar
		RNSReal, RNSImag := bigComplexToRNSScalar(eval.parameters.RingQ().AtLevel(level), &op0.PlaintextScale.Value, bignum.ToComplex(op1, eval.parameters.PlaintextPrecision()))

		// Generic inplace evaluation
		eval.evaluateWithScalar(level, op0.Value[:1], RNSReal, RNSImag, opOut.Value[:1], eval.parameters.RingQ().AtLevel(level).AddDoubleRNSScalar)

		if op0 != opOut {
			for i := 1; i < len(opOut.Value); i++ {
				copy(opOut.Value[i].Buff, op0.Value[i].Buff) // Resize step ensures identical size
			}
		}

		// Copies the metadata on the output
		*opOut.MetaData = *op0.MetaData

	case []complex128, []float64, []*big.Float, []*bignum.Complex:

		// Retrieves minimum level
		level := utils.Min(op0.Level(), opOut.Level())

		// Resizes output to minimum level
		opOut.Resize(op0.Degree(), level)

		// Instantiates new plaintext from buffer
		pt, err := rlwe.NewPlaintextAtLevelFromPoly(level, eval.buffQ[0])
		if err != nil {
			panic(err)
		}
		*pt.MetaData = *op0.MetaData // Sets the metadata, notably matches scales

		// Encodes the vector on the plaintext
		if err := eval.Encoder.Encode(op1, pt); err != nil {
			return err
		}

		// Generic in place evaluation
		eval.evaluateInPlace(level, op0, pt.El(), opOut, eval.parameters.RingQ().AtLevel(level).Add)
	default:
		return fmt.Errorf("invalid op1.(type): must be rlwe.OperandInterface[ring.Poly], complex128, float64, int, int64, uint, uint64, *big.Int, *big.Float, *bignum.Complex, []complex128, []float64, []*big.Float or []*bignum.Complex, but is %T", op1)
	}

	return
}

// AddNew adds op1 to op0 and returns the result in a newly created element opOut.
// The following types are accepted for op1:
// - rlwe.OperandInterface[ring.Poly]
// - complex128, float64, int, int64, uint, uint64, *big.Int, *big.Float, *bignum.Complex
// - []complex128, []float64, []*big.Float or []*bignum.Complex
// Passing an invalid type will return an error.
func (eval Evaluator) AddNew(op0 *rlwe.Ciphertext, op1 interface{}) (opOut *rlwe.Ciphertext, err error) {
	opOut = NewCiphertext(eval.parameters, op0.Degree(), op0.Level())
	return opOut, eval.Add(op0, op1, opOut)
}

// Sub subtracts op1 from op0 and returns the result in opOut.
// The following types are accepted for op1:
// - rlwe.OperandInterface[ring.Poly]
// - complex128, float64, int, int64, uint, uint64, *big.Int, *big.Float, *bignum.Complex
// - []complex128, []float64, []*big.Float or []*bignum.Complex
// Passing an invalid type will return an error.
func (eval Evaluator) Sub(op0 *rlwe.Ciphertext, op1 interface{}, opOut *rlwe.Ciphertext) (err error) {

	switch op1 := op1.(type) {
	case rlwe.OperandInterface[ring.Poly]:

		// Checks operand validity and retrieves minimum level
		_, level, err := eval.InitOutputBinaryOp(op0.El(), op1.El(), utils.Max(op0.Degree(), op1.Degree()), opOut.El())

		if err != nil {
			return err
		}

		// Generic inplace evaluation
		eval.evaluateInPlace(level, op0, op1.El(), opOut, eval.parameters.RingQ().AtLevel(level).Sub)

		// Negates high degree ciphertext coefficients if the degree of the second operand is larger than the first operand
		if op0.Degree() < op1.Degree() {
			for i := op0.Degree() + 1; i < op1.Degree()+1; i++ {
				eval.parameters.RingQ().AtLevel(level).Neg(opOut.Value[i], opOut.Value[i])
			}
		}
	case complex128, float64, int, int64, uint, uint64, *big.Int, *big.Float, *bignum.Complex:

		// Retrieves minimum level
		level := utils.Min(op0.Level(), opOut.Level())

		// Resizes output to minimum level
		opOut.Resize(op0.Degree(), level)

		// Convertes the scalar to a complex RNS scalar
		RNSReal, RNSImag := bigComplexToRNSScalar(eval.parameters.RingQ().AtLevel(level), &op0.PlaintextScale.Value, bignum.ToComplex(op1, eval.parameters.PlaintextPrecision()))

		// Generic inplace evaluation
		eval.evaluateWithScalar(level, op0.Value[:1], RNSReal, RNSImag, opOut.Value[:1], eval.parameters.RingQ().AtLevel(level).SubDoubleRNSScalar)

		if op0 != opOut {
			for i := 1; i < len(opOut.Value); i++ {
				copy(opOut.Value[i].Buff, op0.Value[i].Buff) // Resize step ensures identical size
			}
		}

		// Copies the metadata on the output
		*opOut.MetaData = *op0.MetaData

	case []complex128, []float64, []*big.Float, []*bignum.Complex:

		// Retrieves minimum level
		level := utils.Min(op0.Level(), opOut.Level())

		// Resizes output to minimum level
		opOut.Resize(op0.Degree(), level)

		// Instantiates new plaintext from buffer
		pt, err := rlwe.NewPlaintextAtLevelFromPoly(level, eval.buffQ[0])
		if err != nil {
			panic(err)
		}
		*pt.MetaData = *op0.MetaData

		// Encodes the vector on the plaintext
		if err := eval.Encoder.Encode(op1, pt); err != nil {
			return err
		}

		// Generic inplace evaluation
		eval.evaluateInPlace(level, op0, pt.El(), opOut, eval.parameters.RingQ().AtLevel(level).Sub)

	default:
		return fmt.Errorf("invalid op1.(type): must be rlwe.OperandInterface[ring.Poly], complex128, float64, int, int64, uint, uint64, *big.Int, *big.Float, *bignum.Complex, []complex128, []float64, []*big.Float or []*bignum.Complex, but is %T", op1)
	}

	return
}

// SubNew subtracts op1 from op0 and returns the result in a newly created element opOut.
// The following types are accepted for op1:
// - rlwe.OperandInterface[ring.Poly]
// - complex128, float64, int, int64, uint, uint64, *big.Int, *big.Float, *bignum.Complex
// - []complex128, []float64, []*big.Float or []*bignum.Complex
// Passing an invalid type will return an error.
func (eval Evaluator) SubNew(op0 *rlwe.Ciphertext, op1 interface{}) (opOut *rlwe.Ciphertext, err error) {
	opOut = NewCiphertext(eval.parameters, op0.Degree(), op0.Level())
	return opOut, eval.Sub(op0, op1, opOut)
}

func (eval Evaluator) evaluateInPlace(level int, c0 *rlwe.Ciphertext, c1 *rlwe.Operand[ring.Poly], opOut *rlwe.Ciphertext, evaluate func(ring.Poly, ring.Poly, ring.Poly)) {

	var tmp0, tmp1 *rlwe.Ciphertext

	maxDegree := utils.Max(c0.Degree(), c1.Degree())
	minDegree := utils.Min(c0.Degree(), c1.Degree())

	// Else resizes the receiver element
	opOut.El().Resize(maxDegree, opOut.Level())

	c0Scale := c0.PlaintextScale
	c1Scale := c1.PlaintextScale

	if opOut.Level() > level {
		eval.DropLevel(opOut, opOut.Level()-utils.Min(c0.Level(), c1.Level()))
	}

	cmp := c0.PlaintextScale.Cmp(c1.PlaintextScale)

	var err error

	// Checks whether or not the receiver element is the same as one of the input elements
	// and acts accordingly to avoid unnecessary element creation or element overwriting,
	// and scales properly the element before the evaluation.
	if opOut == c0 {

		if cmp == 1 {

			ratioFlo := c0Scale.Div(c1Scale).Value

			ratioInt, _ := ratioFlo.Int(nil)

			if ratioInt.Cmp(new(big.Int).SetUint64(0)) == 1 {

				tmp1, err = rlwe.NewCiphertextAtLevelFromPoly(level, eval.BuffCt.Value[:c1.Degree()+1])
				if err != nil {
					panic(err)
				}
				*tmp1.MetaData = *opOut.MetaData

				if err = eval.Mul(&rlwe.Ciphertext{Operand: *c1}, ratioInt, tmp1); err != nil {
					return
				}
			}

		} else if cmp == -1 {

			ratioFlo := c1Scale.Div(c0Scale).Value

			ratioInt, _ := ratioFlo.Int(nil)

			if ratioInt.Cmp(new(big.Int).SetUint64(0)) == 1 {

				if err = eval.Mul(c0, ratioInt, c0); err != nil {
					return
				}

				opOut.PlaintextScale = c1.PlaintextScale

				tmp1 = &rlwe.Ciphertext{Operand: *c1}
			}

		} else {
			tmp1 = &rlwe.Ciphertext{Operand: *c1}
		}

		tmp0 = c0

	} else if &opOut.Operand == c1 {

		if cmp == 1 {

			ratioFlo := c0Scale.Div(c1Scale).Value

			ratioInt, _ := ratioFlo.Int(nil)

			if ratioInt.Cmp(new(big.Int).SetUint64(0)) == 1 {
				if err = eval.Mul(&rlwe.Ciphertext{Operand: *c1}, ratioInt, opOut); err != nil {
					return
				}

				opOut.PlaintextScale = c0.PlaintextScale

				tmp0 = c0
			}

		} else if cmp == -1 {

			ratioFlo := c1Scale.Div(c0Scale).Value

			ratioInt, _ := ratioFlo.Int(nil)

			if ratioInt.Cmp(new(big.Int).SetUint64(0)) == 1 {
				// Will avoid resizing on the output
				tmp0, err = rlwe.NewCiphertextAtLevelFromPoly(level, eval.BuffCt.Value[:c0.Degree()+1])
				if err != nil {
					panic(err)
				}
				*tmp0.MetaData = *opOut.MetaData

				if err = eval.Mul(c0, ratioInt, tmp0); err != nil {
					return
				}
			}

		} else {
			tmp0 = c0
		}

		tmp1 = &rlwe.Ciphertext{Operand: *c1}

	} else {

		if cmp == 1 {

			ratioFlo := c0Scale.Div(c1Scale).Value

			ratioInt, _ := ratioFlo.Int(nil)

			if ratioInt.Cmp(new(big.Int).SetUint64(0)) == 1 {
				// Will avoid resizing on the output
				tmp1, err = rlwe.NewCiphertextAtLevelFromPoly(level, eval.BuffCt.Value[:c1.Degree()+1])
				if err != nil {
					panic(err)
				}
				*tmp1.MetaData = *opOut.MetaData

				if err = eval.Mul(&rlwe.Ciphertext{Operand: *c1}, ratioInt, tmp1); err != nil {
					return
				}

				tmp0 = c0
			}

		} else if cmp == -1 {

			ratioFlo := c1Scale.Div(c0Scale).Value

			ratioInt, _ := ratioFlo.Int(nil)

			if ratioInt.Cmp(new(big.Int).SetUint64(0)) == 1 {

				tmp0, err = rlwe.NewCiphertextAtLevelFromPoly(level, eval.BuffCt.Value[:c0.Degree()+1])
				if err != nil {
					panic(err)
				}
				*tmp0.MetaData = *opOut.MetaData

				if err = eval.Mul(c0, ratioInt, tmp0); err != nil {
					return
				}

				tmp1 = &rlwe.Ciphertext{Operand: *c1}

			}

		} else {
			tmp0 = c0
			tmp1 = &rlwe.Ciphertext{Operand: *c1}
		}
	}

	for i := 0; i < minDegree+1; i++ {
		evaluate(tmp0.Value[i], tmp1.Value[i], opOut.El().Value[i])
	}

	scale := c0.PlaintextScale.Max(c1.PlaintextScale)

	*opOut.MetaData = *c0.MetaData
	opOut.PlaintextScale = scale

	// If the inputs degrees differ, it copies the remaining degree on the receiver.
	// Also checks that the receiver is not one of the inputs to avoid unnecessary work.

	if c0.Degree() > c1.Degree() && &tmp0.Operand != opOut.El() {
		for i := minDegree + 1; i < maxDegree+1; i++ {
			ring.Copy(tmp0.Value[i], opOut.El().Value[i])
		}
	} else if c1.Degree() > c0.Degree() && &tmp1.Operand != opOut.El() {
		for i := minDegree + 1; i < maxDegree+1; i++ {
			ring.Copy(tmp1.Value[i], opOut.El().Value[i])
		}
	}
}

func (eval Evaluator) evaluateWithScalar(level int, p0 []ring.Poly, RNSReal, RNSImag ring.RNSScalar, p1 []ring.Poly, evaluate func(ring.Poly, ring.RNSScalar, ring.RNSScalar, ring.Poly)) {

	// Component wise operation with the following vector:
	// [a + b*psi_qi^2, ....., a + b*psi_qi^2, a - b*psi_qi^2, ...., a - b*psi_qi^2] mod Qi
	// [{                  N/2                }{                N/2               }]
	// Which is equivalent outside of the NTT domain to evaluating a to the first coefficient of op0 and b to the N/2-th coefficient of op0.
	for i, s := range eval.parameters.RingQ().SubRings[:level+1] {
		RNSImag[i] = ring.MRed(RNSImag[i], s.RootsForward[1], s.Modulus, s.MRedConstant)
		RNSReal[i], RNSImag[i] = ring.CRed(RNSReal[i]+RNSImag[i], s.Modulus), ring.CRed(RNSReal[i]+s.Modulus-RNSImag[i], s.Modulus)
	}

	for i := range p0 {
		evaluate(p0[i], RNSReal, RNSImag, p1[i])
	}
}

// ScaleUpNew multiplies op0 by scale and sets its scale to its previous scale times scale returns the result in opOut.
func (eval Evaluator) ScaleUpNew(op0 *rlwe.Ciphertext, scale rlwe.Scale) (opOut *rlwe.Ciphertext, err error) {
	opOut = NewCiphertext(eval.parameters, op0.Degree(), op0.Level())
	return opOut, eval.ScaleUp(op0, scale, opOut)
}

// ScaleUp multiplies op0 by scale and sets its scale to its previous scale times scale returns the result in opOut.
func (eval Evaluator) ScaleUp(op0 *rlwe.Ciphertext, scale rlwe.Scale, opOut *rlwe.Ciphertext) (err error) {
	if err = eval.Mul(op0, scale.Uint64(), opOut); err != nil {
		return fmt.Errorf("cannot ScaleUp: %w", err)
	}
	*opOut.MetaData = *op0.MetaData
	opOut.PlaintextScale = op0.PlaintextScale.Mul(scale)

	return
}

// SetScale sets the scale of the ciphertext to the input scale (consumes a level).
func (eval Evaluator) SetScale(ct *rlwe.Ciphertext, scale rlwe.Scale) (err error) {
	ratioFlo := scale.Div(ct.PlaintextScale).Value
	if err = eval.Mul(ct, &ratioFlo, ct); err != nil {
		return fmt.Errorf("cannot SetScale: %w", err)
	}
	if err = eval.Rescale(ct, scale, ct); err != nil {
		return fmt.Errorf("cannot SetScale: %w", err)
	}
	ct.PlaintextScale = scale
	return
}

// DropLevelNew reduces the level of op0 by levels and returns the result in a newly created element.
// No rescaling is applied during this procedure.
func (eval Evaluator) DropLevelNew(op0 *rlwe.Ciphertext, levels int) (opOut *rlwe.Ciphertext) {
	opOut = op0.CopyNew()
	eval.DropLevel(opOut, levels)
	return
}

// DropLevel reduces the level of op0 by levels and returns the result in op0.
// No rescaling is applied during this procedure.
func (eval Evaluator) DropLevel(op0 *rlwe.Ciphertext, levels int) {
	op0.Resize(op0.Degree(), op0.Level()-levels)
}

// RescaleNew divides op0 by the last modulus in the moduli chain, and repeats this
// procedure (consuming one level each time) until the scale reaches the original scale or before it goes below it, and returns the result
// in a newly created element. Since all the moduli in the moduli chain are generated to be close to the
// original scale, this procedure is equivalent to dividing the input element by the scale and adding
// some error.
// Returns an error if "threshold <= 0", ct.PlaintextScale = 0, ct.Level() = 0, ct.IsNTT() != true
func (eval Evaluator) RescaleNew(op0 *rlwe.Ciphertext, minScale rlwe.Scale) (opOut *rlwe.Ciphertext, err error) {
	opOut = NewCiphertext(eval.parameters, op0.Degree(), op0.Level())
	return opOut, eval.Rescale(op0, minScale, opOut)
}

// Rescale divides op0 by the last modulus in the moduli chain, and repeats this
// procedure (consuming one level each time) until the scale reaches the original scale or before it goes below it, and returns the result
// in opOut. Since all the moduli in the moduli chain are generated to be close to the
// original scale, this procedure is equivalent to dividing the input element by the scale and adding
// some error.
// Returns an error if "minScale <= 0", ct.PlaintextScale = 0, ct.Level() = 0, ct.IsNTT() != true or if ct.Leve() != opOut.Level()
func (eval Evaluator) Rescale(op0 *rlwe.Ciphertext, minScale rlwe.Scale, opOut *rlwe.Ciphertext) (err error) {

	if minScale.Cmp(rlwe.NewScale(0)) != 1 {
		return fmt.Errorf("cannot Rescale: minScale is <0")
	}

	minScale = minScale.Div(rlwe.NewScale(2))

	if op0.PlaintextScale.Cmp(rlwe.NewScale(0)) != 1 {
		return fmt.Errorf("cannot Rescale: ciphertext scale is <0")
	}

	if op0.Level() == 0 {
		return fmt.Errorf("cannot Rescale: input Ciphertext already at level 0")
	}

	if opOut.Degree() != op0.Degree() {
		return fmt.Errorf("cannot Rescale: op0.Degree() != opOut.Degree()")
	}

	*opOut.MetaData = *op0.MetaData

	newLevel := op0.Level()

	ringQ := eval.parameters.RingQ().AtLevel(op0.Level())

	// Divides the scale by each moduli of the modulus chain as long as the scale isn't smaller than minScale/2
	// or until the output Level() would be zero
	var nbRescales int
	for newLevel >= 0 {

		scale := opOut.PlaintextScale.Div(rlwe.NewScale(ringQ.SubRings[newLevel].Modulus))

		if scale.Cmp(minScale) == -1 {
			break
		}

		opOut.PlaintextScale = scale

		nbRescales++
		newLevel--
	}

	if nbRescales > 0 {
		for i := range opOut.Value {
			ringQ.DivRoundByLastModulusManyNTT(nbRescales, op0.Value[i], eval.buffQ[0], opOut.Value[i])
		}
		opOut.Resize(opOut.Degree(), newLevel)
	} else {
		if op0 != opOut {
			opOut.Copy(op0)
		}
	}

	return nil
}

// MulNew multiplies op0 with op1 without relinearization and returns the result in a newly created element opOut.
//
// op1.(type) can be rlwe.OperandInterface[ring.Poly], complex128, float64, int, int64, uint64. *big.Float, *big.Int or *ring.Complex.
//
// If op1.(type) == rlwe.OperandInterface[ring.Poly]:
// - The procedure will return an error if either op0.Degree or op1.Degree > 1.
func (eval Evaluator) MulNew(op0 *rlwe.Ciphertext, op1 interface{}) (opOut *rlwe.Ciphertext, err error) {
	opOut = NewCiphertext(eval.parameters, op0.Degree(), op0.Level())
	return opOut, eval.Mul(op0, op1, opOut)
}

// Mul multiplies op0 with op1 without relinearization and returns the result in opOut.
//
// The following types are accepted for op1:
// - rlwe.OperandInterface[ring.Poly]
// - complex128, float64, int, int64, uint, uint64, *big.Int, *big.Float, *bignum.Complex
// - []complex128, []float64, []*big.Float or []*bignum.Complex
// Passing an invalid type will return an error.
//
// If op1.(type) == rlwe.OperandInterface[ring.Poly]:
// - The procedure will return an error if either op0 or op1 are have a degree higher than 1.
// - The procedure will return an error if opOut.Degree != op0.Degree + op1.Degree.
func (eval Evaluator) Mul(op0 *rlwe.Ciphertext, op1 interface{}, opOut *rlwe.Ciphertext) (err error) {
	switch op1 := op1.(type) {
	case rlwe.OperandInterface[ring.Poly]:

		// Generic in place evaluation
		return eval.mulRelin(op0, op1.El(), false, opOut)

	case complex128, float64, int, int64, uint, uint64, *big.Int, *big.Float, *bignum.Complex:

		// Retrieves the minimum level
		level := utils.Min(op0.Level(), opOut.Level())

		// Resizes output to minimum level
		opOut.Resize(op0.Degree(), level)

		// Convertes the scalar to a *bignum.Complex
		cmplxBig := bignum.ToComplex(op1, eval.parameters.PlaintextPrecision())

		// Gets the ring at the target level
		ringQ := eval.parameters.RingQ().AtLevel(level)

		var scale rlwe.Scale
		if cmplxBig.IsInt() {
			scale = rlwe.NewScale(1) // Scalar is a GaussianInteger, thus no scaling required
		} else {
			scale = rlwe.NewScale(ringQ.SubRings[level].Modulus) // Current modulus scaling factor

			// If DefaultScalingFactor > 2^60, then multiple moduli are used per single rescale
			// thus continues multiplying the scale with the appropriate number of moduli
			for i := 1; i < eval.parameters.PlaintextScaleToModuliRatio(); i++ {
				scale = scale.Mul(rlwe.NewScale(ringQ.SubRings[level-i].Modulus))
			}
		}

		// Convertes the *bignum.Complex to a complex RNS scalar
		RNSReal, RNSImag := bigComplexToRNSScalar(ringQ, &scale.Value, cmplxBig)

		// Generic in place evaluation
		eval.evaluateWithScalar(level, op0.Value, RNSReal, RNSImag, opOut.Value, ringQ.MulDoubleRNSScalar)

		// Copies the metadata on the output
		*opOut.MetaData = *op0.MetaData
		opOut.PlaintextScale = op0.PlaintextScale.Mul(scale) // updates the scaling factor

		return nil

	case []complex128, []float64, []*big.Float, []*bignum.Complex:

		// Retrieves minimum level
		level := utils.Min(op0.Level(), opOut.Level())

		// Resizes output to minimum level
		opOut.Resize(op0.Degree(), level)

		// Gets the ring at the target level
		ringQ := eval.parameters.RingQ().AtLevel(level)

		// Instantiates new plaintext from buffer
		pt, err := rlwe.NewPlaintextAtLevelFromPoly(level, eval.buffQ[0])
		if err != nil {
			panic(err)
		}
		*pt.MetaData = *op0.MetaData
		pt.PlaintextScale = rlwe.NewScale(ringQ.SubRings[level].Modulus)

		// If DefaultScalingFactor > 2^60, then multiple moduli are used per single rescale
		// thus continues multiplying the scale with the appropriate number of moduli
		for i := 1; i < eval.parameters.PlaintextScaleToModuliRatio(); i++ {
			pt.PlaintextScale = pt.PlaintextScale.Mul(rlwe.NewScale(ringQ.SubRings[level-i].Modulus))
		}

		// Encodes the vector on the plaintext
		if err = eval.Encoder.Encode(op1, pt); err != nil {
			return err
		}

		// Generic in place evaluation
		return eval.mulRelin(op0, pt.El(), false, opOut)
	default:
		return fmt.Errorf("op1.(type) must be rlwe.OperandInterface[ring.Poly], complex128, float64, int, int64, uint, uint64, *big.Int, *big.Float, *bignum.Complex, []complex128, []float64, []*big.Float or []*bignum.Complex, but is %T", op1)
	}
}

// MulRelinNew multiplies op0 with op1 with relinearization and returns the result in a newly created element.
//
// The following types are accepted for op1:
// - rlwe.OperandInterface[ring.Poly]
// - complex128, float64, int, int64, uint, uint64, *big.Int, *big.Float, *bignum.Complex
// - []complex128, []float64, []*big.Float or []*bignum.Complex
// Passing an invalid type will return an error.
//
// The procedure will return an error if either op0.Degree or op1.Degree > 1.
// The procedure will return an error if the evaluator was not created with an relinearization key.
func (eval Evaluator) MulRelinNew(op0 *rlwe.Ciphertext, op1 interface{}) (opOut *rlwe.Ciphertext, err error) {
	switch op1 := op1.(type) {
	case rlwe.OperandInterface[ring.Poly]:
		opOut = NewCiphertext(eval.parameters, 1, utils.Min(op0.Level(), op1.Level()))
		return opOut, eval.mulRelin(op0, op1.El(), true, opOut)
	default:
		opOut = NewCiphertext(eval.parameters, 1, op0.Level())
		return opOut, eval.Mul(op0, op1, opOut)
	}
}

// MulRelin multiplies op0 with op1 with relinearization and returns the result in opOut.
//
// The following types are accepted for op1:
// - rlwe.OperandInterface[ring.Poly]
// - complex128, float64, int, int64, uint, uint64, *big.Int, *big.Float, *bignum.Complex
// - []complex128, []float64, []*big.Float or []*bignum.Complex
// Passing an invalid type will return an error.
//
// The procedure will return an error if either op0.Degree or op1.Degree > 1.
// The procedure will return an error if opOut.Degree != op0.Degree + op1.Degree.
// The procedure will return an error if the evaluator was not created with an relinearization key.
func (eval Evaluator) MulRelin(op0 *rlwe.Ciphertext, op1 interface{}, opOut *rlwe.Ciphertext) (err error) {
	switch op1 := op1.(type) {
	case rlwe.OperandInterface[ring.Poly]:
		return eval.mulRelin(op0, op1.El(), true, opOut)
	default:
		return eval.Mul(op0, op1, opOut)
	}
}

func (eval Evaluator) mulRelin(op0 *rlwe.Ciphertext, op1 *rlwe.Operand[ring.Poly], relin bool, opOut *rlwe.Ciphertext) (err error) {

	if op0.Degree()+op1.Degree() > 2 {
		return fmt.Errorf("cannot MulRelin: the sum of the input elements' total degree cannot be larger than 2")
	}

	*opOut.MetaData = *op0.MetaData
	opOut.PlaintextScale = op0.PlaintextScale.Mul(op1.PlaintextScale)

	var c00, c01, c0, c1, c2 ring.Poly

	// Case Ciphertext (x) Ciphertext
	if op0.Degree() == 1 && op1.Degree() == 1 {

		_, level, err := eval.InitOutputBinaryOp(op0.El(), op1.El(), opOut.Degree(), opOut.El())

		if err != nil {
			return err
		}

		ringQ := eval.parameters.RingQ().AtLevel(level)

		c00 = eval.buffQ[0]
		c01 = eval.buffQ[1]

		c0 = opOut.Value[0]
		c1 = opOut.Value[1]

		if !relin {
			opOut.El().Resize(2, level)
			c2 = opOut.Value[2]
		} else {
			opOut.El().Resize(1, level)
			c2 = eval.buffQ[2]
		}

		// Avoid overwriting if the second input is the output
		var tmp0, tmp1 *rlwe.Operand[ring.Poly]
		if op1.El() == opOut.El() {
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
				return fmt.Errorf("cannot MulRelin: Relinearize: %w", err)
			}

			tmpCt := &rlwe.Ciphertext{}
			tmpCt.Value = []ring.Poly{eval.BuffQP[1].Q, eval.BuffQP[2].Q}
			tmpCt.MetaData = &rlwe.MetaData{IsNTT: true}

			eval.GadgetProduct(level, c2, &rlk.GadgetCiphertext, tmpCt)
			ringQ.Add(c0, tmpCt.Value[0], opOut.Value[0])
			ringQ.Add(c1, tmpCt.Value[1], opOut.Value[1])
		}

		// Case Plaintext (x) Ciphertext or Ciphertext (x) Plaintext
	} else {

		_, level, err := eval.InitOutputBinaryOp(op0.El(), op1.El(), opOut.Degree(), opOut.El())

		if err != nil {
			return err
		}

		ringQ := eval.parameters.RingQ().AtLevel(level)

		var c0 ring.Poly
		var c1 []ring.Poly
		if op0.Degree() == 0 {
			c0 = eval.buffQ[0]
			ringQ.MForm(op0.Value[0], c0)
			c1 = op1.El().Value

		} else {
			c0 = eval.buffQ[0]
			ringQ.MForm(op1.El().Value[0], c0)
			c1 = op0.Value
		}

		opOut.El().Resize(op0.Degree()+op1.Degree(), level)

		for i := range c1 {
			ringQ.MulCoeffsMontgomery(c0, c1[i], opOut.Value[i])
		}
	}

	return
}

// MulThenAdd evaluate opOut = opOut + op0 * op1.
//
// The following types are accepted for op1:
// - rlwe.OperandInterface[ring.Poly]
// - complex128, float64, int, int64, uint, uint64, *big.Int, *big.Float, *bignum.Complex
// - []complex128, []float64, []*big.Float or []*bignum.Complex
// Passing an invalid type will return an error.
//
// If op1.(type) is complex128, float64, int, int64, uint64. *big.Float, *big.Int or *ring.Complex:
//
// This function will not modify op0 but will multiply opOut by Q[min(op0.Level(), opOut.Level())] if:
// - op0.PlaintextScale == opOut.PlaintextScale
// - constant is not a Gaussian integer.
//
// If op0.PlaintextScale == opOut.PlaintextScale, and constant is not a Gaussian integer, then the constant will be scaled by
// Q[min(op0.Level(), opOut.Level())] else if opOut.PlaintextScale > op0.PlaintextScale, the constant will be scaled by opOut.PlaintextScale/op0.PlaintextScale.
//
// To correctly use this function, make sure that either op0.PlaintextScale == opOut.PlaintextScale or
// opOut.PlaintextScale = op0.PlaintextScale * Q[min(op0.Level(), opOut.Level())].
//
// If op1.(type) is []complex128, []float64, []*big.Float or []*bignum.Complex:
// - If opOut.PlaintextScale == op0.PlaintextScale, op1 will be encoded and scaled by Q[min(op0.Level(), opOut.Level())]
// - If opOut.PlaintextScale > op0.PlaintextScale, op1 will be encoded ans scaled by opOut.PlaintextScale/op1.PlaintextScale.
// Then the method will recurse with op1 given as rlwe.OperandInterface[ring.Poly].
//
// If op1.(type) is rlwe.OperandInterface[ring.Poly], the multiplication is carried outwithout relinearization and:
//
// This function will return an error if op0.PlaintextScale > opOut.PlaintextScale and user must ensure that opOut.PlaintextScale <= op0.PlaintextScale * op1.PlaintextScale.
// If opOut.PlaintextScale < op0.PlaintextScale * op1.PlaintextScale, then scales up opOut before adding the result.
// Additionally, the procedure will return an error if:
// - either op0 or op1 are have a degree higher than 1.
// - opOut.Degree != op0.Degree + op1.Degree.
// - opOut = op0 or op1.
func (eval Evaluator) MulThenAdd(op0 *rlwe.Ciphertext, op1 interface{}, opOut *rlwe.Ciphertext) (err error) {
	switch op1 := op1.(type) {
	case rlwe.OperandInterface[ring.Poly]:

		// Generic in place evaluation
		return eval.mulRelinThenAdd(op0, op1.El(), false, opOut)
	case complex128, float64, int, int64, uint, uint64, *big.Int, *big.Float, *bignum.Complex:

		// Retrieves the minimum level
		level := utils.Min(op0.Level(), opOut.Level())

		// Resizes the output to the minimum level
		opOut.Resize(opOut.Degree(), level)

		// Gets the ring at the minimum level
		ringQ := eval.parameters.RingQ().AtLevel(level)

		// Convertes the scalar to a *bignum.Complex
		cmplxBig := bignum.ToComplex(op1, eval.parameters.PlaintextPrecision())

		var scaleRLWE rlwe.Scale

		// If op0 and opOut scales are identical, but the op1 is not a Gaussian integer then multiplies opOut by scaleRLWE.
		// This ensures noiseless addition with opOut = scaleRLWE * opOut + op0 * round(scalar * scaleRLWE).
		if cmp := op0.PlaintextScale.Cmp(opOut.PlaintextScale); cmp == 0 {

			if cmplxBig.IsInt() {
				scaleRLWE = rlwe.NewScale(1)
			} else {
				scaleRLWE = rlwe.NewScale(ringQ.SubRings[level].Modulus)

				for i := 1; i < eval.parameters.PlaintextScaleToModuliRatio(); i++ {
					scaleRLWE = scaleRLWE.Mul(rlwe.NewScale(ringQ.SubRings[level-i].Modulus))
				}

				scaleInt := new(big.Int)
				scaleRLWE.Value.Int(scaleInt)
				if err = eval.Mul(opOut, scaleInt, opOut); err != nil {
					return fmt.Errorf("cannot MulThenAdd: %w", err)
				}
				opOut.PlaintextScale = opOut.PlaintextScale.Mul(scaleRLWE)
			}

		} else if cmp == -1 { // opOut.PlaintextScale > op0.PlaintextScale then the scaling factor for op1 becomes the quotient between the two scales
			scaleRLWE = opOut.PlaintextScale.Div(op0.PlaintextScale)
		} else {
			return fmt.Errorf("cannot MulThenAdd: op0.PlaintextScale > opOut.PlaintextScale is not supported")
		}

		RNSReal, RNSImag := bigComplexToRNSScalar(ringQ, &scaleRLWE.Value, cmplxBig)

		eval.evaluateWithScalar(level, op0.Value, RNSReal, RNSImag, opOut.Value, ringQ.MulDoubleRNSScalarThenAdd)

		return

	case []complex128, []float64, []*big.Float, []*bignum.Complex:

		// Retrieves minimum level
		level := utils.Min(op0.Level(), opOut.Level())

		// Resizes output to minimum level
		opOut.Resize(opOut.Degree(), level)

		// Gets the ring at the target level
		ringQ := eval.parameters.RingQ().AtLevel(level)

		var scaleRLWE rlwe.Scale
		if cmp := op0.PlaintextScale.Cmp(opOut.PlaintextScale); cmp == 0 { // If op0 and opOut scales are identical then multiplies opOut by scaleRLWE.

			scaleRLWE = rlwe.NewScale(ringQ.SubRings[level].Modulus)

			for i := 1; i < eval.parameters.PlaintextScaleToModuliRatio(); i++ {
				scaleRLWE = scaleRLWE.Mul(rlwe.NewScale(ringQ.SubRings[level-i].Modulus))
			}

			scaleInt := new(big.Int)
			scaleRLWE.Value.Int(scaleInt)
			if err = eval.Mul(opOut, scaleInt, opOut); err != nil {
				return fmt.Errorf("cannot MulThenAdd: %w", err)
			}
			opOut.PlaintextScale = opOut.PlaintextScale.Mul(scaleRLWE)

		} else if cmp == -1 { // opOut.PlaintextScale > op0.PlaintextScale then the scaling factor for op1 becomes the quotient between the two scales
			scaleRLWE = opOut.PlaintextScale.Div(op0.PlaintextScale)
		} else {
			return fmt.Errorf("cannot MulThenAdd: op0.PlaintextScale > opOut.PlaintextScale is not supported")
		}

		// Instantiates new plaintext from buffer
		pt, err := rlwe.NewPlaintextAtLevelFromPoly(level, eval.buffQ[0])
		if err != nil {
			panic(err)
		}
		*pt.MetaData = *op0.MetaData
		pt.PlaintextScale = scaleRLWE

		// Encodes the vector on the plaintext
		if err = eval.Encoder.Encode(op1, pt); err != nil {
			return err
		}

		// Generic in place evaluation
		return eval.mulRelinThenAdd(op0, pt.El(), false, opOut)

	default:
		return fmt.Errorf("op1.(type) must be rlwe.OperandInterface[ring.Poly], complex128, float64, int, int64, uint, uint64, *big.Int, *big.Float, *bignum.Complex, []complex128, []float64, []*big.Float or []*bignum.Complex, but is %T", op1)
	}
}

// MulRelinThenAdd multiplies op0 with op1 with relinearization and adds the result on opOut.
//
// The following types are accepted for op1:
// - rlwe.OperandInterface[ring.Poly]
// - complex128, float64, int, int64, uint, uint64, *big.Int, *big.Float, *bignum.Complex
// - []complex128, []float64, []*big.Float or []*bignum.Complex
// Passing an invalid type will return an error.
//
// User must ensure that opOut.PlaintextScale <= op0.PlaintextScale * op1.PlaintextScale.
//
// If opOut.PlaintextScale < op0.PlaintextScale * op1.PlaintextScale, then scales up opOut before adding the result.
//
// The procedure will return an error if either op0.Degree or op1.Degree > 1.
// The procedure will return an error if opOut.Degree != op0.Degree + op1.Degree.
// The procedure will return an error if the evaluator was not created with an relinearization key.
// The procedure will return an error if opOut = op0 or op1.
func (eval Evaluator) MulRelinThenAdd(op0 *rlwe.Ciphertext, op1 interface{}, opOut *rlwe.Ciphertext) (err error) {
	switch op1 := op1.(type) {
	case rlwe.OperandInterface[ring.Poly]:
		if op1.Degree() == 0 {
			return eval.MulThenAdd(op0, op1, opOut)
		} else {
			return eval.mulRelinThenAdd(op0, op1.El(), true, opOut)
		}
	default:
		return eval.MulThenAdd(op0, op1, opOut)
	}
}

func (eval Evaluator) mulRelinThenAdd(op0 *rlwe.Ciphertext, op1 *rlwe.Operand[ring.Poly], relin bool, opOut *rlwe.Ciphertext) (err error) {

	_, level, err := eval.InitOutputBinaryOp(op0.El(), op1.El(), utils.Max(op0.Degree(), op1.Degree()), opOut.El())
	if err != nil {
		return err
	}

	if op0.Degree()+op1.Degree() > 2 {
		return fmt.Errorf("cannot MulRelinThenAdd: the sum of the input elements' degree cannot be larger than 2")
	}

	if op0.El() == opOut.El() || op1.El() == opOut.El() {
		return fmt.Errorf("cannot MulRelinThenAdd: opOut must be different from op0 and op1")
	}

	resScale := op0.PlaintextScale.Mul(op1.PlaintextScale)

	if opOut.PlaintextScale.Cmp(resScale) == -1 {
		ratio := resScale.Div(opOut.PlaintextScale)
		// Only scales up if int(ratio) >= 2
		if ratio.Float64() >= 2.0 {
			if err = eval.Mul(opOut, &ratio.Value, opOut); err != nil {
				return fmt.Errorf("cannot MulRelinThenAdd: %w", err)
			}
			opOut.PlaintextScale = resScale
		}
	}

	ringQ := eval.parameters.RingQ().AtLevel(level)

	var c00, c01, c0, c1, c2 ring.Poly

	// Case Ciphertext (x) Ciphertext
	if op0.Degree() == 1 && op1.Degree() == 1 {

		c00 = eval.buffQ[0]
		c01 = eval.buffQ[1]

		c0 = opOut.Value[0]
		c1 = opOut.Value[1]

		if !relin {
			opOut.El().Resize(2, level)
			c2 = opOut.Value[2]
		} else {
			// No resize here since we add on opOut
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
				return fmt.Errorf("cannot relinearize: %w", err)
			}

			ringQ.MulCoeffsMontgomery(c01, tmp1.Value[1], c2) // c2 += c[1]*c[1]

			tmpCt := &rlwe.Ciphertext{}
			tmpCt.Value = []ring.Poly{eval.BuffQP[1].Q, eval.BuffQP[2].Q}
			tmpCt.MetaData = &rlwe.MetaData{IsNTT: true}

			eval.GadgetProduct(level, c2, &rlk.GadgetCiphertext, tmpCt)
			ringQ.Add(c0, tmpCt.Value[0], c0)
			ringQ.Add(c1, tmpCt.Value[1], c1)
		} else {
			ringQ.MulCoeffsMontgomeryThenAdd(c01, tmp1.Value[1], c2) // c2 += c[1]*c[1]
		}

		// Case Plaintext (x) Ciphertext or Ciphertext (x) Plaintext
	} else {

		if opOut.Degree() < op0.Degree() {
			opOut.Resize(op0.Degree(), level)
		}

		c00 := eval.buffQ[0]

		ringQ.MForm(op1.El().Value[0], c00)
		for i := range op0.Value {
			ringQ.MulCoeffsMontgomeryThenAdd(op0.Value[i], c00, opOut.Value[i])
		}
	}

	return
}

// RelinearizeNew applies the relinearization procedure on op0 and returns the result in a newly
// created Ciphertext. The input Ciphertext must be of degree two.
func (eval Evaluator) RelinearizeNew(op0 *rlwe.Ciphertext) (opOut *rlwe.Ciphertext, err error) {
	opOut = NewCiphertext(eval.parameters, 1, op0.Level())
	return opOut, eval.Relinearize(op0, opOut)
}

// ApplyEvaluationKeyNew applies the rlwe.EvaluationKey on op0 and returns the result on a new ciphertext opOut.
func (eval Evaluator) ApplyEvaluationKeyNew(op0 *rlwe.Ciphertext, evk *rlwe.EvaluationKey) (opOut *rlwe.Ciphertext, err error) {
	opOut = NewCiphertext(eval.parameters, op0.Degree(), op0.Level())
	return opOut, eval.ApplyEvaluationKey(op0, evk, opOut)
}

// RotateNew rotates the columns of op0 by k positions to the left, and returns the result in a newly created element.
// The method will return an error if the evaluator hasn't been given an evaluation key set with the appropriate GaloisKey.
func (eval Evaluator) RotateNew(op0 *rlwe.Ciphertext, k int) (opOut *rlwe.Ciphertext, err error) {
	opOut = NewCiphertext(eval.parameters, op0.Degree(), op0.Level())
	return opOut, eval.Rotate(op0, k, opOut)
}

// Rotate rotates the columns of op0 by k positions to the left and returns the result in opOut.
// The method will return an error if the evaluator hasn't been given an evaluation key set with the appropriate GaloisKey.
func (eval Evaluator) Rotate(op0 *rlwe.Ciphertext, k int, opOut *rlwe.Ciphertext) (err error) {
	if err = eval.Automorphism(op0, eval.parameters.GaloisElement(k), opOut); err != nil {
		return fmt.Errorf("cannot Rotate: %w", err)
	}
	return
}

// ConjugateNew conjugates op0 (which is equivalent to a row rotation) and returns the result in a newly created element.
// The method will return an error if the evaluator hasn't been given an evaluation key set with the appropriate GaloisKey.
func (eval Evaluator) ConjugateNew(op0 *rlwe.Ciphertext) (opOut *rlwe.Ciphertext, err error) {
	opOut = NewCiphertext(eval.parameters, op0.Degree(), op0.Level())
	return opOut, eval.Conjugate(op0, opOut)
}

// Conjugate conjugates op0 (which is equivalent to a row rotation) and returns the result in opOut.
// The method will return an error if the evaluator hasn't been given an evaluation key set with the appropriate GaloisKey.
func (eval Evaluator) Conjugate(op0 *rlwe.Ciphertext, opOut *rlwe.Ciphertext) (err error) {

	if eval.parameters.RingType() == ring.ConjugateInvariant {
		return fmt.Errorf("cannot Conjugate: method is not supported when parameters.RingType() == ring.ConjugateInvariant")
	}

	if err = eval.Automorphism(op0, eval.parameters.GaloisElementInverse(), opOut); err != nil {
		return fmt.Errorf("cannot Conjugate: %w", err)
	}

	return
}

// RotateHoistedNew takes an input Ciphertext and a list of rotations and returns a map of Ciphertext, where each element of the map is the input Ciphertext
// rotation by one element of the list. It is much faster than sequential calls to Rotate.
func (eval Evaluator) RotateHoistedNew(ctIn *rlwe.Ciphertext, rotations []int) (opOut map[int]*rlwe.Ciphertext, err error) {
	opOut = make(map[int]*rlwe.Ciphertext)
	for _, i := range rotations {
		opOut[i] = NewCiphertext(eval.parameters, 1, ctIn.Level())
	}

	return opOut, eval.RotateHoisted(ctIn, rotations, opOut)
}

// RotateHoisted takes an input Ciphertext and a list of rotations and populates a map of pre-allocated Ciphertexts,
// where each element of the map is the input Ciphertext rotation by one element of the list.
// It is much faster than sequential calls to Rotate.
func (eval Evaluator) RotateHoisted(ctIn *rlwe.Ciphertext, rotations []int, opOut map[int]*rlwe.Ciphertext) (err error) {
	levelQ := ctIn.Level()
	eval.DecomposeNTT(levelQ, eval.parameters.MaxLevelP(), eval.parameters.PCount(), ctIn.Value[1], ctIn.IsNTT, eval.BuffDecompQP)
	for _, i := range rotations {
		if err = eval.AutomorphismHoisted(levelQ, ctIn, eval.BuffDecompQP, eval.parameters.GaloisElement(i), opOut[i]); err != nil {
			return fmt.Errorf("cannot RotateHoisted: %w", err)
		}
	}

	return
}

func (eval Evaluator) RotateHoistedLazyNew(level int, rotations []int, ct *rlwe.Ciphertext, c2DecompQP []ringqp.Poly) (cOut map[int]*rlwe.Operand[ringqp.Poly], err error) {
	cOut = make(map[int]*rlwe.Operand[ringqp.Poly])
	for _, i := range rotations {
		if i != 0 {
			cOut[i] = rlwe.NewOperandQP(eval.parameters.Parameters, 1, level, eval.parameters.MaxLevelP())
			if err = eval.AutomorphismHoistedLazy(level, ct, c2DecompQP, eval.parameters.GaloisElement(i), cOut[i]); err != nil {
				return nil, fmt.Errorf("cannot RotateHoistedLazyNew: %w", err)
			}
		}
	}

	return
}

// Parameters returns the Parametrs of the underlying struct as an rlwe.ParametersInterface.
func (eval Evaluator) Parameters() rlwe.ParametersInterface {
	return eval.parameters
}

// ShallowCopy creates a shallow copy of this evaluator in which all the read-only data-structures are
// shared with the receiver and the temporary buffers are reallocated. The receiver and the returned
// Evaluators can be used concurrently.
func (eval Evaluator) ShallowCopy() *Evaluator {
	return &Evaluator{
		parameters:       eval.parameters,
		Encoder:          NewEncoder(eval.parameters),
		Evaluator:        eval.Evaluator.ShallowCopy(),
		evaluatorBuffers: newEvaluatorBuffers(eval.parameters),
	}
}

// WithKey creates a shallow copy of the receiver Evaluator for which the new EvaluationKey is evaluationKey
// and where the temporary buffers are shared. The receiver and the returned Evaluators cannot be used concurrently.
func (eval Evaluator) WithKey(evk rlwe.EvaluationKeySet) *Evaluator {
	return &Evaluator{
		Evaluator:        eval.Evaluator.WithKey(evk),
		parameters:       eval.parameters,
		Encoder:          eval.Encoder,
		evaluatorBuffers: eval.evaluatorBuffers,
	}
}
