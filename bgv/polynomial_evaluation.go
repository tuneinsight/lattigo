package bgv

import (
	"fmt"
	"math/big"
	"math/bits"

	"github.com/tuneinsight/lattigo/v4/rlwe"
	"github.com/tuneinsight/lattigo/v4/utils"
	"github.com/tuneinsight/lattigo/v4/utils/bignum/polynomial"
)

// NewPowerBasis creates a new PowerBasis from the input ciphertext.
// The input ciphertext is treated as the base monomial X used to 
// generate the other powers X^{n}.
func NewPowerBasis(ct *rlwe.Ciphertext) rlwe.PowerBasis {
	return rlwe.NewPowerBasis(ct, polynomial.Monomial)
}

func (eval Evaluator) Polynomial(input interface{}, p interface{}, InvariantTensoring bool, targetScale rlwe.Scale) (opOut *rlwe.Ciphertext, err error) {

	var polyVec rlwe.PolynomialVector
	switch p := p.(type) {
	case polynomial.Polynomial:
		polyVec = rlwe.PolynomialVector{Value: []rlwe.Polynomial{{Polynomial: p, MaxDeg: p.Degree(), Lead: true, Lazy: false}}}
	case rlwe.Polynomial:
		polyVec = rlwe.PolynomialVector{Value: []rlwe.Polynomial{p}}
	case rlwe.PolynomialVector:
		polyVec = p
	default:
		return nil, fmt.Errorf("cannot Polynomial: invalid polynomial type: %T", p)
	}

	polyEval := PolynomialEvaluator{
		Evaluator:          &eval,
		InvariantTensoring: InvariantTensoring,
	}

	var powerbasis rlwe.PowerBasis
	switch input := input.(type) {
	case *rlwe.Ciphertext:

		if level, depth := input.Level(), polyVec.Value[0].Depth(); level < depth {
			return nil, fmt.Errorf("%d levels < %d log(d) -> cannot evaluate poly", level, depth)
		}

		powerbasis = rlwe.NewPowerBasis(input, polynomial.Monomial)

	case rlwe.PowerBasis:
		if input.Value[1] == nil {
			return nil, fmt.Errorf("cannot evaluatePolyVector: given PowerBasis[1] is empty")
		}
		powerbasis = input
	default:
		return nil, fmt.Errorf("cannot evaluatePolyVector: invalid input, must be either *rlwe.Ciphertext or *PowerBasis")
	}

	logDegree := bits.Len64(uint64(polyVec.Value[0].Degree()))
	logSplit := polynomial.OptimalSplit(logDegree)

	var odd, even bool
	for _, p := range polyVec.Value {
		odd, even = odd || p.IsOdd, even || p.IsEven
	}

	// Computes all the powers of two with relinearization
	// This will recursively compute and store all powers of two up to 2^logDegree
	if err = powerbasis.GenPower(1<<(logDegree-1), false, polyEval); err != nil {
		return nil, err
	}

	// Computes the intermediate powers, starting from the largest, without relinearization if possible
	for i := (1 << logSplit) - 1; i > 2; i-- {
		if !(even || odd) || (i&1 == 0 && even) || (i&1 == 1 && odd) {
			if err = powerbasis.GenPower(i, polyVec.Value[0].Lazy, polyEval); err != nil {
				return nil, err
			}
		}
	}

	PS := polyVec.GetPatersonStockmeyerPolynomial(eval.Parameters(), powerbasis.Value[1].Level(), powerbasis.Value[1].PlaintextScale, targetScale, &dummyEvaluator{eval.Parameters().(Parameters), InvariantTensoring})

	if opOut, err = rlwe.EvaluatePatersonStockmeyerPolynomialVector(PS, powerbasis, polyEval); err != nil {
		return nil, err
	}

	return opOut, err
}

type dummyEvaluator struct {
	params             Parameters
	InvariantTensoring bool
}

func (d dummyEvaluator) PolynomialDepth(degree int) int {
	if d.InvariantTensoring {
		return 0
	}
	return bits.Len64(uint64(degree)) - 1
}

// Rescale rescales the target DummyOperand n times and returns it.
func (d dummyEvaluator) Rescale(op0 *rlwe.DummyOperand) {
	if !d.InvariantTensoring {
		op0.PlaintextScale = op0.PlaintextScale.Div(rlwe.NewScale(d.params.Q()[op0.Level]))
		op0.Level--
	}
}

// Mul multiplies two DummyOperand, stores the result the taret DummyOperand and returns the result.
func (d dummyEvaluator) MulNew(op0, op1 *rlwe.DummyOperand) (op2 *rlwe.DummyOperand) {
	op2 = new(rlwe.DummyOperand)
	op2.Level = utils.Min(op0.Level, op1.Level)
	op2.PlaintextScale = op0.PlaintextScale.Mul(op1.PlaintextScale)
	if d.InvariantTensoring {
		params := d.params
		qModTNeg := new(big.Int).Mod(params.RingQ().ModulusAtLevel[op2.Level], new(big.Int).SetUint64(params.T())).Uint64()
		qModTNeg = params.T() - qModTNeg
		op2.PlaintextScale = op2.PlaintextScale.Div(params.NewScale(qModTNeg))
	}

	return
}

func (d dummyEvaluator) UpdateLevelAndScaleBabyStep(lead bool, tLevelOld int, tScaleOld rlwe.Scale) (tLevelNew int, tScaleNew rlwe.Scale) {
	tLevelNew = tLevelOld
	tScaleNew = tScaleOld
	if !d.InvariantTensoring && lead {
		tScaleNew = tScaleOld.Mul(d.params.NewScale(d.params.Q()[tLevelOld]))
	}
	return
}

func (d dummyEvaluator) UpdateLevelAndScaleGiantStep(lead bool, tLevelOld int, tScaleOld, xPowScale rlwe.Scale) (tLevelNew int, tScaleNew rlwe.Scale) {

	Q := d.params.Q()

	tLevelNew = tLevelOld
	tScaleNew = tScaleOld.Div(xPowScale)

	// tScaleNew = targetScale*currentQi/XPow.PlaintextScale
	if !d.InvariantTensoring {

		var currentQi uint64
		if lead {
			currentQi = Q[tLevelNew]
		} else {
			currentQi = Q[tLevelNew+1]
		}

		tScaleNew = tScaleNew.Mul(d.params.NewScale(currentQi))

	} else {

		T := d.params.T()

		// -Q mod T
		qModTNeg := new(big.Int).Mod(d.params.RingQ().ModulusAtLevel[tLevelNew], new(big.Int).SetUint64(T)).Uint64()
		qModTNeg = T - qModTNeg
		tScaleNew = tScaleNew.Mul(d.params.NewScale(qModTNeg))
	}

	if !d.InvariantTensoring {
		tLevelNew++
	}

	return
}

func NewPolynomialEvaluator(eval *Evaluator, InvariantTensoring bool) *PolynomialEvaluator {
	return &PolynomialEvaluator{Evaluator: eval, InvariantTensoring: InvariantTensoring}
}

type PolynomialEvaluator struct {
	*Evaluator
	InvariantTensoring bool
}

func (polyEval PolynomialEvaluator) Parameters() rlwe.ParametersInterface {
	return polyEval.Evaluator.Parameters()
}

func (polyEval PolynomialEvaluator) Mul(op0 *rlwe.Ciphertext, op1 interface{}, op2 *rlwe.Ciphertext) {
	if !polyEval.InvariantTensoring {
		polyEval.Evaluator.Mul(op0, op1, op2)
	} else {
		polyEval.Evaluator.MulInvariant(op0, op1, op2)
	}
}

func (polyEval PolynomialEvaluator) MulRelin(op0 *rlwe.Ciphertext, op1 interface{}, op2 *rlwe.Ciphertext) {
	if !polyEval.InvariantTensoring {
		polyEval.Evaluator.MulRelin(op0, op1, op2)
	} else {
		polyEval.Evaluator.MulRelinInvariant(op0, op1, op2)
	}
}

func (polyEval PolynomialEvaluator) MulNew(op0 *rlwe.Ciphertext, op1 interface{}) (op2 *rlwe.Ciphertext) {
	if !polyEval.InvariantTensoring {
		return polyEval.Evaluator.MulNew(op0, op1)
	} else {
		return polyEval.Evaluator.MulInvariantNew(op0, op1)
	}
}

func (polyEval PolynomialEvaluator) MulRelinNew(op0 *rlwe.Ciphertext, op1 interface{}) (op2 *rlwe.Ciphertext) {
	if !polyEval.InvariantTensoring {
		return polyEval.Evaluator.MulRelinNew(op0, op1)
	} else {
		return polyEval.Evaluator.MulRelinInvariantNew(op0, op1)
	}
}

func (polyEval PolynomialEvaluator) Rescale(op0, op1 *rlwe.Ciphertext) (err error) {
	if !polyEval.InvariantTensoring {
		return polyEval.Evaluator.Rescale(op0, op1)
	}
	return
}

func (polyEval PolynomialEvaluator) EvaluatePolynomialVectorFromPowerBasis(targetLevel int, pol rlwe.PolynomialVector, pb rlwe.PowerBasis, targetScale rlwe.Scale) (res *rlwe.Ciphertext, err error) {

	X := pb.Value

	params := polyEval.Evaluator.Parameters().(Parameters)
	slotsIndex := pol.SlotsIndex
	slots := params.RingT().N()
	even := pol.IsEven()
	odd := pol.IsOdd()

	// Retrieve the degree of the highest degree non-zero coefficient
	// TODO: optimize for nil/zero coefficients
	minimumDegreeNonZeroCoefficient := len(pol.Value[0].Coeffs) - 1
	if even && !odd {
		minimumDegreeNonZeroCoefficient--
	}

	// Get the minimum non-zero degree coefficient
	maximumCiphertextDegree := 0
	for i := pol.Value[0].Degree(); i > 0; i-- {
		if x, ok := X[i]; ok {
			maximumCiphertextDegree = utils.Max(maximumCiphertextDegree, x.Degree())
		}
	}

	// If an index slot is given (either multiply polynomials or masking)
	if slotsIndex != nil {

		var toEncode bool

		// Allocates temporary buffer for coefficients encoding
		values := make([]uint64, slots)

		// If the degree of the poly is zero
		if minimumDegreeNonZeroCoefficient == 0 {

			// Allocates the output ciphertext
			res = rlwe.NewCiphertext(params, 1, targetLevel)
			res.PlaintextScale = targetScale

			// Looks for non-zero coefficients among the degree 0 coefficients of the polynomials
			for i, p := range pol.Value {
				if c := p.Coeffs[0].Uint64(); c != 0 {
					toEncode = true
					for _, j := range slotsIndex[i] {
						values[j] = c
					}
				}
			}

			// If a non-zero coefficient was found, encode the values, adds on the ciphertext, and returns
			if toEncode {
				pt := rlwe.NewPlaintextAtLevelFromPoly(targetLevel, &res.Value[0])
				pt.PlaintextScale = res.PlaintextScale
				pt.IsNTT = NTTFlag
				if err = polyEval.Encode(values, pt); err != nil {
					return
				}
			}

			return
		}

		// Allocates the output ciphertext
		res = rlwe.NewCiphertext(params, maximumCiphertextDegree, targetLevel)
		res.PlaintextScale = targetScale

		// Allocates a temporary plaintext to encode the values
		buffq := polyEval.Evaluator.BuffQ()
		pt := rlwe.NewPlaintextAtLevelFromPoly(targetLevel, &buffq[0]) // buffQ[0] is safe in this case
		pt.PlaintextScale = targetScale
		pt.IsNTT = NTTFlag

		// Looks for a non-zero coefficient among the degree zero coefficient of the polynomials
		for i, p := range pol.Value {
			if c := p.Coeffs[0].Uint64(); c != 0 {
				toEncode = true
				for _, j := range slotsIndex[i] {
					values[j] = c
				}
			}
		}

		// If a non-zero degree coefficient was found, encode and adds the values on the output
		// ciphertext
		if toEncode {
			// Add would actually scale the plaintext accordingly,
			// but encoding with the correct scale is slightly faster
			polyEval.Add(res, values, res)
			toEncode = false
		}

		// Loops starting from the highest degree coefficient
		for key := pol.Value[0].Degree(); key > 0; key-- {

			var reset bool
			// Loops over the polynomials
			for i, p := range pol.Value {

				// Looks for a non-zero coefficient
				if c := p.Coeffs[key].Uint64(); c != 0 {
					toEncode = true

					// Resets the temporary array to zero
					// is needed if a zero coefficient
					// is at the place of a previous non-zero
					// coefficient
					if !reset {
						for j := range values {
							values[j] = 0
						}
						reset = true
					}

					// Copies the coefficient on the temporary array
					// according to the slot map index
					for _, j := range slotsIndex[i] {
						values[j] = c
					}
				}
			}

			// If a non-zero degree coefficient was found, encode and adds the values on the output
			// ciphertext
			if toEncode {

				// MulAndAdd would actually scale the plaintext accordingly,
				// but encoding with the correct scale is slightly faster
				polyEval.MulThenAdd(X[key], values, res)
				toEncode = false
			}
		}

	} else {

		c := pol.Value[0].Coeffs[0].Uint64()

		if minimumDegreeNonZeroCoefficient == 0 {

			res = rlwe.NewCiphertext(params, 1, targetLevel)
			res.PlaintextScale = targetScale

			if c != 0 {
				polyEval.Add(res, c, res)
			}

			return
		}

		res = rlwe.NewCiphertext(params, maximumCiphertextDegree, targetLevel)
		res.PlaintextScale = targetScale

		if c != 0 {
			polyEval.Add(res, c, res)
		}

		for key := pol.Value[0].Degree(); key > 0; key-- {
			if c = pol.Value[0].Coeffs[key].Uint64(); key != 0 && c != 0 {
				// MulScalarAndAdd automatically scales c to match the scale of res.
				polyEval.MulThenAdd(X[key], c, res)
			}
		}
	}

	return
}
