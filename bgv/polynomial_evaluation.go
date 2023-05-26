package bgv

import (
	"fmt"
	"math/big"
	"math/bits"
	"runtime"

	"github.com/tuneinsight/lattigo/v4/rlwe"
	"github.com/tuneinsight/lattigo/v4/utils"
	"github.com/tuneinsight/lattigo/v4/utils/bignum/polynomial"
)

func (eval *Evaluator) Polynomial(input interface{}, p interface{}, invariantTensoring bool, targetScale rlwe.Scale) (opOut *rlwe.Ciphertext, err error) {

	var polyVec *rlwe.PolynomialVector
	switch p := p.(type) {
	case *polynomial.Polynomial:
		polyVec = &rlwe.PolynomialVector{Value: []*rlwe.Polynomial{&rlwe.Polynomial{Polynomial: p, MaxDeg: p.Degree(), Lead: true, Lazy: false}}}
	case *rlwe.Polynomial:
		polyVec = &rlwe.PolynomialVector{Value: []*rlwe.Polynomial{p}}
	case *rlwe.PolynomialVector:
		polyVec = p
	default:
		return nil, fmt.Errorf("cannot Polynomial: invalid polynomial type: %T", p)
	}

	var powerbasis *PowerBasis
	switch input := input.(type) {
	case *rlwe.Ciphertext:

		if level, depth := input.Level(), polyVec.Value[0].Depth(); level < depth {
			return nil, fmt.Errorf("%d levels < %d log(d) -> cannot evaluate poly", level, depth)
		}

		powerbasis = NewPowerBasis(input)

	case *PowerBasis:
		if input.Value[1] == nil {
			return nil, fmt.Errorf("cannot evaluatePolyVector: given PowerBasis[1] is empty")
		}
		powerbasis = input
	default:
		return nil, fmt.Errorf("cannot evaluatePolyVector: invalid input, must be either *rlwe.Ciphertext or *PowerBasis")
	}

	logDegree := bits.Len64(uint64(polyVec.Value[0].Degree()))
	logSplit := polynomial.OptimalSplit(logDegree)

	var odd, even bool = false, false
	for _, p := range polyVec.Value {
		odd, even = odd || p.IsOdd, even || p.IsEven
	}

	// Computes all the powers of two with relinearization
	// This will recursively compute and store all powers of two up to 2^logDegree
	if err = powerbasis.GenPower(1<<(logDegree-1), false, invariantTensoring, eval); err != nil {
		return nil, err
	}

	// Computes the intermediate powers, starting from the largest, without relinearization if possible
	for i := (1 << logSplit) - 1; i > 2; i-- {
		if !(even || odd) || (i&1 == 0 && even) || (i&1 == 1 && odd) {
			if err = powerbasis.GenPower(i, polyVec.Value[0].Lazy, invariantTensoring, eval); err != nil {
				return nil, err
			}
		}
	}

	PS := polyVec.GetPatersonStockmeyerPolynomial(eval.params.Parameters, powerbasis.Value[1].Level(), powerbasis.Value[1].Scale, targetScale, &dummyEvaluator{eval.params, invariantTensoring})

	polyEval := &polynomialEvaluator{
		Evaluator:          eval,
		Encoder:            NewEncoder(eval.params),
		invariantTensoring: invariantTensoring,
	}

	if opOut, err = rlwe.EvaluatePatersonStockmeyerPolynomialVector(PS, powerbasis.PowerBasis, polyEval); err != nil {
		return nil, err
	}

	powerbasis = nil

	runtime.GC()
	return opOut, err
}

type dummyEvaluator struct {
	params             Parameters
	invariantTensoring bool
}

func (d *dummyEvaluator) PolynomialDepth(degree int) int {
	if d.invariantTensoring {
		return 0
	}
	return bits.Len64(uint64(degree)) - 1
}

// Rescale rescales the target DummyOperand n times and returns it.
func (d *dummyEvaluator) Rescale(op0 *rlwe.DummyOperand) {
	if !d.invariantTensoring {
		op0.Scale = op0.Scale.Div(rlwe.NewScale(d.params.Q()[op0.Level]))
		op0.Level--
	}
}

// Mul multiplies two DummyOperand, stores the result the taret DummyOperand and returns the result.
func (d *dummyEvaluator) MulNew(op0, op1 *rlwe.DummyOperand) (op2 *rlwe.DummyOperand) {
	op2 = new(rlwe.DummyOperand)
	op2.Level = utils.Min(op0.Level, op1.Level)
	op2.Scale = op0.Scale.Mul(op1.Scale)
	if d.invariantTensoring {
		params := d.params
		qModTNeg := new(big.Int).Mod(params.RingQ().ModulusAtLevel[op2.Level], new(big.Int).SetUint64(params.T())).Uint64()
		qModTNeg = params.T() - qModTNeg
		op2.Scale = op2.Scale.Div(params.NewScale(qModTNeg))
	}
	return
}

func (d *dummyEvaluator) UpdateLevelAndScaleBabyStep(lead bool, tLevelOld int, tScaleOld rlwe.Scale) (tLevelNew int, tScaleNew rlwe.Scale) {
	tLevelNew = tLevelOld
	tScaleNew = tScaleOld
	if !d.invariantTensoring && lead {
		tScaleNew = tScaleOld.Mul(d.params.NewScale(d.params.Q()[tLevelOld]))
	}
	return
}

func (d *dummyEvaluator) UpdateLevelAndScaleGiantStep(lead bool, tLevelOld int, tScaleOld, xPowScale rlwe.Scale) (tLevelNew int, tScaleNew rlwe.Scale) {

	Q := d.params.Q()

	tLevelNew = tLevelOld
	tScaleNew = tScaleOld.Div(xPowScale)

	// tScaleNew = targetScale*currentQi/XPow.Scale
	if !d.invariantTensoring {

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

	if !d.invariantTensoring {
		tLevelNew++
	}

	return
}

type polynomialEvaluator struct {
	*Evaluator
	*Encoder
	invariantTensoring bool
}

func (polyEval *polynomialEvaluator) Mul(op0 *rlwe.Ciphertext, op1 interface{}, op2 *rlwe.Ciphertext) {
	if !polyEval.invariantTensoring {
		polyEval.Evaluator.Mul(op0, op1, op2)
	} else {
		polyEval.Evaluator.MulInvariant(op0, op1, op2)
	}
}

func (polyEval *polynomialEvaluator) Rescale(op0, op1 *rlwe.Ciphertext) (err error) {
	if !polyEval.invariantTensoring {
		return polyEval.Evaluator.Rescale(op0, op1)
	}
	return
}

func (polyEval *polynomialEvaluator) EvaluatePolynomialVectorFromPowerBasis(targetLevel int, pol *rlwe.PolynomialVector, pb *rlwe.PowerBasis, targetScale rlwe.Scale) (res *rlwe.Ciphertext, err error) {

	X := pb.Value

	params := polyEval.Evaluator.params
	slotsIndex := pol.SlotsIndex
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
		values := make([]uint64, params.N())

		// If the degree of the poly is zero
		if minimumDegreeNonZeroCoefficient == 0 {

			// Allocates the output ciphertext
			res = rlwe.NewCiphertext(params.Parameters, 1, targetLevel)
			res.Scale = targetScale

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
				pt := rlwe.NewPlaintextAtLevelFromPoly(targetLevel, res.Value[0])
				pt.Scale = res.Scale
				pt.IsNTT = true
				polyEval.Encode(values, pt)
			}

			return
		}

		// Allocates the output ciphertext
		res = rlwe.NewCiphertext(params.Parameters, maximumCiphertextDegree, targetLevel)
		res.Scale = targetScale

		// Allocates a temporary plaintext to encode the values
		pt := rlwe.NewPlaintextAtLevelFromPoly(targetLevel, polyEval.Evaluator.BuffQ()[0]) // buffQ[0] is safe in this case
		pt.Scale = targetScale
		pt.IsNTT = true

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
			pt.Scale = res.Scale
			polyEval.Encode(values, pt)
			polyEval.Add(res, pt, res)
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
				pt.Scale = targetScale.Div(X[key].Scale)
				polyEval.Encode(values, pt)
				polyEval.MulThenAdd(X[key], pt, res)
				toEncode = false
			}
		}

	} else {

		c := pol.Value[0].Coeffs[0].Uint64()

		if minimumDegreeNonZeroCoefficient == 0 {

			res = rlwe.NewCiphertext(params.Parameters, 1, targetLevel)
			res.Scale = targetScale

			if c != 0 {
				polyEval.Add(res, c, res)
			}

			return
		}

		res = rlwe.NewCiphertext(params.Parameters, maximumCiphertextDegree, targetLevel)
		res.Scale = targetScale

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

func isOddOrEvenPolynomial(coeffs []uint64) (odd, even bool) {
	even = true
	odd = true
	for i, c := range coeffs {
		isnotzero := c != 0
		odd = odd && !(i&1 == 0 && isnotzero)
		even = even && !(i&1 == 1 && isnotzero)
		if !odd && !even {
			break
		}
	}

	return
}
