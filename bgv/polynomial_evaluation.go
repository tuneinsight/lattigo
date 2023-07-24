package bgv

import (
	"fmt"
	"math/big"
	"math/bits"

	"github.com/tuneinsight/lattigo/v4/hebase"
	"github.com/tuneinsight/lattigo/v4/rlwe"
	"github.com/tuneinsight/lattigo/v4/utils"
	"github.com/tuneinsight/lattigo/v4/utils/bignum"
)

func (eval Evaluator) Polynomial(input interface{}, p interface{}, InvariantTensoring bool, targetScale rlwe.Scale) (opOut *rlwe.Ciphertext, err error) {

	var polyVec hebase.PolynomialVector
	switch p := p.(type) {
	case bignum.Polynomial:
		polyVec = hebase.PolynomialVector{Value: []hebase.Polynomial{{Polynomial: p, MaxDeg: p.Degree(), Lead: true, Lazy: false}}}
	case hebase.Polynomial:
		polyVec = hebase.PolynomialVector{Value: []hebase.Polynomial{p}}
	case hebase.PolynomialVector:
		polyVec = p
	default:
		return nil, fmt.Errorf("cannot Polynomial: invalid polynomial type: %T", p)
	}

	polyEval := PolynomialEvaluator{
		Evaluator:          &eval,
		InvariantTensoring: InvariantTensoring,
	}

	var powerbasis hebase.PowerBasis
	switch input := input.(type) {
	case *rlwe.Ciphertext:

		if level, depth := input.Level(), polyVec.Value[0].Depth(); level < depth {
			return nil, fmt.Errorf("%d levels < %d log(d) -> cannot evaluate poly", level, depth)
		}

		powerbasis = hebase.NewPowerBasis(input, bignum.Monomial)

	case hebase.PowerBasis:
		if input.Value[1] == nil {
			return nil, fmt.Errorf("cannot evaluatePolyVector: given PowerBasis[1] is empty")
		}
		powerbasis = input
	default:
		return nil, fmt.Errorf("cannot evaluatePolyVector: invalid input, must be either *hebase.Ciphertext or *PowerBasis")
	}

	logDegree := bits.Len64(uint64(polyVec.Value[0].Degree()))
	logSplit := bignum.OptimalSplit(logDegree)

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

	PS := polyVec.GetPatersonStockmeyerPolynomial(*eval.GetParameters(), powerbasis.Value[1].Level(), powerbasis.Value[1].Scale, targetScale, &dummyEvaluator{*eval.GetParameters(), InvariantTensoring})

	if opOut, err = hebase.EvaluatePatersonStockmeyerPolynomialVector(PS, powerbasis, polyEval); err != nil {
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
func (d dummyEvaluator) Rescale(op0 *hebase.DummyOperand) {
	if !d.InvariantTensoring {
		op0.Scale = op0.Scale.Div(rlwe.NewScale(d.params.Q()[op0.Level]))
		op0.Level--
	}
}

// Mul multiplies two DummyOperand, stores the result the taret DummyOperand and returns the result.
func (d dummyEvaluator) MulNew(op0, op1 *hebase.DummyOperand) (opOut *hebase.DummyOperand) {
	opOut = new(hebase.DummyOperand)
	opOut.Level = utils.Min(op0.Level, op1.Level)
	opOut.Scale = op0.Scale.Mul(op1.Scale)
	if d.InvariantTensoring {
		params := d.params
		qModTNeg := new(big.Int).Mod(params.RingQ().ModulusAtLevel[opOut.Level], new(big.Int).SetUint64(params.T())).Uint64()
		qModTNeg = params.T() - qModTNeg
		opOut.Scale = opOut.Scale.Div(params.NewScale(qModTNeg))
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

	// tScaleNew = targetScale*currentQi/XPow.Scale
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

func (polyEval PolynomialEvaluator) Mul(op0 *rlwe.Ciphertext, op1 interface{}, opOut *rlwe.Ciphertext) (err error) {
	if !polyEval.InvariantTensoring {
		return polyEval.Evaluator.Mul(op0, op1, opOut)
	} else {
		return polyEval.Evaluator.MulScaleInvariant(op0, op1, opOut)
	}
}

func (polyEval PolynomialEvaluator) MulRelin(op0 *rlwe.Ciphertext, op1 interface{}, opOut *rlwe.Ciphertext) (err error) {
	if !polyEval.InvariantTensoring {
		return polyEval.Evaluator.MulRelin(op0, op1, opOut)
	} else {
		return polyEval.Evaluator.MulRelinScaleInvariant(op0, op1, opOut)
	}
}

func (polyEval PolynomialEvaluator) MulNew(op0 *rlwe.Ciphertext, op1 interface{}) (opOut *rlwe.Ciphertext, err error) {
	if !polyEval.InvariantTensoring {
		return polyEval.Evaluator.MulNew(op0, op1)
	} else {
		return polyEval.Evaluator.MulScaleInvariantNew(op0, op1)
	}
}

func (polyEval PolynomialEvaluator) MulRelinNew(op0 *rlwe.Ciphertext, op1 interface{}) (opOut *rlwe.Ciphertext, err error) {
	if !polyEval.InvariantTensoring {
		return polyEval.Evaluator.MulRelinNew(op0, op1)
	} else {
		return polyEval.Evaluator.MulRelinScaleInvariantNew(op0, op1)
	}
}

func (polyEval PolynomialEvaluator) Rescale(op0, op1 *rlwe.Ciphertext) (err error) {
	if !polyEval.InvariantTensoring {
		return polyEval.Evaluator.Rescale(op0, op1)
	}
	return
}

func (polyEval PolynomialEvaluator) EvaluatePolynomialVectorFromPowerBasis(targetLevel int, pol hebase.PolynomialVector, pb hebase.PowerBasis, targetScale rlwe.Scale) (res *rlwe.Ciphertext, err error) {

	X := pb.Value

	params := *polyEval.Evaluator.GetParameters()
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
			res = NewCiphertext(params, 1, targetLevel)
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
				pt, err := rlwe.NewPlaintextAtLevelFromPoly(targetLevel, res.Value[0])
				if err != nil {
					panic(err)
				}
				pt.Scale = res.Scale
				pt.IsNTT = NTTFlag
				pt.IsBatched = true
				if err = polyEval.Encode(values, pt); err != nil {
					return nil, err
				}
			}

			return
		}

		// Allocates the output ciphertext
		res = NewCiphertext(params, maximumCiphertextDegree, targetLevel)
		res.Scale = targetScale

		// Allocates a temporary plaintext to encode the values
		buffq := polyEval.Evaluator.BuffQ()
		pt, err := rlwe.NewPlaintextAtLevelFromPoly(targetLevel, buffq[0]) // buffQ[0] is safe in this case
		if err != nil {
			panic(err)
		}
		pt.IsBatched = true
		pt.Scale = targetScale
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
			if err := polyEval.Add(res, values, res); err != nil {
				return nil, err
			}

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
				if err = polyEval.MulThenAdd(X[key], values, res); err != nil {
					return nil, err
				}
				toEncode = false
			}
		}

	} else {

		c := pol.Value[0].Coeffs[0].Uint64()

		if minimumDegreeNonZeroCoefficient == 0 {

			res = NewCiphertext(params, 1, targetLevel)
			res.Scale = targetScale

			if c != 0 {
				if err := polyEval.Add(res, c, res); err != nil {
					return nil, err
				}
			}

			return
		}

		res = NewCiphertext(params, maximumCiphertextDegree, targetLevel)
		res.Scale = targetScale

		if c != 0 {
			if err := polyEval.Add(res, c, res); err != nil {
				return nil, err
			}
		}

		for key := pol.Value[0].Degree(); key > 0; key-- {
			if c = pol.Value[0].Coeffs[key].Uint64(); key != 0 && c != 0 {
				// MulScalarAndAdd automatically scales c to match the scale of res.
				if err := polyEval.MulThenAdd(X[key], c, res); err != nil {
					return nil, err
				}
			}
		}
	}

	return
}
