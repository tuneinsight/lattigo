package circuits

import (
	"fmt"
	"math/big"
	"math/bits"

	"github.com/tuneinsight/lattigo/v4/bgv"
	"github.com/tuneinsight/lattigo/v4/rlwe"
	"github.com/tuneinsight/lattigo/v4/utils"
	"github.com/tuneinsight/lattigo/v4/utils/bignum"
)

type BGVPolyEvaluator struct {
	*PolynomialEvaluator
	bgv.Parameters
	InvariantTensoring bool
}

// NewBGVPowerBasis is a wrapper of NewPolynomialBasis.
// This function creates a new powerBasis from the input ciphertext.
// The input ciphertext is treated as the base monomial X used to
// generate the other powers X^{n}.
func NewBGVPowerBasis(ct *rlwe.Ciphertext) PowerBasis {
	return NewPowerBasis(ct, bignum.Monomial)
}

// NewBGVPolynomial is a wrapper of NewPolynomial.
// This function creates a new polynomial from the input coefficients.
// This polynomial can be evaluated on a ciphertext.
func NewBGVPolynomial[T int64 | uint64](coeffs []T) Polynomial {
	return NewPolynomial(bignum.NewPolynomial(bignum.Monomial, coeffs, nil))
}

// NewBGVPolynomialVector is a wrapper of NewPolynomialVector.
// This function creates a new PolynomialVector from the input polynomials and the desired function mapping.
// This polynomial vector can be evaluated on a ciphertext.
func NewBGVPolynomialVector(polys []Polynomial, mapping map[int][]int) (PolynomialVector, error) {
	return NewPolynomialVector(polys, mapping)
}

func NewBGVPolynomialEvaluator(params bgv.Parameters, eval *bgv.Evaluator, InvariantTensoring bool) *BGVPolyEvaluator {
	e := new(BGVPolyEvaluator)

	if InvariantTensoring {
		e.PolynomialEvaluator = &PolynomialEvaluator{BGVScaleInvariantEvaluator{eval}, eval.GetEvaluatorBuffer()}
	} else {
		e.PolynomialEvaluator = &PolynomialEvaluator{eval, eval.GetEvaluatorBuffer()}
	}

	e.InvariantTensoring = InvariantTensoring
	e.Parameters = params
	return e
}

func (eval *BGVPolyEvaluator) Polynomial(input interface{}, p interface{}, targetScale rlwe.Scale) (opOut *rlwe.Ciphertext, err error) {

	var polyVec PolynomialVector
	switch p := p.(type) {
	case bignum.Polynomial:
		polyVec = PolynomialVector{Value: []Polynomial{{Polynomial: p, MaxDeg: p.Degree(), Lead: true, Lazy: false}}}
	case Polynomial:
		polyVec = PolynomialVector{Value: []Polynomial{p}}
	case PolynomialVector:
		polyVec = p
	default:
		return nil, fmt.Errorf("cannot Polynomial: invalid polynomial type: %T", p)
	}

	var powerbasis PowerBasis
	switch input := input.(type) {
	case *rlwe.Ciphertext:
		if level, depth := input.Level(), polyVec.Value[0].Depth(); level < depth {
			return nil, fmt.Errorf("%d levels < %d log(d) -> cannot evaluate poly", level, depth)
		}
		powerbasis = NewPowerBasis(input, bignum.Monomial)
	case PowerBasis:
		if input.Value[1] == nil {
			return nil, fmt.Errorf("cannot evaluatePolyVector: given PowerBasis[1] is empty")
		}
		powerbasis = input
	default:
		return nil, fmt.Errorf("cannot evaluatePolyVector: invalid input, must be either *Ciphertext or *PowerBasis")
	}

	logDegree := bits.Len64(uint64(polyVec.Value[0].Degree()))
	logSplit := bignum.OptimalSplit(logDegree)

	var odd, even bool
	for _, p := range polyVec.Value {
		odd, even = odd || p.IsOdd, even || p.IsEven
	}

	// Computes all the powers of two with relinearization
	// This will recursively compute and store all powers of two up to 2^logDegree
	if err = powerbasis.GenPower(1<<(logDegree-1), false, eval); err != nil {
		return nil, err
	}

	// Computes the intermediate powers, starting from the largest, without relinearization if possible
	for i := (1 << logSplit) - 1; i > 2; i-- {
		if !(even || odd) || (i&1 == 0 && even) || (i&1 == 1 && odd) {
			if err = powerbasis.GenPower(i, polyVec.Value[0].Lazy, eval); err != nil {
				return nil, err
			}
		}
	}

	PS := polyVec.GetPatersonStockmeyerPolynomial(eval.Parameters.Parameters, powerbasis.Value[1].Level(), powerbasis.Value[1].Scale, targetScale, &dummyBGVPolyEvaluator{eval.Parameters, eval.InvariantTensoring})

	if opOut, err = eval.EvaluatePatersonStockmeyerPolynomialVector(eval, PS, powerbasis); err != nil {
		return nil, err
	}

	return opOut, err
}

type dummyBGVPolyEvaluator struct {
	params             bgv.Parameters
	InvariantTensoring bool
}

func (d dummyBGVPolyEvaluator) PolynomialDepth(degree int) int {
	if d.InvariantTensoring {
		return 0
	}
	return bits.Len64(uint64(degree)) - 1
}

// Rescale rescales the target DummyOperand n times and returns it.
func (d dummyBGVPolyEvaluator) Rescale(op0 *DummyOperand) {
	if !d.InvariantTensoring {
		op0.Scale = op0.Scale.Div(rlwe.NewScale(d.params.Q()[op0.Level]))
		op0.Level--
	}
}

// Mul multiplies two DummyOperand, stores the result the taret DummyOperand and returns the result.
func (d dummyBGVPolyEvaluator) MulNew(op0, op1 *DummyOperand) (opOut *DummyOperand) {
	opOut = new(DummyOperand)
	opOut.Level = utils.Min(op0.Level, op1.Level)

	if d.InvariantTensoring {
		opOut.Scale = bgv.MulScaleInvariant(d.params, op0.Scale, op1.Scale, opOut.Level)
	} else {
		opOut.Scale = op0.Scale.Mul(op1.Scale)
	}

	return
}

func (d dummyBGVPolyEvaluator) UpdateLevelAndScaleBabyStep(lead bool, tLevelOld int, tScaleOld rlwe.Scale) (tLevelNew int, tScaleNew rlwe.Scale) {
	tLevelNew = tLevelOld
	tScaleNew = tScaleOld
	if !d.InvariantTensoring && lead {
		tScaleNew = tScaleOld.Mul(d.params.NewScale(d.params.Q()[tLevelOld]))
	}
	return
}

func (d dummyBGVPolyEvaluator) UpdateLevelAndScaleGiantStep(lead bool, tLevelOld int, tScaleOld, xPowScale rlwe.Scale) (tLevelNew int, tScaleNew rlwe.Scale) {

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

		T := d.params.PlaintextModulus()

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

type BGVScaleInvariantEvaluator struct {
	*bgv.Evaluator
}

func (polyEval BGVScaleInvariantEvaluator) Mul(op0 *rlwe.Ciphertext, op1 interface{}, opOut *rlwe.Ciphertext) (err error) {
	return polyEval.MulScaleInvariant(op0, op1, opOut)
}

func (polyEval BGVScaleInvariantEvaluator) MulRelin(op0 *rlwe.Ciphertext, op1 interface{}, opOut *rlwe.Ciphertext) (err error) {
	return polyEval.Evaluator.MulRelinScaleInvariant(op0, op1, opOut)
}

func (polyEval BGVScaleInvariantEvaluator) MulNew(op0 *rlwe.Ciphertext, op1 interface{}) (opOut *rlwe.Ciphertext, err error) {
	return polyEval.Evaluator.MulScaleInvariantNew(op0, op1)
}

func (polyEval BGVScaleInvariantEvaluator) MulRelinNew(op0 *rlwe.Ciphertext, op1 interface{}) (opOut *rlwe.Ciphertext, err error) {
	return polyEval.Evaluator.MulRelinScaleInvariantNew(op0, op1)
}

func (polyEval BGVScaleInvariantEvaluator) Rescale(op0, op1 *rlwe.Ciphertext) (err error) {
	return nil
}

func (eval BGVPolyEvaluator) EvaluatePolynomialVectorFromPowerBasis(targetLevel int, pol PolynomialVector, pb PowerBasis, targetScale rlwe.Scale) (res *rlwe.Ciphertext, err error) {

	X := pb.Value

	params := eval.Parameters
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
			res = bgv.NewCiphertext(params, 1, targetLevel)
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
				pt.IsNTT = bgv.NTTFlag
				pt.IsBatched = true
				if err = eval.Encode(values, pt); err != nil {
					return nil, err
				}
			}

			return
		}

		// Allocates the output ciphertext
		res = bgv.NewCiphertext(params, maximumCiphertextDegree, targetLevel)
		res.Scale = targetScale

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
			if err := eval.Add(res, values, res); err != nil {
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
				if err = eval.MulThenAdd(X[key], values, res); err != nil {
					return nil, err
				}
				toEncode = false
			}
		}

	} else {

		c := pol.Value[0].Coeffs[0].Uint64()

		if minimumDegreeNonZeroCoefficient == 0 {

			res = bgv.NewCiphertext(params, 1, targetLevel)
			res.Scale = targetScale

			if c != 0 {
				if err := eval.Add(res, c, res); err != nil {
					return nil, err
				}
			}

			return
		}

		res = bgv.NewCiphertext(params, maximumCiphertextDegree, targetLevel)
		res.Scale = targetScale

		if c != 0 {
			if err := eval.Add(res, c, res); err != nil {
				return nil, err
			}
		}

		for key := pol.Value[0].Degree(); key > 0; key-- {
			if c = pol.Value[0].Coeffs[key].Uint64(); key != 0 && c != 0 {
				// MulScalarAndAdd automatically scales c to match the scale of res.
				if err := eval.MulThenAdd(X[key], c, res); err != nil {
					return nil, err
				}
			}
		}
	}

	return
}
