package integer

import (
	"github.com/tuneinsight/lattigo/v4/bgv"
	"github.com/tuneinsight/lattigo/v4/circuits"
	"github.com/tuneinsight/lattigo/v4/rlwe"
	"github.com/tuneinsight/lattigo/v4/utils"
	"github.com/tuneinsight/lattigo/v4/utils/bignum"
)

type PolynomialEvaluator struct {
	circuits.PolynomialEvaluator
	bgv.Parameters
	InvariantTensoring bool
}

// NewIntegerPowerBasis is a wrapper of NewPolynomialBasis.
// This function creates a new powerBasis from the input ciphertext.
// The input ciphertext is treated as the base monomial X used to
// generate the other powers X^{n}.
func NewIntegerPowerBasis(ct *rlwe.Ciphertext) circuits.PowerBasis {
	return circuits.NewPowerBasis(ct, bignum.Monomial)
}

// NewIntegerPolynomial is a wrapper of NewPolynomial.
// This function creates a new polynomial from the input coefficients.
// This polynomial can be evaluated on a ciphertext.
func NewIntegerPolynomial[T Integer](coeffs []T) circuits.Polynomial {
	return circuits.NewPolynomial(bignum.NewPolynomial(bignum.Monomial, coeffs, nil))
}

func NewPolynomialEvaluator(params bgv.Parameters, eval *bgv.Evaluator, InvariantTensoring bool) *PolynomialEvaluator {
	e := new(PolynomialEvaluator)

	if InvariantTensoring {
		e.PolynomialEvaluator = circuits.PolynomialEvaluator{EvaluatorForPolyEval: scaleInvariantEvaluator{eval}, EvaluatorBuffers: eval.GetEvaluatorBuffer()}
	} else {
		e.PolynomialEvaluator = circuits.PolynomialEvaluator{EvaluatorForPolyEval: eval, EvaluatorBuffers: eval.GetEvaluatorBuffer()}
	}

	e.InvariantTensoring = InvariantTensoring
	e.Parameters = params
	return e
}

func (eval PolynomialEvaluator) Polynomial(input interface{}, p interface{}, targetScale rlwe.Scale) (opOut *rlwe.Ciphertext, err error) {
	return circuits.EvaluatePolynomial(eval.PolynomialEvaluator, eval, input, p, targetScale, 1, &simIntegerPolynomialEvaluator{eval.Parameters, eval.InvariantTensoring})
}

type scaleInvariantEvaluator struct {
	*bgv.Evaluator
}

func (polyEval scaleInvariantEvaluator) Mul(op0 *rlwe.Ciphertext, op1 interface{}, opOut *rlwe.Ciphertext) (err error) {
	return polyEval.MulScaleInvariant(op0, op1, opOut)
}

func (polyEval scaleInvariantEvaluator) MulRelin(op0 *rlwe.Ciphertext, op1 interface{}, opOut *rlwe.Ciphertext) (err error) {
	return polyEval.Evaluator.MulRelinScaleInvariant(op0, op1, opOut)
}

func (polyEval scaleInvariantEvaluator) MulNew(op0 *rlwe.Ciphertext, op1 interface{}) (opOut *rlwe.Ciphertext, err error) {
	return polyEval.Evaluator.MulScaleInvariantNew(op0, op1)
}

func (polyEval scaleInvariantEvaluator) MulRelinNew(op0 *rlwe.Ciphertext, op1 interface{}) (opOut *rlwe.Ciphertext, err error) {
	return polyEval.Evaluator.MulRelinScaleInvariantNew(op0, op1)
}

func (polyEval scaleInvariantEvaluator) Rescale(op0, op1 *rlwe.Ciphertext) (err error) {
	return nil
}

func (eval PolynomialEvaluator) EvaluatePolynomialVectorFromPowerBasis(targetLevel int, pol circuits.PolynomialVector, pb circuits.PowerBasis, targetScale rlwe.Scale) (res *rlwe.Ciphertext, err error) {

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
