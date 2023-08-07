package float

import (
	"math/big"

	"github.com/tuneinsight/lattigo/v4/circuits"
	"github.com/tuneinsight/lattigo/v4/ckks"
	"github.com/tuneinsight/lattigo/v4/rlwe"
	"github.com/tuneinsight/lattigo/v4/utils"
	"github.com/tuneinsight/lattigo/v4/utils/bignum"
)

type PolynomialEvaluator struct {
	circuits.PolynomialEvaluator
	Parameters ckks.Parameters
}

// NewFloatPowerBasis is a wrapper of NewPolynomialBasis.
// This function creates a new powerBasis from the input ciphertext.
// The input ciphertext is treated as the base monomial X used to
// generate the other powers X^{n}.
func NewFloatPowerBasis(ct *rlwe.Ciphertext, basis bignum.Basis) circuits.PowerBasis {
	return circuits.NewPowerBasis(ct, basis)
}

func NewPolynomialEvaluator(params ckks.Parameters, eval circuits.EvaluatorForPolyEval) *PolynomialEvaluator {
	e := new(PolynomialEvaluator)
	e.PolynomialEvaluator = circuits.PolynomialEvaluator{EvaluatorForPolyEval: eval, EvaluatorBuffers: eval.GetEvaluatorBuffer()}
	e.Parameters = params
	return e
}

// Polynomial evaluates a polynomial in standard basis on the input Ciphertext in ceil(log2(deg+1)) levels.
// Returns an error if the input ciphertext does not have enough level to carry out the full polynomial evaluation.
// Returns an error if something is wrong with the scale.
// If the polynomial is given in Chebyshev basis, then a change of basis ct' = (2/(b-a)) * (ct + (-a-b)/(b-a))
// is necessary before the polynomial evaluation to ensure correctness.
// input must be either *rlwe.Ciphertext or *PolynomialBasis.
// pol: a *bignum.Polynomial, *Polynomial or *PolynomialVector
// targetScale: the desired output scale. This value shouldn't differ too much from the original ciphertext scale. It can
// for example be used to correct small deviations in the ciphertext scale and reset it to the default scale.
func (eval PolynomialEvaluator) Polynomial(input interface{}, p interface{}, targetScale rlwe.Scale) (opOut *rlwe.Ciphertext, err error) {
	levelsConsummedPerRescaling := eval.Parameters.LevelsConsummedPerRescaling()
	return circuits.EvaluatePolynomial(eval.PolynomialEvaluator, eval, input, p, targetScale, levelsConsummedPerRescaling, &simEvaluator{eval.Parameters, levelsConsummedPerRescaling})
}

func (eval PolynomialEvaluator) EvaluatePolynomialVectorFromPowerBasis(targetLevel int, pol circuits.PolynomialVector, pb circuits.PowerBasis, targetScale rlwe.Scale) (res *rlwe.Ciphertext, err error) {

	// Map[int] of the powers [X^{0}, X^{1}, X^{2}, ...]
	X := pb.Value

	// Retrieve the number of slots
	logSlots := X[1].LogDimensions
	slots := 1 << logSlots.Cols

	params := eval.Parameters
	slotsIndex := pol.SlotsIndex
	even := pol.IsEven()
	odd := pol.IsOdd()

	// Retrieve the degree of the highest degree non-zero coefficient
	// TODO: optimize for nil/zero coefficients
	minimumDegreeNonZeroCoefficient := len(pol.Value[0].Coeffs) - 1
	if even && !odd {
		minimumDegreeNonZeroCoefficient--
	}

	// Gets the maximum degree of the ciphertexts among the power basis
	// TODO: optimize for nil/zero coefficients, odd/even polynomial
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
		values := make([]*bignum.Complex, slots)

		// If the degree of the poly is zero
		if minimumDegreeNonZeroCoefficient == 0 {

			// Allocates the output ciphertext
			res = ckks.NewCiphertext(params, 1, targetLevel)
			res.Scale = targetScale
			res.LogDimensions = logSlots

			// Looks for non-zero coefficients among the degree 0 coefficients of the polynomials
			if even {
				for i, p := range pol.Value {
					if !isZero(p.Coeffs[0]) {
						toEncode = true
						for _, j := range slotsIndex[i] {
							values[j] = p.Coeffs[0]
						}
					}
				}
			}

			// If a non-zero coefficient was found, encode the values, adds on the ciphertext, and returns
			if toEncode {
				pt := &rlwe.Plaintext{}
				pt.Value = res.Value[0]
				pt.MetaData = res.MetaData
				if err = eval.Encode(values, pt); err != nil {
					return nil, err
				}
			}

			return
		}

		// Allocates the output ciphertext
		res = ckks.NewCiphertext(params, maximumCiphertextDegree, targetLevel)
		res.Scale = targetScale
		res.LogDimensions = logSlots

		// Looks for a non-zero coefficient among the degree zero coefficient of the polynomials
		if even {
			for i, p := range pol.Value {
				if !isZero(p.Coeffs[0]) {
					toEncode = true
					for _, j := range slotsIndex[i] {
						values[j] = p.Coeffs[0]
					}
				}
			}
		}

		// If a non-zero degre coefficient was found, encode and adds the values on the output
		// ciphertext
		if toEncode {
			if err = eval.Add(res, values, res); err != nil {
				return
			}
			toEncode = false
		}

		// Loops starting from the highest degree coefficient
		for key := pol.Value[0].Degree(); key > 0; key-- {

			var reset bool

			if !(even || odd) || (key&1 == 0 && even) || (key&1 == 1 && odd) {

				// Loops over the polynomials
				for i, p := range pol.Value {

					// Looks for a non-zero coefficient
					if !isZero(p.Coeffs[key]) {
						toEncode = true

						// Resets the temporary array to zero
						// is needed if a zero coefficient
						// is at the place of a previous non-zero
						// coefficient
						if !reset {
							for j := range values {
								if values[j] != nil {
									values[j][0].SetFloat64(0)
									values[j][1].SetFloat64(0)
								}
							}
							reset = true
						}

						// Copies the coefficient on the temporary array
						// according to the slot map index
						for _, j := range slotsIndex[i] {
							values[j] = p.Coeffs[key]
						}
					}
				}
			}

			// If a non-zero degre coefficient was found, encode and adds the values on the output
			// ciphertext
			if toEncode {
				if err = eval.MulThenAdd(X[key], values, res); err != nil {
					return
				}
				toEncode = false
			}
		}

	} else {

		var c *bignum.Complex
		if even && !isZero(pol.Value[0].Coeffs[0]) {
			c = pol.Value[0].Coeffs[0]
		}

		if minimumDegreeNonZeroCoefficient == 0 {

			res = ckks.NewCiphertext(params, 1, targetLevel)
			res.Scale = targetScale
			res.LogDimensions = logSlots

			if !isZero(c) {
				if err = eval.Add(res, c, res); err != nil {
					return
				}
			}

			return
		}

		res = ckks.NewCiphertext(params, maximumCiphertextDegree, targetLevel)
		res.Scale = targetScale
		res.LogDimensions = logSlots

		if c != nil {
			if err = eval.Add(res, c, res); err != nil {
				return
			}
		}

		for key := pol.Value[0].Degree(); key > 0; key-- {
			if c = pol.Value[0].Coeffs[key]; key != 0 && !isZero(c) && (!(even || odd) || (key&1 == 0 && even) || (key&1 == 1 && odd)) {
				if err = eval.MulThenAdd(X[key], c, res); err != nil {
					return
				}
			}
		}
	}

	return
}

func isZero(c *bignum.Complex) bool {
	zero := new(big.Float)
	return c == nil || (c[0].Cmp(zero) == 0 && c[1].Cmp(zero) == 0)
}
