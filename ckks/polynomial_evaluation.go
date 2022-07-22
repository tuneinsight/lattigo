package ckks

import (
	"fmt"
	"math/big"
	"math/bits"
	"runtime"

	"github.com/tuneinsight/lattigo/v3/ring"
	"github.com/tuneinsight/lattigo/v3/rlwe"
	"github.com/tuneinsight/lattigo/v3/utils"
)

// EvaluatePoly evaluates a polynomial in standard basis on the input Ciphertext in ceil(log2(deg+1)) levels.
// Returns an error if the input ciphertext does not have enough level to carry out the full polynomial evaluation.
// Returns an error if something is wrong with the scale.
// If the polynomial is given in Chebyshev basis, then a change of basis ct' = (2/(b-a)) * (ct + (-a-b)/(b-a))
// is necessary before the polynomial evaluation to ensure correctness.
// Coefficients of the polynomial with an absolute value smaller than "IsNegligbleThreshold" will automatically be set to zero
// if the polynomial is "even" or "odd" (to ensure that the even or odd property remains valid
// after the "splitCoeffs" polynomial decomposition).
// input must be either *Ciphertext or *PolynomialBasis.
// pol: a Polynomial
// targetScale: the desired output scale. This value shouldn't differ too much from the original ciphertext scale. It can
// for example be used to correct small deviations in the ciphertext scale and reset it to the default scale.
func (eval *evaluator) EvaluatePoly(input interface{}, pol Polynomial, targetScale float64) (opOut *Ciphertext, err error) {

	var monomialBasis *PolynomialBasis
	if monomialBasis, err = eval.genPolynomialBasis(input, pol); err != nil {
		return
	}

	polyEval := &polynomialEvaluator{}
	polyEval.Evaluator = eval
	polyEval.PolynomialBasis = *monomialBasis
	polyEval.giant, polyEval.baby = pol.BSGSSplit()

	if opOut, err = polyEval.recurse(monomialBasis.Value[1].Level()-polyEval.giant+1, targetScale, evalPoly{pol, true, pol.Degree()}); err != nil {
		return nil, err
	}

	polyEval.Relinearize(opOut, opOut)

	if err = polyEval.Rescale(opOut, targetScale, opOut); err != nil {
		return nil, err
	}

	opOut.Scale = targetScale // solves float64 precision issues

	polyEval = nil
	runtime.GC()
	return opOut, err
}

// checkEnoughLevels checks that enough levels are available to evaluate the polynomial.
// Also checks if c is a Gaussian integer or not. If not, then one more level is needed
// to evaluate the polynomial.
func checkEnoughLevels(levels, depth int, c complex128) (err error) {

	if real(c) != float64(int64(real(c))) || imag(c) != float64(int64(imag(c))) {
		depth++
	}

	if levels < depth {
		return fmt.Errorf("%d levels < %d log(d) -> cannot evaluate", levels, depth)
	}

	return nil
}

type polynomialEvaluator struct {
	Evaluator
	Encoder
	PolynomialBasis
	giant int
	baby  int
}

type evalPoly struct {
	coefficients
	lead   bool
	maxDeg int
}

func (p *evalPoly) splitBSGS(split int) (polyq, polyr evalPoly) {

	polyq = evalPoly{}
	polyr = evalPoly{}

	polyq.coefficients, polyr.coefficients = p.coefficients.SplitBSGS(split)

	polyq.lead = p.lead
	polyq.maxDeg = p.maxDeg

	if p.maxDeg == p.Degree() {
		polyr.maxDeg = split - 1
	} else {
		polyr.maxDeg = p.maxDeg - (p.Degree() - split + 1)
	}

	return
}

func (polyEval *polynomialEvaluator) recurse(targetLevel int, targetScale float64, pol evalPoly) (res *Ciphertext, err error) {

	params := polyEval.Evaluator.(*evaluator).params

	baby := polyEval.baby

	// Recursively computes the evaluation of the Chebyshev polynomial using a baby-set giant-step algorithm.
	if pol.Degree() < (1 << baby) {

		if pol.lead && polyEval.baby > 1 && pol.maxDeg%(1<<(baby+1)) > (1<<(baby-1)) {

			giant := int(bits.Len64(uint64(pol.Degree())))
			baby := giant >> 1

			polyEvalBis := new(polynomialEvaluator)
			polyEvalBis.Evaluator = polyEval.Evaluator
			polyEvalBis.Encoder = polyEval.Encoder
			polyEvalBis.giant = giant
			polyEvalBis.baby = baby
			polyEvalBis.PolynomialBasis = polyEval.PolynomialBasis

			return polyEvalBis.recurse(targetLevel, targetScale, pol)
		}

		if pol.lead {
			targetScale *= params.QiFloat64(targetLevel)
		}

		switch poly := pol.coefficients.(type) {
		case *coefficientsComplex128:
			return polyEval.evaluatePolyFromPolynomialBasisComplex128(targetScale, targetLevel, poly)
		case *coefficientsBSGSComplex128:
			//fmt.Println(poly)
			//return polyEval.evaluatePolyFromPolynomialBasisComplex128(targetScale, targetLevel, poly)
		case *coefficientsBSGSPlaintext:
			return polyEval.evaluatePolyFromPlaintext(poly.coeffs[0])
		case *coefficientsBSGSCiphertext:
			return polyEval.evaluatePolyFromCiphertext(poly.coeffs[0])
		}
	}

	var nextPower = 1 << polyEval.baby
	for nextPower < (pol.Degree()>>1)+1 {
		nextPower <<= 1
	}

	coeffsq, coeffsr := pol.splitBSGS(nextPower)

	XPow := polyEval.PolynomialBasis.Value[nextPower]

	var currentQi float64
	if pol.lead {
		currentQi = params.QiFloat64(targetLevel)
	} else {
		currentQi = params.QiFloat64(targetLevel + 1)
	}

	if res, err = polyEval.recurse(targetLevel+1, targetScale*currentQi/XPow.Scale, coeffsq); err != nil {
		return nil, err
	}

	if res != nil {
		if res.Degree() == 2 {
			polyEval.Relinearize(res, res)
		}

		if err = polyEval.Rescale(res, params.DefaultScale(), res); err != nil {
			return nil, err
		}

		polyEval.Mul(res, XPow, res)
	}

	var tmp *Ciphertext
	if tmp, err = polyEval.recurse(res.Level(), res.Scale, coeffsr); err != nil {
		return nil, err
	}

	if res != nil {
		polyEval.Add(res, tmp, res)
	} else {
		res = tmp
	}

	return
}

func (polyEval *polynomialEvaluator) evaluatePolyFromPlaintext(pt []*Plaintext) (res *Ciphertext, err error) {

	X := polyEval.PolynomialBasis.Value

	if pt[0] != nil {
		res = &Ciphertext{Ciphertext: &rlwe.Ciphertext{Value: []*ring.Poly{pt[0].Value.CopyNew()}}, Scale: pt[0].Scale}
	}

	for i := 1; i < len(pt); i++ {
		if pt[i] != nil {
			if res == nil {
				res = polyEval.MulNew(X[i], pt[i])
			} else {
				polyEval.MulAndAdd(X[i], pt[i], res)
			}
		}
	}

	return
}

func (polyEval *polynomialEvaluator) evaluatePolyFromCiphertext(ct []*Ciphertext) (res *Ciphertext, err error) {

	X := polyEval.PolynomialBasis.Value

	if ct[0] != nil {
		res = ct[0].CopyNew()
	}

	for i := 1; i < len(ct); i++ {
		if ct[i] != nil {
			if res == nil {
				res = polyEval.MulNew(X[i], ct[i])
			} else {
				polyEval.MulAndAdd(X[i], ct[i], res)
			}
		}
	}

	return
}

func (polyEval *polynomialEvaluator) evaluatePolyFromPolynomialBasisComplex128(targetScale float64, level int, pol *coefficientsComplex128) (res *Ciphertext, err error) {

	X := polyEval.PolynomialBasis.Value

	params := polyEval.Evaluator.(*evaluator).params
	slotsIndex := pol.slotsIndex

	minimumDegreeNonZeroCoefficient := pol.Degree()

	if pol.even {
		minimumDegreeNonZeroCoefficient--
	}

	// Get the minimum non-zero degree coefficient
	maximumCiphertextDegree := 0
	for i := pol.Degree(); i > 0; i-- {
		if x, ok := X[i]; ok {
			maximumCiphertextDegree = utils.MaxInt(maximumCiphertextDegree, x.Degree())
		}
	}

	// If an index slot is given (either multiply polynomials or masking)
	if slotsIndex != nil {

		if polyEval.Encoder == nil {
			polyEval.Encoder = NewEncoder(params)
		}

		var toEncode bool

		// Allocates temporary buffer for coefficients encoding
		values := make([]complex128, params.Slots())

		// If the degree of the poly is zero
		if minimumDegreeNonZeroCoefficient == 0 {

			// Allocates the output ciphertext
			res = NewCiphertext(params, 1, level, targetScale)

			// Looks for non-zero coefficients among the degree 0 coefficients of the polynomials
			for i, c := range pol.coeffs {
				if c[0] != 0 {
					toEncode = true
					for _, j := range slotsIndex[i] {
						values[j] = c[0]
					}
				}
			}

			// If a non-zero coefficient was found, encode the values, adds on the ciphertext, and returns
			if toEncode {
				pt := NewPlaintextAtLevelFromPoly(level, res.Value[0])
				pt.Scale = res.Scale
				polyEval.EncodeSlots(values, pt, params.LogSlots())
			}

			return
		}

		// Allocates the output ciphertext
		res = NewCiphertext(params, maximumCiphertextDegree, level, targetScale)

		// Allocates a temporary plaintext to encode the values
		pt := NewPlaintextAtLevelFromPoly(level, polyEval.Evaluator.BuffCt().Value[0])

		// Looks for a non-zero coefficient among the degree zero coefficient of the polynomials
		for i, c := range pol.coeffs {
			if c[0] != 0 {
				toEncode = true
				for _, j := range slotsIndex[i] {
					values[j] = c[0]
				}
			}
		}

		// If a non-zero degre coefficient was found, encode and adds the values on the output
		// ciphertext
		if toEncode {
			pt.Scale = res.Scale
			polyEval.EncodeSlots(values, pt, params.LogSlots())
			polyEval.Add(res, pt, res)
			toEncode = false
		}

		// Loops starting from the highest degree coefficient
		for key := pol.Degree(); key > 0; key-- {

			var reset bool
			// Loops over the polynomials
			for i, c := range pol.coeffs {

				// Looks for a non-zero coefficient
				if c[key] != 0 {
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
						values[j] = c[key]
					}
				}
			}

			// If a non-zero degre coefficient was found, encode and adds the values on the output
			// ciphertext
			if toEncode {
				pt.Scale = targetScale / X[key].Scale
				polyEval.EncodeSlots(values, pt, params.LogSlots())
				polyEval.MulAndAdd(X[key], pt, res)
				toEncode = false
			}
		}

	} else {

		c := pol.coeffs[0][0]

		if minimumDegreeNonZeroCoefficient == 0 {

			res = NewCiphertext(params, 1, level, targetScale)

			if isNotNegligible(c) {
				polyEval.AddConst(res, c, res)
			}

			return
		}

		res = NewCiphertext(params, maximumCiphertextDegree, level, targetScale)

		if c != 0 {
			polyEval.AddConst(res, c, res)
		}

		cRealFlo, cImagFlo, constScale := ring.NewFloat(0, 128), ring.NewFloat(0, 128), ring.NewFloat(0, 128)
		cRealBig, cImagBig := ring.NewUint(0), ring.NewUint(0)

		for key := pol.Degree(); key > 0; key-- {

			c = pol.coeffs[0][key]

			if key != 0 && c != 0 {

				cRealFlo.SetFloat64(real(c))
				cImagFlo.SetFloat64(imag(c))
				constScale.SetFloat64(targetScale / X[key].Scale)

				// Target scale * rescale-scale / power basis scale
				cRealFlo.Mul(cRealFlo, constScale)
				cImagFlo.Mul(cImagFlo, constScale)

				if cRealFlo.Sign() < 0 {
					cRealFlo.Sub(cRealFlo, new(big.Float).SetFloat64(0.5))
				} else {
					cRealFlo.Add(cRealFlo, new(big.Float).SetFloat64(0.5))
				}

				if cImagFlo.Sign() < 0 {
					cImagFlo.Sub(cImagFlo, new(big.Float).SetFloat64(0.5))
				} else {
					cImagFlo.Add(cImagFlo, new(big.Float).SetFloat64(0.5))
				}

				cRealFlo.Int(cRealBig)
				cImagFlo.Int(cImagBig)

				polyEval.MultByGaussianIntegerAndAdd(X[key], cRealBig, cImagBig, res)
			}
		}
	}

	return
}
