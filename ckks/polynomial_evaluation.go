package ckks

import (
	"fmt"
	"math/big"
	"math/bits"
	"runtime"

	"github.com/tuneinsight/lattigo/v4/rlwe"
	"github.com/tuneinsight/lattigo/v4/utils"
	"github.com/tuneinsight/lattigo/v4/utils/bignum"
	"github.com/tuneinsight/lattigo/v4/utils/bignum/polynomial"
)

type poly struct {
	*polynomial.Polynomial
	MaxDeg int  // Always set to len(Coeffs)-1
	Lead   bool // Always set to true
	Lazy   bool // Flag for lazy-relinearization
}

func newPolynomial(p *polynomial.Polynomial) *poly {
	return &poly{
		Polynomial: p,
		MaxDeg:     p.Degree(),
		Lead:       true,
	}
}

type polynomialVector struct {
	Value      []*poly
	SlotsIndex map[int][]int
}

// checkEnoughLevels checks that enough levels are available to evaluate the polynomial.
// Also checks if c is a Gaussian integer or not. If not, then one more level is needed
// to evaluate the polynomial.
func checkEnoughLevels(levels, depth int) (err error) {

	if levels < depth {
		return fmt.Errorf("%d levels < %d log(d) -> cannot evaluate", levels, depth)
	}

	return nil
}

type polynomialEvaluator struct {
	*Evaluator
	PowerBasis
	slotsIndex map[int][]int
	logDegree  int
	logSplit   int
	isOdd      bool
	isEven     bool
}

// EvaluatePoly evaluates a polynomial in standard basis on the input Ciphertext in ceil(log2(deg+1)) levels.
// Returns an error if the input ciphertext does not have enough level to carry out the full polynomial evaluation.
// Returns an error if something is wrong with the scale.
// If the polynomial is given in Chebyshev basis, then a change of basis ct' = (2/(b-a)) * (ct + (-a-b)/(b-a))
// is necessary before the polynomial evaluation to ensure correctness.
// input must be either *rlwe.Ciphertext or *PolynomialBasis.
// pol: a *Polynomial
// targetScale: the desired output scale. This value shouldn't differ too much from the original ciphertext scale. It can
// for example be used to correct small deviations in the ciphertext scale and reset it to the default scale.
func (eval *Evaluator) EvaluatePoly(input interface{}, p *polynomial.Polynomial, targetScale rlwe.Scale) (opOut *rlwe.Ciphertext, err error) {
	return eval.evaluatePolyVector(input, polynomialVector{Value: []*poly{newPolynomial(p)}}, targetScale)
}

// EvaluatePolyVector evaluates a vector of Polynomials on the input Ciphertext in ceil(log2(deg+1)) levels.
// Returns an error if the input Ciphertext does not have enough level to carry out the full polynomial evaluation.
// Returns an error if something is wrong with the scale.
// Returns an error if polynomials are not all in the same basis.
// Returns an error if polynomials do not all have the same degree.
// If the polynomials are given in Chebyshev basis, then a change of basis ct' = (2/(b-a)) * (ct + (-a-b)/(b-a))
// is necessary before the polynomial evaluation to ensure correctness.
// input: must be either *rlwe.Ciphertext or *PolynomialBasis.
// pols: a slice of up to 'n' *Polynomial ('n' being the maximum number of slots), indexed from 0 to n-1.
// slotsIndex: a map[int][]int indexing as key the polynomial to evaluate and as value the index of the slots on which to evaluate the polynomial indexed by the key.
// targetScale: the desired output scale. This value shouldn't differ too much from the original ciphertext scale. It can
// for example be used to correct small deviations in the ciphertext scale and reset it to the default scale.
//
// Example: if pols = []*Polynomial{pol0, pol1} and slotsIndex = map[int][]int:{0:[1, 2, 4, 5, 7], 1:[0, 3]},
// then pol0 will be applied to slots [1, 2, 4, 5, 7], pol1 to slots [0, 3] and the slot 6 will be zero-ed.
func (eval *Evaluator) EvaluatePolyVector(input interface{}, polys []*polynomial.Polynomial, slotsIndex map[int][]int, targetScale rlwe.Scale) (opOut *rlwe.Ciphertext, err error) {
	var maxDeg int
	var basis polynomial.Basis
	for i := range polys {
		maxDeg = utils.Max(maxDeg, polys[i].Degree())
		basis = polys[i].Basis
	}

	for i := range polys {
		if basis != polys[i].Basis {
			return nil, fmt.Errorf("polynomial basis must be the same for all polynomials in a polynomial vector")
		}

		if maxDeg != polys[i].Degree() {
			return nil, fmt.Errorf("polynomial degree must all be the same")
		}
	}

	polyvec := make([]*poly, len(polys))

	for i := range polys {
		polyvec[i] = newPolynomial(polys[i])
	}

	return eval.evaluatePolyVector(input, polynomialVector{Value: polyvec, SlotsIndex: slotsIndex}, targetScale)
}

func optimalSplit(logDegree int) (logSplit int) {
	logSplit = logDegree >> 1
	a := (1 << logSplit) + (1 << (logDegree - logSplit)) + logDegree - logSplit - 3
	b := (1 << (logSplit + 1)) + (1 << (logDegree - logSplit - 1)) + logDegree - logSplit - 4
	if a > b {
		logSplit++
	}

	return
}

func (eval *Evaluator) evaluatePolyVector(input interface{}, pol polynomialVector, targetScale rlwe.Scale) (opOut *rlwe.Ciphertext, err error) {

	var powerbasis *PowerBasis
	switch input := input.(type) {
	case *rlwe.Ciphertext:
		powerbasis = NewPowerBasis(input, pol.Value[0].Basis)
	case *PowerBasis:
		if input.Value[1] == nil {
			return nil, fmt.Errorf("cannot evaluatePolyVector: given PowerBasis.Value[1] is empty")
		}
		powerbasis = input
	default:
		return nil, fmt.Errorf("cannot evaluatePolyVector: invalid input, must be either *rlwe.Ciphertext or *PowerBasis")
	}

	nbModuliPerRescale := eval.params.DefaultScaleModuliRatio()

	if err := checkEnoughLevels(powerbasis.Value[1].Level(), nbModuliPerRescale*pol.Value[0].Depth()); err != nil {
		return nil, err
	}

	logDegree := bits.Len64(uint64(pol.Value[0].Degree()))
	logSplit := optimalSplit(logDegree)

	var odd, even bool = false, false
	for _, p := range pol.Value {
		odd, even = odd || p.IsOdd, even || p.IsEven
	}

	// Computes all the powers of two with relinearization
	// This will recursively compute and store all powers of two up to 2^logDegree
	if err = powerbasis.GenPower(1<<logDegree, false, targetScale, eval); err != nil {
		return nil, err
	}

	// Computes the intermediate powers, starting from the largest, without relinearization if possible
	for i := (1 << logSplit) - 1; i > 2; i-- {
		if !(even || odd) || (i&1 == 0 && even) || (i&1 == 1 && odd) {
			if err = powerbasis.GenPower(i, pol.Value[0].Lazy, targetScale, eval); err != nil {
				return nil, err
			}
		}
	}

	polyEval := &polynomialEvaluator{}
	polyEval.slotsIndex = pol.SlotsIndex
	polyEval.Evaluator = eval
	polyEval.PowerBasis = *powerbasis
	polyEval.logDegree = logDegree
	polyEval.logSplit = logSplit
	polyEval.isOdd = odd
	polyEval.isEven = even

	if opOut, err = polyEval.recurse(powerbasis.Value[1].Level()-nbModuliPerRescale*(logDegree-1), targetScale, pol); err != nil {
		return nil, err
	}

	if opOut.Degree() == 2 {
		polyEval.Relinearize(opOut, opOut)
	}

	if err = polyEval.Rescale(opOut, targetScale, opOut); err != nil {
		return nil, err
	}

	opOut.Scale = targetScale

	polyEval = nil
	runtime.GC()
	return opOut, err
}

func (p *poly) factorize(n int) (pq, pr *poly) {

	ppq, ppr := p.Polynomial.Factorize(n)

	pq = &poly{Polynomial: ppq}
	pr = &poly{Polynomial: ppr}

	pq.MaxDeg = p.MaxDeg

	if p.MaxDeg == p.Degree() {
		pr.MaxDeg = n - 1
	} else {
		pr.MaxDeg = p.MaxDeg - (p.Degree() - n + 1)
	}

	if p.Lead {
		pq.Lead = true
	}

	return
}

func (p *polynomialVector) factorize(n int) (polyq, polyr polynomialVector) {

	coeffsq := make([]*poly, len(p.Value))
	coeffsr := make([]*poly, len(p.Value))

	for i, p := range p.Value {
		coeffsq[i], coeffsr[i] = p.factorize(n)
	}

	return polynomialVector{Value: coeffsq}, polynomialVector{Value: coeffsr}
}

func (polyEval *polynomialEvaluator) recurse(targetLevel int, targetScale rlwe.Scale, pol polynomialVector) (res *rlwe.Ciphertext, err error) {

	params := polyEval.Evaluator.params

	logSplit := polyEval.logSplit

	nbModuliPerRescale := params.DefaultScaleModuliRatio()

	// Recursively computes the evaluation of the Chebyshev polynomial using a baby-set giant-step algorithm.
	if pol.Value[0].Degree() < (1 << logSplit) {

		if pol.Value[0].Lead && polyEval.logSplit > 1 && pol.Value[0].MaxDeg%(1<<(logSplit+1)) > (1<<(logSplit-1)) {

			logDegree := int(bits.Len64(uint64(pol.Value[0].Degree())))
			logSplit := logDegree >> 1

			polyEvalBis := new(polynomialEvaluator)
			polyEvalBis.Evaluator = polyEval.Evaluator
			polyEvalBis.slotsIndex = polyEval.slotsIndex
			polyEvalBis.logDegree = logDegree
			polyEvalBis.logSplit = logSplit
			polyEvalBis.PowerBasis = polyEval.PowerBasis
			polyEvalBis.isOdd = polyEval.isOdd
			polyEvalBis.isEven = polyEval.isEven

			return polyEvalBis.recurse(targetLevel, targetScale, pol)
		}

		if pol.Value[0].Lead {

			targetScale = targetScale.Mul(rlwe.NewScale(params.Q()[targetLevel]))

			for i := 1; i < nbModuliPerRescale; i++ {
				targetScale = targetScale.Mul(rlwe.NewScale(params.Q()[targetLevel-i]))
			}
		}

		return polyEval.evaluatePolyFromPowerBasis(targetScale, targetLevel, pol)
	}

	var nextPower = 1 << polyEval.logSplit
	for nextPower < (pol.Value[0].Degree()>>1)+1 {
		nextPower <<= 1
	}

	coeffsq, coeffsr := pol.factorize(nextPower)

	XPow := polyEval.PowerBasis.Value[nextPower]

	level := targetLevel

	var qi *big.Int
	if pol.Value[0].Lead {
		qi = bignum.NewInt(params.Q()[level])
		for i := 1; i < nbModuliPerRescale; i++ {
			qi.Mul(qi, bignum.NewInt(params.Q()[level-i]))
		}
	} else {
		qi = bignum.NewInt(params.Q()[level+nbModuliPerRescale])
		for i := 1; i < nbModuliPerRescale; i++ {
			qi.Mul(qi, bignum.NewInt(params.Q()[level+nbModuliPerRescale-i]))
		}
	}

	targetScale = targetScale.Mul(rlwe.NewScale(qi))
	targetScale = targetScale.Div(XPow.Scale)

	if res, err = polyEval.recurse(targetLevel+nbModuliPerRescale, targetScale, coeffsq); err != nil {
		return nil, err
	}
	if res.Degree() == 2 {
		polyEval.Relinearize(res, res)
	}

	if err = polyEval.Rescale(res, params.DefaultScale(), res); err != nil {
		return nil, err
	}

	polyEval.Mul(res, XPow, res)

	var tmp *rlwe.Ciphertext
	if tmp, err = polyEval.recurse(res.Level(), res.Scale, coeffsr); err != nil {
		return nil, err
	}

	polyEval.Add(res, tmp, res)

	tmp = nil

	return
}

func (polyEval *polynomialEvaluator) evaluatePolyFromPowerBasis(targetScale rlwe.Scale, level int, pol polynomialVector) (res *rlwe.Ciphertext, err error) {

	// Map[int] of the powers [X^{0}, X^{1}, X^{2}, ...]
	X := polyEval.PowerBasis.Value

	// Retrieve the number of slots
	logSlots := X[1].LogSlots
	slots := 1 << X[1].LogSlots

	params := polyEval.Evaluator.params
	slotsIndex := polyEval.slotsIndex

	// Retrieve the degree of the highest degree non-zero coefficient
	// TODO: optimize for nil/zero coefficients
	minimumDegreeNonZeroCoefficient := len(pol.Value[0].Coeffs) - 1
	if polyEval.isEven && !polyEval.isOdd {
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

	// Retrieve flags for even/odd
	even := polyEval.isEven
	odd := polyEval.isOdd

	// If an index slot is given (either multiply polynomials or masking)
	if slotsIndex != nil {

		var toEncode bool

		// Allocates temporary buffer for coefficients encoding
		values := make([]*bignum.Complex, slots)

		// If the degree of the poly is zero
		if minimumDegreeNonZeroCoefficient == 0 {

			// Allocates the output ciphertext
			res = NewCiphertext(params, 1, level)
			res.Scale = targetScale
			res.LogSlots = logSlots

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
				if err = polyEval.Evaluator.Encode(values, pt); err != nil {
					return nil, err
				}
			}

			return
		}

		// Allocates the output ciphertext
		res = NewCiphertext(params, maximumCiphertextDegree, level)
		res.Scale = targetScale
		res.LogSlots = logSlots

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
			polyEval.Add(res, values, res)
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
				polyEval.MulThenAdd(X[key], values, res)
				toEncode = false
			}
		}

	} else {

		var c *bignum.Complex
		if even && !isZero(pol.Value[0].Coeffs[0]) {
			c = pol.Value[0].Coeffs[0]
		}

		if minimumDegreeNonZeroCoefficient == 0 {

			res = NewCiphertext(params, 1, level)
			res.Scale = targetScale
			res.LogSlots = logSlots

			if !isZero(c) {
				polyEval.Add(res, c, res)
			}

			return
		}

		res = NewCiphertext(params, maximumCiphertextDegree, level)
		res.Scale = targetScale
		res.LogSlots = logSlots

		if c != nil {
			polyEval.Add(res, c, res)
		}

		for key := pol.Value[0].Degree(); key > 0; key-- {
			if c = pol.Value[0].Coeffs[key]; key != 0 && !isZero(c) && (!(even || odd) || (key&1 == 0 && even) || (key&1 == 1 && odd)) {
				polyEval.Evaluator.MulThenAdd(X[key], c, res)
			}
		}
	}

	return
}

func isZero(c *bignum.Complex) bool {
	zero := new(big.Float)
	return c == nil || (c[0].Cmp(zero) == 0 && c[1].Cmp(zero) == 0)
}
