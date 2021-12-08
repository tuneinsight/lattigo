package ckks

import (
	"fmt"
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/utils"
	"math"
	"math/big"
	"math/bits"
	"runtime"
)

// Polynomial is a struct storing the coefficients of a polynomial
// that then can be evaluated on the ciphertext
type Polynomial struct {
	MaxDeg int
	Coeffs []complex128
	Lead   bool
	A      float64
	B      float64
	Basis  PolynomialBasis
}

// PolynomialVector is a struct storing a vector of polynomials
// that can be evaluated on a ciphertext.
type PolynomialVector struct {
	Encoder    Encoder
	Value      []*Polynomial
	SlotsIndex map[int][]int
}

// NewPolynomialVector creates a new PolynomialVector, whose degree
// is the maximum degree among all the polynomials.
// pols : a list of polynomials (indexed)
// slotsIndex : a map with as input the index of the polynomial to apply and output the index of the slots
// to apply the polynomial to.
// endoder : a ckks.Encoder
// The method will panic if the input polynmials are not all of the same Basis or do not all have the
// same degree (padd with 0 coeffs if needed).
// Example: if pols = []*Polynomial{pol0, pol1} and slotsIndex = map[int][]int:{0:[1, 2, 4, 5, 7],
// 1:[0, 3]}, then pol0 will be applied to slots [1, 2, 4, 5, 7], pol1 to slots [0, 3] and the slot
// 6 will be multiplied by zero.
func NewPolynomialVector(pols []*Polynomial, slotsIndex map[int][]int, encoder Encoder) PolynomialVector {
	var maxDeg int
	var basis PolynomialBasis
	for i := range pols {
		maxDeg = utils.MaxInt(maxDeg, pols[i].MaxDeg)
		basis = pols[i].Basis
	}

	for i := range pols {
		if basis != pols[i].Basis {
			panic("polynomial basis must be the same for all polynomials in a polynomial vector")
		}

		if maxDeg != pols[i].MaxDeg {
			panic("polynomial degree must all be the same")
		}
	}

	return PolynomialVector{Encoder: encoder, Value: pols, SlotsIndex: slotsIndex}
}

// PolynomialBasis is a type for the polynomials basis
type PolynomialBasis int

const (
	// StandardBasis : x^(a+b) = x^a * x^b
	StandardBasis = PolynomialBasis(0)
	// ChebyshevBasis : x^(a+b) = 2 * x^a *x^b - T^|a-b|
	ChebyshevBasis = PolynomialBasis(1)
)

// IsNegligbleThreshold : threshold under which a coefficient
// of a polynomial is ignored.
const IsNegligbleThreshold float64 = 1e-14

// Depth returns the number of levels needed to evaluate the polynomial.
func (p *Polynomial) Depth() int {
	return int(math.Ceil(math.Log2(float64(len(p.Coeffs)))))
}

// Degree returns the degree of the polynomial
func (p *Polynomial) Degree() int {
	return len(p.Coeffs) - 1
}

// NewPoly creates a new Poly from the input coefficients
func NewPoly(coeffs []complex128) (p *Polynomial) {
	c := make([]complex128, len(coeffs))
	copy(c, coeffs)
	return &Polynomial{Coeffs: c, MaxDeg: len(c) - 1, Lead: true}
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
	slotsIndex map[int][]int
	powerBasis map[int]*Ciphertext
	logDegree  int
	logSplit   int
}

// EvaluatePoly evaluates a Polynomial or a PolynomialVector in on the input Ciphertext in ceil(log2(deg+1)) levels.
// Returns an error if the input ciphertext does not have enough level to carry out the full polynomial evaluation.
// Returns an error if something is wrong with the scale.
// If the polynomial is given in Chebyshev basis, then a change of basis ct' = (2/(b-a)) * (ct + (-a-b)/(b-a))
// is necessary before the polynomial evaluation to ensure correctness.
// Coefficients of the polynomial with an absolute value smaller than "IsNegligbleThreshold" will automatically be set to zero
// if the polynomial is "even" or "odd" (to ensure that the even or odd property remains valid
// after the "splitCoeffs" polynomial decomposition).
func (eval *evaluator) EvaluatePoly(ct0 *Ciphertext, pol interface{}, targetScale float64) (opOut *Ciphertext, err error) {
	switch el := pol.(type) {
	case *Polynomial:
		return eval.evaluatePolyVector(ct0, PolynomialVector{Value: []*Polynomial{el}}, targetScale)
	case PolynomialVector:
		return eval.evaluatePolyVector(ct0, el, targetScale)
	default:
		return nil, fmt.Errorf("input polynomial must either be *ckks.Polynomial or *ckks.PolynomialVector")
	}
}

func (eval *evaluator) evaluatePolyVector(ct0 *Ciphertext, pol PolynomialVector, targetScale float64) (opOut *Ciphertext, err error) {

	if pol.SlotsIndex != nil && pol.Encoder == nil {
		return nil, fmt.Errorf("cannot EvaluatePolyVector, missing Encoder input")
	}

	if err := checkEnoughLevels(ct0.Level(), pol.Value[0].Depth(), 1); err != nil {
		return ct0, err
	}

	poweBasis := make(map[int]*Ciphertext)

	poweBasis[1] = ct0.CopyNew()

	logDegree := bits.Len64(uint64(pol.Value[0].Degree()))
	logSplit := (logDegree >> 1) //optimalSplit(logDegree) //

	var odd, even bool
	for _, p := range pol.Value {
		tmp0, tmp1 := isOddOrEvenPolynomial(p.Coeffs)
		odd, even = odd && tmp0, even && tmp1
	}

	for i := 2; i < (1 << logSplit); i++ {
		if !(even || odd) || (i&1 == 0 && even) || (i&1 == 1 && odd) {
			if err = computePowerBasis(i, poweBasis, targetScale, pol.Value[0].Basis, eval); err != nil {
				return nil, err
			}
		}
	}

	for i := logSplit; i < logDegree; i++ {
		if err = computePowerBasis(1<<i, poweBasis, targetScale, pol.Value[0].Basis, eval); err != nil {
			return nil, err
		}
	}

	polyEval := &polynomialEvaluator{}
	polyEval.slotsIndex = pol.SlotsIndex
	polyEval.Evaluator = eval
	polyEval.Encoder = pol.Encoder
	polyEval.powerBasis = poweBasis
	polyEval.logDegree = logDegree
	polyEval.logSplit = logSplit

	opOut, err = polyEval.recurse(targetScale, pol)

	opOut.Scale = targetScale // solves float64 precision issues

	polyEval = nil
	runtime.GC()
	return opOut, err
}

func computePowerBasis(n int, C map[int]*Ciphertext, scale float64, basis PolynomialBasis, evaluator *evaluator) (err error) {

	if C[n] == nil {

		// Computes the index required to compute the asked ring evaluation
		var a, b, c int
		if n&(n-1) == 0 {
			a, b = n/2, n/2 //Necessary for optimal depth
		} else {
			// [Lee et al. 2020] : High-Precision and Low-Complexity Approximate Homomorphic Encryption by Error Variance Minimization
			// Maximize the number of odd terms of Chebyshev basis
			k := int(math.Ceil(math.Log2(float64(n)))) - 1
			a = (1 << k) - 1
			b = n + 1 - (1 << k)

			if basis == ChebyshevBasis {
				c = int(math.Abs(float64(a) - float64(b))) // Cn = 2*Ca*Cb - Cc, n = a+b and c = abs(a-b)
			}
		}
		// Recurses on the given indexes
		if err = computePowerBasis(a, C, scale, basis, evaluator); err != nil {
			return err
		}
		if err = computePowerBasis(b, C, scale, basis, evaluator); err != nil {
			return err
		}

		// Computes C[n] = C[a]*C[b]
		C[n] = evaluator.MulRelinNew(C[a], C[b])

		if err = evaluator.Rescale(C[n], scale, C[n]); err != nil {
			return err
		}

		if basis == ChebyshevBasis {

			// Computes C[n] = 2*C[a]*C[b]
			evaluator.Add(C[n], C[n], C[n])

			// Computes C[n] = 2*C[a]*C[b] - C[c]
			if c == 0 {
				evaluator.AddConst(C[n], -1, C[n])
			} else {
				// Since C[0] is not stored (but rather seen as the constant 1), only recurses on c if c!= 0
				if err = computePowerBasis(c, C, scale, basis, evaluator); err != nil {
					return err
				}
				evaluator.Sub(C[n], C[c], C[n])
			}
		}
	}

	return nil
}

func splitCoeffs(coeffs *Polynomial, split int) (coeffsq, coeffsr *Polynomial) {

	// Splits a polynomial p such that p = q*C^degree + r.
	coeffsr = &Polynomial{}
	coeffsr.Coeffs = make([]complex128, split)
	if coeffs.MaxDeg == coeffs.Degree() {
		coeffsr.MaxDeg = split - 1
	} else {
		coeffsr.MaxDeg = coeffs.MaxDeg - (coeffs.Degree() - split + 1)
	}

	for i := 0; i < split; i++ {
		coeffsr.Coeffs[i] = coeffs.Coeffs[i]
	}

	coeffsq = &Polynomial{}
	coeffsq.Coeffs = make([]complex128, coeffs.Degree()-split+1)
	coeffsq.MaxDeg = coeffs.MaxDeg

	coeffsq.Coeffs[0] = coeffs.Coeffs[split]

	if coeffs.Basis == StandardBasis {
		for i := split + 1; i < coeffs.Degree()+1; i++ {
			coeffsq.Coeffs[i-split] = coeffs.Coeffs[i]
		}
	} else if coeffs.Basis == ChebyshevBasis {
		for i, j := split+1, 1; i < coeffs.Degree()+1; i, j = i+1, j+1 {
			coeffsq.Coeffs[i-split] = 2 * coeffs.Coeffs[i]
			coeffsr.Coeffs[split-j] -= coeffs.Coeffs[i]
		}
	}

	if coeffs.Lead {
		coeffsq.Lead = true
	}

	coeffsq.Basis, coeffsr.Basis = coeffs.Basis, coeffs.Basis

	return
}

func splitCoeffsPolyVector(poly PolynomialVector, split int) (polyq, polyr PolynomialVector) {
	coeffsq := make([]*Polynomial, len(poly.Value))
	coeffsr := make([]*Polynomial, len(poly.Value))
	for i, p := range poly.Value {
		coeffsq[i], coeffsr[i] = splitCoeffs(p, split)
	}

	return PolynomialVector{Value: coeffsq}, PolynomialVector{Value: coeffsr}
}

func (polyEval *polynomialEvaluator) recurse(targetScale float64, pol PolynomialVector) (res *Ciphertext, err error) {

	logSplit := polyEval.logSplit

	// Recursively computes the evaluation of the Chebyshev polynomial using a baby-set giant-step algorithm.
	if pol.Value[0].Degree() < (1 << logSplit) {

		if pol.Value[0].Lead && polyEval.logSplit > 1 && pol.Value[0].MaxDeg%(1<<(logSplit+1)) > (1<<(logSplit-1)) {

			logDegree := int(bits.Len64(uint64(pol.Value[0].Degree())))
			logSplit := logDegree >> 1

			polyEvalBis := new(polynomialEvaluator)
			polyEvalBis.Evaluator = polyEval.Evaluator
			polyEvalBis.Encoder = polyEval.Encoder
			polyEvalBis.logDegree = logDegree
			polyEvalBis.logSplit = logSplit
			polyEvalBis.powerBasis = polyEval.powerBasis

			return polyEvalBis.recurse(targetScale, pol)
		}

		return polyEval.evaluatePolyFromPowerBasis(targetScale, pol)
	}

	var nextPower = 1 << polyEval.logSplit
	for nextPower < (pol.Value[0].Degree()>>1)+1 {
		nextPower <<= 1
	}

	coeffsq, coeffsr := splitCoeffsPolyVector(pol, nextPower)

	XPow := polyEval.powerBasis[nextPower]

	level := XPow.Level() - 1

	if coeffsq.Value[0].MaxDeg >= 1<<(polyEval.logDegree-1) && coeffsq.Value[0].Lead {
		level++
	}

	currentQi := polyEval.Evaluator.(*evaluator).params.QiFloat64(level)

	if res, err = polyEval.recurse(targetScale*currentQi/XPow.Scale, coeffsq); err != nil {
		return nil, err
	}

	var tmp *Ciphertext
	if tmp, err = polyEval.recurse(targetScale, coeffsr); err != nil {
		return nil, err
	}

	if res.Level() > tmp.Level() {
		for res.Level() != tmp.Level()+1 {
			polyEval.DropLevel(res, 1)
		}
	}

	polyEval.MulRelin(res, XPow, res)

	if res.Level() > tmp.Level() {
		if err = polyEval.Rescale(res, targetScale, res); err != nil {
			return nil, err
		}
		polyEval.Add(res, tmp, res)
	} else {
		polyEval.Add(res, tmp, res)
		if err = polyEval.Rescale(res, targetScale, res); err != nil {
			return nil, err
		}
	}

	tmp = nil

	return
}

func (polyEval *polynomialEvaluator) evaluatePolyFromPowerBasis(targetScale float64, pol PolynomialVector) (res *Ciphertext, err error) {

	X := polyEval.powerBasis

	params := polyEval.Evaluator.(*evaluator).params
	slotsIndex := polyEval.slotsIndex

	minimumDegreeNonZeroCoefficient := 0

	// Get the minimum non-zero degree coefficient
	for i := pol.Value[0].Degree(); i > 0; i-- {
		for _, p := range pol.Value {
			if isNotNegligible(p.Coeffs[i]) {
				minimumDegreeNonZeroCoefficient = utils.MaxInt(minimumDegreeNonZeroCoefficient, i)
				break
			}
		}
	}

	// If an index slot is given (either multiply polynomials or masking)
	if slotsIndex != nil {

		var toEncode bool

		// Allocates temporary buffer for coefficients encoding
		values := make([]complex128, params.Slots())

		// If the degree of the poly is zero
		if minimumDegreeNonZeroCoefficient == 0 {

			// Allocates the output ciphertext
			res = NewCiphertext(params, 1, X[1].Level(), targetScale)

			// Looks for non-zero coefficients among the degree 0 coefficients of the polynomials
			for i, p := range pol.Value {
				if isNotNegligible(p.Coeffs[0]) {
					toEncode = true
					for _, j := range slotsIndex[i] {
						values[j] = p.Coeffs[0]
					}
				}
			}

			// If a non-zero coefficient was found, encode the values, adds on the ciphertext, and returns
			if toEncode {
				pt := NewPlaintextAtLevelFromPoly(X[1].Level(), res.Value[0])
				pt.Scale = res.Scale
				polyEval.EncodeSlots(values, pt, params.LogSlots())
			}

			return
		}

		// Current modulus qi for the next rescale
		currentQi := params.QiFloat64(X[(minimumDegreeNonZeroCoefficient)].Level())

		// Target output scale before the rescaling, such that the rescaling does not change the ciphertext scale
		ctScale := targetScale * currentQi

		// Allocates the output ciphertext
		res = NewCiphertext(params, 1, X[minimumDegreeNonZeroCoefficient].Level(), ctScale)

		// Allocates a temporary plaintext to encode the values
		pt := NewPlaintextAtLevelFromPoly(X[minimumDegreeNonZeroCoefficient].Level(), polyEval.Evaluator.CtxPool().Value[0])

		// Looks for a non-zero coefficient among the degree zero coefficient of the polynomials
		for i, p := range pol.Value {
			if isNotNegligible(p.Coeffs[0]) {
				toEncode = true
				for _, j := range slotsIndex[i] {
					values[j] = p.Coeffs[0]
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
		for key := pol.Value[0].Degree(); key > 0; key-- {

			var reset bool
			// Loops over the polynomials
			for i, p := range pol.Value {

				// Looks for a non-zero coefficient
				if isNotNegligible(p.Coeffs[key]) {
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
						values[j] = p.Coeffs[key]
					}
				}
			}

			// If a non-zero degre coefficient was found, encode and adds the values on the output
			// ciphertext
			if toEncode {
				pt.Scale = targetScale * currentQi / X[key].Scale
				polyEval.EncodeSlots(values, pt, params.LogSlots())
				polyEval.MulAndAdd(X[key], pt, res)
				toEncode = false
			}
		}

	} else {

		c := pol.Value[0].Coeffs[0]

		if minimumDegreeNonZeroCoefficient == 0 {

			res = NewCiphertext(params, 1, X[1].Level(), targetScale)

			if isNotNegligible(c) {
				polyEval.AddConst(res, c, res)
			}

			return
		}

		currentQi := params.QiFloat64(X[(minimumDegreeNonZeroCoefficient)].Level())
		ctScale := targetScale * currentQi
		res = NewCiphertext(params, 1, X[minimumDegreeNonZeroCoefficient].Level(), ctScale)

		if isNotNegligible(c) {
			polyEval.AddConst(res, c, res)
		}

		cRealFlo, cImagFlo, constScale := ring.NewFloat(0, 128), ring.NewFloat(0, 128), ring.NewFloat(0, 128)
		cRealBig, cImagBig := ring.NewUint(0), ring.NewUint(0)

		for key := pol.Value[0].Degree(); key > 0; key-- {

			c = pol.Value[0].Coeffs[key]

			if key != 0 && isNotNegligible(c) {

				cRealFlo.SetFloat64(real(c))
				cImagFlo.SetFloat64(imag(c))
				constScale.SetFloat64(targetScale * currentQi / X[key].Scale)

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

	if err = polyEval.Rescale(res, targetScale, res); err != nil {
		return nil, err
	}

	return
}

func isNotNegligible(c complex128) bool {
	return (math.Abs(real(c)) > IsNegligbleThreshold || math.Abs(imag(c)) > IsNegligbleThreshold)
}

func isOddOrEvenPolynomial(coeffs []complex128) (odd, even bool) {
	even = true
	odd = true
	for i, c := range coeffs {
		isnotnegligible := isNotNegligible(c)
		odd = odd && !(i&1 == 0 && isnotnegligible)
		even = even && !(i&1 == 1 && isnotnegligible)
		if !odd && !even {
			break
		}
	}

	// If even or odd, then sets the expected zero coefficients to zero
	if even || odd {
		var start int
		if even {
			start = 1
		}
		for i := start; i < len(coeffs); i += 2 {
			coeffs[i] = complex(0, 0)
		}
	}

	return
}
