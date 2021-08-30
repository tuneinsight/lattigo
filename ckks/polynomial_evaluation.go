package ckks

import (
	"fmt"
	"github.com/ldsec/lattigo/v2/ring"
	"math"
	"math/bits"
)

// Poly is a struct storing the coeffients of a polynomial
// that then can be evaluated on the ciphertext
type Poly struct {
	MaxDeg int
	Coeffs []complex128
	Lead   bool
	A      complex128
	B      complex128
	Basis  PolynomialBasis
}

// PolynomialBasis is a type for the polynomials basis
type PolynomialBasis int

const (
	// StandardBasis : x^(a+b) = x^a * x^b
	StandardBasis = PolynomialBasis(0)
	// ChebyshevBasis : x^(a+b) = 2 * x^a *x^b - T^|a-b|
	ChebyshevBasis = PolynomialBasis(1)
)

// Depth returns the number of levels needed to evaluate the polynomial.
func (p *Poly) Depth() int {
	return int(math.Ceil(math.Log2(float64(len(p.Coeffs)))))
}

// Degree returns the degree of the polynomial
func (p *Poly) Degree() int {
	return len(p.Coeffs) - 1
}

// NewPoly creates a new Poly from the input coefficients
func NewPoly(coeffs []complex128) (p *Poly) {

	p = new(Poly)
	p.Coeffs = make([]complex128, len(coeffs))
	copy(p.Coeffs, coeffs)
	p.MaxDeg = len(coeffs) - 1
	p.Lead = true

	return
}

// checkEnoughLevels checks that enough levels are available to evaluate the polynomial.
// Also checks if c is a gaussian integer or not. If not, then one more level is needed
// to evaluate the polynomial.
func checkEnoughLevels(levels int, pol *Poly, c complex128) (err error) {

	depth := pol.Depth()

	if real(c) != float64(int64(real(c))) || imag(c) != float64(int64(imag(c))) {
		depth++
	}

	if levels < depth {
		return fmt.Errorf("%d levels < %d log(d) -> cannot evaluate", levels, depth)
	}

	return nil
}

// EvaluatePoly evaluates a polynomial in standard basis on the input Ciphertext in ceil(log2(deg+1)) levels.
// Returns an error if the input ciphertext does not have enough level to carry out the full polynomial evaluation.
// Returns an error if something is wrong with the scale.
// If the polynomial is given in Chebyshev basis, then a change of basis ct' = (2/(b-a)) * (ct + (-a-b)/(b-a))
// is necessary before the polynomial evaluation to ensure correctness.
func (eval *evaluator) EvaluatePoly(ct0 *Ciphertext, pol *Poly, targetScale float64) (opOut *Ciphertext, err error) {

	if err := checkEnoughLevels(ct0.Level(), pol, 1); err != nil {
		return ct0, err
	}

	C := make(map[int]*Ciphertext)

	C[1] = ct0.CopyNew()

	logDegree := bits.Len64(uint64(pol.Degree()))
	logSplit := (logDegree >> 1) //optimalSplit(logDegree) //

	odd, even := isOddOrEvenPolynomial(pol.Coeffs)

	for i := 2; i < (1 << logSplit); i++ {
		if i&1 == 0 && even {
			if err = computePowerBasis(i, C, targetScale, pol.Basis, eval); err != nil {
				return nil, err
			}
		} else if i&1 == 1 && odd {
			if err = computePowerBasis(i, C, targetScale, pol.Basis, eval); err != nil {
				return nil, err
			}
		} else if !even && !odd {
			if err = computePowerBasis(i, C, targetScale, pol.Basis, eval); err != nil {
				return nil, err
			}
		}
	}

	for i := logSplit; i < logDegree; i++ {
		if err = computePowerBasis(1<<i, C, targetScale, pol.Basis, eval); err != nil {
			return nil, err
		}
	}

	opOut, err = recurse(targetScale, logSplit, logDegree, pol, C, eval)

	opOut.Scale = targetScale // solves float64 precision issues

	C = nil
	return opOut, err
}

func computePowerBasis(n int, C map[int]*Ciphertext, scale float64, basis PolynomialBasis, evaluator *evaluator) (err error) {

	if C[n] == nil {

		// Computes the index required to compute the asked ring evaluation
		var a, b, c int
		if n&(n-1) == 0 {
			a, b = n/2, n/2 //Necessary for depth optimality
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

func splitCoeffs(coeffs *Poly, split int) (coeffsq, coeffsr *Poly) {

	// Splits a polynomial p such that p = q*C^degree + r.

	coeffsr = new(Poly)
	coeffsr.Coeffs = make([]complex128, split)
	if coeffs.MaxDeg == coeffs.Degree() {
		coeffsr.MaxDeg = split - 1
	} else {
		coeffsr.MaxDeg = coeffs.MaxDeg - (coeffs.Degree() - split + 1)
	}

	for i := 0; i < split; i++ {
		coeffsr.Coeffs[i] = coeffs.Coeffs[i]
	}

	coeffsq = new(Poly)
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

	return coeffsq, coeffsr
}

func recurse(targetScale float64, logSplit, logDegree int, coeffs *Poly, C map[int]*Ciphertext, evaluator *evaluator) (res *Ciphertext, err error) {

	// Recursively computes the evalution of the Chebyshev polynomial using a baby-set giant-step algorithm.
	if coeffs.Degree() < (1 << logSplit) {

		if coeffs.Lead && logSplit > 1 && coeffs.MaxDeg%(1<<(logSplit+1)) > (1<<(logSplit-1)) {

			logDegree = int(bits.Len64(uint64(coeffs.Degree())))
			logSplit = logDegree >> 1

			return recurse(targetScale, logSplit, logDegree, coeffs, C, evaluator)
		}

		return evaluatePolyFromPowerBasis(targetScale, coeffs, logSplit, C, evaluator)
	}

	var nextPower = 1 << logSplit
	for nextPower < (coeffs.Degree()>>1)+1 {
		nextPower <<= 1
	}

	coeffsq, coeffsr := splitCoeffs(coeffs, nextPower)

	level := C[nextPower].Level() - 1

	if coeffsq.MaxDeg >= 1<<(logDegree-1) && coeffsq.Lead {
		level++
	}

	currentQi := evaluator.params.QiFloat64(level)

	if res, err = recurse(targetScale*currentQi/C[nextPower].Scale, logSplit, logDegree, coeffsq, C, evaluator); err != nil {
		return nil, err
	}

	var tmp *Ciphertext
	if tmp, err = recurse(targetScale, logSplit, logDegree, coeffsr, C, evaluator); err != nil {
		return nil, err
	}

	if res.Level() > tmp.Level() {
		for res.Level() != tmp.Level()+1 {
			evaluator.DropLevel(res, 1)
		}
	}

	evaluator.MulRelin(res, C[nextPower], res)

	if res.Level() > tmp.Level() {
		if err = evaluator.Rescale(res, targetScale, res); err != nil {
			return nil, err
		}
		evaluator.Add(res, tmp, res)
	} else {
		evaluator.Add(res, tmp, res)
		if err = evaluator.Rescale(res, targetScale, res); err != nil {
			return nil, err
		}
	}

	tmp = nil

	return
}

func evaluatePolyFromPowerBasis(targetScale float64, coeffs *Poly, logSplit int, C map[int]*Ciphertext, evaluator *evaluator) (res *Ciphertext, err error) {

	minimumDegreeNonZeroCoefficient := 0

	for i := coeffs.Degree(); i > 0; i-- {
		if isNotNegligible(coeffs.Coeffs[i]) {
			minimumDegreeNonZeroCoefficient = i
			break
		}
	}

	c := coeffs.Coeffs[0]

	if minimumDegreeNonZeroCoefficient == 0 {

		res = NewCiphertext(evaluator.params, 1, C[1].Level(), targetScale)

		if isNotNegligible(c) {
			evaluator.AddConst(res, c, res)
		}

		return
	}

	minimumDegreeNonZeroCoefficient = coeffs.Degree()

	currentQi := evaluator.params.QiFloat64(C[(minimumDegreeNonZeroCoefficient)].Level())

	ctScale := targetScale * currentQi

	res = NewCiphertext(evaluator.params, 1, C[minimumDegreeNonZeroCoefficient].Level(), ctScale)

	if isNotNegligible(c) {
		evaluator.AddConst(res, c, res)
	}

	cRealFlo, cImagFlo, constScale := ring.NewFloat(0, 128), ring.NewFloat(0, 128), ring.NewFloat(0, 128)
	cRealBig, cImagBig := ring.NewUint(0), ring.NewUint(0)

	for key := coeffs.Degree(); key > 0; key-- {

		c = coeffs.Coeffs[key]

		if key != 0 && isNotNegligible(c) {

			cRealFlo.SetFloat64(real(c))
			cImagFlo.SetFloat64(imag(c))
			constScale.SetFloat64(targetScale * currentQi / C[key].Scale)

			// Target scale * rescale-scale / power basis scale
			cRealFlo.Mul(cRealFlo, constScale).Int(cRealBig)
			cImagFlo.Mul(cImagFlo, constScale).Int(cImagBig)

			evaluator.MultByGaussianIntegerAndAdd(C[key], cRealBig, cImagBig, res)
		}
	}

	if err = evaluator.Rescale(res, targetScale, res); err != nil {
		return nil, err
	}

	return
}

func isNotNegligible(c complex128) bool {
	return (math.Abs(real(c)) > 1e-14 || math.Abs(imag(c)) > 1e-14)
}

func isOddOrEvenPolynomial(coeffs []complex128) (odd, even bool) {
	even = true
	odd = true
	for i, c := range coeffs {
		isnotnegligible := isNotNegligible(c)
		if i&1 == 0 && isnotnegligible {
			odd = false
		}

		if i&1 == 1 && isnotnegligible {
			even = false
		}

		if !odd && !even {
			break
		}
	}

	if even {
		for i := 1; i < len(coeffs); i += 2 {
			coeffs[i] = complex(0, 0)
		}
	} else if odd {
		for i := 0; i < len(coeffs); i += 2 {
			coeffs[i] = complex(0, 0)
		}
	}

	return
}
