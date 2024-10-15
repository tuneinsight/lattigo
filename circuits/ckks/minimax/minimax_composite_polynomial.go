// Package minimax implements a homomorphic minimax circuit for the CKKS scheme.
package minimax

import (
	"fmt"
	"math"
	"math/big"

	"github.com/tuneinsight/lattigo/v6/utils"
	"github.com/tuneinsight/lattigo/v6/utils/bignum"
)

// Polynomial is a struct storing P(x) = pk(x) o pk-1(x) o ... o p1(x) o p0(x).
type Polynomial []bignum.Polynomial

// NewPolynomial creates a new Polynomial from a list of coefficients.
// Coefficients are expected to be given in the Chebyshev basis.
func NewPolynomial(coeffsStr [][]string) Polynomial {
	polys := make([]bignum.Polynomial, len(coeffsStr))

	for i := range coeffsStr {

		coeffs := parseCoeffs(coeffsStr[i])

		poly := bignum.NewPolynomial(
			bignum.Chebyshev,
			coeffs,
			&bignum.Interval{
				A: *bignum.NewFloat(-1, coeffs[0].Prec()),
				B: *bignum.NewFloat(1, coeffs[0].Prec()),
			},
		)

		polys[i] = poly
	}

	return Polynomial(polys)
}

func (mcp Polynomial) MaxDepth() (depth int) {
	for i := range mcp {
		depth = utils.Max(depth, mcp[i].Depth())
	}
	return
}

func (mcp Polynomial) Evaluate(x interface{}) (y *bignum.Complex) {
	y = mcp[0].Evaluate(x)

	for _, p := range mcp[1:] {
		y = p.Evaluate(y)
	}

	return
}

// CoeffsSignX2Cheby (from https://eprint.iacr.org/2019/1234.pdf) are the coefficients
// of 1.5*x - 0.5*x^3 in Chebyshev basis.
// Evaluating this polynomial on values already close to -1, or 1 ~doubles the number of
// of correct digits.
// For example, if x = -0.9993209 then p(x) = -0.999999308
// This polynomial can be composed after the minimax composite polynomial to double the
// output precision (up to the scheme precision) each time it is evaluated.
var CoeffsSignX2Cheby = []string{"0", "1.125", "0", "-0.125"}

// CoeffsSignX4Cheby (from https://eprint.iacr.org/2019/1234.pdf) are the coefficients
// of 35/16 * x - 35/16 * x^3 + 21/16 * x^5 - 5/16 * x^7 in Chebyshev basis.
// Evaluating this polynomial on values already close to -1, or 1 ~quadruples the number of
// of correct digits.
// For example, if x = -0.9993209 then p(x) = -0.9999999999990705
// This polynomial can be composed after the minimax composite polynomial to quadruple the
// output precision (up to the scheme precision) each time it is evaluated.
var CoeffsSignX4Cheby = []string{"0", "1.1962890625", "0", "-0.2392578125", "0", "0.0478515625", "0", "-0.0048828125"}

// GenMinimaxCompositePolynomialForSign generates the minimax composite polynomial
// P(x) = pk(x) o pk-1(x) o ... o p1(x) o p0(x) of the sign function in their interval
// [min-err, -2^{-alpha}] U [2^{-alpha}, max+err] where alpha is the desired distinguishing
// precision between two values and err an upperbound on the scheme error.
//
// The sign function is defined as: -1 if -1 <= x < 0, 0 if x = 0, 1 if 0 < x <= 1.
//
// See [GenMinimaxCompositePolynomial] for information about how to instantiate and
// parameterize each input value of the algorithm.
func GenMinimaxCompositePolynomialForSign(prec uint, logalpha, logerr int, deg []int) {

	coeffs := GenMinimaxCompositePolynomial(prec, logalpha, logerr, deg, bignum.Sign)

	decimals := int(float64(logalpha)/math.Log2(10)+0.5) + 10

	fmt.Println("COEFFICIENTS:")
	fmt.Printf("{\n")
	for i := range coeffs {
		PrettyPrintCoefficients(decimals, coeffs[i], true, false, false)
	}
	fmt.Printf("},\n")
}

// GenMinimaxCompositePolynomial generates the minimax composite polynomial
// P(x) = pk(x) o pk-1(x) o ... o p1(x) o p0(x) for the provided function in the interval
// in their interval [min-err, -2^{-alpha}] U [2^{-alpha}, max+err] where alpha is
// the desired distinguishing precision between two values and err an upperbound on
// the scheme error.
//
// The user must provide the following inputs:
//   - prec: the bit precision of the big.Float values used by the algorithm to compute the polynomials.
//     This will impact the speed of the algorithm.
//     A too low precision can prevent convergence or induce a slope zero during the zero finding.
//     A sign that the precision is too low is when the iteration continue without the error getting smaller.
//   - logalpha: log2(alpha)
//   - logerr: log2(err), the upperbound on the scheme precision. Usually this value should be smaller or equal to logalpha.
//     Correctly setting this value is mandatory for correctness, because if x is outside of the interval
//     (i.e. smaller than -1-e or greater than 1+e), then the values will explode during the evaluation.
//     Note that it is not required to apply change of interval [-1, 1] -> [-1-e, 1+e] because the function to evaluate
//     is the sign (i.e. it will evaluate to the same value).
//   - deg: the degree of each polynomial, ordered as follow [deg(p0(x)), deg(p1(x)), ..., deg(pk(x))].
//     It is highly recommended that deg(p0) <= deg(p1) <= ... <= deg(pk) for optimal approximation.
//
// The polynomials are returned in the Chebyshev basis and pre-scaled for
// the interval [-1, 1] (no further scaling is required on the ciphertext).
//
// Be aware that finding the minimax polynomials can take a while (in the order of minutes for high precision when using large degree polynomials).
//
// The function will print information about each step of the computation in real time so that it can be monitored.
//
// The underlying algorithm use the multi-interval Remez algorithm of https://eprint.iacr.org/2020/834.pdf.
func GenMinimaxCompositePolynomial(prec uint, logalpha, logerr int, deg []int, f func(*big.Float) *big.Float) (coeffs [][]*big.Float) {
	decimals := int(float64(logalpha)/math.Log2(10)+0.5) + 10

	// Precision of the output value of the sign polynomial
	alpha := math.Exp2(-float64(logalpha))

	// Expected upperbound scheme error
	e := bignum.NewFloat(math.Exp2(-float64(logerr)), prec)

	// Maximum number of iterations
	maxIters := 50

	// Scan step for finding zeroes of the error function
	scanStep := bignum.NewFloat(1e-3, prec)

	// Interval [-1, alpha] U [alpha, 1]
	intervals := []bignum.Interval{
		{A: *bignum.NewFloat(-1, prec), B: *bignum.NewFloat(-alpha, prec), Nodes: 1 + ((deg[0] + 1) >> 1)},
		{A: *bignum.NewFloat(alpha, prec), B: *bignum.NewFloat(1, prec), Nodes: 1 + ((deg[0] + 1) >> 1)},
	}

	// Adds the error to the interval
	// [A, -alpha] U [alpha, B] becomes [A-e, -alpha] U [alpha, B+e]
	intervals[0].A.Sub(&intervals[0].A, e)
	intervals[1].B.Add(&intervals[1].B, e)

	// Parameters of the minimax approximation
	params := bignum.RemezParameters{
		Function:        f,
		Basis:           bignum.Chebyshev,
		Intervals:       intervals,
		ScanStep:        scanStep,
		Prec:            prec,
		OptimalScanStep: true,
	}

	fmt.Printf("P[0]\n")
	fmt.Printf("Interval: [%.*f, %.*f] U [%.*f, %.*f]\n", decimals, &intervals[0].A, decimals, &intervals[0].B, decimals, &intervals[1].A, decimals, &intervals[1].B)
	r := bignum.NewRemez(params)
	r.Approximate(maxIters, alpha)
	//r.ShowCoeffs(decimals)
	r.ShowError(decimals)
	fmt.Println()

	coeffs = make([][]*big.Float, len(deg))

	for i := 1; i < len(deg); i++ {

		// New interval as [-(1+max_err), -(1-min_err)] U [1-min_err, 1+max_err]
		maxInterval := bignum.NewFloat(1, prec)
		maxInterval.Add(maxInterval, r.MaxErr)

		minInterval := bignum.NewFloat(1, prec)
		minInterval.Sub(minInterval, r.MinErr)

		// Extends the new interval by the scheme error
		// [-(1+max_err), -(1-min_err)] U [1-min_err, 1 + max_err] becomes [-(1+max_err+e), -(1-min_err-e)] U [1-min_err-e, 1+max_err+e]
		maxInterval.Add(maxInterval, e)
		minInterval.Sub(minInterval, e)

		intervals = []bignum.Interval{
			{A: *new(big.Float).Neg(maxInterval), B: *new(big.Float).Neg(minInterval), Nodes: 1 + ((deg[i] + 1) >> 1)},
			{A: *minInterval, B: *maxInterval, Nodes: 1 + ((deg[i] + 1) >> 1)},
		}

		coeffs[i-1] = make([]*big.Float, deg[i-1]+1)
		for j := range coeffs[i-1] {
			coeffs[i-1][j] = new(big.Float).Set(r.Coeffs[j])
			coeffs[i-1][j].Quo(coeffs[i-1][j], maxInterval) // Interval normalization
		}

		params := bignum.RemezParameters{
			Function:        f,
			Basis:           bignum.Chebyshev,
			Intervals:       intervals,
			ScanStep:        scanStep,
			Prec:            prec,
			OptimalScanStep: true,
		}

		fmt.Printf("P[%d]\n", i)
		fmt.Printf("Interval: [%.*f, %.*f] U [%.*f, %.*f]\n", decimals, &intervals[0].A, decimals, &intervals[0].B, decimals, &intervals[1].A, decimals, &intervals[1].B)
		r = bignum.NewRemez(params)
		r.Approximate(maxIters, alpha)
		//r.ShowCoeffs(decimals)
		r.ShowError(decimals)
		fmt.Println()
	}

	// Since this is the last polynomial, we can skip the interval scaling.
	coeffs[len(deg)-1] = make([]*big.Float, deg[len(deg)-1]+1)
	for j := range coeffs[len(deg)-1] {
		coeffs[len(deg)-1][j] = new(big.Float).Set(r.Coeffs[j])
	}

	f64, _ := r.MaxErr.Float64()
	fmt.Printf("Output Precision: %f\n", math.Log2(f64))
	fmt.Println()

	return coeffs
}

// PrettyPrintCoefficients prints the coefficients formatted.
// If odd = true, even coefficients are zeroed.
// If even = true, odd coefficients are zeroed.
func PrettyPrintCoefficients(decimals int, coeffs []*big.Float, odd, even, first bool) {
	fmt.Printf("{")
	for i, c := range coeffs {
		if (i&1 == 1 && odd) || (i&1 == 0 && even) || (i == 0 && first) {
			fmt.Printf("\"%.*f\", ", decimals, c)
		} else {
			fmt.Printf("\"0\", ")
		}

	}
	fmt.Printf("},\n")
}

func parseCoeffs(coeffsStr []string) (coeffs []*big.Float) {

	var prec uint
	for _, c := range coeffsStr {
		prec = utils.Max(prec, uint(len(c)))
	}

	prec = uint(float64(prec)*3.3219280948873626 + 0.5) // max(float64, digits * log2(10))

	coeffs = make([]*big.Float, len(coeffsStr))
	for i := range coeffsStr {
		coeffs[i], _ = new(big.Float).SetPrec(prec).SetString(coeffsStr[i])
	}

	return
}
