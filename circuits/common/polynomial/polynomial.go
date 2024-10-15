// Package polynomial bundles generic parts of the homomorphic polynomial evaluation circuit.
package polynomial

import (
	"fmt"
	"math/bits"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/utils"
	"github.com/tuneinsight/lattigo/v6/utils/bignum"
)

// Polynomial is a struct for representing plaintext polynomials
// for their homomorphic evaluation in an encrypted point. The
// type wraps a [bignum.Polynomial] along with several evaluation-
// related parameters.
type Polynomial struct {
	bignum.Polynomial
	MaxDeg int        // Always set to len(Coeffs)-1
	Lead   bool       // Always set to true
	Lazy   bool       // Flag for lazy-relinearization
	Level  int        // Metadata for BSGS polynomial evaluation
	Scale  rlwe.Scale // Metadata for BSGS polynomial evaluation
}

// NewPolynomial returns an instantiated Polynomial for the
// provided [bignum.Polynomial].
func NewPolynomial(poly bignum.Polynomial) Polynomial {
	return Polynomial{
		Polynomial: poly,
		MaxDeg:     len(poly.Coeffs) - 1,
		Lead:       true,
		Lazy:       false,
	}
}

// Factorize factorizes p as X^{n} * pq + pr.
func (p Polynomial) Factorize(n int) (pq, pr Polynomial) {

	pq = Polynomial{}
	pr = Polynomial{}

	pq.Polynomial, pr.Polynomial = p.Polynomial.Factorize(n)

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

// PatersonStockmeyerPolynomial is a struct that stores
// the Paterson Stockmeyer decomposition of a polynomial.
// The decomposition of P(X) is given as sum pi(X) * X^{2^{n}}
// where degree(pi(X)) =~ sqrt(degree(P(X)))
type PatersonStockmeyerPolynomial struct {
	Degree int
	Base   int
	Level  int
	Scale  rlwe.Scale
	Value  []Polynomial
}

// PatersonStockmeyerPolynomial returns the Paterson Stockmeyer polynomial decomposition of the target polynomial.
// The decomposition is done with the power of two basis.
func (p Polynomial) PatersonStockmeyerPolynomial(params rlwe.ParameterProvider, inputLevel int, inputScale, outputScale rlwe.Scale, eval SimEvaluator) PatersonStockmeyerPolynomial {

	// ceil(log2(degree))
	logDegree := bits.Len64(uint64(p.Degree()))

	// optimal ratio between degree(pi(X)) et degree(P(X))
	logSplit := bignum.OptimalSplit(logDegree)

	// Initializes the simulated polynomial evaluation
	pb := SimPowerBasis{}
	pb[1] = &SimOperand{
		Level: inputLevel,
		Scale: inputScale,
	}

	// Generates the simulated powers (to get the scaling factors)
	pb.GenPower(params, 1<<logDegree, eval)
	for i := (1 << logSplit) - 1; i > 2; i-- {
		pb.GenPower(params, i, eval)
	}

	// Simulates the homomorphic evaluation with levels and scaling factors to retrieve the scaling factor of each pi(X).
	PSPoly, _ := recursePS(params, logSplit, inputLevel-eval.PolynomialDepth(p.Degree()), p, pb, outputScale, eval)

	return PatersonStockmeyerPolynomial{
		Degree: p.Degree(),
		Base:   1 << logSplit,
		Level:  inputLevel,
		Scale:  outputScale,
		Value:  PSPoly,
	}
}

// recursePS is a recursive implementation of a polynomial evaluation via the Paterson Stockmeyer algorithm with a power of two decomposition.
func recursePS(params rlwe.ParameterProvider, logSplit, targetLevel int, p Polynomial, pb SimPowerBasis, outputScale rlwe.Scale, eval SimEvaluator) ([]Polynomial, *SimOperand) {

	if p.Degree() < (1 << logSplit) {

		if p.Lead && logSplit > 1 && p.MaxDeg > (1<<bits.Len64(uint64(p.MaxDeg)))-(1<<(logSplit-1)) {

			logDegree := int(bits.Len64(uint64(p.Degree())))
			logSplit := bignum.OptimalSplit(logDegree)

			return recursePS(params, logSplit, targetLevel, p, pb, outputScale, eval)
		}

		p.Level, p.Scale = eval.UpdateLevelAndScaleBabyStep(p.Lead, targetLevel, outputScale)

		return []Polynomial{p}, &SimOperand{Level: p.Level, Scale: p.Scale}
	}

	var nextPower = 1 << logSplit
	for nextPower < (p.Degree()>>1)+1 {
		nextPower <<= 1
	}

	XPow := pb[nextPower]

	coeffsq, coeffsr := p.Factorize(nextPower)

	tLevelNew, tScaleNew := eval.UpdateLevelAndScaleGiantStep(p.Lead, targetLevel, outputScale, XPow.Scale)

	bsgsQ, res := recursePS(params, logSplit, tLevelNew, coeffsq, pb, tScaleNew, eval)

	eval.Rescale(res)
	res = eval.MulNew(res, XPow)

	bsgsR, tmp := recursePS(params, logSplit, targetLevel, coeffsr, pb, res.Scale, eval)

	// This checks that the underlying algorithm behaves as expected, which will always be
	// the case, unless the user provides an incorrect custom implementation.
	if !tmp.Scale.InDelta(res.Scale, float64(rlwe.ScalePrecision-12)) {
		panic(fmt.Errorf("recursePS: res.Scale != tmp.Scale: %v != %v", &res.Scale.Value, &tmp.Scale.Value))
	}

	return append(bsgsQ, bsgsR...), res
}

// PolynomialVector is a struct storing a set of polynomials and a mapping that
// indicates on which slot each polynomial has to be independently evaluated.
// For example, if we are given two polynomials P0(X) and P1(X) and the folling mapping: map[int][]int{0:[0, 1, 2], 1:[3, 4, 5]},
// then the polynomial evaluation on a vector [a, b, c, d, e, f, g, h] will evaluate to [P0(a), P0(b), P0(c), P1(d), P1(e), P1(f), 0, 0]
type PolynomialVector struct {
	Value   []Polynomial
	Mapping map[int][]int
}

// NewPolynomialVector instantiates a new [PolynomialVector] from a set of [bignum.Polynomial] and a mapping indicating
// which polynomial has to be evaluated on which slot.
// For example, if we are given two polynomials P0(X) and P1(X) and the following mapping: map[int][]int{0:[0, 1, 2], 1:[3, 4, 5]},
// then the polynomial evaluation on a vector [a, b, c, d, e, f, g, h] will evaluate to [P0(a), P0(b), P0(c), P1(d), P1(e), P1(f), 0, 0]
func NewPolynomialVector(polys []bignum.Polynomial, mapping map[int][]int) (PolynomialVector, error) {
	var maxDeg int
	var basis bignum.Basis
	for i := range polys {
		maxDeg = utils.Max(maxDeg, polys[i].Degree())
		basis = polys[i].Basis
	}

	for i := range polys {
		if basis != polys[i].Basis {
			return PolynomialVector{}, fmt.Errorf("polynomial basis must be the same for all polynomials in a polynomial vector")
		}

		if maxDeg != polys[i].Degree() {
			return PolynomialVector{}, fmt.Errorf("polynomial degree must all be the same")
		}
	}

	polyvec := make([]Polynomial, len(polys))

	for i := range polyvec {
		polyvec[i] = NewPolynomial(polys[i])
	}

	return PolynomialVector{
		Value:   polyvec,
		Mapping: mapping,
	}, nil
}

// IsEven returns true if all underlying polynomials are even,
// i.e. all odd powers are zero.
func (p PolynomialVector) IsEven() (even bool) {
	even = true
	for _, poly := range p.Value {
		even = even && poly.IsEven
	}
	return
}

// IsOdd returns true if all underlying polynomials are odd,
// i.e. all even powers are zero.
func (p PolynomialVector) IsOdd() (odd bool) {
	odd = true
	for _, poly := range p.Value {
		odd = odd && poly.IsOdd
	}
	return
}

// Factorize factorizes the underlying Polynomial vector p into p = polyq * X^{n} + polyr.
func (p PolynomialVector) Factorize(n int) (polyq, polyr PolynomialVector) {

	coeffsq := make([]Polynomial, len(p.Value))
	coeffsr := make([]Polynomial, len(p.Value))

	for i, p := range p.Value {
		coeffsq[i], coeffsr[i] = p.Factorize(n)
	}

	return PolynomialVector{Value: coeffsq, Mapping: p.Mapping}, PolynomialVector{Value: coeffsr, Mapping: p.Mapping}
}

// PatersonStockmeyerPolynomialVector is a struct implementing the
// Paterson Stockmeyer decomposition of a PolynomialVector.
// See [PatersonStockmeyerPolynomial] for additional information.
type PatersonStockmeyerPolynomialVector struct {
	Value   []PatersonStockmeyerPolynomial
	Mapping map[int][]int
}

// PatersonStockmeyerPolynomial returns the Paterson Stockmeyer polynomial decomposition of the target PolynomialVector.
// The decomposition is done with the power of two basis
func (p PolynomialVector) PatersonStockmeyerPolynomial(params rlwe.Parameters, inputLevel int, inputScale, outputScale rlwe.Scale, eval SimEvaluator) PatersonStockmeyerPolynomialVector {
	Value := make([]PatersonStockmeyerPolynomial, len(p.Value))
	for i := range Value {
		Value[i] = p.Value[i].PatersonStockmeyerPolynomial(params, inputLevel, inputScale, outputScale, eval)
	}

	return PatersonStockmeyerPolynomialVector{
		Value:   Value,
		Mapping: p.Mapping,
	}
}
