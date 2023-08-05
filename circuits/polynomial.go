package circuits

import (
	"fmt"
	"math/bits"

	"github.com/tuneinsight/lattigo/v4/rlwe"
	"github.com/tuneinsight/lattigo/v4/utils"
	"github.com/tuneinsight/lattigo/v4/utils/bignum"
)

// Polynomial is a struct for representing plaintext polynomials
// for their homomorphic evaluation in an encrypted point. The
// type wraps a bignum.Polynomial along with several evaluation-
// related parameters.
type Polynomial struct {
	bignum.Polynomial
	MaxDeg int        // Always set to len(Coeffs)-1
	Lead   bool       // Always set to true
	Lazy   bool       // Flag for lazy-relinearization
	Level  int        // Metadata for BSGS polynomial evaluation
	Scale  rlwe.Scale // Metatata for BSGS polynomial evaluation
}

// NewPolynomial returns an instantiated Polynomial for the
// provided bignum.Polynomial.
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

type PatersonStockmeyerPolynomial struct {
	Degree int
	Base   int
	Level  int
	Scale  rlwe.Scale
	Value  []Polynomial
}

func (p Polynomial) GetPatersonStockmeyerPolynomial(params rlwe.ParameterProvider, inputLevel int, inputScale, outputScale rlwe.Scale, eval SimEvaluator) PatersonStockmeyerPolynomial {

	logDegree := bits.Len64(uint64(p.Degree()))
	logSplit := bignum.OptimalSplit(logDegree)

	pb := SimPowerBasis{}
	pb[1] = &SimOperand{
		Level: inputLevel,
		Scale: inputScale,
	}

	pb.GenPower(params, 1<<logDegree, eval)
	for i := (1 << logSplit) - 1; i > 2; i-- {
		pb.GenPower(params, i, eval)
	}

	PSPoly, _ := recursePS(params, logSplit, inputLevel-eval.PolynomialDepth(p.Degree()), p, pb, outputScale, eval)

	return PatersonStockmeyerPolynomial{
		Degree: p.Degree(),
		Base:   1 << logSplit,
		Level:  inputLevel,
		Scale:  outputScale,
		Value:  PSPoly,
	}
}

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

	if !tmp.Scale.InDelta(res.Scale, float64(rlwe.ScalePrecision-12)) {
		panic(fmt.Errorf("recursePS: res.Scale != tmp.Scale: %v != %v", &res.Scale.Value, &tmp.Scale.Value))
	}

	return append(bsgsQ, bsgsR...), res
}

type PolynomialVector struct {
	Value      []Polynomial
	SlotsIndex map[int][]int
}

func NewPolynomialVector(polys []Polynomial, slotsIndex map[int][]int) (PolynomialVector, error) {
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

	copy(polyvec, polys)

	return PolynomialVector{
		Value:      polyvec,
		SlotsIndex: slotsIndex,
	}, nil
}

func (p PolynomialVector) IsEven() (even bool) {
	even = true
	for _, poly := range p.Value {
		even = even && poly.IsEven
	}
	return
}

func (p PolynomialVector) IsOdd() (odd bool) {
	odd = true
	for _, poly := range p.Value {
		odd = odd && poly.IsOdd
	}
	return
}

func (p PolynomialVector) Factorize(n int) (polyq, polyr PolynomialVector) {

	coeffsq := make([]Polynomial, len(p.Value))
	coeffsr := make([]Polynomial, len(p.Value))

	for i, p := range p.Value {
		coeffsq[i], coeffsr[i] = p.Factorize(n)
	}

	return PolynomialVector{Value: coeffsq, SlotsIndex: p.SlotsIndex}, PolynomialVector{Value: coeffsr, SlotsIndex: p.SlotsIndex}
}

type PatersonStockmeyerPolynomialVector struct {
	Value      []PatersonStockmeyerPolynomial
	SlotsIndex map[int][]int
}

// GetPatersonStockmeyerPolynomial returns
func (p PolynomialVector) GetPatersonStockmeyerPolynomial(params rlwe.Parameters, inputLevel int, inputScale, outputScale rlwe.Scale, eval SimEvaluator) PatersonStockmeyerPolynomialVector {
	Value := make([]PatersonStockmeyerPolynomial, len(p.Value))
	for i := range Value {
		Value[i] = p.Value[i].GetPatersonStockmeyerPolynomial(params, inputLevel, inputScale, outputScale, eval)
	}

	return PatersonStockmeyerPolynomialVector{
		Value:      Value,
		SlotsIndex: p.SlotsIndex,
	}
}
