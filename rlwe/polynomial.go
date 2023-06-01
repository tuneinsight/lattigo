package rlwe

import (
	"fmt"
	"math/bits"

	"github.com/tuneinsight/lattigo/v4/utils"
	"github.com/tuneinsight/lattigo/v4/utils/bignum/polynomial"
)

type Polynomial struct {
	*polynomial.Polynomial
	MaxDeg int   // Always set to len(Coeffs)-1
	Lead   bool  // Always set to true
	Lazy   bool  // Flag for lazy-relinearization
	Level  int   // Metadata for BSGS polynomial evaluation
	Scale  Scale // Metatata for BSGS polynomial evaluation
}

func NewPolynomial(poly *polynomial.Polynomial) *Polynomial {
	return &Polynomial{
		Polynomial: poly,
		MaxDeg:     len(poly.Coeffs) - 1,
		Lead:       true,
		Lazy:       false,
	}
}

func (p *Polynomial) Factorize(n int) (pq, pr *Polynomial) {

	pq = &Polynomial{}
	pr = &Polynomial{}

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
	Scale  Scale
	Value  []*Polynomial
}

func (p *Polynomial) GetPatersonStockmeyerPolynomial(params ParametersInterface, inputLevel int, inputScale, outputScale Scale, eval DummyEvaluator) *PatersonStockmeyerPolynomial {

	logDegree := bits.Len64(uint64(p.Degree()))
	logSplit := polynomial.OptimalSplit(logDegree)

	pb := DummyPowerBasis{}
	pb[1] = &DummyOperand{
		Level:          inputLevel,
		PlaintextScale: inputScale,
	}

	pb.GenPower(params, 1<<logDegree, eval)
	for i := (1 << logSplit) - 1; i > 2; i-- {
		pb.GenPower(params, i, eval)
	}

	PSPoly, _ := recursePS(params, logSplit, inputLevel-eval.PolynomialDepth(p.Degree()), p, pb, outputScale, eval)

	return &PatersonStockmeyerPolynomial{
		Degree: p.Degree(),
		Base:   1 << logSplit,
		Level:  inputLevel,
		Scale:  outputScale,
		Value:  PSPoly,
	}
}

func recursePS(params ParametersInterface, logSplit, targetLevel int, p *Polynomial, pb DummyPowerBasis, outputScale Scale, eval DummyEvaluator) ([]*Polynomial, *DummyOperand) {

	if p.Degree() < (1 << logSplit) {

		if p.Lead && logSplit > 1 && p.MaxDeg > (1<<bits.Len64(uint64(p.MaxDeg)))-(1<<(logSplit-1)) {

			logDegree := int(bits.Len64(uint64(p.Degree())))
			logSplit := polynomial.OptimalSplit(logDegree)

			return recursePS(params, logSplit, targetLevel, p, pb, outputScale, eval)
		}

		p.Level, p.Scale = eval.UpdateLevelAndScaleBabyStep(p.Lead, targetLevel, outputScale)

		return []*Polynomial{p}, &DummyOperand{Level: p.Level, PlaintextScale: p.Scale}
	}

	var nextPower = 1 << logSplit
	for nextPower < (p.Degree()>>1)+1 {
		nextPower <<= 1
	}

	XPow := pb[nextPower]

	coeffsq, coeffsr := p.Factorize(nextPower)

	tLevelNew, tScaleNew := eval.UpdateLevelAndScaleGiantStep(p.Lead, targetLevel, outputScale, XPow.PlaintextScale)

	bsgsQ, res := recursePS(params, logSplit, tLevelNew, coeffsq, pb, tScaleNew, eval)

	eval.Rescale(res)
	res = eval.MulNew(res, XPow)

	bsgsR, tmp := recursePS(params, logSplit, targetLevel, coeffsr, pb, res.PlaintextScale, eval)

	if !tmp.PlaintextScale.InDelta(res.PlaintextScale, float64(ScalePrecision-12)) {
		panic(fmt.Errorf("recursePS: res.PlaintextScale != tmp.PlaintextScale: %v != %v", &res.PlaintextScale.Value, &tmp.PlaintextScale.Value))
	}

	return append(bsgsQ, bsgsR...), res
}

type PolynomialVector struct {
	Value      []*Polynomial
	SlotsIndex map[int][]int
}

func NewPolynomialVector(polys []*Polynomial, slotsIndex map[int][]int) *PolynomialVector {
	var maxDeg int
	var basis polynomial.Basis
	for i := range polys {
		maxDeg = utils.Max(maxDeg, polys[i].Degree())
		basis = polys[i].Basis
	}

	for i := range polys {
		if basis != polys[i].Basis {
			panic(fmt.Errorf("polynomial basis must be the same for all polynomials in a polynomial vector"))
		}

		if maxDeg != polys[i].Degree() {
			panic(fmt.Errorf("polynomial degree must all be the same"))
		}
	}

	polyvec := make([]*Polynomial, len(polys))

	copy(polyvec, polys)

	return &PolynomialVector{
		Value:      polyvec,
		SlotsIndex: slotsIndex,
	}
}

func (p *PolynomialVector) IsEven() (even bool) {
	even = true
	for _, poly := range p.Value {
		even = even && poly.IsEven
	}
	return
}

func (p *PolynomialVector) IsOdd() (odd bool) {
	odd = true
	for _, poly := range p.Value {
		odd = odd && poly.IsOdd
	}
	return
}

func (p *PolynomialVector) Factorize(n int) (polyq, polyr *PolynomialVector) {

	coeffsq := make([]*Polynomial, len(p.Value))
	coeffsr := make([]*Polynomial, len(p.Value))

	for i, p := range p.Value {
		coeffsq[i], coeffsr[i] = p.Factorize(n)
	}

	return &PolynomialVector{Value: coeffsq, SlotsIndex: p.SlotsIndex}, &PolynomialVector{Value: coeffsr, SlotsIndex: p.SlotsIndex}
}

type PatersonStockmeyerPolynomialVector struct {
	Value      []*PatersonStockmeyerPolynomial
	SlotsIndex map[int][]int
}

// GetPatersonStockmeyerPolynomial returns
func (p *PolynomialVector) GetPatersonStockmeyerPolynomial(params ParametersInterface, inputLevel int, inputScale, outputScale Scale, eval DummyEvaluator) *PatersonStockmeyerPolynomialVector {
	Value := make([]*PatersonStockmeyerPolynomial, len(p.Value))
	for i := range Value {
		Value[i] = p.Value[i].GetPatersonStockmeyerPolynomial(params, inputLevel, inputScale, outputScale, eval)
	}

	return &PatersonStockmeyerPolynomialVector{
		Value:      Value,
		SlotsIndex: p.SlotsIndex,
	}
}
