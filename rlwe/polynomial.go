package rlwe

import (
	"fmt"
	"math/big"
	"math/bits"

	"github.com/tuneinsight/lattigo/v4/utils"

	"github.com/tuneinsight/lattigo/v4/utils/bignum"
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

// GetPatersonStockmeyerPolynomial
func (p *Polynomial) GetPatersonStockmeyerPolynomial(params Parameters, inputLevel int, inputScale, outputScale Scale) *PatersonStockmeyerPolynomial {

	logDegree := bits.Len64(uint64(p.Degree()))
	logSplit := polynomial.OptimalSplit(logDegree)

	nbModuliPerRescale := params.DefaultScaleModuliRatio()

	targetLevel := inputLevel - nbModuliPerRescale*(logDegree-1)

	pb := DummyPowerBasis{}
	pb[1] = &DummyOperand{
		Level: inputLevel,
		Scale: inputScale,
	}

	pb.GenPower(params, 1<<logDegree, nbModuliPerRescale)
	for i := (1 << logSplit) - 1; i > 2; i-- {
		pb.GenPower(params, i, nbModuliPerRescale)
	}

	PSPoly, _ := recursePS(params, logSplit, targetLevel, nbModuliPerRescale, p, pb, outputScale)

	return &PatersonStockmeyerPolynomial{
		Degree: p.Degree(),
		Base:   1 << logSplit,
		Level:  inputLevel,
		Scale:  outputScale,
		Value:  PSPoly,
	}
}

func recursePS(params Parameters, logSplit, targetLevel, nbModuliPerRescale int, p *Polynomial, pb DummyPowerBasis, outputScale Scale) ([]*Polynomial, *DummyOperand) {

	if p.Degree() < (1 << logSplit) {

		if p.Lead && logSplit > 1 && p.MaxDeg > (1<<bits.Len64(uint64(p.MaxDeg)))-(1<<(logSplit-1)) {

			logDegree := int(bits.Len64(uint64(p.Degree())))
			logSplit := polynomial.OptimalSplit(logDegree)

			return recursePS(params, logSplit, targetLevel, nbModuliPerRescale, p, pb, outputScale)
		}

		if p.Lead {
			for i := 0; i < nbModuliPerRescale; i++ {
				outputScale = outputScale.Mul(NewScale(params.Q()[targetLevel-i]))
			}
		}

		p.Level = targetLevel
		p.Scale = outputScale

		return []*Polynomial{p}, &DummyOperand{Level: targetLevel, Scale: outputScale}
	}

	var nextPower = 1 << logSplit
	for nextPower < (p.Degree()>>1)+1 {
		nextPower <<= 1
	}

	XPow := pb[nextPower]

	coeffsq, coeffsr := p.Factorize(nextPower)

	var qi *big.Int
	if p.Lead {
		qi = bignum.NewInt(params.Q()[targetLevel])
		for i := 1; i < nbModuliPerRescale; i++ {
			qi.Mul(qi, bignum.NewInt(params.Q()[targetLevel-i]))
		}
	} else {
		qi = bignum.NewInt(params.Q()[targetLevel+nbModuliPerRescale])
		for i := 1; i < nbModuliPerRescale; i++ {
			qi.Mul(qi, bignum.NewInt(params.Q()[targetLevel+nbModuliPerRescale-i]))
		}
	}

	tScaleNew := outputScale.Mul(NewScale(qi))
	tScaleNew = tScaleNew.Div(XPow.Scale)

	bsgsQ, res := recursePS(params, logSplit, targetLevel+nbModuliPerRescale, nbModuliPerRescale, coeffsq, pb, tScaleNew)

	res.Rescale(params, nbModuliPerRescale)
	res.Mul(res, XPow)

	bsgsR, tmp := recursePS(params, logSplit, targetLevel, nbModuliPerRescale, coeffsr, pb, res.Scale)

	if !tmp.Scale.InDelta(res.Scale, float64(ScalePrecision-12)) {
		panic(fmt.Errorf("recursePS: res.Scale != tmp.Scale: %v != %v", &res.Scale.Value, &tmp.Scale.Value))
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
func (p *PolynomialVector) GetPatersonStockmeyerPolynomial(params Parameters, inputLevel int, inputScale, outputScale Scale) *PatersonStockmeyerPolynomialVector {
	Value := make([]*PatersonStockmeyerPolynomial, len(p.Value))
	for i := range Value {
		Value[i] = p.Value[i].GetPatersonStockmeyerPolynomial(params, inputLevel, inputScale, outputScale)
	}

	return &PatersonStockmeyerPolynomialVector{
		Value:      Value,
		SlotsIndex: p.SlotsIndex,
	}
}
