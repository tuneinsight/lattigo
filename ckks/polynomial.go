package ckks

import (
	"fmt"
	"math"
	"math/bits"

	"github.com/tuneinsight/lattigo/v3/ring"
	"github.com/tuneinsight/lattigo/v3/utils"
)

// BasisType is a type for the polynomials basis
type BasisType int

const (
	// Monomial : x^(a+b) = x^a * x^b
	Monomial = BasisType(0)
	// Chebyshev : T_(a+b) = 2 * T_a * T_b - T_(|a-b|)
	Chebyshev = BasisType(1)
)

// IsNegligbleThreshold : threshold under which a coefficient
// of a polynomial is ignored.
const IsNegligbleThreshold float64 = 1e-14

// Polynomial is a struct storing the coefficients of a polynomial
// that then can be evaluated on the ciphertext
type Polynomial struct {
	coefficients
}

func NewPolynomial(basisType BasisType, coeffs interface{}, slotsIndex map[int][]int) (poly Polynomial, err error) {

	var coeffsInterface coefficients
	var odd, even bool = true, true

	switch coeffs := coeffs.(type) {
	case []complex128:

		c := make([]complex128, len(coeffs))
		copy(c, coeffs)

		tmp0, tmp1 := isOddOrEvenPolynomial(c)
		odd, even = odd && tmp0, even && tmp1

		coeffsInterface = &coefficientsComplex128{
			coeffs:     [][]complex128{c},
			slotsIndex: slotsIndex,
			maxDeg:     len(coeffs) - 1,
			lead:       true,
			odd:        odd,
			even:       even,
			basis:      basisType,
		}

	case []float64:

		c := make([]complex128, len(coeffs))
		for i := range coeffs {
			c[i] = complex(coeffs[i], 0)
		}

		tmp0, tmp1 := isOddOrEvenPolynomial(c)
		odd, even = odd && tmp0, even && tmp1

		coeffsInterface = &coefficientsComplex128{
			coeffs:     [][]complex128{c},
			slotsIndex: slotsIndex,
			maxDeg:     len(coeffs) - 1,
			lead:       true,
			odd:        odd,
			even:       even,
			basis:      basisType,
		}

	case [][]complex128:
		{

			if slotsIndex == nil {
				return poly, fmt.Errorf("NewPolynomial: cannot have slotsIndex == nil if coeffs == [][]complex128")
			}

			var maxDeg int
			for i := range coeffs {
				maxDeg = utils.MaxInt(maxDeg, len(coeffs[i]))
			}

			c := make([][]complex128, len(coeffs))
			for i := range c {
				c[i] = make([]complex128, maxDeg)
				copy(c[i], coeffs[i])

				tmp0, tmp1 := isOddOrEvenPolynomial(c[i])
				odd, even = odd && tmp0, even && tmp1
			}

			coeffsInterface = &coefficientsComplex128{
				coeffs:     c,
				slotsIndex: slotsIndex,
				maxDeg:     maxDeg - 1,
				lead:       true,
				odd:        odd,
				even:       even,
				basis:      basisType,
			}
		}

	case [][]float64:
		{
			if slotsIndex == nil {
				return poly, fmt.Errorf("NewPolynomial: cannot have slotsIndex == nil if coeffs == [][]float64")
			}

			var maxDeg int
			for i := range coeffs {
				maxDeg = utils.MaxInt(maxDeg, len(coeffs[i]))
			}

			c := make([][]complex128, len(coeffs))
			for i := range c {

				c[i] = make([]complex128, maxDeg)

				for j := range coeffs[i] {
					c[i][j] = complex(coeffs[i][j], 0)
				}

				tmp0, tmp1 := isOddOrEvenPolynomial(c[i])
				odd, even = odd && tmp0, even && tmp1
			}

			coeffsInterface = &coefficientsComplex128{
				coeffs:     c,
				slotsIndex: slotsIndex,
				maxDeg:     maxDeg - 1,
				lead:       true,
				odd:        odd,
				even:       even,
				basis:      basisType,
			}
		}

	default:
		return poly, fmt.Errorf("NewPolynomial: invalid coeffs.(type)")
	}

	return Polynomial{
		coefficients: coeffsInterface,
	}, nil
}

type coefficients interface {
	Basis() BasisType
	Depth() int
	Degree() int
}

type coefficientsComplex128 struct {
	basis      BasisType
	coeffs     [][]complex128
	slotsIndex map[int][]int
	odd, even  bool
	lead       bool
	maxDeg     int
}

func (c *coefficientsComplex128) Basis() BasisType {
	return c.basis
}

func (c *coefficientsComplex128) Depth() int {
	return bits.Len64(uint64(c.Degree()))
}

func (c *coefficientsComplex128) Degree() int {
	return len(c.coeffs[0]) - 1
}

func (c *coefficientsComplex128) splitBSGS(split int) (polyq, polyr coefficientsComplex128) {

	polyq = coefficientsComplex128{}
	polyr = coefficientsComplex128{}

	polyq.coeffs = make([][]complex128, len(c.coeffs))
	polyr.coeffs = make([][]complex128, len(c.coeffs))

	for i, p := range c.coeffs {
		polyq.coeffs[i], polyr.coeffs[i] = splitCoeffsBSGS(p, split, c.basis)
	}

	polyq.lead = c.lead
	polyq.maxDeg = c.maxDeg

	if c.maxDeg == c.Degree() {
		polyr.maxDeg = split - 1
	} else {
		polyr.maxDeg = c.maxDeg - (c.Degree() - split + 1)
	}

	polyq.basis = c.basis
	polyr.basis = c.basis

	polyq.slotsIndex = c.slotsIndex
	polyr.slotsIndex = c.slotsIndex

	return
}

func splitCoeffsBSGS(coeffs []complex128, split int, basis BasisType) (coeffsq, coeffsr []complex128) {

	// Splits a polynomial p such that p = q*C^degree + r.
	coeffsr = make([]complex128, split)
	coeffsq = make([]complex128, len(coeffs)-split)

	coeffsq[0] = coeffs[split]
	for i := 0; i < split; i++ {
		coeffsr[i] = coeffs[i]
	}

	if basis == Monomial {
		for i := split + 1; i < len(coeffs); i++ {
			coeffsq[i-split] = coeffs[i]
		}
	} else if basis == Chebyshev {
		for i, j := split+1, 1; i < len(coeffs); i, j = i+1, j+1 {
			coeffsq[i-split] = 2 * coeffs[i]
			coeffsr[split-j] -= coeffs[i]
		}
	}

	return
}

type coefficientsBSGSComplex128 struct {
	basis      BasisType
	coeffs     [][][]complex128
	slotsIndex map[int][]int
}

func (c *coefficientsBSGSComplex128) Basis() BasisType {
	return c.basis
}

func (c *coefficientsBSGSComplex128) Depth() int {
	return bits.Len64(uint64(c.Degree()))
}

func (c *coefficientsBSGSComplex128) Degree() int {
	var deg int
	for _, ci := range c.coeffs {
		deg += len(ci)
	}
	return deg - 1
}

type coefficientsBSGSPlaintext struct {
	basis  BasisType
	coeffs [][]*Plaintext
}

func (c *coefficientsBSGSPlaintext) Basis() BasisType {
	return c.basis
}

func (c *coefficientsBSGSPlaintext) Depth() int {
	return bits.Len64(uint64(c.Degree()))
}

func (c *coefficientsBSGSPlaintext) Degree() int {
	var deg int
	for _, ci := range c.coeffs {
		deg += len(ci)
	}
	return deg - 1
}

type coefficientsBSGSCiphertext struct {
	basis  BasisType
	coeffs [][]*Ciphertext
}

func (c *coefficientsBSGSCiphertext) Basis() BasisType {
	return c.basis
}

func (c *coefficientsBSGSCiphertext) Depth() int {
	return bits.Len64(uint64(c.Degree()))
}

func (c *coefficientsBSGSCiphertext) Degree() int {
	var deg int
	for _, ci := range c.coeffs {
		deg += len(ci)
	}
	return deg - 1
}

func (p *Polynomial) Encode(level int, ecd Encoder, inputScale, outputScale float64) (err error) {

	params := ecd.(*encoderComplex128).params

	var coeffs [][][]complex128
	switch coeffInterface := p.coefficients.(type) {
	case *coefficientsComplex128:
		coeffs = getScaledBSGSCoefficients(params, level, inputScale, coeffInterface, outputScale)
	default:
		return fmt.Errorf("Polynomial.Encode(*): underlying polynomial coefficient must be *coefficientsComplex128")
	}

	_ = coeffs

	return
}

func getScaledBSGSCoefficients(params Parameters, level int, scale float64, pol *coefficientsComplex128, targetScale float64) (coeffs [][][]complex128) {

	dummbpb := newDummyPolynomialBasis(params, &dummyCiphertext{level, scale})

	logDegree := bits.Len64(uint64(pol.Degree()))
	logSplit := optimalSplit(logDegree)

	isRingStandard := params.RingType() == ring.Standard

	odd, even := pol.odd, pol.even

	for i := (1 << logSplit) - 1; i > 1; i-- {
		if !(even || odd) || (i&1 == 0 && even) || (i&1 == 1 && odd) {
			dummbpb.GenPower(i, isRingStandard, targetScale)
		}
	}

	for i := logSplit; i < logDegree; i++ {
		dummbpb.GenPower(1<<i, false, targetScale)
	}

	polyEval := &dummyPolynomialEvaluator{}
	polyEval.params = params
	polyEval.dummyPolynomialBasis = *dummbpb
	polyEval.logDegree = logDegree
	polyEval.logSplit = logSplit
	polyEval.isOdd = odd
	polyEval.isEven = even
	polyEval.bsgsCoeffs = make([][][]complex128, len(pol.coeffs))

	polyEval.recurse(dummbpb.Value[1].Level-logDegree+1, targetScale, *pol)

	return polyEval.bsgsCoeffs
}

func (polyEval *dummyPolynomialEvaluator) recurse(targetLevel int, targetScale float64, pol coefficientsComplex128) (res *dummyCiphertext) {

	params := polyEval.params

	logSplit := polyEval.logSplit

	// Recursively computes the evaluation of the Chebyshev polynomial using a baby-set giant-step algorithm.
	if pol.Degree() < (1 << logSplit) {

		if pol.lead && polyEval.logSplit > 1 && pol.maxDeg%(1<<(logSplit+1)) > (1<<(logSplit-1)) {

			logDegree := int(bits.Len64(uint64(pol.Degree())))
			logSplit := logDegree >> 1

			polyEvalBis := new(dummyPolynomialEvaluator)
			polyEvalBis.params = params
			polyEvalBis.logDegree = logDegree
			polyEvalBis.logSplit = logSplit
			polyEvalBis.dummyPolynomialBasis = polyEval.dummyPolynomialBasis
			polyEvalBis.isOdd = polyEval.isOdd
			polyEvalBis.isEven = polyEval.isEven
			polyEvalBis.bsgsCoeffs = polyEval.bsgsCoeffs

			return polyEvalBis.recurse(targetLevel, targetScale, pol)
		}

		if pol.lead {
			targetScale *= params.QiFloat64(targetLevel)
		}

		return polyEval.evaluatePolyFromPolynomialBasis(targetScale, targetLevel, pol)
	}

	var nextPower = 1 << polyEval.logSplit
	for nextPower < (pol.Degree()>>1)+1 {
		nextPower <<= 1
	}

	coeffsq, coeffsr := pol.splitBSGS(nextPower)

	XPow := polyEval.dummyPolynomialBasis.Value[nextPower]

	level := targetLevel

	var currentQi float64
	if pol.lead {
		currentQi = params.QiFloat64(level)
	} else {
		currentQi = params.QiFloat64(level + 1)
	}

	res = polyEval.recurse(targetLevel+1, targetScale*currentQi/XPow.Scale, coeffsq)

	res.rescale(params, params.DefaultScale())

	res.Scale *= XPow.Scale

	tmp := polyEval.recurse(res.Level, res.Scale, coeffsr)

	res.Level = utils.MinInt(tmp.Level, res.Level)

	return
}

func (polyEval *dummyPolynomialEvaluator) evaluatePolyFromPolynomialBasis(targetScale float64, level int, pol coefficientsComplex128) (res *dummyCiphertext) {

	X := polyEval.dummyPolynomialBasis.Value

	values := make([][]complex128, pol.Degree()+1)
	for i := range values {
		values[i] = make([]complex128, len(pol.coeffs))
	}

	for i, c := range pol.coeffs {
		if isNotNegligible(c[0]) {
			values[0][i] = c[0] * complex(targetScale, 0)
		}
	}

	for i := 1; i < pol.Degree(); i++ {
		for j, c := range pol.coeffs {
			if isNotNegligible(c[i]) {
				values[i][j] = c[i] * complex(targetScale/X[i].Scale, 0)
			}
		}
	}

	polyEval.bsgsCoeffs = append([][][]complex128{values}, polyEval.bsgsCoeffs...)

	return &dummyCiphertext{level, targetScale}
}

type dummyCiphertext struct {
	Level int
	Scale float64
}

func (d *dummyCiphertext) rescale(params Parameters, minScale float64) {
	var nbRescales int
	for d.Level-nbRescales >= 0 && d.Scale/float64(params.Q()[d.Level-nbRescales]) >= minScale/2 {
		d.Scale /= (float64(params.Q()[d.Level-nbRescales]))
		nbRescales++
	}

	d.Level -= nbRescales
}

type dummyPolynomialBasis struct {
	Value  map[int]*dummyCiphertext
	params Parameters
}

func newDummyPolynomialBasis(params Parameters, ct *dummyCiphertext) (p *dummyPolynomialBasis) {
	p = new(dummyPolynomialBasis)
	p.params = params
	p.Value = make(map[int]*dummyCiphertext)
	p.Value[1] = ct
	return
}

func (p *dummyPolynomialBasis) GenPower(n int, lazy bool, scale float64) {
	if p.Value[n] == nil {
		p.genPower(n, lazy, scale)
		p.Value[n].rescale(p.params, scale)
	}
}

func (p *dummyPolynomialBasis) genPower(n int, lazy bool, scale float64) (err error) {
	if p.Value[n] == nil {

		isPow2 := n&(n-1) == 0

		// Computes the index required to compute the asked ring evaluation
		var a, b int
		if isPow2 {
			a, b = n/2, n/2 //Necessary for optimal depth
		} else {
			// [Lee et al. 2020] : High-Precision and Low-Complexity Approximate Homomorphic Encryption by Error Variance Minimization
			// Maximize the number of odd terms of Chebyshev basis
			k := int(math.Ceil(math.Log2(float64(n)))) - 1
			a = (1 << k) - 1
			b = n + 1 - (1 << k)
		}

		// Recurses on the given indexes
		p.genPower(a, lazy, scale)
		p.genPower(b, lazy, scale)

		p.Value[a].rescale(p.params, scale)
		p.Value[b].rescale(p.params, scale)

		// Computes C[n] = C[a]*C[b]
		p.Value[n] = new(dummyCiphertext)
		p.Value[n].Level = utils.MinInt(p.Value[a].Level, p.Value[b].Level)
		if lazy && !isPow2 {
			p.Value[n].Scale = p.Value[a].Scale * p.Value[b].Scale
		} else {
			p.Value[n].Scale = p.Value[a].Scale * p.Value[b].Scale
			p.Value[n].rescale(p.params, scale)
		}
	}
	return
}

type dummyPolynomialEvaluator struct {
	dummyPolynomialBasis
	params     Parameters
	logDegree  int
	logSplit   int
	isOdd      bool
	isEven     bool
	bsgsCoeffs [][][]complex128
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
