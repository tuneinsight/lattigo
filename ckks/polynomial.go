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

// NewPolynomial creates a new polynomial from:
//
// basisType: ckks.Monomial or ckks.Chebyshev
//
// coeffs:
//		- singe polynomial    : [coeffs]complex128 or [coeffs]float64
//      - multiple polynomials: [poly][coeffs]complex128 or [poly][coeffs]float64
//
// slotsIndex: in the case of multiple polynomials a map[poly]slot must be provided so that the evaluator will
//             know which polynomial to evaluate on each slot (which can also be no polynomial). A slotIndex map
//             can also be provided in the context of single polynomial, to only evaluate the polynomial on
//             specific slots (instead of all), resulting in zero values in slots where the polynomial isn't evaluated.
//
// The struct Polynomial can then be given to the ckks.Evaluator to evaluate the polynomial on a *ckks.Ciphertext.
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

func (p *Polynomial) Encode(ecd Encoder, level int, inputScale, outputScale float64) (ptPoly *Polynomial, err error) {

	params := ecd.(*encoderComplex128).params

	ptCoeffs := &coefficientsBSGSPlaintext{}

	switch coeffInterface := p.coefficients.(type) {
	case *coefficientsComplex128:

		ptCoeffs.basis = coeffInterface.basis
		ptCoeffs.odd = coeffInterface.odd
		ptCoeffs.even = coeffInterface.even

		getScaledBSGSCoefficients(params, ecd, nil, level, inputScale, coeffInterface, outputScale, ptCoeffs)

	default:
		return nil, fmt.Errorf("Polynomial.Encode(*): underlying polynomial is already encoded or encrypted")
	}

	return &Polynomial{ptCoeffs}, nil
}

func (p *Polynomial) Encrypt(ecd Encoder, enc Encryptor, level int, inputScale, outputScale float64) (ctPoly *Polynomial, err error) {

	params := ecd.(*encoderComplex128).params

	ctCoeffs := &coefficientsBSGSCiphertext{}

	switch coeffInterface := p.coefficients.(type) {
	case *coefficientsComplex128:

		ctCoeffs.basis = coeffInterface.basis
		ctCoeffs.odd = coeffInterface.odd
		ctCoeffs.even = coeffInterface.even

		getScaledBSGSCoefficients(params, ecd, enc, level, inputScale, coeffInterface, outputScale, ctCoeffs)

	case *coefficientsBSGSPlaintext:

		ctCoeffs.basis = coeffInterface.basis
		ctCoeffs.odd = coeffInterface.odd
		ctCoeffs.even = coeffInterface.even

		pt := coeffInterface.coeffs
		ct := make([][]*Ciphertext, len(pt))
		for i := range ct {
			ct[i] = make([]*Ciphertext, len(pt[i]))

			for j := range ct[i] {
				if pt[i][j] != nil {
					ct[i][j] = enc.EncryptNew(pt[i][j])
				}
			}
		}

	default:
		return nil, fmt.Errorf("Polynomial.Encrypt(*): underlying polynomial is already encrypted")
	}

	return &Polynomial{ctCoeffs}, nil
}

// coefficients is an interface to manage different types
// of coefficient representation.
//
//There four types currently:
// coefficientsComplex128
// coefficientsBSGSComplex128
// coefficientsBSGSPlaintext
// coefficientsBSGSCiphertext
type coefficients interface {
	Basis() BasisType
	Depth() (depth int)
	Degree() (degree int)
	OddEven() (odd, even bool)
	BSGSSplit() (giant, baby int)
}

func optimalSplit(giant int) (baby int) {
	baby = giant >> 1
	a := (1 << baby) + (1 << (giant - baby)) + giant - baby - 3
	b := (1 << (baby + 1)) + (1 << (giant - baby - 1)) + giant - baby - 4
	if a > b {
		baby++
	}

	return
}

// coefficientsComplex128: regular coefficients in complex128
// can store multiple polynomials: [#poly][coefficients]
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

func (c *coefficientsComplex128) Depth() (depth int) {
	return bits.Len64(uint64(c.Degree()))
}

func (c *coefficientsComplex128) Degree() (degree int) {
	return len(c.coeffs[0]) - 1
}

func (c *coefficientsComplex128) OddEven() (odd, even bool) {
	return c.odd, c.even
}

func (c *coefficientsComplex128) BSGSSplit() (giant, baby int) {
	return c.Depth(), optimalSplit(c.Depth())
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

// coefficientsBSGSComplex128: coefficients in baby-step giant-step format
// can store multiple polynomials: [giant-step][baby-step][#poly]
type coefficientsBSGSComplex128 struct {
	basis      BasisType
	coeffs     [][][]complex128
	slotsIndex map[int][]int
	odd, even  bool
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

func (c *coefficientsBSGSComplex128) OddEven() (odd, even bool) {
	return c.odd, c.even
}

func (c *coefficientsBSGSComplex128) BSGSSplit() (giant, baby int) {
	return len(c.coeffs[0]), len(c.coeffs[0][0])
}

// coefficientsBSGSPlaintext: coefficients in baby-step giant-step format
// and encoded in plaintext.
// can store multiple polynomials: [giant-step][baby-step]*Plaintext{#poly}
type coefficientsBSGSPlaintext struct {
	basis     BasisType
	coeffs    [][]*Plaintext
	odd, even bool
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

func (c *coefficientsBSGSPlaintext) OddEven() (odd, even bool) {
	return c.odd, c.even
}

func (c *coefficientsBSGSPlaintext) BSGSSplit() (giant, baby int) {
	return len(c.coeffs), len(c.coeffs[0])
}

// coefficientsBSGSPlaintext: coefficients in baby-step giant-step format
// and encrypted.
// can store multiple polynomials: [giant-step][baby-step]*Ciphertext{#poly}
type coefficientsBSGSCiphertext struct {
	basis     BasisType
	coeffs    [][]*Ciphertext
	odd, even bool
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

func (c *coefficientsBSGSCiphertext) OddEven() (odd, even bool) {
	return c.odd, c.even
}

func (c *coefficientsBSGSCiphertext) BSGSSplit() (giant, baby int) {
	return len(c.coeffs), len(c.coeffs[0])
}

type dummyPolynomialEvaluator struct {
	dummyPolynomialBasis
	params          Parameters
	ecd             Encoder
	enc             Encryptor
	giant           int
	baby            int
	coeffsInterface coefficients
}

func getScaledBSGSCoefficients(params Parameters, ecd Encoder, enc Encryptor, level int, scale float64, polIn *coefficientsComplex128, targetScale float64, polOut coefficients) {

	dummbpb := newDummyPolynomialBasis(params, &dummyCiphertext{level, scale})

	giant, baby := polIn.BSGSSplit()

	isRingStandard := params.RingType() == ring.Standard

	odd, even := polIn.OddEven()

	for i := (1 << baby) - 1; i > 1; i-- {
		if !(even || odd) || (i&1 == 0 && even) || (i&1 == 1 && odd) {
			dummbpb.GenPower(i, isRingStandard, targetScale)
		}
	}

	for i := baby; i < giant; i++ {
		dummbpb.GenPower(1<<i, false, targetScale)
	}

	polyEval := &dummyPolynomialEvaluator{}
	polyEval.params = params
	polyEval.ecd = ecd
	polyEval.enc = enc
	polyEval.dummyPolynomialBasis = *dummbpb
	polyEval.giant = giant
	polyEval.baby = baby

	switch c := polOut.(type) {
	case *coefficientsBSGSComplex128:

		c.basis = polIn.basis
		c.even, c.odd = polIn.OddEven()
		c.slotsIndex = polIn.slotsIndex

		polyEval.coeffsInterface = c

	case *coefficientsBSGSPlaintext:

		c.basis = polIn.basis
		c.even, c.odd = polIn.OddEven()

		polyEval.coeffsInterface = c

	case *coefficientsBSGSCiphertext:

		c.basis = polIn.basis
		c.even, c.odd = polIn.OddEven()

		polyEval.coeffsInterface = c
	}

	polyEval.recurse(dummbpb.Value[1].Level-giant+1, targetScale, *polIn)
}

func (polyEval *dummyPolynomialEvaluator) recurse(targetLevel int, targetScale float64, pol coefficientsComplex128) (res *dummyCiphertext) {

	params := polyEval.params

	baby := polyEval.baby

	// Recursively computes the evaluation of the Chebyshev polynomial using a baby-set giant-step algorithm.
	if pol.Degree() < (1 << baby) {

		if pol.lead && polyEval.baby > 1 && pol.maxDeg%(1<<(baby+1)) > (1<<(baby-1)) {

			giant := int(bits.Len64(uint64(pol.Degree())))
			baby := giant >> 1

			polyEvalBis := new(dummyPolynomialEvaluator)
			polyEvalBis.params = polyEval.params
			polyEvalBis.ecd = polyEval.ecd
			polyEvalBis.enc = polyEval.enc
			polyEvalBis.giant = giant
			polyEvalBis.baby = baby
			polyEvalBis.dummyPolynomialBasis = polyEval.dummyPolynomialBasis
			polyEvalBis.coeffsInterface = polyEval.coeffsInterface

			return polyEvalBis.recurse(targetLevel, targetScale, pol)
		}

		if pol.lead {
			targetScale *= params.QiFloat64(targetLevel)
		}

		switch coeffsInterface := polyEval.coeffsInterface.(type) {
		case *coefficientsBSGSComplex128:
			return polyEval.evaluatePolyFromPolynomialBasisComplex128(targetScale, targetLevel, pol, coeffsInterface)
		case *coefficientsBSGSPlaintext:
			return polyEval.evaluatePolyFromPolynomialBasisPlaintext(targetScale, targetLevel, pol, coeffsInterface)
		case *coefficientsBSGSCiphertext:
			return polyEval.evaluatePolyFromPolynomialBasisCiphertext(targetScale, targetLevel, pol, coeffsInterface)
		}
	}

	var nextPower = 1 << polyEval.baby
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

func (polyEval *dummyPolynomialEvaluator) evaluatePolyFromPolynomialBasisComplex128(targetScale float64, level int, polIn coefficientsComplex128, polOut *coefficientsBSGSComplex128) (res *dummyCiphertext) {

	X := polyEval.dummyPolynomialBasis.Value

	values := make([][]complex128, polIn.Degree()+1)
	for i := range values {
		values[i] = make([]complex128, len(polIn.coeffs))
	}

	for i, c := range polIn.coeffs {
		if isNotNegligible(c[0]) {
			values[0][i] = c[0] * complex(targetScale, 0)
		}
	}

	for i := 1; i < polIn.Degree(); i++ {
		for j, c := range polIn.coeffs {
			if isNotNegligible(c[i]) {
				values[i][j] = c[i] * complex(targetScale/X[i].Scale, 0)
			}
		}
	}

	return &dummyCiphertext{level, targetScale}
}

func (polyEval *dummyPolynomialEvaluator) evaluatePolyFromPolynomialBasisPlaintext(targetScale float64, level int, polIn coefficientsComplex128, polOut *coefficientsBSGSPlaintext) (res *dummyCiphertext) {

	params := polyEval.params
	ecd := polyEval.ecd

	slotsIndex := polIn.slotsIndex

	X := polyEval.dummyPolynomialBasis.Value

	polys := polIn.coeffs

	degree := polIn.Degree()
	nbPoly := len(polys)

	//[#poly][coefficients]
	pt := make([]*Plaintext, degree+1)

	values := make([]complex128, params.Slots())

	var toEncode bool

	for i := 0; i < nbPoly; i++ {

		c := polys[i][0] // [poly][coeff]

		if c != 0 {

			toEncode = true

			for _, j := range slotsIndex[i] {
				values[j] = c
			}
		}
	}

	if toEncode {

		pt[0] = ecd.EncodeNew(values, level, targetScale, params.LogSlots())

		for i := range values {
			values[i] = 0
		}

		toEncode = false
	}

	for i := 1; i < degree+1; i++ {

		for j := 0; j < nbPoly; j++ {

			c := polys[j][i] // [poly][coeff]

			if c != 0 {

				toEncode = true

				for _, k := range slotsIndex[j] {
					values[k] = c
				}
			}
		}

		if toEncode {

			pt[i] = ecd.EncodeNew(values, level, targetScale/X[i].Scale, params.LogSlots())

			for i := range values {
				values[i] = 0
			}

			toEncode = false
		}
	}

	return &dummyCiphertext{level, targetScale}
}

func (polyEval *dummyPolynomialEvaluator) evaluatePolyFromPolynomialBasisCiphertext(targetScale float64, level int, polIn coefficientsComplex128, polOut *coefficientsBSGSCiphertext) (res *dummyCiphertext) {

	params := polyEval.params
	ecd := polyEval.ecd
	enc := polyEval.enc

	slotsIndex := polIn.slotsIndex

	X := polyEval.dummyPolynomialBasis.Value

	polys := polIn.coeffs

	degree := polIn.Degree()
	nbPoly := len(polys)

	//[#poly][coefficients]
	ct := make([]*Ciphertext, degree+1)

	values := make([]complex128, params.Slots())
	pt := NewPlaintext(params, level, 0)

	var toEncrypt bool

	for i := 0; i < nbPoly; i++ {

		c := polys[i][0] // [poly][coeff]

		if c != 0 {

			toEncrypt = true

			for _, j := range slotsIndex[i] {
				values[j] = c
			}
		}
	}

	if toEncrypt {

		pt.Scale = targetScale
		ecd.Encode(values, pt, params.LogSlots())
		ct[0] = enc.EncryptNew(pt)

		for i := range values {
			values[i] = 0
		}

		toEncrypt = false
	}

	for i := 1; i < degree+1; i++ {

		for j := 0; j < nbPoly; j++ {

			c := polys[j][i] // [poly][coeff]

			if c != 0 {

				toEncrypt = true

				for _, k := range slotsIndex[j] {
					values[k] = c
				}
			}
		}

		if toEncrypt {

			pt.Scale = targetScale / X[i].Scale
			ecd.Encode(values, pt, params.LogSlots())
			ct[i] = enc.EncryptNew(pt)

			for i := range values {
				values[i] = 0
			}

			toEncrypt = false
		}
	}

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
