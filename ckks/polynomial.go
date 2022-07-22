package ckks

import (
	"encoding/binary"
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

func (p *Polynomial) MarshalBinary() (data []byte, err error) {

	if data, err = p.coefficients.MarshalBinary(); err != nil {
		return nil, err
	}

	data = append([]byte{0}, data...)

	switch p.coefficients.(type) {
	case *coefficientsComplex128:
		data[0] = 0
	case *coefficientsBSGSComplex128:
		data[0] = 1
	case *coefficientsBSGSPlaintext:
		data[0] = 2
	case *coefficientsBSGSCiphertext:
		data[0] = 3
	}

	return
}

func (p *Polynomial) UnmarshalBinary(data []byte) (err error) {

	switch data[0] {
	case 0:
		p.coefficients = new(coefficientsComplex128)
	case 1:
		p.coefficients = new(coefficientsBSGSComplex128)
	case 2:
		p.coefficients = new(coefficientsBSGSPlaintext)
	case 3:
		p.coefficients = new(coefficientsBSGSCiphertext)
	}

	return p.coefficients.UnmarshalBinary(data[1:])
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
func NewPolynomial(basisType BasisType, coeffs interface{}, slotsIndex [][]int) (poly Polynomial, err error) {

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

func (p *Polynomial) Encode(ecd Encoder, level int, inputScale, outputScale float64) (ptPoly Polynomial, err error) {

	params := ecd.(*encoderComplex128).params

	ptCoeffs := &coefficientsBSGSPlaintext{}

	switch coeffInterface := p.coefficients.(type) {
	case *coefficientsComplex128:

		ptCoeffs.basis = coeffInterface.basis
		ptCoeffs.odd = coeffInterface.odd
		ptCoeffs.even = coeffInterface.even

		getScaledBSGSCoefficients(params, ecd, nil, level, inputScale, *coeffInterface, outputScale, ptCoeffs)

	default:
		return ptPoly, fmt.Errorf("Polynomial.Encode(*): underlying polynomial is already encoded or encrypted")
	}

	return Polynomial{ptCoeffs}, nil
}

func (p *Polynomial) Encrypt(ecd Encoder, enc Encryptor, level int, inputScale, outputScale float64) (ctPoly Polynomial, err error) {

	params := ecd.(*encoderComplex128).params

	ctCoeffs := &coefficientsBSGSCiphertext{}

	switch coeffInterface := p.coefficients.(type) {
	case *coefficientsComplex128:

		ctCoeffs.basis = coeffInterface.basis
		ctCoeffs.odd = coeffInterface.odd
		ctCoeffs.even = coeffInterface.even

		getScaledBSGSCoefficients(params, ecd, enc, level, inputScale, *coeffInterface, outputScale, ctCoeffs)

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
		return ctPoly, fmt.Errorf("Polynomial.Encrypt(*): underlying polynomial is already encrypted")
	}

	return Polynomial{ctCoeffs}, nil
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
	SplitBSGS(split int) (coeffq, coeffsr coefficients)
	MarshalBinary() (data []byte, err error)
	UnmarshalBinary(data []byte) (err error)
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
	slotsIndex [][]int
	odd, even  bool
}

func (c *coefficientsComplex128) MarshalBinary() (data []byte, err error) {
	data = make([]byte, 3+4)

	data[0] = uint8(c.basis)

	if c.odd {
		data[1] = 1
	}

	if c.even {
		data[2] = 1
	}

	if c.slotsIndex != nil {

		var ptr int
		binary.LittleEndian.PutUint32(data[3:], uint32(len(c.slotsIndex)))

		for _, slots := range c.slotsIndex {

			dataSlots := make([]byte, 4+len(slots)<<2)

			binary.LittleEndian.PutUint32(dataSlots, uint32(len(slots)))

			ptr += 4
			for _, j := range slots {
				binary.LittleEndian.PutUint32(dataSlots[ptr:], uint32(j))
				ptr += 4
			}

			data = append(data, dataSlots...)
		}
	}

	data = append(data, make([]byte, 4)...)
	binary.LittleEndian.PutUint32(data[len(data)-4:], uint32(len(c.coeffs)))

	var ptr int
	for _, coeffs := range c.coeffs {

		dataCoeffs := make([]byte, 4+len(coeffs)<<5)

		binary.LittleEndian.PutUint32(dataCoeffs, uint32(len(coeffs)))
		ptr += 4

		for _, coeff := range coeffs {
			binary.LittleEndian.PutUint64(dataCoeffs[ptr:], math.Float64bits(real(coeff)))
			ptr += 8
			binary.LittleEndian.PutUint64(dataCoeffs[ptr:], math.Float64bits(imag(coeff)))
			ptr += 8
		}

		data = append(data, dataCoeffs...)
	}

	return
}

//
func (c *coefficientsComplex128) UnmarshalBinary(data []byte) (err error) {

	c.basis = BasisType(int(data[0]))

	c.odd = data[1] == 1
	c.even = data[2] == 1

	nbSlotIndex := int(binary.LittleEndian.Uint32(data[3:]))

	var ptr int
	ptr += 7

	if nbSlotIndex != 0 {

		c.slotsIndex = make([][]int, nbSlotIndex)

		for i := 0; i < nbSlotIndex; i++ {

			nbSlots := binary.LittleEndian.Uint32(data[ptr:])
			ptr += 4

			slots := make([]int, nbSlots)

			for i := range slots {
				slots[i] = int(binary.LittleEndian.Uint32(data[ptr:]))
				ptr += 4
			}

			c.slotsIndex[i] = slots
		}
	}

	c.coeffs = make([][]complex128, int(binary.LittleEndian.Uint32(data[ptr:])))
	ptr += 4

	for i := range c.coeffs {

		c.coeffs[i] = make([]complex128, int(binary.LittleEndian.Uint32(data[ptr:])))
		ptr += 4

		for j := range c.coeffs[i] {
			real := math.Float64frombits(binary.LittleEndian.Uint64(data[ptr:]))
			ptr += 8
			imag := math.Float64frombits(binary.LittleEndian.Uint64(data[ptr:]))
			ptr += 8

			c.coeffs[i][j] = complex(real, imag)
		}
	}

	return
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
	depth := c.Depth()
	return depth, optimalSplit(depth)
}

func (c *coefficientsComplex128) SplitBSGS(split int) (q, r coefficients) {

	polyq := &coefficientsComplex128{}
	polyr := &coefficientsComplex128{}

	polyq.coeffs = make([][]complex128, len(c.coeffs))
	polyr.coeffs = make([][]complex128, len(c.coeffs))

	for i, p := range c.coeffs {
		polyq.coeffs[i], polyr.coeffs[i] = splitCoeffsBSGS(p, split, c.basis)
	}

	polyq.basis = c.basis
	polyr.basis = c.basis

	polyq.slotsIndex = c.slotsIndex
	polyr.slotsIndex = c.slotsIndex

	return polyq, polyr
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
	slotsIndex [][]int
	odd, even  bool
}

func (c *coefficientsBSGSComplex128) MarshalBinary() (data []byte, err error) {
	return nil, fmt.Errorf("unsupported method")
}

func (c *coefficientsBSGSComplex128) UnmarshalBinary(data []byte) (err error) {
	return fmt.Errorf("unsupported method")
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
	depth := c.Depth()
	return depth, optimalSplit(depth)
}

func (c *coefficientsBSGSComplex128) SplitBSGS(split int) (q, r coefficients) {

	polyq := &coefficientsBSGSComplex128{}
	polyr := &coefficientsBSGSComplex128{}

	polyq.basis = c.basis
	polyr.basis = c.basis

	polyq.slotsIndex = c.slotsIndex
	polyr.slotsIndex = c.slotsIndex

	var deg, n int
	for i := range c.coeffs {

		if deg >= split {
			break
		}

		deg += len(c.coeffs[i])
		n++
	}

	polyq.coeffs = c.coeffs[n:]
	polyr.coeffs = c.coeffs[:n]

	return polyq, polyr
}

// coefficientsBSGSPlaintext: coefficients in baby-step giant-step format
// and encoded in plaintext.
// can store multiple polynomials: [giant-step][baby-step]*Plaintext{#poly}
type coefficientsBSGSPlaintext struct {
	basis     BasisType
	coeffs    [][]*Plaintext
	odd, even bool
}

func (c *coefficientsBSGSPlaintext) MarshalBinary() (data []byte, err error) {
	return nil, fmt.Errorf("unsupported method: there is no reason to do that (instead marshal non-encoded coefficients)")

}

func (c *coefficientsBSGSPlaintext) UnmarshalBinary(data []byte) (err error) {
	return fmt.Errorf("unsupported method: there is no reason to do that (instead unmarshal non-encoded coefficients)")
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
	depth := c.Depth()
	return depth, optimalSplit(depth)
}

func (c *coefficientsBSGSPlaintext) SplitBSGS(split int) (q, r coefficients) {

	polyq := &coefficientsBSGSPlaintext{}
	polyr := &coefficientsBSGSPlaintext{}

	polyq.basis = c.basis
	polyr.basis = c.basis

	var deg, n int
	for i := range c.coeffs {

		if deg >= split {
			break
		}

		deg += len(c.coeffs[i])
		n++
	}

	polyq.coeffs = c.coeffs[n:]
	polyr.coeffs = c.coeffs[:n]

	return polyq, polyr
}

// coefficientsBSGSPlaintext: coefficients in baby-step giant-step format
// and encrypted.
// can store multiple polynomials: [giant-step][baby-step]*Ciphertext{#poly}
type coefficientsBSGSCiphertext struct {
	basis     BasisType
	coeffs    [][]*Ciphertext
	odd, even bool
}

func (c *coefficientsBSGSCiphertext) MarshalBinary() (data []byte, err error) {
	data = make([]byte, 3+4)

	data[0] = uint8(c.basis)
	if c.odd {
		data[1] = 1
	}
	if c.even {
		data[1] = 1
	}

	ct := c.coeffs

	binary.LittleEndian.PutUint32(data[3:], uint32(len(ct)))

	for i := range ct {

		deg := make([]byte, 4)

		binary.LittleEndian.PutUint32(deg, uint32(len(ct[i])))

		data = append(data, deg...)

		for j := range c.coeffs[i] {

			if ct[i][j] != nil {

				var dataCt []byte
				if dataCt, err = ct[i][j].MarshalBinary(); err != nil {
					return nil, err
				}

				data = append(data, 1)

				ctLen := make([]byte, 4)
				binary.LittleEndian.PutUint32(ctLen, uint32(len(dataCt)))
				data = append(data, ctLen...)
				data = append(data, dataCt...)

			} else {
				data = append(data, 0)
			}
		}
	}

	return
}

func (c *coefficientsBSGSCiphertext) UnmarshalBinary(data []byte) (err error) {

	c.basis = BasisType(data[0])
	c.odd = data[1] == 1
	c.even = data[2] == 1

	ct := make([][]*Ciphertext, binary.LittleEndian.Uint32(data[3:]))

	var ptr = 7
	for i := range ct {

		ct[i] = make([]*Ciphertext, binary.LittleEndian.Uint32(data[ptr:]))
		ptr += 4

		for j := range ct[i] {

			ptr++

			if data[ptr-1] == 1 {

				ctLen := int(binary.LittleEndian.Uint32(data[ptr:]))
				ptr += 4

				ct[i][j] = new(Ciphertext)
				if err = ct[i][j].UnmarshalBinary(data[ptr : ptr+ctLen]); err != nil {
					return
				}
				ptr += ctLen
			}
		}
	}

	c.coeffs = ct
	return
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
	depth := c.Depth()
	return depth, optimalSplit(depth)
}

func (c *coefficientsBSGSCiphertext) SplitBSGS(split int) (q, r coefficients) {

	polyq := &coefficientsBSGSCiphertext{}
	polyr := &coefficientsBSGSCiphertext{}

	polyq.basis = c.basis
	polyr.basis = c.basis

	var deg, n int
	for i := range c.coeffs {

		if deg >= split {
			break
		}

		deg += len(c.coeffs[i])
		n++
	}

	polyq.coeffs = c.coeffs[n:]
	polyr.coeffs = c.coeffs[:n]

	return polyq, polyr
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

func getScaledBSGSCoefficients(params Parameters, ecd Encoder, enc Encryptor, level int, scale float64, polIn coefficientsComplex128, targetScale float64, polOut coefficients) {

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

		c.coeffs = [][]*Plaintext{}

		polyEval.coeffsInterface = c

	case *coefficientsBSGSCiphertext:

		c.basis = polIn.basis
		c.even, c.odd = polIn.OddEven()

		c.coeffs = [][]*Ciphertext{}

		polyEval.coeffsInterface = c
	}

	polyEval.recurse(dummbpb.Value[1].Level-giant+1, targetScale, evalPoly{&polIn, true, polIn.Degree()})
}

func (polyEval *dummyPolynomialEvaluator) recurse(targetLevel int, targetScale float64, pol evalPoly) (res *dummyCiphertext) {

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
			return polyEval.evaluatePolyFromPolynomialBasisComplex128(targetScale, targetLevel, pol.coefficients.(*coefficientsComplex128), coeffsInterface)
		case *coefficientsBSGSPlaintext:
			return polyEval.evaluatePolyFromPolynomialBasisPlaintext(targetScale, targetLevel, pol.coefficients.(*coefficientsComplex128), coeffsInterface)
		case *coefficientsBSGSCiphertext:
			return polyEval.evaluatePolyFromPolynomialBasisCiphertext(targetScale, targetLevel, pol.coefficients.(*coefficientsComplex128), coeffsInterface)
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

func (polyEval *dummyPolynomialEvaluator) evaluatePolyFromPolynomialBasisComplex128(targetScale float64, level int, polIn *coefficientsComplex128, polOut *coefficientsBSGSComplex128) (res *dummyCiphertext) {

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

func (polyEval *dummyPolynomialEvaluator) evaluatePolyFromPolynomialBasisPlaintext(targetScale float64, level int, polIn *coefficientsComplex128, polOut *coefficientsBSGSPlaintext) (res *dummyCiphertext) {

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

	polOut.coeffs = append([][]*Plaintext{pt}, polOut.coeffs...)

	return &dummyCiphertext{level, targetScale}
}

func (polyEval *dummyPolynomialEvaluator) evaluatePolyFromPolynomialBasisCiphertext(targetScale float64, level int, polIn *coefficientsComplex128, polOut *coefficientsBSGSCiphertext) (res *dummyCiphertext) {

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

	var values []complex128
	var pt *Plaintext
	if slotsIndex != nil {
		values = make([]complex128, params.Slots())
		pt = NewPlaintext(params, level, 0)
	} else {
		values = make([]complex128, 1)
	}

	var toEncrypt bool

	for i := 0; i < nbPoly; i++ {

		c := polys[i][0] // [poly][coeff]

		if c != 0 {

			toEncrypt = true

			if slotsIndex != nil {
				for _, j := range slotsIndex[i] {
					values[j] = c
				}
			} else {
				values[0] = c
			}
		}
	}

	if toEncrypt {
		if slotsIndex != nil {
			pt.Scale = targetScale
			ecd.Encode(values, pt, params.LogSlots())
			ct[0] = enc.EncryptNew(pt)
		} else {
			ct[0] = enc.EncryptZeroNew(level, targetScale)
			addConst(params, ct[0], values[0], ct[0])
		}

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

				if slotsIndex != nil {
					for _, k := range slotsIndex[j] {
						values[k] = c
					}
				} else {
					values[0] = c
				}

			}
		}

		if toEncrypt {

			if slotsIndex != nil {
				pt.Scale = targetScale / X[i].Scale
				ecd.Encode(values, pt, params.LogSlots())
				ct[i] = enc.EncryptNew(pt)
			} else {
				ct[i] = enc.EncryptZeroNew(level, targetScale/X[i].Scale)
				addConst(params, ct[i], values[0], ct[i])
			}

			for i := range values {
				values[i] = 0
			}

			toEncrypt = false
		}
	}

	polOut.coeffs = append([][]*Ciphertext{ct}, polOut.coeffs...)

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
