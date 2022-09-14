package rlwe

import (
	"encoding/binary"
	"fmt"
	"math"
	"math/bits"

	"github.com/tuneinsight/lattigo/v3/utils"
)

// PolynomialBasisType is a type for the polynomials basis
type PolynomialBasisType int

const (
	// Monomial : x^(a+b) = x^a * x^b
	Monomial = PolynomialBasisType(0)
	// Chebyshev : T_(a+b) = 2 * T_a * T_b - T_(|a-b|)
	Chebyshev = PolynomialBasisType(1)
)

// Polynomial stores the coefficient in baby-step giant-step format.
// can store multiple polynomials: [#poly][coefficients]
type Polynomial struct {
	Coefficients
	Basis      PolynomialBasisType
	SlotsIndex [][]int
	Odd, Even  bool
}

// NewPolynomial creates a new polynomial from:
//
// PolynomialBasis: ckks.Monomial or ckks.Chebyshev
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
func NewPolynomial(polynomialBasis PolynomialBasisType, coeffs interface{}, slotsIndex [][]int) (Polynomial, error) {

	var odd, even bool

	var poly Polynomial

	switch coeffs := coeffs.(type) {
	case []complex128:

		c := make([]complex128, len(coeffs))
		copy(c, coeffs)

		tmp0, tmp1 := isOddOrEvenFloat(c)
		odd, even = odd && tmp0, even && tmp1

		poly = Polynomial{
			Coefficients: Coefficients{Value: [][]complex128{c}},
			SlotsIndex:   slotsIndex,
			Odd:          odd,
			Even:         even,
			Basis:        polynomialBasis,
		}

	case []float64:

		c := make([]complex128, len(coeffs))
		for i := range coeffs {
			c[i] = complex(coeffs[i], 0)
		}

		tmp0, tmp1 := isOddOrEvenFloat(c)
		odd, even = odd && tmp0, even && tmp1

		poly = Polynomial{
			Coefficients: Coefficients{Value: [][]complex128{c}},
			SlotsIndex:   slotsIndex,
			Odd:          odd,
			Even:         even,
			Basis:        polynomialBasis,
		}

	case [][]complex128:
		{

			if slotsIndex == nil {
				return Polynomial{}, fmt.Errorf("NewPolynomial: cannot have slotsIndex == nil if coeffs == [][]complex128")
			}

			var maxDeg int
			for i := range coeffs {
				maxDeg = utils.MaxInt(maxDeg, len(coeffs[i]))
			}

			c := make([][]complex128, len(coeffs))
			for i := range c {
				c[i] = make([]complex128, maxDeg)
				copy(c[i], coeffs[i])

				tmp0, tmp1 := isOddOrEvenFloat(c[i])
				odd, even = odd && tmp0, even && tmp1
			}

			poly = Polynomial{
				Coefficients: Coefficients{Value: c},
				SlotsIndex:   slotsIndex,
				Odd:          odd,
				Even:         even,
				Basis:        polynomialBasis,
			}
		}

	case [][]float64:
		{
			if slotsIndex == nil {
				return Polynomial{}, fmt.Errorf("NewPolynomial: cannot have slotsIndex == nil if coeffs == [][]float64")
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

				tmp0, tmp1 := isOddOrEvenFloat(c[i])
				odd, even = odd && tmp0, even && tmp1
			}

			poly = Polynomial{
				Coefficients: Coefficients{Value: c},
				SlotsIndex:   slotsIndex,
				Odd:          odd,
				Even:         even,
				Basis:        polynomialBasis,
			}
		}

	default:
		return Polynomial{}, fmt.Errorf("NewPolynomial: invalid coeffs.(type)")
	}

	return poly, nil
}

func (p *Polynomial) IsEncrypted() bool {
	switch coeffs := p.Coefficients.Value.(type) {
	case [][]Operand:
		for i := range coeffs {
			for j := range coeffs[i] {
				if coeffs[i][j] != nil {

					if coeffs[i][j].El().Degree() == 1 {
						return true
					} else {
						return false
					}
				}
			}
		}
		return false
	default:
		return false
	}
}

func (p *Polynomial) BasisType() PolynomialBasisType {
	return p.Basis
}

func (p *Polynomial) Depth() int {
	return bits.Len64(uint64(p.Degree()))
}

func (p *Polynomial) Degree() int {
	var deg int
	switch coeffs := p.Coefficients.Value.(type) {
	case [][]float64:
		deg = len(coeffs[0])
	case [][]complex128:
		deg = len(coeffs[0])
	case [][]uint64:
		deg = len(coeffs[0])
	case [][]Operand:
		for _, ci := range coeffs {
			deg += len(ci)
		}
	default:
		panic("invalid type")
	}

	return deg - 1
}

func (p *Polynomial) OddEven() (odd, even bool) {
	return p.Odd, p.Even
}

func (p *Polynomial) BSGSSplit() (giant, baby int) {
	depth := p.Depth()
	return depth, OptimalSplit(depth)
}

func (p *Polynomial) SplitBSGS(split int) (coeffsq, coeffsr Polynomial) {

	coeffsq = Polynomial{}
	coeffsq.Basis = p.Basis
	coeffsq.Odd = p.Odd
	coeffsq.Even = p.Even
	coeffsq.SlotsIndex = p.SlotsIndex

	coeffsr = Polynomial{}
	coeffsr.Basis = p.Basis
	coeffsr.Odd = p.Odd
	coeffsr.Even = p.Even
	coeffsr.SlotsIndex = p.SlotsIndex

	switch coeffs := p.Coefficients.Value.(type) {
	case [][]uint64:
		coeffsq.Coefficients.Value, coeffsr.Coefficients.Value = splitBSGS(coeffs, split, p.Basis)
	case [][]float64:
		coeffsq.Coefficients.Value, coeffsr.Coefficients.Value = splitBSGS(coeffs, split, p.Basis)
	case [][]complex128:
		coeffsq.Coefficients.Value, coeffsr.Coefficients.Value = splitBSGS(coeffs, split, p.Basis)
	case [][]Operand:
		coeffsq.Coefficients.Value, coeffsr.Coefficients.Value = splitBSGSRLWE(coeffs, split)
	}

	return
}

func splitBSGSRLWE(coeffs [][]Operand, split int) (coeffsq, coeffsr [][]Operand) {

	var deg, n int
	for i := range coeffs {

		if deg >= split {
			break
		}

		deg += len(coeffs[i])
		n++
	}

	return coeffs[n:], coeffs[:n]
}

func splitBSGS[T uint64 | float64 | complex128](coeffs [][]T, split int, basis PolynomialBasisType) (coeffsq, coeffsr [][]T) {

	coeffsq = make([][]T, len(coeffs))
	coeffsr = make([][]T, len(coeffs))

	for i, c := range coeffs {

		// Splits a polynomial p such that p = q*C^degree + r.
		r := make([]T, split)
		q := make([]T, len(c)-split)

		q[0] = c[split]
		for j := 0; j < split; j++ {
			r[j] = c[j]
		}

		if basis == Monomial {
			for j := split + 1; j < len(c); j++ {
				q[j-split] = c[j]
			}
		} else if basis == Chebyshev {
			for j, k := split+1, 1; j < len(c); j, k = j+1, k+1 {
				q[j-split] = 2 * c[j]
				r[split-k] -= c[j]
			}
		}

		coeffsq[i] = q
		coeffsr[i] = r
	}

	return
}

// IsNegligbleThreshold : threshold under which a coefficient
// of a polynomial is ignored.
const IsNegligbleThreshold float64 = 1e-14

func IsNotNegligible(c complex128) bool {
	return (math.Abs(real(c)) > IsNegligbleThreshold || math.Abs(imag(c)) > IsNegligbleThreshold)
}

func isOddOrEvenFloat(coeffs []complex128) (odd, even bool) {
	even = true
	odd = true
	for i, c := range coeffs {
		isNotNegligible := IsNotNegligible(c)
		odd = odd && !(i&1 == 0 && isNotNegligible)
		even = even && !(i&1 == 1 && isNotNegligible)
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

type Coefficients struct {
	Scale
	Value interface{}
}

func (c *Coefficients) MarshalBinary() (data []byte, err error) {

	data = make([]byte, 5)

	switch Value := c.Value.(type) {
	case [][]uint64:

		data[0] = 0

		binary.LittleEndian.PutUint32(data[1:], uint32(len(Value)))

		for _, coeffs := range Value {

			var ptr int

			tmp := make([]byte, 4+len(coeffs)<<3)
			binary.LittleEndian.PutUint32(tmp, uint32(len(coeffs)))
			ptr += 4

			for _, c := range coeffs {
				binary.LittleEndian.PutUint64(tmp[ptr:], c)
				ptr += 8
			}

			data = append(data, tmp...)
		}

		return

	case [][]complex128:

		data[0] = 1

		binary.LittleEndian.PutUint32(data[1:], uint32(len(Value)))

		for _, coeffs := range Value {

			var ptr int

			tmp := make([]byte, 4+len(coeffs)<<4)
			binary.LittleEndian.PutUint32(tmp, uint32(len(coeffs)))
			ptr += 4

			for _, c := range coeffs {
				binary.LittleEndian.PutUint64(tmp[ptr:], math.Float64bits(real(c)))
				ptr += 8
				binary.LittleEndian.PutUint64(tmp[ptr:], math.Float64bits(imag(c)))
				ptr += 8
			}

			data = append(data, tmp...)
		}

		return

	case [][]Operand:

		switch Value[0][0].(type) {
		case *Plaintext:
			return nil, fmt.Errorf("unsupported method: there is no reason to do that (instead marshal non-encoded coefficients)")
		case *Ciphertext:

			data[0] = 2

			ct := Value

			binary.LittleEndian.PutUint32(data[1:], uint32(len(ct)))

			for i := range ct {

				deg := make([]byte, 4)

				binary.LittleEndian.PutUint32(deg, uint32(len(ct[i])))

				data = append(data, deg...)

				for j := range ct[i] {

					if ct[i][j] != nil {

						var dataCt []byte
						if dataCt, err = ct[i][j].El().MarshalBinary(); err != nil {
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
		}

		return

	default:
		return
	}
}

func (c *Coefficients) UnmarshalBinary(data []byte) (err error) {

	ptr := 1

	switch data[0] {
	case 0:

		Value := make([][]uint64, int(binary.LittleEndian.Uint32(data[ptr:])))
		ptr += 4

		for i := range Value {

			coeffs := make([]uint64, binary.LittleEndian.Uint32(data[ptr:]))
			ptr += 4

			for j := range coeffs {
				coeffs[j] = binary.LittleEndian.Uint64(data[ptr:])
				ptr += 8
			}

			Value[i] = coeffs
		}

		c.Value = Value

		return
	case 1:

		Value := make([][]complex128, int(binary.LittleEndian.Uint32(data[ptr:])))
		ptr += 4

		for i := range Value {

			coeffs := make([]complex128, binary.LittleEndian.Uint32(data[ptr:]))
			ptr += 4

			for j := range coeffs {
				coeffs[j] = complex(math.Float64frombits(binary.LittleEndian.Uint64(data[ptr:])), math.Float64frombits(binary.LittleEndian.Uint64(data[ptr+8:])))
				ptr += 16
			}

			Value[i] = coeffs
		}

		c.Value = Value

		return
	case 2:

		cts := make([][]Operand, binary.LittleEndian.Uint32(data[ptr:]))

		ptr += 4
		for i := range cts {

			cts[i] = make([]Operand, binary.LittleEndian.Uint32(data[ptr:]))
			ptr += 4

			for j := range cts[i] {

				ptr++

				if data[ptr-1] == 1 {

					ctLen := int(binary.LittleEndian.Uint32(data[ptr:]))
					ptr += 4

					ct := new(Ciphertext)
					ct.Scale = c.Scale.CopyNew()
					if err = ct.UnmarshalBinary(data[ptr : ptr+ctLen]); err != nil {
						return
					}

					cts[i][j] = ct

					ptr += ctLen
				}
			}
		}

		c.Value = cts

	default:
		return
	}

	return
}

func OptimalSplit(giant int) (baby int) {
	baby = giant >> 1
	a := (1 << baby) + (1 << (giant - baby)) + giant - baby - 3
	b := (1 << (baby + 1)) + (1 << (giant - baby - 1)) + giant - baby - 4
	if a > b {
		baby++
	}

	return
}

func (p *Polynomial) MarshalBinary() (data []byte, err error) {

	data = make([]byte, 3)

	data[0] = uint8(p.Basis)

	if p.Odd {
		data[1] = 1
	}

	if p.Even {
		data[2] = 1
	}

	var dataCoeffs []byte
	if dataCoeffs, err = p.Coefficients.MarshalBinary(); err != nil {
		return
	}

	data = append(data, make([]byte, 8)...)
	binary.LittleEndian.PutUint64(data[len(data)-8:], uint64(len(dataCoeffs)))
	data = append(data, dataCoeffs...)

	data = append(data, make([]byte, 4)...)
	binary.LittleEndian.PutUint32(data[len(data)-4:], uint32(len(p.SlotsIndex)))

	var dataSlotsIndex []byte
	if p.SlotsIndex != nil {

		dataSlotsIndex = make([]byte, 4)

		var ptr int
		binary.LittleEndian.PutUint32(dataSlotsIndex, uint32(len(p.SlotsIndex)))

		for _, slots := range p.SlotsIndex {

			tmp := make([]byte, 4+len(slots)<<2)

			binary.LittleEndian.PutUint32(tmp, uint32(len(slots)))

			ptr += 4
			for _, j := range slots {
				binary.LittleEndian.PutUint32(tmp[ptr:], uint32(j))
				ptr += 4
			}

			dataSlotsIndex = append(dataSlotsIndex, tmp...)
		}

		data = append(data, dataSlotsIndex...)
	}

	return
}

func (p *Polynomial) UnmarshalBinary(data []byte) (err error) {

	p.Basis = PolynomialBasisType(int(data[0]))
	p.Odd = data[1] == 1
	p.Even = data[2] == 1

	lenDataCoeffs := binary.LittleEndian.Uint64(data[3:])

	if err = p.Coefficients.UnmarshalBinary(data[11:]); err != nil {
		return
	}

	lenSlotsIndex := int(binary.LittleEndian.Uint32(data[11+lenDataCoeffs:]))

	if lenSlotsIndex != 0 {

		ptr := 8 + 8 + 3 + lenDataCoeffs

		p.SlotsIndex = make([][]int, lenSlotsIndex)

		for i := 0; i < lenSlotsIndex; i++ {

			nbSlots := binary.LittleEndian.Uint32(data[ptr:])
			ptr += 4

			slots := make([]int, nbSlots)

			for i := range slots {
				slots[i] = int(binary.LittleEndian.Uint32(data[ptr:]))
				ptr += 4
			}

			p.SlotsIndex[i] = slots
		}
	}

	return
}
