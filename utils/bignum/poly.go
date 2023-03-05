package bignum

import (
	"fmt"
	"math"
	"math/big"
)

// BasisType is a type for the polynomials basis
type BasisType int

const (
	// Monomial : x^(a+b) = x^a * x^b
	Monomial = BasisType(0)
	// Chebyshev : T_(a+b) = 2 * T_a * T_b - T_(|a-b|)
	Chebyshev = BasisType(1)
)

type Interval struct {
	A, B *big.Float
}

type Polynomial struct {
	BasisType
	Interval
	Coeffs []*Complex
	IsOdd  bool
	IsEven bool
}

// NewPolynomial creates a new polynomial from the input parameters:
// basis: either `Monomial` or `Chebyshev`
// coeffs: []complex128, []float64, []*Complex or []*big.Float
// interval: [2]float64{a, b} or *Interval
func NewPolynomial(basis BasisType, coeffs interface{}, interval interface{}) *Polynomial {
	var coefficients []*Complex

	switch coeffs := coeffs.(type) {
	case []complex128:
		coefficients = make([]*Complex, len(coeffs))
		for i := range coeffs {
			if c := coeffs[i]; c != 0 {
				coefficients[i] = &Complex{
					new(big.Float).SetFloat64(real(c)),
					new(big.Float).SetFloat64(imag(c)),
				}
			}
		}
	case []float64:
		coefficients = make([]*Complex, len(coeffs))
		for i := range coeffs {
			if c := coeffs[i]; c != 0 {
				coefficients[i] = &Complex{
					new(big.Float).SetFloat64(c),
					new(big.Float),
				}
			}
		}
	case []*Complex:
		coefficients = make([]*Complex, len(coeffs))
		copy(coefficients, coeffs)
	case []*big.Float:
		coefficients = make([]*Complex, len(coeffs))
		for i := range coeffs {
			if coeffs[i] != nil {
				coefficients[i] = &Complex{
					new(big.Float).Set(coeffs[i]),
					new(big.Float),
				}
			}
		}
	default:
		panic(fmt.Sprintf("invalid coefficient type, allowed types are []{complex128, float64, *Complex, *big.Float} but is %T", coeffs))
	}

	inter := Interval{}
	switch interval := interval.(type) {
	case [2]float64:
		inter.A = new(big.Float).SetFloat64(interval[0])
		inter.B = new(big.Float).SetFloat64(interval[1])
	case *Interval:
		inter.A = new(big.Float).Set(interval.A)
		inter.B = new(big.Float).Set(interval.B)
	case nil:

	default:
		panic(fmt.Sprintf("invalid interval type, allowed types are [2]float64 or *Interval, but is %T", interval))
	}

	return &Polynomial{
		BasisType: basis,
		Interval:  inter,
		Coeffs:    coefficients,
		IsOdd:     true,
		IsEven:    true,
	}
}

// ChangeOfBasis returns change of basis required to evaluate the polynomial
// Change of basis is defined as follow:
// - Monomial: scalar=1, constant=0.
// - Chebyshev: scalar=2/(b-a), constant = (-a-b)/(b-a).
func (p *Polynomial) ChangeOfBasis() (scalar, constant *big.Float) {

	switch p.BasisType {
	case Monomial:
		scalar = new(big.Float).SetInt64(1)
		constant = new(big.Float)
	case Chebyshev:
		num := new(big.Float).Sub(p.B, p.A)

		// 2 / (b-a)
		scalar = new(big.Float).Quo(new(big.Float).SetInt64(2), num)

		// (-b-a)/(b-a)
		constant = new(big.Float).Set(p.B)
		constant.Neg(constant)
		constant.Sub(constant, p.A)
		constant.Quo(constant, num)
	default:
		panic(fmt.Sprintf("invalid basis type, allowed types are `Monomial` or `Chebyshev` but is %T", p.BasisType))
	}

	return
}

// Depth returns the number of sequential multiplications needed to evaluate the polynomial.
func (p *Polynomial) Depth() int {
	return int(math.Ceil(math.Log2(float64(p.Degree()))))
}

// Degree returns the degree of the polynomial.
func (p *Polynomial) Degree() int {
	return len(p.Coeffs) - 1
}

// Evaluate takes x a *big.Float or *big.Complex and returns y = P(x).
// The precision of x is used as reference precision for y.
func (p *Polynomial) Evaluate(x interface{}) (y *Complex) {

	var xcmplx *Complex
	switch x := x.(type) {
	case *big.Float:
		xcmplx = ToComplex(x, x.Prec())
	case *Complex:
		xcmplx = ToComplex(x, x.Prec())
	default:
		panic(fmt.Errorf("cannot Evaluate: accepted x.(type) are *big.Float and *Complex but x is %T", x))
	}

	coeffs := p.Coeffs

	n := len(coeffs)

	mul := NewComplexMultiplier()

	switch p.BasisType {
	case Monomial:
		y = coeffs[n-1].Copy()
		y.SetPrec(xcmplx.Prec())
		for i := n - 2; i >= 0; i-- {
			mul.Mul(y, xcmplx, y)
			if coeffs[i] != nil {
				y.Add(y, coeffs[i])
			}
		}

	case Chebyshev:

		tmp := &Complex{new(big.Float), new(big.Float)}

		scalar, constant := p.ChangeOfBasis()

		xcmplx[0].Mul(xcmplx[0], scalar)
		xcmplx[1].Mul(xcmplx[1], scalar)

		xcmplx[0].Add(xcmplx[0], constant)
		xcmplx[1].Add(xcmplx[1], constant)

		TPrev := &Complex{new(big.Float).SetInt64(1), new(big.Float)}

		T := xcmplx
		if coeffs[0] != nil {
			y = coeffs[0].Copy()
		} else {
			y = &Complex{new(big.Float), new(big.Float)}
		}

		y.SetPrec(xcmplx.Prec())

		two := new(big.Float).SetInt64(2)
		for i := 1; i < n; i++ {

			if coeffs[i] != nil {
				mul.Mul(T, coeffs[i], tmp)
				y.Add(y, tmp)
			}

			tmp[0].Mul(xcmplx[0], two)
			tmp[1].Mul(xcmplx[1], two)

			mul.Mul(tmp, T, tmp)
			tmp.Sub(tmp, TPrev)

			TPrev = T.Copy()
			T = tmp.Copy()
		}

	default:
		panic(fmt.Sprintf("invalid basis type, allowed types are `Monomial` or `Chebyshev` but is %T", p.BasisType))
	}

	return
}
