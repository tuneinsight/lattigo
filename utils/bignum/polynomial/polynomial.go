// Package polynomial provides helper for polynomials, approximation of functions using polynomials and their evaluation.
package polynomial

import (
	"fmt"
	"math"
	"math/big"

	"github.com/tuneinsight/lattigo/v4/utils/bignum"
)

// Basis is a type for the polynomials basis
type Basis int

const (
	// Monomial : x^(a+b) = x^a * x^b
	Monomial = Basis(0)
	// Chebyshev : T_(a+b) = 2 * T_a * T_b - T_(|a-b|)
	Chebyshev = Basis(1)
)

type Interval struct {
	A, B *big.Float
}

type Polynomial struct {
	Basis
	Interval
	Coeffs []*bignum.Complex
	IsOdd  bool
	IsEven bool
}

// NewPolynomial creates a new polynomial from the input parameters:
// basis: either `Monomial` or `Chebyshev`
// coeffs: []bignum.Complex128, []float64, []*bignum.Complex or []*big.Float
// interval: [2]float64{a, b} or *Interval
func NewPolynomial(basis Basis, coeffs interface{}, interval interface{}) *Polynomial {
	var coefficients []*bignum.Complex

	switch coeffs := coeffs.(type) {
	case []complex128:
		coefficients = make([]*bignum.Complex, len(coeffs))
		for i := range coeffs {
			if c := coeffs[i]; c != 0 {
				coefficients[i] = &bignum.Complex{
					new(big.Float).SetFloat64(real(c)),
					new(big.Float).SetFloat64(imag(c)),
				}
			}
		}
	case []float64:
		coefficients = make([]*bignum.Complex, len(coeffs))
		for i := range coeffs {
			if c := coeffs[i]; c != 0 {
				coefficients[i] = &bignum.Complex{
					new(big.Float).SetFloat64(c),
					new(big.Float),
				}
			}
		}
	case []*bignum.Complex:
		coefficients = make([]*bignum.Complex, len(coeffs))
		copy(coefficients, coeffs)
	case []*big.Float:
		coefficients = make([]*bignum.Complex, len(coeffs))
		for i := range coeffs {
			if coeffs[i] != nil {
				coefficients[i] = &bignum.Complex{
					new(big.Float).Set(coeffs[i]),
					new(big.Float),
				}
			}
		}
	default:
		panic(fmt.Sprintf("invalid coefficient type, allowed types are []{bignum.Complex128, float64, *bignum.Complex, *big.Float} but is %T", coeffs))
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
		Basis:    basis,
		Interval: inter,
		Coeffs:   coefficients,
		IsOdd:    true,
		IsEven:   true,
	}
}

// ChangeOfBasis returns change of basis required to evaluate the polynomial
// Change of basis is defined as follow:
// - Monomial: scalar=1, constant=0.
// - Chebyshev: scalar=2/(b-a), constant = (-a-b)/(b-a).
func (p *Polynomial) ChangeOfBasis() (scalar, constant *big.Float) {

	switch p.Basis {
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
		panic(fmt.Sprintf("invalid basis type, allowed types are `Monomial` or `Chebyshev` but is %T", p.Basis))
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

// Evaluate takes x a *big.Float or *big.bignum.Complex and returns y = P(x).
// The precision of x is used as reference precision for y.
func (p *Polynomial) Evaluate(x interface{}) (y *bignum.Complex) {

	var xcmplx *bignum.Complex
	switch x := x.(type) {
	case *big.Float:
		xcmplx = bignum.ToComplex(x, x.Prec())
	case *bignum.Complex:
		xcmplx = bignum.ToComplex(x, x.Prec())
	default:
		panic(fmt.Errorf("cannot Evaluate: accepted x.(type) are *big.Float and *bignum.Complex but x is %T", x))
	}

	coeffs := p.Coeffs

	n := len(coeffs)

	mul := bignum.NewComplexMultiplier()

	switch p.Basis {
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

		tmp := &bignum.Complex{new(big.Float), new(big.Float)}

		scalar, constant := p.ChangeOfBasis()

		xcmplx[0].Mul(xcmplx[0], scalar)
		xcmplx[1].Mul(xcmplx[1], scalar)

		xcmplx[0].Add(xcmplx[0], constant)
		xcmplx[1].Add(xcmplx[1], constant)

		TPrev := &bignum.Complex{new(big.Float).SetInt64(1), new(big.Float)}

		T := xcmplx
		if coeffs[0] != nil {
			y = coeffs[0].Copy()
		} else {
			y = &bignum.Complex{new(big.Float), new(big.Float)}
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
		panic(fmt.Sprintf("invalid basis type, allowed types are `Monomial` or `Chebyshev` but is %T", p.Basis))
	}

	return
}

// Factorize factorizes p as X^{n} * pq + pr.
func (p *Polynomial) Factorize(n int) (pq, pr *Polynomial) {

	// ns a polynomial p such that p = q*C^degree + r.
	pr = &Polynomial{}
	pr.Coeffs = make([]*bignum.Complex, n)
	for i := 0; i < n; i++ {
		if p.Coeffs[i] != nil {
			pr.Coeffs[i] = p.Coeffs[i].Copy()
		}
	}

	pq = &Polynomial{}
	pq.Coeffs = make([]*bignum.Complex, p.Degree()-n+1)

	if p.Coeffs[n] != nil {
		pq.Coeffs[0] = p.Coeffs[n].Copy()
	}

	odd := p.IsOdd
	even := p.IsEven

	switch p.Basis {
	case Monomial:
		for i := n + 1; i < p.Degree()+1; i++ {
			if p.Coeffs[i] != nil && (!(even || odd) || (i&1 == 0 && even) || (i&1 == 1 && odd)) {
				pq.Coeffs[i-n] = p.Coeffs[i].Copy()
			}
		}
	case Chebyshev:

		for i, j := n+1, 1; i < p.Degree()+1; i, j = i+1, j+1 {
			if p.Coeffs[i] != nil && (!(even || odd) || (i&1 == 0 && even) || (i&1 == 1 && odd)) {
				pq.Coeffs[i-n] = p.Coeffs[i].Copy()
				pq.Coeffs[i-n].Add(pq.Coeffs[i-n], pq.Coeffs[i-n])

				if pr.Coeffs[n-j] != nil {
					pr.Coeffs[n-j].Sub(pr.Coeffs[n-j], p.Coeffs[i])
				} else {
					pr.Coeffs[n-j] = p.Coeffs[i].Copy()
					pr.Coeffs[n-j][0].Neg(pr.Coeffs[n-j][0])
					pr.Coeffs[n-j][1].Neg(pr.Coeffs[n-j][1])
				}
			}
		}
	}

	pq.Basis, pr.Basis = p.Basis, p.Basis
	pq.IsOdd, pr.IsOdd = p.IsOdd, p.IsOdd
	pq.IsEven, pr.IsEven = p.IsEven, p.IsEven
	pq.Interval, pr.Interval = p.Interval, p.Interval

	return
}
