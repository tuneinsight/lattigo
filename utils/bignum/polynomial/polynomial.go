// Package polynomial provides helper for polynomials, approximation of functions using polynomials and their evaluation.
package polynomial

import (
	"fmt"
	"math"
	"math/big"

	"github.com/tuneinsight/lattigo/v4/utils/bignum"
)

type Polynomial struct {
	MetaData
	Coeffs []*bignum.Complex
}

func (p *Polynomial) Clone() *Polynomial {
	Coeffs := make([]*bignum.Complex, len(p.Coeffs))
	for i := range Coeffs {
		Coeffs[i] = p.Coeffs[i].Clone()
	}

	return &Polynomial{
		MetaData: p.MetaData,
		Coeffs:   Coeffs,
	}
}

// NewPolynomial creates a new polynomial from the input parameters:
// basis: either `Monomial` or `Chebyshev`
// coeffs: []bignum.Complex128, []float64, []*bignum.Complex or []*big.Float
// interval: [2]float64{a, b} or *Interval
func NewPolynomial(basis Basis, coeffs interface{}, interval interface{}) *Polynomial {
	var coefficients []*bignum.Complex

	switch coeffs := coeffs.(type) {
	case []uint64:
		coefficients = make([]*bignum.Complex, len(coeffs))
		for i, c := range coeffs {
			coefficients[i] = &bignum.Complex{
				new(big.Float).SetUint64(c),
				new(big.Float),
			}
		}
	case []complex128:
		coefficients = make([]*bignum.Complex, len(coeffs))
		for i, c := range coeffs {
			coefficients[i] = &bignum.Complex{
				new(big.Float).SetFloat64(real(c)),
				new(big.Float).SetFloat64(imag(c)),
			}
		}
	case []float64:
		coefficients = make([]*bignum.Complex, len(coeffs))
		for i, c := range coeffs {
			coefficients[i] = &bignum.Complex{
				new(big.Float).SetFloat64(c),
				new(big.Float),
			}
		}
	case []*bignum.Complex:
		coefficients = make([]*bignum.Complex, len(coeffs))
		copy(coefficients, coeffs)
	case []*big.Float:
		coefficients = make([]*bignum.Complex, len(coeffs))
		for i, c := range coeffs {
			coefficients[i] = &bignum.Complex{
				new(big.Float).Set(c),
				new(big.Float),
			}
		}
	default:
		panic(fmt.Sprintf("invalid coefficient type, allowed types are []{bignum.Complex128, float64, *bignum.Complex, *big.Float} but is %T", coeffs))
	}

	inter := Interval{}
	switch interval := interval.(type) {
	case [2]float64:
		inter.A = *new(big.Float).SetFloat64(interval[0])
		inter.B = *new(big.Float).SetFloat64(interval[1])
	case *Interval:
		inter.A = interval.A
		inter.B = interval.B
	case nil:

	default:
		panic(fmt.Sprintf("invalid interval type, allowed types are [2]float64 or *Interval, but is %T", interval))
	}

	return &Polynomial{
		MetaData: MetaData{
			Basis:    basis,
			Interval: inter,
			IsOdd:    true,
			IsEven:   true,
		},
		Coeffs: coefficients,
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
		num := new(big.Float).Sub(&p.B, &p.A)

		// 2 / (b-a)
		scalar = new(big.Float).Quo(new(big.Float).SetInt64(2), num)

		// (-b-a)/(b-a)
		constant = new(big.Float).Set(&p.B)
		constant.Neg(constant)
		constant.Sub(constant, &p.A)
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

// EvaluateModP evalutes the polynomial modulo p, treating each coefficient as
// integer variables and returning the result as *big.Int in the interval [0, P-1].
func (p *Polynomial) EvaluateModP(xInt, PInt *big.Int) (yInt *big.Int) {

	degree := p.Degree()

	yInt = p.Coeffs[degree].Int()

	for i := degree - 1; i >= 0; i-- {
		yInt.Mul(yInt, xInt)
		yInt.Mod(yInt, PInt)
		yInt.Add(yInt, p.Coeffs[i].Int())
	}

	if yInt.Cmp(new(big.Int)) == -1 {
		yInt.Add(yInt, PInt)
	}

	return
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
		y = coeffs[n-1].Clone()
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
			y = coeffs[0].Clone()
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

			TPrev = T.Clone()
			T = tmp.Clone()
		}

	default:
		panic(fmt.Sprintf("invalid basis type, allowed types are `Monomial` or `Chebyshev` but is %T", p.Basis))
	}

	return
}

// Factorize factorizes p as X^{n} * pq + pr.
func (p *Polynomial) Factorize(n int) (pq, pr *Polynomial) {

	if n < p.Degree()>>1 {
		panic("cannot Factorize: n < p.Degree()/2")
	}

	// ns a polynomial p such that p = q*C^degree + r.
	pr = &Polynomial{}
	pr.Coeffs = make([]*bignum.Complex, n)
	for i := 0; i < n; i++ {
		if p.Coeffs[i] != nil {
			pr.Coeffs[i] = p.Coeffs[i].Clone()
		}
	}

	pq = &Polynomial{}
	pq.Coeffs = make([]*bignum.Complex, p.Degree()-n+1)

	if p.Coeffs[n] != nil {
		pq.Coeffs[0] = p.Coeffs[n].Clone()
	}

	odd := p.IsOdd
	even := p.IsEven

	switch p.Basis {
	case Monomial:
		for i := n + 1; i < p.Degree()+1; i++ {
			if p.Coeffs[i] != nil && (!(even || odd) || (i&1 == 0 && even) || (i&1 == 1 && odd)) {
				pq.Coeffs[i-n] = p.Coeffs[i].Clone()
			}
		}
	case Chebyshev:

		for i, j := n+1, 1; i < p.Degree()+1; i, j = i+1, j+1 {
			if p.Coeffs[i] != nil && (!(even || odd) || (i&1 == 0 && even) || (i&1 == 1 && odd)) {
				pq.Coeffs[i-n] = p.Coeffs[i].Clone()
				pq.Coeffs[i-n].Add(pq.Coeffs[i-n], pq.Coeffs[i-n])

				if pr.Coeffs[n-j] != nil {
					pr.Coeffs[n-j].Sub(pr.Coeffs[n-j], p.Coeffs[i])
				} else {
					pr.Coeffs[n-j] = p.Coeffs[i].Clone()
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
