package bignum

import (
	"math/big"
)

// MonomialEval evaluates y = sum x^i * poly[i].
func MonomialEval(x *big.Float, poly []*big.Float) (y *big.Float) {
	n := len(poly) - 1
	y = new(big.Float).Set(poly[n-1])
	for i := n - 2; i >= 0; i-- {
		y.Mul(y, x)
		y.Add(y, poly[i])
	}
	return
}

// ChebyshevEval evaluates y = sum Ti(x) * poly[i], where T0(x) = 1, T1(x) = (2x-a-b)/(b-a) and T{i+j}(x) = 2TiTj(x)- T|i-j|(x).
func ChebyshevEval(x *big.Float, poly []*big.Float, inter Interval) (y *big.Float) {

	precision := x.Prec()

	two := NewFloat(2, precision)
	var tmp, u = new(big.Float), new(big.Float)
	var T, Tprev, Tnext = new(big.Float), new(big.Float), new(big.Float)

	// u = (2*x - (a+b))/(b-a)
	u.Set(x)
	u.Mul(u, two)
	u.Sub(u, &inter.A)
	u.Sub(u, &inter.B)
	tmp.Set(&inter.B)
	tmp.Sub(tmp, &inter.A)
	u.Quo(u, tmp)

	Tprev.SetPrec(precision)
	Tprev.SetFloat64(1)
	T.Set(u)
	y = new(big.Float).Set(poly[0])

	for i := 1; i < len(poly); i++ {
		y.Add(y, tmp.Mul(T, poly[i]))
		Tnext.Mul(two, u)
		Tnext.Mul(Tnext, T)
		Tnext.Sub(Tnext, Tprev)
		Tprev.Set(T)
		T.Set(Tnext)
	}

	return
}
