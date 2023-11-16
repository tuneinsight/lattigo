package bignum

import (
	"math/big"
)

// ChebyshevApproximation computes a Chebyshev approximation of the input function, for the range [-a, b] of degree degree.
// f.(type) can be either :
//   - func(Complex128)Complex128
//   - func(float64)float64
//   - func(*big.Float)*big.Float
//   - func(*Complex)*Complex
//
// The reference precision is taken from the values stored in the Interval struct.
func ChebyshevApproximation(f interface{}, interval Interval) (pol Polynomial) {

	var fCmplx func(*Complex) *Complex

	switch f := f.(type) {
	case func(x complex128) (y complex128):
		fCmplx = func(x *Complex) (y *Complex) {
			yCmplx := f(x.Complex128())
			return &Complex{new(big.Float).SetFloat64(real(yCmplx)), new(big.Float).SetFloat64(imag(yCmplx))}
		}
	case func(x float64) (y float64):
		fCmplx = func(x *Complex) (y *Complex) {
			xf64, _ := x[0].Float64()
			return &Complex{new(big.Float).SetFloat64(f(xf64)), new(big.Float)}
		}
	case func(x *big.Float) (y *big.Float):
		fCmplx = func(x *Complex) (y *Complex) {
			return &Complex{f(x[0]), new(big.Float)}
		}
	case func(x *Complex) *Complex:
		fCmplx = f
	}

	nodes := chebyshevNodes(interval.Nodes+1, interval)

	fi := make([]*Complex, len(nodes))

	x := NewComplex()
	x.SetPrec(interval.A.Prec())

	for i := range nodes {
		x[0].Set(nodes[i])
		fi[i] = fCmplx(x)
	}

	return NewPolynomial(Chebyshev, chebyCoeffs(nodes, fi, interval), &interval)
}

func chebyshevNodes(n int, interval Interval) (nodes []*big.Float) {

	prec := interval.A.Prec()

	nodes = make([]*big.Float, n)

	half := new(big.Float).SetPrec(prec).SetFloat64(0.5)

	x := new(big.Float).Add(&interval.A, &interval.B)
	x.Mul(x, half)
	y := new(big.Float).Sub(&interval.B, &interval.A)
	y.Mul(y, half)

	PiOverN := Pi(prec)
	PiOverN.Quo(PiOverN, new(big.Float).SetInt64(int64(n)))

	for k := 1; k < n+1; k++ {
		up := new(big.Float).SetPrec(prec).SetFloat64(float64(k) - 0.5)
		up.Mul(up, PiOverN)
		up = Cos(up)
		up.Mul(up, y)
		up.Add(up, x)
		nodes[n-k] = up
	}

	return
}

func chebyCoeffs(nodes []*big.Float, fi []*Complex, interval Interval) (coeffs []*Complex) {

	prec := interval.A.Prec()

	n := len(nodes)

	coeffs = make([]*Complex, n)
	for i := range coeffs {
		coeffs[i] = NewComplex().SetPrec(prec)
	}

	u := NewComplex().SetPrec(prec)

	mul := NewComplexMultiplier()

	tmp := NewComplex().SetPrec(prec)

	two := new(big.Float).SetPrec(prec).SetInt64(2)

	minusab := new(big.Float).Set(&interval.A)
	minusab.Neg(minusab)
	minusab.Sub(minusab, &interval.B)

	bminusa := new(big.Float).Set(&interval.B)
	bminusa.Sub(bminusa, &interval.A)

	Tnext := NewComplex().SetPrec(prec)

	for i := 0; i < n; i++ {

		u[0].Mul(nodes[i], two)
		u[0].Add(u[0], minusab)
		u[0].Quo(u[0], bminusa)

		Tprev := NewComplex().SetPrec(prec)
		Tprev[0].SetFloat64(1)

		T := u.Clone()

		for j := 0; j < n; j++ {

			mul.Mul(fi[i], Tprev, tmp)
			coeffs[j].Add(coeffs[j], tmp)

			mul.Mul(u, T, Tnext)
			Tnext[0].Mul(Tnext[0], two)
			Tnext[1].Mul(Tnext[1], two)
			Tnext.Sub(Tnext, Tprev)

			Tprev.Set(T)
			T.Set(Tnext)
		}
	}

	NHalf := new(big.Float).SetInt64(int64(n))

	coeffs[0][0].Quo(coeffs[0][0], NHalf)
	coeffs[0][1].Quo(coeffs[0][1], NHalf)

	NHalf.Quo(NHalf, two)

	for i := 1; i < n; i++ {
		coeffs[i][0].Quo(coeffs[i][0], NHalf)
		coeffs[i][1].Quo(coeffs[i][1], NHalf)
	}

	return
}

func chebyshevBasisInPlace(deg int, x *big.Float, inter Interval, poly []*big.Float) {

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
	poly[0].Set(Tprev)

	for i := 1; i < deg; i++ {
		Tnext.Mul(two, u)
		Tnext.Mul(Tnext, T)
		Tnext.Sub(Tnext, Tprev)
		Tprev.Set(T)
		T.Set(Tnext)
		poly[i].Set(Tprev)
	}
}
