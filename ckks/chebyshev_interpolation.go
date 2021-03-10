package ckks

import (
	"math"
)

// ChebyshevInterpolation is a struct storing the coefficients, degree and range of a Chebyshev interpolation polynomial.
type ChebyshevInterpolation struct {
	Poly
	a complex128
	b complex128
}

// A returns the start of the approximation interval.
func (c *ChebyshevInterpolation) A() complex128 {
	return c.a
}

// B returns the end of the approximation interval.
func (c *ChebyshevInterpolation) B() complex128 {
	return c.b
}

// Approximate computes a Chebyshev approximation of the input function, for the range [-a, b] of degree degree.
// To be used in conjunction with the function EvaluateCheby.
func Approximate(function func(complex128) complex128, a, b complex128, degree int) (cheby *ChebyshevInterpolation) {

	cheby = new(ChebyshevInterpolation)
	cheby.a = a
	cheby.b = b
	cheby.maxDeg = degree
	cheby.lead = true

	nodes := chebyshevNodes(degree+1, a, b)

	fi := make([]complex128, len(nodes))
	for i := range nodes {
		fi[i] = function(nodes[i])
	}

	cheby.coeffs = chebyCoeffs(nodes, fi, a, b)

	return
}

func chebyshevNodes(n int, a, b complex128) (u []complex128) {
	u = make([]complex128, n)
	var x, y complex128
	for k := 1; k < n+1; k++ {
		x = 0.5 * (a + b)
		y = 0.5 * (b - a)
		u[k-1] = x + y*complex(math.Cos((float64(k)-0.5)*(3.141592653589793/float64(n))), 0)
	}
	return
}

func chebyCoeffs(nodes, fi []complex128, a, b complex128) (coeffs []complex128) {

	var u, Tprev, T, Tnext complex128

	n := len(nodes)

	coeffs = make([]complex128, n)

	for i := 0; i < n; i++ {

		u = (2*nodes[i] - a - b) / (b - a)
		Tprev = 1
		T = u

		for j := 0; j < n; j++ {
			coeffs[j] += fi[i] * Tprev
			Tnext = 2*u*T - Tprev
			Tprev = T
			T = Tnext
		}
	}

	coeffs[0] /= complex(float64(n), 0)
	for i := 1; i < n; i++ {
		coeffs[i] *= (2.0 / complex(float64(n), 0))
	}

	return
}
