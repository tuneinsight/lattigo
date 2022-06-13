package ckks

import (
	"math"
)

// Approximate computes a Chebyshev approximation of the input function, for the range [-a, b] of degree degree.
// function.(type) can be either func(complex128)complex128 or func(float64)float64
// To be used in conjunction with the function EvaluateCheby.
func Approximate(function interface{}, a, b float64, degree int) (pol *Polynomial) {

	nodes := chebyshevNodes(degree+1, a, b)

	fi := make([]complex128, len(nodes))

	switch f := function.(type) {
	case func(complex128) complex128:
		for i := range nodes {
			fi[i] = f(complex(nodes[i], 0))
		}
	case func(float64) float64:
		for i := range nodes {
			fi[i] = complex(f(nodes[i]), 0)
		}
	default:
		panic("function must be either func(complex128)complex128 or func(float64)float64")
	}

	pol = NewPoly(chebyCoeffs(nodes, fi, a, b))
	pol.A = a
	pol.B = b
	pol.MaxDeg = degree
	pol.Lead = true
	pol.BasisType = Chebyshev

	return
}

func chebyshevNodes(n int, a, b float64) (u []float64) {
	u = make([]float64, n)
	x, y := 0.5*(a+b), 0.5*(b-a)
	for k := 1; k < n+1; k++ {
		u[k-1] = x + y*math.Cos((float64(k)-0.5)*(3.141592653589793/float64(n)))
	}
	return
}

func chebyCoeffs(nodes []float64, fi []complex128, a, b float64) (coeffs []complex128) {

	var u, Tprev, T, Tnext complex128

	n := len(nodes)

	coeffs = make([]complex128, n)

	for i := 0; i < n; i++ {

		u = complex((2*nodes[i]-a-b)/(b-a), 0)
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
