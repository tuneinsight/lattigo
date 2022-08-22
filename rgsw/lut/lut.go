// Package lut implements look-up tables evaluation for R-LWE schemes.
package lut

import (
	"github.com/tuneinsight/lattigo/v3/ring"
)

// InitLUT takes a function g, and creates a LUT polynomial for the function in the interval [a, b].
// Inputs to the LUT evaluation are assumed to have been normalized with the change of basis (2*x - a - b)/(b-a).
// Interval [a, b] should take into account the "drift" of the value x, caused by the change of modulus from Q to 2N.
func InitLUT(g func(x float64) (y float64), scale float64, ringQ *ring.Ring, a, b float64) (F *ring.Poly) {
	F = ringQ.NewPoly()
	Q := ringQ.Modulus

	// Discretization interval
	interval := 2.0 / float64(ringQ.N)

	for j, qi := range Q {

		// Interval [-1, 0] of g(x)
		for i := 0; i < (ringQ.N>>1)+1; i++ {
			F.Coeffs[j][i] = scaleUp(g(normalizeInv(-interval*float64(i), a, b)), scale, qi)
		}

		// Interval ]0, 1[ of g(x)
		for i := (ringQ.N >> 1) + 1; i < ringQ.N; i++ {
			F.Coeffs[j][i] = scaleUp(-g(normalizeInv(interval*float64(ringQ.N-i), a, b)), scale, qi)
		}
	}

	ringQ.NTT(F, F)

	return
}
