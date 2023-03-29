// Package polynomial provides helper for polynomials, approximation of functions using polynomials and their evaluation.
package polynomial

// Basis is a type for the polynomials basis
type Basis int

const (
	// Monomial : x^(a+b) = x^a * x^b
	Monomial = Basis(0)
	// Chebyshev : T_(a+b) = 2 * T_a * T_b - T_(|a-b|)
	Chebyshev = Basis(1)
)
