package heint

import (
	"github.com/tuneinsight/lattigo/v5/he"
	"github.com/tuneinsight/lattigo/v5/utils/bignum"
)

// Polynomial is a type wrapping the type he.Polynomial.
type Polynomial he.Polynomial

// NewPolynomial creates a new Polynomial from a list of coefficients []T.
func NewPolynomial[T Integer](coeffs []T) Polynomial {
	return Polynomial(he.NewPolynomial(bignum.NewPolynomial(bignum.Monomial, coeffs, nil)))
}

// PolynomialVector is a type wrapping the type he.PolynomialVector.
type PolynomialVector he.PolynomialVector

// Depth returns the depth of the target PolynomialVector.
func (p PolynomialVector) Depth() int {
	return p.Value[0].Depth()
}

func NewPolynomialVector[T Integer](polys [][]T, mapping map[int][]int) (PolynomialVector, error) {

	ps := make([]bignum.Polynomial, len(polys))

	for i := range ps {
		ps[i] = bignum.NewPolynomial(bignum.Monomial, polys[i], nil)
	}

	p, err := he.NewPolynomialVector(ps, mapping)
	return PolynomialVector(p), err
}
