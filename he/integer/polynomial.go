package integer

import (
	"github.com/tuneinsight/lattigo/v4/circuits"
	"github.com/tuneinsight/lattigo/v4/utils/bignum"
)

// Polynomial is a type wrapping the type circuits.Polynomial.
type Polynomial circuits.Polynomial

// NewPolynomial creates a new Polynomial from a list of coefficients []T.
func NewPolynomial[T Integer](coeffs []T) Polynomial {
	return Polynomial(circuits.NewPolynomial(bignum.NewPolynomial(bignum.Monomial, coeffs, nil)))
}

// PolynomialVector is a type wrapping the type circuits.PolynomialVector.
type PolynomialVector circuits.PolynomialVector

// Depth returns the depth of the target PolynomialVector.
func (p PolynomialVector) Depth() int {
	return p.Value[0].Depth()
}

func NewPolynomialVector[T Integer](polys [][]T, mapping map[int][]int) (PolynomialVector, error) {

	ps := make([]bignum.Polynomial, len(polys))

	for i := range ps {
		ps[i] = bignum.NewPolynomial(bignum.Monomial, polys[i], nil)
	}

	p, err := circuits.NewPolynomialVector(ps, mapping)
	return PolynomialVector(p), err
}
