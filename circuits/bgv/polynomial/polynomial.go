package polynomial

import (
	"github.com/tuneinsight/lattigo/v5/circuits/common/polynomial"
	"github.com/tuneinsight/lattigo/v5/schemes/bgv"
	"github.com/tuneinsight/lattigo/v5/utils/bignum"
)

// Polynomial is a type wrapping the type he.Polynomial.
type Polynomial polynomial.Polynomial

// NewPolynomial creates a new Polynomial from a list of coefficients []T.
func NewPolynomial[T bgv.Integer](coeffs []T) Polynomial {
	return Polynomial(polynomial.NewPolynomial(bignum.NewPolynomial(bignum.Monomial, coeffs, nil)))
}

// PolynomialVector is a type wrapping the type he.PolynomialVector.
type PolynomialVector polynomial.PolynomialVector

// Depth returns the depth of the target PolynomialVector.
func (p PolynomialVector) Depth() int {
	return p.Value[0].Depth()
}

func NewPolynomialVector[T bgv.Integer](polys [][]T, mapping map[int][]int) (PolynomialVector, error) {

	ps := make([]bignum.Polynomial, len(polys))

	for i := range ps {
		ps[i] = bignum.NewPolynomial(bignum.Monomial, polys[i], nil)
	}

	p, err := polynomial.NewPolynomialVector(ps, mapping)
	return PolynomialVector(p), err
}
