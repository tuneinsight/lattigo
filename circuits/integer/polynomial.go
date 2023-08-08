package integer

import (
	"github.com/tuneinsight/lattigo/v4/circuits"
	"github.com/tuneinsight/lattigo/v4/utils/bignum"
)

type Polynomial circuits.Polynomial

func NewPolynomial[T Integer](coeffs []T) Polynomial {
	return Polynomial(circuits.NewPolynomial(bignum.NewPolynomial(bignum.Monomial, coeffs, nil)))
}

type PolynomialVector circuits.PolynomialVector

func NewPolynomialVector[T Integer](polys [][]T, mapping map[int][]int) (PolynomialVector, error) {

	ps := make([]bignum.Polynomial, len(polys))

	for i := range ps {
		ps[i] = bignum.NewPolynomial(bignum.Monomial, polys[i], nil)
	}

	p, err := circuits.NewPolynomialVector(ps, mapping)
	return PolynomialVector(p), err
}
