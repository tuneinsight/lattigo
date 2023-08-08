package float

import (
	"github.com/tuneinsight/lattigo/v4/circuits"
	"github.com/tuneinsight/lattigo/v4/utils/bignum"
)

type Polynomial circuits.Polynomial

func NewPolynomial(poly bignum.Polynomial) Polynomial {
	return Polynomial(circuits.NewPolynomial(poly))
}

type PolynomialVector circuits.PolynomialVector

func NewPolynomialVector(polys []bignum.Polynomial, mapping map[int][]int) (PolynomialVector, error) {
	p, err := circuits.NewPolynomialVector(polys, mapping)
	return PolynomialVector(p), err
}
