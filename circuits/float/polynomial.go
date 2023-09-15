package float

import (
	"github.com/tuneinsight/lattigo/v4/circuits"
	"github.com/tuneinsight/lattigo/v4/utils/bignum"
)

// Polynomial is a type wrapping the type circuits.Polynomial.
type Polynomial circuits.Polynomial

// NewPolynomial creates a new Polynomial from a bignum.Polynomial.
func NewPolynomial(poly bignum.Polynomial) Polynomial {
	return Polynomial(circuits.NewPolynomial(poly))
}

// PolynomialVector is a type wrapping the type circuits.PolynomialVector.
type PolynomialVector circuits.PolynomialVector

// Depth returns the depth of the target PolynomialVector.
func (p PolynomialVector) Depth() int {
	return p.Value[0].Depth()
}

// NewPolynomialVector creates a new PolynomialVector from a list of bignum.Polynomial and a mapping
// map[poly_index][slots_index] which stores which polynomial has to be evaluated on which slot.
// Slots that are not referenced in this mapping will be evaluated to zero.
// User must ensure that a same slot is not referenced twice.
func NewPolynomialVector(polys []bignum.Polynomial, mapping map[int][]int) (PolynomialVector, error) {
	p, err := circuits.NewPolynomialVector(polys, mapping)
	return PolynomialVector(p), err
}
