package hefloat

import (
	"math/big"

	"github.com/tuneinsight/lattigo/v5/he"
	"github.com/tuneinsight/lattigo/v5/utils/bignum"
)

// Polynomial is a type wrapping the type he.Polynomial.
type Polynomial he.Polynomial

// NewPolynomial creates a new Polynomial from a bignum.Polynomial.
func NewPolynomial(poly bignum.Polynomial) Polynomial {
	return Polynomial(he.NewPolynomial(poly))
}

// PolynomialVector is a type wrapping the type he.PolynomialVector.
type PolynomialVector he.PolynomialVector

// Depth returns the depth of the target PolynomialVector.
func (p PolynomialVector) Depth() int {
	return p.Value[0].Depth()
}

// NewPolynomialVector creates a new PolynomialVector from a list of bignum.Polynomial and a mapping
// map[poly_index][slots_index] which stores which polynomial has to be evaluated on which slot.
// Slots that are not referenced in this mapping will be evaluated to zero.
// User must ensure that a same slot is not referenced twice.
func NewPolynomialVector(polys []bignum.Polynomial, mapping map[int][]int) (PolynomialVector, error) {
	p, err := he.NewPolynomialVector(polys, mapping)
	return PolynomialVector(p), err
}

func (p PolynomialVector) ChangeOfBasis(slots int) (scalar, constant []*big.Float) {

	scalar = make([]*big.Float, slots)
	constant = make([]*big.Float, slots)

	for i := 0; i < slots; i++ {
		scalar[i] = new(big.Float)
		constant[i] = new(big.Float)
	}

	for i := range p.Mapping {
		m := p.Mapping[i]
		s, c := p.Value[i].ChangeOfBasis()
		for _, j := range m {
			scalar[j] = s
			constant[j] = c
		}
	}

	return
}
