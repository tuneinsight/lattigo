package ckks

import (
	"github.com/tuneinsight/lattigo/v4/hebase"
	"github.com/tuneinsight/lattigo/v4/rlwe"
	"github.com/tuneinsight/lattigo/v4/utils/bignum"
)

// NewPowerBasis is a wrapper of hebase.NewPolynomialBasis.
// This function creates a new powerBasis from the input ciphertext.
// The input ciphertext is treated as the base monomial X used to
// generate the other powers X^{n}.
func NewPowerBasis(ct *rlwe.Ciphertext, basis bignum.Basis) hebase.PowerBasis {
	return hebase.NewPowerBasis(ct, basis)
}

// NewPolynomial is a wrapper of hebase.NewPolynomial.
// This function creates a new polynomial from the input coefficients.
// This polynomial can be evaluated on a ciphertext.
func NewPolynomial(poly bignum.Polynomial) hebase.Polynomial {
	return hebase.NewPolynomial(poly)
}

// NewPolynomialVector is a wrapper of hebase.NewPolynomialVector.
// This function creates a new PolynomialVector from the input polynomials and the desired function mapping.
// This polynomial vector can be evaluated on a ciphertext.
func NewPolynomialVector(polys []hebase.Polynomial, mapping map[int][]int) (hebase.PolynomialVector, error) {
	return hebase.NewPolynomialVector(polys, mapping)
}

type PolynomialEvaluator struct {
	*Evaluator
}

func NewPolynomialEvaluator(eval *Evaluator) *PolynomialEvaluator {
	return &PolynomialEvaluator{eval}
}
