package bgv

import (
	"github.com/tuneinsight/lattigo/v4/hebase"
	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/rlwe"
	"github.com/tuneinsight/lattigo/v4/rlwe/ringqp"
	"github.com/tuneinsight/lattigo/v4/utils/bignum"
)

// NewPowerBasis is a wrapper of hebase.NewPolynomialBasis.
// This function creates a new powerBasis from the input ciphertext.
// The input ciphertext is treated as the base monomial X used to
// generate the other powers X^{n}.
func NewPowerBasis(ct *rlwe.Ciphertext) hebase.PowerBasis {
	return hebase.NewPowerBasis(ct, bignum.Monomial)
}

// NewPolynomial is a wrapper of hebase.NewPolynomial.
// This function creates a new polynomial from the input coefficients.
// This polynomial can be evaluated on a ciphertext.
func NewPolynomial[T int64 | uint64](coeffs []T) hebase.Polynomial {
	return hebase.NewPolynomial(bignum.NewPolynomial(bignum.Monomial, coeffs, nil))
}

// NewPolynomialVector is a wrapper of hebase.NewPolynomialVector.
// This function creates a new PolynomialVector from the input polynomials and the desired function mapping.
// This polynomial vector can be evaluated on a ciphertext.
func NewPolynomialVector(polys []hebase.Polynomial, mapping map[int][]int) (hebase.PolynomialVector, error) {
	return hebase.NewPolynomialVector(polys, mapping)
}

// LinearTransformationParametersLiteral is a struct defining the parameterization of a linear transformation.
// See hebase.LinearTranfromationParameters for additional informations about each fields.
type LinearTransformationParametersLiteral[T int64 | uint64] struct {
	Diagonals                map[int][]T
	Level                    int
	Scale                    rlwe.Scale
	LogDimensions            ring.Dimensions
	LogBabyStepGianStepRatio int
}

// NewLinearTransformationParmeters creates a new hebase.LinearTransformationParameters from the provided LinearTransformationParametersLiteral.
func NewLinearTransformationParmeters[T int64 | uint64](params LinearTransformationParametersLiteral[T]) hebase.LinearTranfromationParameters[T] {
	return hebase.MemLinearTransformationParameters[T]{
		Diagonals:                params.Diagonals,
		Level:                    params.Level,
		Scale:                    params.Scale,
		LogDimensions:            params.LogDimensions,
		LogBabyStepGianStepRatio: params.LogBabyStepGianStepRatio,
	}
}

// GaloisElementsForLinearTransformation returns the list of Galois elements required to evaluate the linear transformation.
func GaloisElementsForLinearTransformation[T int64 | uint64](params rlwe.ParametersInterface, lt hebase.LinearTranfromationParameters[T]) (galEls []uint64) {
	return hebase.GaloisElementsForLinearTransformation(params, lt)
}

// NewLinearTransformation allocates a new LinearTransformation with zero values and according to the provided parameters.
func NewLinearTransformation[T int64 | uint64](params rlwe.ParametersInterface, lt hebase.LinearTranfromationParameters[T]) hebase.LinearTransformation {
	return hebase.NewLinearTransformation(params, lt)
}

// EncodeLinearTransformation encodes a linear transformation on a pre-allocated linear transformation.
// The method will return an error if the non-zero diagonals between the pre-allocated linear transformation and the parameters of the linear transformation to encode do not match.
func EncodeLinearTransformation[T int64 | uint64](allocated hebase.LinearTransformation, params hebase.LinearTranfromationParameters[T], ecd *Encoder) (err error) {
	return hebase.EncodeLinearTransformation[T](allocated, params, &encoder[T, ringqp.Poly]{ecd})
}
