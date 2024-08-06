package polynomial

import (
	"fmt"

	"github.com/tuneinsight/lattigo/v6/circuits/common/polynomial"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
	"github.com/tuneinsight/lattigo/v6/utils/bignum"
)

// Evaluator is a wrapper of the [polynomial.Evaluator].
// All fields of this struct are public, enabling custom instantiations.
type Evaluator struct {
	ckks.Parameters
	polynomial.Evaluator[*bignum.Complex]
}

// NewEvaluator instantiates a new [Evaluator] from a [ckks.Evaluator].
// This method is allocation free.
func NewEvaluator(params ckks.Parameters, eval *ckks.Evaluator) *Evaluator {
	return &Evaluator{
		Parameters: params,
		Evaluator: polynomial.Evaluator[*bignum.Complex]{
			Evaluator:         eval,
			CoefficientGetter: CoefficientGetter{values: make([]*bignum.Complex, params.MaxSlots())},
		},
	}
}

// Evaluate evaluates a polynomial on the input Ciphertext in ceil(log2(deg+1)) levels.
// Returns an error if the input ciphertext does not have enough levels to carry out the full polynomial evaluation.
// Returns an error if something is wrong with the scale.
//
// If the polynomial is given in Chebyshev basis, then the user must apply change of basis
// ct' = scale * ct + offset before the polynomial evaluation to ensure correctness.
// The values `scale` and `offet` can be obtained from the polynomial with the method .ChangeOfBasis().
//
// pol: a *[bignum.Polynomial], *[Polynomial] or *[PolynomialVector]
// targetScale: the desired output scale. This value shouldn't differ too much from the original ciphertext scale. It can
// for example be used to correct small deviations in the ciphertext scale and reset it to the default scale.
func (eval Evaluator) Evaluate(ct *rlwe.Ciphertext, p interface{}, targetScale rlwe.Scale) (opOut *rlwe.Ciphertext, err error) {

	var phe interface{}
	switch p := p.(type) {
	case Polynomial:
		phe = polynomial.Polynomial(p)
	case PolynomialVector:
		phe = polynomial.PolynomialVector(p)
	default:
		phe = p
	}

	levelsConsumedPerRescaling := eval.Parameters.LevelsConsumedPerRescaling()

	return eval.Evaluator.Evaluate(ct, phe, targetScale, levelsConsumedPerRescaling, &simEvaluator{eval.Parameters, levelsConsumedPerRescaling})
}

// EvaluateFromPowerBasis evaluates a polynomial using the provided [polynomial.PowerBasis], holding pre-computed powers of X.
// This method is the same as [Evaluate] except that the encrypted input is a [polynomial.PowerBasis].
// See [Evaluate] for additional information.
func (eval Evaluator) EvaluateFromPowerBasis(pb polynomial.PowerBasis, p interface{}, targetScale rlwe.Scale) (opOut *rlwe.Ciphertext, err error) {

	var phe interface{}
	switch p := p.(type) {
	case Polynomial:
		phe = polynomial.Polynomial(p)
	case PolynomialVector:
		phe = polynomial.PolynomialVector(p)
	default:
		phe = p
	}

	levelsConsumedPerRescaling := eval.Parameters.LevelsConsumedPerRescaling()

	if _, ok := pb.Value[1]; !ok {
		return nil, fmt.Errorf("cannot EvaluateFromPowerBasis: X^{1} is nil")
	}

	return eval.Evaluator.Evaluate(pb, phe, targetScale, levelsConsumedPerRescaling, &simEvaluator{eval.Parameters, levelsConsumedPerRescaling})
}

// CoefficientGetter is a struct that implements the
// [polynomial.CoefficientGetter][*bignum.Complex] interface.
type CoefficientGetter struct {
	values []*bignum.Complex
}

// GetVectorCoefficient return a slice []*[bignum.Complex] containing the k-th coefficient
// of each polynomial of [polynomial.PolynomialVector] indexed by its Mapping.
// See [polynomial.PolynomialVector] for additional information about the Mapping.
func (c CoefficientGetter) GetVectorCoefficient(pol polynomial.PolynomialVector, k int) (values []*bignum.Complex) {

	values = c.values

	for j := range values {
		values[j] = nil
	}

	mapping := pol.Mapping

	for i, p := range pol.Value {
		for _, j := range mapping[i] {
			values[j] = p.Coeffs[k]
		}
	}

	return
}

// GetSingleCoefficient returns the k-th coefficient of Polynomial as the type *[bignum.Complex].
func (c CoefficientGetter) GetSingleCoefficient(pol polynomial.Polynomial, k int) (value *bignum.Complex) {
	return pol.Coeffs[k]
}
