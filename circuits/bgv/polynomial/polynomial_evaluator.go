package polynomial

import (
	"fmt"

	"github.com/tuneinsight/lattigo/v6/circuits/common/polynomial"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/schemes/bgv"
)

type Evaluator struct {
	bgv.Parameters
	polynomial.Evaluator[uint64]
	InvariantTensoring bool
}

// NewEvaluator instantiates a new PolynomialEvaluator from a [bgv.Evaluator].
func NewEvaluator(params bgv.Parameters, eval *bgv.Evaluator) *Evaluator {

	return &Evaluator{
		Parameters: params,
		Evaluator: polynomial.Evaluator[uint64]{
			Evaluator:         eval,
			CoefficientGetter: CoefficientGetter{values: make([]uint64, params.MaxSlots())},
		},
		InvariantTensoring: eval.ScaleInvariant,
	}
}

// Evaluate evaluates a polynomial on the input Ciphertext in ceil(log2(deg+1)) levels.
// Returns an error if the input ciphertext does not have enough level to carry out the full polynomial evaluation.
// Returns an error if something is wrong with the scale.
// If the polynomial is given in Chebyshev basis, then a change of basis ct' = (2/(b-a)) * (ct + (-a-b)/(b-a))
// is necessary before the polynomial evaluation to ensure correctness.
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

	return eval.Evaluator.Evaluate(ct, phe, targetScale, 1, &simEvaluator{eval.Parameters, eval.InvariantTensoring})
}

// EvaluateFromPowerBasis evaluates a polynomial using the provided [polynomial.PowerBasis], holding pre-computed powers of X.
// This method is the same as Evaluate except that the encrypted input is a PowerBasis.
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

	if _, ok := pb.Value[1]; !ok {
		return nil, fmt.Errorf("cannot EvaluateFromPowerBasis: X^{1} is nil")
	}

	return eval.Evaluator.Evaluate(pb, phe, targetScale, 1, &simEvaluator{eval.Parameters, eval.InvariantTensoring})
}

// CoefficientGetter is a struct that implements the
// [polynomial.CoefficientGetter][uint64] interface.
type CoefficientGetter struct {
	values []uint64
}

// GetVectorCoefficient return a slice []uint64 containing the k-th coefficient
// of each polynomial of PolynomialVector indexed by its Mapping.
// See [PolynomialVector] for additional information about the Mapping.
func (c CoefficientGetter) GetVectorCoefficient(pol polynomial.PolynomialVector, k int) (values []uint64) {

	values = c.values

	for j := range values {
		values[j] = 0
	}

	mapping := pol.Mapping

	for i, p := range pol.Value {
		for _, j := range mapping[i] {
			values[j] = p.Coeffs[k].Uint64()
		}
	}

	return
}

// GetSingleCoefficient should return the k-th coefficient of Polynomial as the type uint64.
func (c CoefficientGetter) GetSingleCoefficient(pol polynomial.Polynomial, k int) (value uint64) {
	return pol.Coeffs[k].Uint64()
}
