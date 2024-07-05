package polyint

import (
	"fmt"

	"github.com/tuneinsight/lattigo/v5/circuits/polynomial"
	"github.com/tuneinsight/lattigo/v5/core/rlwe"
	"github.com/tuneinsight/lattigo/v5/schemes"
	"github.com/tuneinsight/lattigo/v5/schemes/bfv"
	"github.com/tuneinsight/lattigo/v5/schemes/bgv"
)

type PolynomialEvaluator struct {
	bgv.Parameters
	polynomial.PolynomialEvaluator[uint64]
	InvariantTensoring bool
}

// NewPolynomialEvaluator instantiates a new PolynomialEvaluator from a circuit.Evaluator.
// The default heint.Evaluator is compliant to the circuit.Evaluator interface.
// InvariantTensoring is a boolean that specifies if the evaluator performed the invariant tensoring (BFV-style) or
// the regular tensoring (BGB-style).
func NewPolynomialEvaluator(params bgv.Parameters, eval schemes.Evaluator, InvariantTensoring bool) *PolynomialEvaluator {

	var evalForPoly schemes.Evaluator

	switch eval := eval.(type) {
	case *bgv.Evaluator:
		if InvariantTensoring {
			evalForPoly = scaleInvariantEvaluator{eval}
		} else {
			evalForPoly = eval
		}
	case *bfv.Evaluator:
		if InvariantTensoring {
			evalForPoly = eval
		} else {
			evalForPoly = eval.Evaluator
		}
	default:
		evalForPoly = eval
	}

	return &PolynomialEvaluator{
		Parameters: params,
		PolynomialEvaluator: polynomial.PolynomialEvaluator[uint64]{
			Evaluator:         evalForPoly,
			CoefficientGetter: CoefficientGetter{values: make([]uint64, params.MaxSlots())},
		},
		InvariantTensoring: InvariantTensoring,
	}
}

// Evaluate evaluates a polynomial on the input Ciphertext in ceil(log2(deg+1)) levels.
// Returns an error if the input ciphertext does not have enough level to carry out the full polynomial evaluation.
// Returns an error if something is wrong with the scale.
// If the polynomial is given in Chebyshev basis, then a change of basis ct' = (2/(b-a)) * (ct + (-a-b)/(b-a))
// is necessary before the polynomial evaluation to ensure correctness.
// pol: a *bignum.Polynomial, *Polynomial or *PolynomialVector
// targetScale: the desired output scale. This value shouldn't differ too much from the original ciphertext scale. It can
// for example be used to correct small deviations in the ciphertext scale and reset it to the default scale.
func (eval PolynomialEvaluator) Evaluate(ct *rlwe.Ciphertext, p interface{}, targetScale rlwe.Scale) (opOut *rlwe.Ciphertext, err error) {

	var phe interface{}
	switch p := p.(type) {
	case Polynomial:
		phe = polynomial.Polynomial(p)
	case PolynomialVector:
		phe = polynomial.PolynomialVector(p)
	default:
		phe = p
	}

	return eval.EvaluatePolynomial(ct, phe, targetScale, 1, &simEvaluator{eval.Parameters, eval.InvariantTensoring})
}

// EvaluateFromPowerBasis evaluates a polynomial using the provided PowerBasis, holding pre-computed powers of X.
// This method is the same as Evaluate except that the encrypted input is a PowerBasis.
// See Evaluate for additional information.
func (eval PolynomialEvaluator) EvaluateFromPowerBasis(pb polynomial.PowerBasis, p interface{}, targetScale rlwe.Scale) (opOut *rlwe.Ciphertext, err error) {

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

	return eval.EvaluatePolynomial(pb, phe, targetScale, 1, &simEvaluator{eval.Parameters, eval.InvariantTensoring})
}

// scaleInvariantEvaluator is a struct implementing the interface he.Evaluator with
// scale invariant tensoring (BFV-style).
type scaleInvariantEvaluator struct {
	*bgv.Evaluator
}

func (polyEval scaleInvariantEvaluator) Mul(op0 *rlwe.Ciphertext, op1 rlwe.Operand, opOut *rlwe.Ciphertext) (err error) {
	return polyEval.MulScaleInvariant(op0, op1, opOut)
}

func (polyEval scaleInvariantEvaluator) MulRelin(op0 *rlwe.Ciphertext, op1 rlwe.Operand, opOut *rlwe.Ciphertext) (err error) {
	return polyEval.Evaluator.MulRelinScaleInvariant(op0, op1, opOut)
}

func (polyEval scaleInvariantEvaluator) MulNew(op0 *rlwe.Ciphertext, op1 rlwe.Operand) (opOut *rlwe.Ciphertext, err error) {
	return polyEval.Evaluator.MulScaleInvariantNew(op0, op1)
}

func (polyEval scaleInvariantEvaluator) MulRelinNew(op0 *rlwe.Ciphertext, op1 rlwe.Operand) (opOut *rlwe.Ciphertext, err error) {
	return polyEval.Evaluator.MulRelinScaleInvariantNew(op0, op1)
}

func (polyEval scaleInvariantEvaluator) Rescale(op0, op1 *rlwe.Ciphertext) (err error) {
	return nil
}

// CoefficientGetter is a struct that implements the
// he.CoefficientGetter[uint64] interface.
type CoefficientGetter struct {
	values []uint64
}

// GetVectorCoefficient return a slice []uint64 containing the k-th coefficient
// of each polynomial of PolynomialVector indexed by its Mapping.
// See PolynomialVector for additional information about the Mapping.
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
