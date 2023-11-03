package integer

import (
	"fmt"

	"github.com/tuneinsight/lattigo/v4/bfv"
	"github.com/tuneinsight/lattigo/v4/bgv"
	"github.com/tuneinsight/lattigo/v4/circuits"
	"github.com/tuneinsight/lattigo/v4/rlwe"
	"github.com/tuneinsight/lattigo/v4/utils/bignum"
)

// PolynomialEvaluator is a wrapper of the circuits.PolynomialEvaluator.
// All fields of this struct are public, enabling custom instantiations.
type PolynomialEvaluator struct {
	circuits.EvaluatorForPolynomial
	bgv.Parameters
	InvariantTensoring bool
}

// NewPowerBasis is a wrapper of NewPolynomialBasis.
// This function creates a new powerBasis from the input ciphertext.
// The input ciphertext is treated as the base monomial X used to
// generate the other powers X^{n}.
func NewPowerBasis(ct *rlwe.Ciphertext) circuits.PowerBasis {
	return circuits.NewPowerBasis(ct, bignum.Monomial)
}

// NewPolynomialEvaluator instantiates a new PolynomialEvaluator from a circuit.Evaluator.
// The default *bgv.Evaluator is compliant to the circuit.Evaluator interface.
// InvariantTensoring is a boolean that specifies if the evaluator performed the invariant tensoring (BFV-style) or
// the regular tensoring (BGB-style).
func NewPolynomialEvaluator(params bgv.Parameters, eval circuits.Evaluator, InvariantTensoring bool) *PolynomialEvaluator {

	var evalForPoly circuits.EvaluatorForPolynomial

	switch eval := eval.(type) {
	case *bgv.Evaluator:
		if InvariantTensoring {
			evalForPoly = &defaultCircuitEvaluatorForPolynomial{Evaluator: &scaleInvariantEvaluator{eval}}
		} else {
			evalForPoly = &defaultCircuitEvaluatorForPolynomial{Evaluator: eval}
		}
	case *bfv.Evaluator:
		if InvariantTensoring {
			evalForPoly = &defaultCircuitEvaluatorForPolynomial{Evaluator: eval}
		} else {
			evalForPoly = &defaultCircuitEvaluatorForPolynomial{Evaluator: eval.Evaluator}
		}
	default:
		evalForPoly = &defaultCircuitEvaluatorForPolynomial{Evaluator: eval}
	}

	return &PolynomialEvaluator{
		Parameters:             params,
		EvaluatorForPolynomial: evalForPoly,
		InvariantTensoring:     InvariantTensoring,
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

	var pcircuits interface{}
	switch p := p.(type) {
	case Polynomial:
		pcircuits = circuits.Polynomial(p)
	case PolynomialVector:
		pcircuits = circuits.PolynomialVector(p)
	default:
		pcircuits = p
	}

	return circuits.EvaluatePolynomial(eval.EvaluatorForPolynomial, ct, pcircuits, targetScale, 1, &simEvaluator{eval.Parameters, eval.InvariantTensoring})
}

// EvaluateFromPowerBasis evaluates a polynomial using the provided PowerBasis, holding pre-computed powers of X.
// This method is the same as Evaluate except that the encrypted input is a PowerBasis.
// See Evaluate for additional information.
func (eval PolynomialEvaluator) EvaluateFromPowerBasis(pb circuits.PowerBasis, p interface{}, targetScale rlwe.Scale) (opOut *rlwe.Ciphertext, err error) {

	var pcircuits interface{}
	switch p := p.(type) {
	case Polynomial:
		pcircuits = circuits.Polynomial(p)
	case PolynomialVector:
		pcircuits = circuits.PolynomialVector(p)
	default:
		pcircuits = p
	}

	if _, ok := pb.Value[1]; !ok {
		return nil, fmt.Errorf("cannot EvaluateFromPowerBasis: X^{1} is nil")
	}

	return circuits.EvaluatePolynomial(eval.EvaluatorForPolynomial, pb, pcircuits, targetScale, 1, &simEvaluator{eval.Parameters, eval.InvariantTensoring})
}

// scaleInvariantEvaluator is a struct implementing the interface circuits.Evaluator with
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
// circuits.CoefficientGetter[uint64] interface.
type CoefficientGetter struct {
	Values []uint64
}

// GetVectorCoefficient return a slice []uint64 containing the k-th coefficient
// of each polynomial of PolynomialVector indexed by its Mapping.
// See PolynomialVector for additional information about the Mapping.
func (c CoefficientGetter) GetVectorCoefficient(pol circuits.PolynomialVector, k int) (values []uint64) {

	values = c.Values

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
func (c CoefficientGetter) GetSingleCoefficient(pol circuits.Polynomial, k int) (value uint64) {
	return pol.Coeffs[k].Uint64()
}

// ShallowCopy returns a thread-safe copy of the original CoefficientGetter.
func (c CoefficientGetter) ShallowCopy() circuits.CoefficientGetter[uint64] {
	return &CoefficientGetter{Values: make([]uint64, len(c.Values))}
}

// defaultCircuitEvaluatorForPolynomial is a struct implementing the interface circuits.EvaluatorForPolynomial.
type defaultCircuitEvaluatorForPolynomial struct {
	circuits.Evaluator
}

// EvaluatePatersonStockmeyerPolynomialVector evaluates a pre-decomposed PatersonStockmeyerPolynomialVector on a pre-computed power basis [1, X^{1}, X^{2}, ..., X^{2^{n}}, X^{2^{n+1}}, ..., X^{2^{m}}]
func (eval defaultCircuitEvaluatorForPolynomial) EvaluatePatersonStockmeyerPolynomialVector(poly circuits.PatersonStockmeyerPolynomialVector, pb circuits.PowerBasis) (res *rlwe.Ciphertext, err error) {
	coeffGetter := circuits.CoefficientGetter[uint64](&CoefficientGetter{Values: make([]uint64, pb.Value[1].Slots())})
	return circuits.EvaluatePatersonStockmeyerPolynomialVector(eval, poly, coeffGetter, pb)
}
