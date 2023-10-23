package float

import (
	"fmt"

	"github.com/tuneinsight/lattigo/v4/circuits"
	"github.com/tuneinsight/lattigo/v4/ckks"
	"github.com/tuneinsight/lattigo/v4/rlwe"
	"github.com/tuneinsight/lattigo/v4/utils/bignum"
)

// PolynomialEvaluator is a wrapper of the circuits.PolynomialEvaluator.
// All fields of this struct are public, enabling custom instantiations.
type PolynomialEvaluator struct {
	Parameters ckks.Parameters
	circuits.EvaluatorForPolynomial
}

// NewPowerBasis is a wrapper of NewPolynomialBasis.
// This function creates a new powerBasis from the input ciphertext.
// The input ciphertext is treated as the base monomial X used to
// generate the other powers X^{n}.
func NewPowerBasis(ct *rlwe.Ciphertext, basis bignum.Basis) circuits.PowerBasis {
	return circuits.NewPowerBasis(ct, basis)
}

// NewPolynomialEvaluator instantiates a new PolynomialEvaluator from a circuit.Evaluator.
// The default *ckks.Evaluator is compliant to the circuit.Evaluator interface.
// This method is allocation free.
func NewPolynomialEvaluator(params ckks.Parameters, eval circuits.Evaluator) *PolynomialEvaluator {
	return &PolynomialEvaluator{
		Parameters:             params,
		EvaluatorForPolynomial: &defaultCircuitEvaluatorForPolynomial{Evaluator: eval},
	}
}

// Evaluate evaluates a polynomial on the input Ciphertext in ceil(log2(deg+1)) levels.
// Returns an error if the input ciphertext does not have enough levels to carry out the full polynomial evaluation.
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

	levelsConsumedPerRescaling := eval.Parameters.LevelsConsumedPerRescaling()

	return circuits.EvaluatePolynomial(eval, ct, pcircuits, targetScale, levelsConsumedPerRescaling, &simEvaluator{eval.Parameters, levelsConsumedPerRescaling})
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

	levelsConsumedPerRescaling := eval.Parameters.LevelsConsumedPerRescaling()

	if _, ok := pb.Value[1]; !ok {
		return nil, fmt.Errorf("cannot EvaluateFromPowerBasis: X^{1} is nil")
	}

	return circuits.EvaluatePolynomial(eval, pb, pcircuits, targetScale, levelsConsumedPerRescaling, &simEvaluator{eval.Parameters, levelsConsumedPerRescaling})
}

type CoefficientGetter struct {
	Values []*bignum.Complex
}

func (c CoefficientGetter) Clone() *CoefficientGetter {
	return &CoefficientGetter{Values: make([]*bignum.Complex, len(c.Values))}
}

func (c *CoefficientGetter) GetVectorCoefficient(pol circuits.PolynomialVector, k int) (values []*bignum.Complex) {

	values = c.Values

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

func (c *CoefficientGetter) GetSingleCoefficient(pol circuits.Polynomial, k int) (value *bignum.Complex) {
	return pol.Coeffs[k]
}

func (c CoefficientGetter) ShallowCopy() circuits.CoefficientGetter[*bignum.Complex] {
	return &CoefficientGetter{Values: make([]*bignum.Complex, len(c.Values))}
}

type defaultCircuitEvaluatorForPolynomial struct {
	circuits.Evaluator
}

func (eval defaultCircuitEvaluatorForPolynomial) EvaluatePatersonStockmeyerPolynomialVector(poly circuits.PatersonStockmeyerPolynomialVector, pb circuits.PowerBasis) (res *rlwe.Ciphertext, err error) {
	coeffGetter := circuits.CoefficientGetter[*bignum.Complex](&CoefficientGetter{Values: make([]*bignum.Complex, pb.Value[1].Slots())})
	return circuits.EvaluatePatersonStockmeyerPolynomialVector(eval, poly, coeffGetter, pb)
}
