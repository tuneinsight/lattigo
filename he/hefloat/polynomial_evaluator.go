package hefloat

import (
	"fmt"

	"github.com/tuneinsight/lattigo/v5/core/rlwe"
	"github.com/tuneinsight/lattigo/v5/he"
	"github.com/tuneinsight/lattigo/v5/utils/bignum"
)

// PolynomialEvaluator is a wrapper of the he.PolynomialEvaluator.
// All fields of this struct are public, enabling custom instantiations.
type PolynomialEvaluator struct {
	Parameters Parameters
	he.EvaluatorForPolynomial
}

// NewPowerBasis is a wrapper of NewPolynomialBasis.
// This function creates a new powerBasis from the input ciphertext.
// The input ciphertext is treated as the base monomial X used to
// generate the other powers X^{n}.
func NewPowerBasis(ct *rlwe.Ciphertext, basis bignum.Basis) he.PowerBasis {
	return he.NewPowerBasis(ct, basis)
}

// NewPolynomialEvaluator instantiates a new PolynomialEvaluator from a circuit.Evaluator.
// The default hefloat.Evaluator is compliant to the circuit.Evaluator interface.
// This method is allocation free.
func NewPolynomialEvaluator(params Parameters, eval he.Evaluator) *PolynomialEvaluator {
	return &PolynomialEvaluator{
		Parameters:             params,
		EvaluatorForPolynomial: &defaultCircuitEvaluatorForPolynomial{Evaluator: eval},
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
// pol: a *bignum.Polynomial, *Polynomial or *PolynomialVector
// targetScale: the desired output scale. This value shouldn't differ too much from the original ciphertext scale. It can
// for example be used to correct small deviations in the ciphertext scale and reset it to the default scale.
func (eval PolynomialEvaluator) Evaluate(ct *rlwe.Ciphertext, p interface{}, targetScale rlwe.Scale) (opOut *rlwe.Ciphertext, err error) {

	var phe interface{}
	switch p := p.(type) {
	case Polynomial:
		phe = he.Polynomial(p)
	case PolynomialVector:
		phe = he.PolynomialVector(p)
	default:
		phe = p
	}

	levelsConsumedPerRescaling := eval.Parameters.LevelsConsumedPerRescaling()

	return he.EvaluatePolynomial(eval, ct, phe, targetScale, levelsConsumedPerRescaling, &simEvaluator{eval.Parameters, levelsConsumedPerRescaling})
}

// EvaluateFromPowerBasis evaluates a polynomial using the provided PowerBasis, holding pre-computed powers of X.
// This method is the same as Evaluate except that the encrypted input is a PowerBasis.
// See Evaluate for additional information.
func (eval PolynomialEvaluator) EvaluateFromPowerBasis(pb he.PowerBasis, p interface{}, targetScale rlwe.Scale) (opOut *rlwe.Ciphertext, err error) {

	var phe interface{}
	switch p := p.(type) {
	case Polynomial:
		phe = he.Polynomial(p)
	case PolynomialVector:
		phe = he.PolynomialVector(p)
	default:
		phe = p
	}

	levelsConsumedPerRescaling := eval.Parameters.LevelsConsumedPerRescaling()

	if _, ok := pb.Value[1]; !ok {
		return nil, fmt.Errorf("cannot EvaluateFromPowerBasis: X^{1} is nil")
	}

	return he.EvaluatePolynomial(eval, pb, phe, targetScale, levelsConsumedPerRescaling, &simEvaluator{eval.Parameters, levelsConsumedPerRescaling})
}

// CoefficientGetter is a struct that implements the
// he.CoefficientGetter[*bignum.Complex] interface.
type CoefficientGetter struct {
	Values []*bignum.Complex
}

// GetVectorCoefficient return a slice []*bignum.Complex containing the k-th coefficient
// of each polynomial of PolynomialVector indexed by its Mapping.
// See PolynomialVector for additional information about the Mapping.
func (c CoefficientGetter) GetVectorCoefficient(pol he.PolynomialVector, k int) (values []*bignum.Complex) {

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

// GetSingleCoefficient returns the k-th coefficient of Polynomial as the type *bignum.Complex.
func (c CoefficientGetter) GetSingleCoefficient(pol he.Polynomial, k int) (value *bignum.Complex) {
	return pol.Coeffs[k]
}

// ShallowCopy returns a thread-safe copy of the original CoefficientGetter.
func (c CoefficientGetter) ShallowCopy() he.CoefficientGetter[*bignum.Complex] {
	return &CoefficientGetter{Values: make([]*bignum.Complex, len(c.Values))}
}

// defaultCircuitEvaluatorForPolynomial is a struct implementing the interface he.EvaluatorForPolynomial.
type defaultCircuitEvaluatorForPolynomial struct {
	he.Evaluator
}

// EvaluatePatersonStockmeyerPolynomialVector evaluates a pre-decomposed PatersonStockmeyerPolynomialVector on a pre-computed power basis [1, X^{1}, X^{2}, ..., X^{2^{n}}, X^{2^{n+1}}, ..., X^{2^{m}}]
func (eval defaultCircuitEvaluatorForPolynomial) EvaluatePatersonStockmeyerPolynomialVector(poly he.PatersonStockmeyerPolynomialVector, pb he.PowerBasis) (res *rlwe.Ciphertext, err error) {
	coeffGetter := he.CoefficientGetter[*bignum.Complex](&CoefficientGetter{Values: make([]*bignum.Complex, pb.Value[1].Slots())})
	return he.EvaluatePatersonStockmeyerPolynomialVector(eval, poly, coeffGetter, pb)
}
