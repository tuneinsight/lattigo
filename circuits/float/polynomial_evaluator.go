package float

import (
	"fmt"

	"github.com/tuneinsight/lattigo/v4/circuits"
	"github.com/tuneinsight/lattigo/v4/ckks"
	"github.com/tuneinsight/lattigo/v4/rlwe"
	"github.com/tuneinsight/lattigo/v4/utils/bignum"
)

// PolynomialEvaluator is a wrapper of the circuits.PolynomialEvaluator.
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

// NewPolynomialEvaluator instantiates a new PolynomialEvaluator.
// eval can be a circuit.Evaluator, in which case it will use the default circuit.[...] polynomial
// evaluation function, or it can be an interface implementing circuits.EvaluatorForPolynomial, in
// which case it will use this interface to evaluate the polynomial.
func NewPolynomialEvaluator(params ckks.Parameters, eval interface{}) *PolynomialEvaluator {
	e := new(PolynomialEvaluator)

	switch eval := eval.(type) {
	case *ckks.Evaluator:
		e.EvaluatorForPolynomial = &defaultCircuitEvaluatorForPolynomial{Evaluator: eval}
	case circuits.EvaluatorForPolynomial:
		e.EvaluatorForPolynomial = eval
	default:
		panic(fmt.Sprintf("invalid eval type: must be circuits.Evaluator or circuits.EvaluatorForPolynomial but is %T", eval))
	}

	e.Parameters = params
	return e
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

	levelsConsummedPerRescaling := eval.Parameters.LevelsConsummedPerRescaling()

	return circuits.EvaluatePolynomial(eval, ct, pcircuits, targetScale, levelsConsummedPerRescaling, &simEvaluator{eval.Parameters, levelsConsummedPerRescaling})
}

// EvaluateFromPowerBasis evaluates a polynomial using the provided PowerBasis, holding pre-computed powers of X.
// This method is the same as Evaluate except that the encrypted input is a PowerBasis.
// See Evaluate for additional informations.
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

	levelsConsummedPerRescaling := eval.Parameters.LevelsConsummedPerRescaling()

	if _, ok := pb.Value[1]; !ok {
		return nil, fmt.Errorf("cannot EvaluateFromPowerBasis: X^{1} is nil")
	}

	return circuits.EvaluatePolynomial(eval, pb, pcircuits, targetScale, levelsConsummedPerRescaling, &simEvaluator{eval.Parameters, levelsConsummedPerRescaling})
}

type CoefficientGetter struct {
	Values []*bignum.Complex
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

type defaultCircuitEvaluatorForPolynomial struct {
	circuits.Evaluator
}

func (eval defaultCircuitEvaluatorForPolynomial) EvaluatePatersonStockmeyerPolynomialVector(poly circuits.PatersonStockmeyerPolynomialVector, pb circuits.PowerBasis) (res *rlwe.Ciphertext, err error) {
	coeffGetter := circuits.CoefficientGetter[*bignum.Complex](&CoefficientGetter{Values: make([]*bignum.Complex, pb.Value[1].Slots())})
	return circuits.EvaluatePatersonStockmeyerPolynomialVector(eval, poly, coeffGetter, pb)
}
