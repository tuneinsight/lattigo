package circuits

import (
	"fmt"
	"math/bits"

	"github.com/tuneinsight/lattigo/v4/rlwe"
	"github.com/tuneinsight/lattigo/v4/utils/bignum"
)

// EvaluatorForPolynomial defines a set of common and scheme agnostic method that are necessary to instantiate a PolynomialVectorEvaluator.
type EvaluatorForPolynomial interface {
	rlwe.ParameterProvider
	Evaluator
	Encode(values interface{}, pt *rlwe.Plaintext) (err error)
	GetEvaluatorBuffer() *rlwe.EvaluatorBuffers // TODO extract
}

// PolynomialVectorEvaluator defines a scheme agnostic method to evaluate P(X) = sum ci * X^{i}.
type PolynomialVectorEvaluator interface {
	EvaluatePolynomialVectorFromPowerBasis(targetLevel int, pol PolynomialVector, pb PowerBasis, targetScale rlwe.Scale) (res *rlwe.Ciphertext, err error)
}

// PolynomialEvaluator is an evaluator used to evaluate polynomials on ciphertexts.
type PolynomialEvaluator struct {
	EvaluatorForPolynomial
	*rlwe.EvaluatorBuffers
}

// EvaluatePolynomial is a generic and scheme agnostic method to evaluate polynomials on rlwe.Ciphertexts.
func EvaluatePolynomial(eval PolynomialEvaluator, evalp PolynomialVectorEvaluator, input interface{}, p interface{}, targetScale rlwe.Scale, levelsConsummedPerRescaling int, SimEval SimEvaluator) (opOut *rlwe.Ciphertext, err error) {

	var polyVec PolynomialVector
	switch p := p.(type) {
	case bignum.Polynomial:
		polyVec = PolynomialVector{Value: []Polynomial{{Polynomial: p, MaxDeg: p.Degree(), Lead: true, Lazy: false}}}
	case Polynomial:
		polyVec = PolynomialVector{Value: []Polynomial{p}}
	case PolynomialVector:
		polyVec = p
	default:
		return nil, fmt.Errorf("cannot Polynomial: invalid polynomial type, must be either bignum.Polynomial, circuits.Polynomial or circuits.PolynomialVector, but is %T", p)
	}

	var powerbasis PowerBasis
	switch input := input.(type) {
	case *rlwe.Ciphertext:
		powerbasis = NewPowerBasis(input, polyVec.Value[0].Basis)
	case PowerBasis:
		if input.Value[1] == nil {
			return nil, fmt.Errorf("cannot evaluatePolyVector: given PowerBasis.Value[1] is empty")
		}
		powerbasis = input
	default:
		return nil, fmt.Errorf("cannot evaluatePolyVector: invalid input, must be either *rlwe.Ciphertext or *PowerBasis")
	}

	if level, depth := powerbasis.Value[1].Level(), levelsConsummedPerRescaling*polyVec.Value[0].Depth(); level < depth {
		return nil, fmt.Errorf("%d levels < %d log(d) -> cannot evaluate poly", level, depth)
	}

	logDegree := bits.Len64(uint64(polyVec.Value[0].Degree()))
	logSplit := bignum.OptimalSplit(logDegree)

	var odd, even = false, false
	for _, p := range polyVec.Value {
		odd, even = odd || p.IsOdd, even || p.IsEven
	}

	// Computes all the powers of two with relinearization
	// This will recursively compute and store all powers of two up to 2^logDegree
	if err = powerbasis.GenPower(1<<(logDegree-1), false, eval); err != nil {
		return nil, err
	}

	// Computes the intermediate powers, starting from the largest, without relinearization if possible
	for i := (1 << logSplit) - 1; i > 2; i-- {
		if !(even || odd) || (i&1 == 0 && even) || (i&1 == 1 && odd) {
			if err = powerbasis.GenPower(i, polyVec.Value[0].Lazy, eval); err != nil {
				return nil, err
			}
		}
	}

	PS := polyVec.GetPatersonStockmeyerPolynomial(*eval.GetRLWEParameters(), powerbasis.Value[1].Level(), powerbasis.Value[1].Scale, targetScale, SimEval)

	if opOut, err = eval.EvaluatePatersonStockmeyerPolynomialVector(evalp, PS, powerbasis); err != nil {
		return nil, err
	}

	return opOut, err
}

// EvaluatePatersonStockmeyerPolynomialVector evaluates a pre-decomposed PatersonStockmeyerPolynomialVector on a pre-computed power basis [1, X^{1}, X^{2}, ..., X^{2^{n}}, X^{2^{n+1}}, ..., X^{2^{m}}]
func (eval PolynomialEvaluator) EvaluatePatersonStockmeyerPolynomialVector(pvEval PolynomialVectorEvaluator, poly PatersonStockmeyerPolynomialVector, pb PowerBasis) (res *rlwe.Ciphertext, err error) {

	type Poly struct {
		Degree int
		Value  *rlwe.Ciphertext
	}

	split := len(poly.Value[0].Value)

	tmp := make([]*Poly, split)

	nbPoly := len(poly.Value)

	// Small steps
	for i := range tmp {

		polyVec := PolynomialVector{
			Value:   make([]Polynomial, nbPoly),
			Mapping: poly.Mapping,
		}

		// Transposes the polynomial matrix
		for j := 0; j < nbPoly; j++ {
			polyVec.Value[j] = poly.Value[j].Value[i]
		}

		level := poly.Value[0].Value[i].Level
		scale := poly.Value[0].Value[i].Scale

		idx := split - i - 1
		tmp[idx] = new(Poly)
		tmp[idx].Degree = poly.Value[0].Value[i].Degree()
		if tmp[idx].Value, err = pvEval.EvaluatePolynomialVectorFromPowerBasis(level, polyVec, pb, scale); err != nil {
			return nil, fmt.Errorf("cannot EvaluatePolynomialVectorFromPowerBasis: polynomial[%d]: %w", i, err)
		}
	}

	// Loops as long as there is more than one sub-polynomial
	for len(tmp) != 1 {

		for i := 0; i < len(tmp); i++ {

			// If we reach the end of the list it means we weren't able to combine
			// the last two sub-polynomials which necessarily implies that that the
			// last one has degree smaller than the previous one and that there is
			// no next polynomial to combine it with.
			// Therefore we update it's degree to the one of the previous one.
			if i == len(tmp)-1 {
				tmp[i].Degree = tmp[i-1].Degree

				// If two consecutive sub-polynomials, from ascending degree order, have the
				// same degree, we combine them.
			} else if tmp[i].Degree == tmp[i+1].Degree {

				even, odd := tmp[i], tmp[i+1]

				deg := 1 << bits.Len64(uint64(tmp[i].Degree))

				if err = eval.EvaluateMonomial(even.Value, odd.Value, pb.Value[deg]); err != nil {
					return nil, err
				}

				odd.Degree = 2*deg - 1
				tmp[i] = nil

				i++
			}
		}

		// Discards processed sub-polynomials
		var idx int
		for i := range tmp {
			if tmp[i] != nil {
				tmp[idx] = tmp[i]
				idx++
			}
		}

		tmp = tmp[:idx]
	}

	if tmp[0].Value.Degree() == 2 {
		if err = eval.Relinearize(tmp[0].Value, tmp[0].Value); err != nil {
			return nil, fmt.Errorf("cannot EvaluatePatersonStockmeyerPolynomial: %w", err)
		}
	}

	if err = eval.Rescale(tmp[0].Value, tmp[0].Value); err != nil {
		return nil, err
	}

	return tmp[0].Value, nil
}

// EvaluateMonomial evaluates a monomial of the form a + b * X^{pow} and writes the results in b.
func (eval PolynomialEvaluator) EvaluateMonomial(a, b, xpow *rlwe.Ciphertext) (err error) {

	if b.Degree() == 2 {
		if err = eval.Relinearize(b, b); err != nil {
			return fmt.Errorf("evalMonomial: %w", err)
		}
	}

	if err = eval.Rescale(b, b); err != nil {
		return fmt.Errorf("evalMonomial: %w", err)
	}

	if err = eval.Mul(b, xpow, b); err != nil {
		return fmt.Errorf("evalMonomial: %w", err)
	}

	if !a.Scale.InDelta(b.Scale, float64(rlwe.ScalePrecision-12)) {
		return fmt.Errorf("evalMonomial: scale discrepency: (rescale(b) * X^{n}).Scale = %v != a.Scale = %v", &a.Scale.Value, &b.Scale.Value)
	}

	if err = eval.Add(b, a, b); err != nil {
		return fmt.Errorf("evalMonomial: %w", err)
	}

	return
}
