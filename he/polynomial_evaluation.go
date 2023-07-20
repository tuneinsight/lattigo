package he

import (
	"fmt"
	"math/bits"

	"github.com/tuneinsight/lattigo/v4/rlwe"
)

// PolynomialEvaluatorInterface defines the set of common and scheme agnostic homomorphic operations
// that are required for the encrypted evaluation of plaintext polynomial.
type PolynomialEvaluatorInterface interface {
	EvaluatorInterface
	EvaluatePolynomialVectorFromPowerBasis(targetLevel int, pol PolynomialVector, pb PowerBasis, targetScale rlwe.Scale) (res *rlwe.Ciphertext, err error)
}

func EvaluatePatersonStockmeyerPolynomialVector(poly PatersonStockmeyerPolynomialVector, pb PowerBasis, eval PolynomialEvaluatorInterface) (res *rlwe.Ciphertext, err error) {

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
			Value:      make([]Polynomial, nbPoly),
			SlotsIndex: poly.SlotsIndex,
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
		if tmp[idx].Value, err = eval.EvaluatePolynomialVectorFromPowerBasis(level, polyVec, pb, scale); err != nil {
			return nil, fmt.Errorf("cannot EvaluatePatersonStockmeyerPolynomial: polynomial[%d]: %w", i, err)
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

				if err = evalMonomial(even.Value, odd.Value, pb.Value[deg], eval); err != nil {
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

// Evaluates a = a + b * xpow
func evalMonomial(a, b, xpow *rlwe.Ciphertext, eval PolynomialEvaluatorInterface) (err error) {

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

	if !a.PlaintextScale.InDelta(b.PlaintextScale, float64(rlwe.ScalePrecision-12)) {
		return fmt.Errorf("evalMonomial: scale discrepency: (rescale(b) * X^{n}).Scale = %v != a.Scale = %v", &a.PlaintextScale.Value, &b.PlaintextScale.Value)
	}

	if err = eval.Add(b, a, b); err != nil {
		return fmt.Errorf("evalMonomial: %w", err)
	}

	return
}
