package hefloat

import (
	"fmt"

	"github.com/tuneinsight/lattigo/v5/core/rlwe"
	"github.com/tuneinsight/lattigo/v5/he"
	"github.com/tuneinsight/lattigo/v5/ring"
)

// EvaluatorForMinimaxCompositePolynomial defines a set of common and scheme agnostic method that are necessary to instantiate a MinimaxCompositePolynomialEvaluator.
type EvaluatorForMinimaxCompositePolynomial interface {
	he.Evaluator
	ConjugateNew(ct *rlwe.Ciphertext) (ctConj *rlwe.Ciphertext, err error)
}

// MinimaxCompositePolynomialEvaluator is an evaluator used to evaluate composite polynomials on ciphertexts.
// All fields of this struct are publics, enabling custom instantiations.
type MinimaxCompositePolynomialEvaluator struct {
	EvaluatorForMinimaxCompositePolynomial
	PolynomialEvaluator
	he.Bootstrapper[rlwe.Ciphertext]
	Parameters Parameters
}

// NewMinimaxCompositePolynomialEvaluator instantiates a new MinimaxCompositePolynomialEvaluator.
// The default hefloat.Evaluator is compliant to the EvaluatorForMinimaxCompositePolynomial interface.
// This method is allocation free.
func NewMinimaxCompositePolynomialEvaluator(params Parameters, eval EvaluatorForMinimaxCompositePolynomial, bootstrapper he.Bootstrapper[rlwe.Ciphertext]) *MinimaxCompositePolynomialEvaluator {
	return &MinimaxCompositePolynomialEvaluator{eval, *NewPolynomialEvaluator(params, eval), bootstrapper, params}
}

// Evaluate evaluates the provided MinimaxCompositePolynomial on the input ciphertext.
func (eval MinimaxCompositePolynomialEvaluator) Evaluate(ct *rlwe.Ciphertext, mcp MinimaxCompositePolynomial) (res *rlwe.Ciphertext, err error) {

	params := eval.Parameters

	btp := eval.Bootstrapper

	levelsConsumedPerRescaling := params.LevelsConsumedPerRescaling()

	// Checks that the number of levels available after the bootstrapping is enough to evaluate all polynomials
	if maxDepth := mcp.MaxDepth() * levelsConsumedPerRescaling; params.MaxLevel() < maxDepth+btp.MinimumInputLevel() {
		return nil, fmt.Errorf("parameters do not enable the evaluation of the minimax composite polynomial, required levels is %d but parameters only provide %d levels", maxDepth+btp.MinimumInputLevel(), params.MaxLevel())
	}

	res = ct.CopyNew()

	for _, poly := range mcp {

		// Checks that res has enough level to evaluate the next polynomial, else bootstrap
		if res.Level() < poly.Depth()*params.LevelsConsumedPerRescaling()+btp.MinimumInputLevel() {
			if res, err = btp.Bootstrap(res); err != nil {
				return
			}
		}

		// Define the scale that res must have after the polynomial evaluation.
		// If we use the regular CKKS (with complex values), we chose a scale to be
		// half of the desired scale, so that (x + conj(x)/2) has the correct scale.
		var targetScale rlwe.Scale
		if params.RingType() == ring.Standard {
			targetScale = params.DefaultScale().Div(rlwe.NewScale(2))
		} else {
			targetScale = params.DefaultScale()
		}

		// Evaluate the polynomial
		if res, err = eval.PolynomialEvaluator.Evaluate(res, poly, targetScale); err != nil {
			return nil, fmt.Errorf("evaluate polynomial: %w", err)
		}

		// Clean the imaginary part (else it tends to explode)
		if params.RingType() == ring.Standard {

			// Reassigns the scale back to the original one
			res.Scale = res.Scale.Mul(rlwe.NewScale(2))

			var resConj *rlwe.Ciphertext
			if resConj, err = eval.ConjugateNew(res); err != nil {
				return
			}

			if err = eval.Add(res, resConj, res); err != nil {
				return
			}
		}
	}

	// Avoids float errors
	res.Scale = ct.Scale

	return
}
