package float

import (
	"fmt"
	//"math/big"

	"github.com/tuneinsight/lattigo/v4/circuits"
	"github.com/tuneinsight/lattigo/v4/ckks"
	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/rlwe"
)

// EvaluatorForMinimaxCompositePolynomial defines a set of common and scheme agnostic method that are necessary to instantiate a MinimaxCompositePolynomialEvaluator.
type EvaluatorForMinimaxCompositePolynomial interface {
	circuits.EvaluatorForPolynomialEvaluation
	circuits.Evaluator
	ConjugateNew(ct *rlwe.Ciphertext) (ctConj *rlwe.Ciphertext, err error)
}

// MinimaxCompositePolynomialEvaluator is an evaluator used to evaluate composite polynomials on ciphertexts.
type MinimaxCompositePolynomialEvaluator struct {
	EvaluatorForMinimaxCompositePolynomial
	*PolynomialEvaluator
	rlwe.Bootstrapper
	Parameters ckks.Parameters
}

// NewMinimaxCompositePolynomialEvaluator instantiates a new MinimaxCompositePolynomialEvaluator from an EvaluatorForMinimaxCompositePolynomial.
// This method is allocation free.
func NewMinimaxCompositePolynomialEvaluator(params ckks.Parameters, eval EvaluatorForMinimaxCompositePolynomial, bootstrapper rlwe.Bootstrapper) *MinimaxCompositePolynomialEvaluator {
	return &MinimaxCompositePolynomialEvaluator{eval, NewPolynomialEvaluator(params, eval), bootstrapper, params}
}

func (eval MinimaxCompositePolynomialEvaluator) Evaluate(ct *rlwe.Ciphertext, mcp MinimaxCompositePolynomial) (res *rlwe.Ciphertext, err error) {

	params := eval.Parameters

	btp := eval.Bootstrapper

	levelsConsummedPerRescaling := params.LevelsConsummedPerRescaling()

	// Checks that the number of levels available after the bootstrapping is enough to evaluate all polynomials
	if maxDepth := mcp.MaxDepth() * levelsConsummedPerRescaling; params.MaxLevel() < maxDepth+btp.MinimumInputLevel() {
		return nil, fmt.Errorf("parameters do not enable the evaluation of the minimax composite polynomial, required levels is %d but parameters only provide %d levels", maxDepth+btp.MinimumInputLevel(), params.MaxLevel())
	}

	res = ct.CopyNew()

	for _, poly := range mcp {

		// Checks that res has enough level to evaluate the next polynomial, else bootstrap
		if res.Level() < poly.Depth()*params.LevelsConsummedPerRescaling()+btp.MinimumInputLevel() {
			if res, err = btp.Bootstrap(res); err != nil {
				return
			}
		}

		// Define the scale that res must have after the polynomial evaluation.
		// If we use the regular CKKS (with complex values), we chose a scale to be
		// half of the desired scale, so that (x + conj(x)/2) has the correct scale.
		var targetScale rlwe.Scale
		if params.RingType() == ring.Standard {
			targetScale = res.Scale.Div(rlwe.NewScale(2))
		} else {
			targetScale = res.Scale
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

	return
}
