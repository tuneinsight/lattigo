package minimax

import (
	"fmt"

	"github.com/tuneinsight/lattigo/v6/circuits/ckks/bootstrapping"
	"github.com/tuneinsight/lattigo/v6/circuits/ckks/polynomial"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
)

// Evaluator is an evaluator used to evaluate composite polynomials on ciphertexts.
// All fields of this struct are publics, enabling custom instantiations.
type Evaluator struct {
	*ckks.Evaluator
	PolyEval   *polynomial.Evaluator
	BtsEval    bootstrapping.Bootstrapper
	Parameters ckks.Parameters
}

// NewEvaluator instantiates a new Evaluator.
// This method is allocation free.
func NewEvaluator(params ckks.Parameters, eval *ckks.Evaluator, btsEval bootstrapping.Bootstrapper) *Evaluator {
	return &Evaluator{eval, polynomial.NewEvaluator(params, eval), btsEval, params}
}

// Evaluate evaluates the provided MinimaxCompositePolynomial on the input ciphertext.
func (eval Evaluator) Evaluate(ct *rlwe.Ciphertext, mcp Polynomial) (res *rlwe.Ciphertext, err error) {

	params := eval.Parameters

	btp := eval.BtsEval

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
		if res, err = eval.PolyEval.Evaluate(res, poly, targetScale); err != nil {
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
