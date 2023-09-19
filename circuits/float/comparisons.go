package float

import (
	"math/big"

	"github.com/tuneinsight/lattigo/v4/rlwe"
	"github.com/tuneinsight/lattigo/v4/utils"
	"github.com/tuneinsight/lattigo/v4/utils/bignum"
)

// ComparisonEvaluator is an evaluator providing an API for homomorphic comparisons.
type ComparisonEvaluator struct {
	MinimaxCompositePolynomialEvaluator
	MinimaxCompositeSignPolynomial MinimaxCompositePolynomial
}

// NewComparisonEvaluator instantiates a new ComparisonEvaluator from a MinimaxCompositePolynomialEvaluator and a MinimaxCompositePolynomial.
// The MinimaxCompositePolynomial must be a composite minimax approximation of the sign function: f(x) = 1 if x > 0, -1 if x < 0, else 0.
// This polynomial will define the internal precision of all computation performed by this evaluator and it can be obtained with the function
// GenMinimaxCompositePolynomialForSign.
func NewComparisonEvaluator(eval *MinimaxCompositePolynomialEvaluator, signPoly MinimaxCompositePolynomial) *ComparisonEvaluator {
	return &ComparisonEvaluator{*eval, signPoly}
}

// Sign evaluates f(x) = 1 if x > 0, -1 if x < 0, else 0.
// This will ensure that sign.Scale = params.DefaultScale().
func (eval ComparisonEvaluator) Sign(op0 *rlwe.Ciphertext) (sign *rlwe.Ciphertext, err error) {
	return eval.Evaluate(op0, eval.MinimaxCompositeSignPolynomial)
}

// Step evaluates f(x) = 1 if x > 0, 0 if x < 0, else 0.5 (i.e. (sign+1)/2).
// This will ensure that step.Scale = params.DefaultScale().
func (eval ComparisonEvaluator) Step(op0 *rlwe.Ciphertext) (step *rlwe.Ciphertext, err error) {

	n := len(eval.MinimaxCompositeSignPolynomial)

	stepPoly := make([]bignum.Polynomial, n)

	for i := 0; i < n; i++ {
		stepPoly[i] = eval.MinimaxCompositeSignPolynomial[i]
	}

	half := new(big.Float).SetFloat64(0.5)

	// (x+1)/2
	lastPoly := eval.MinimaxCompositeSignPolynomial[n-1].Clone()
	for i := range lastPoly.Coeffs {
		lastPoly.Coeffs[i][0].Mul(lastPoly.Coeffs[i][0], half)
	}
	lastPoly.Coeffs[0][0].Add(lastPoly.Coeffs[0][0], half)

	stepPoly[n-1] = lastPoly

	return eval.Evaluate(op0, stepPoly)
}

// Max returns the smooth maximum of op0 and op1, which is defined as: op0 * x + op1 * (1-x) where x = step(diff = op0-op1).
// Use must ensure that:
//   - op0 + op1 is in the interval [-1, 1].
//   - op0.Scale = op1.Scale.
//
// This method ensures that max.Scale = params.DefaultScale.
func (eval ComparisonEvaluator) Max(op0, op1 *rlwe.Ciphertext) (max *rlwe.Ciphertext, err error) {

	// step * diff
	var stepdiff *rlwe.Ciphertext
	if stepdiff, err = eval.stepdiff(op0, op1); err != nil {
		return
	}

	// max = step * diff + op1
	if err = eval.Add(stepdiff, op1, stepdiff); err != nil {
		return
	}

	return stepdiff, nil
}

// Min returns the smooth min of op0 and op1, which is defined as: op0 * (1-x) + op1 * x where x = step(diff = op0-op1)
// Use must ensure that:
//   - op0 + op1 is in the interval [-1, 1].
//   - op0.Scale = op1.Scale.
//
// This method ensures that min.Scale = params.DefaultScale.
func (eval ComparisonEvaluator) Min(op0, op1 *rlwe.Ciphertext) (min *rlwe.Ciphertext, err error) {

	// step * diff
	var stepdiff *rlwe.Ciphertext
	if stepdiff, err = eval.stepdiff(op0, op1); err != nil {
		return
	}

	// min = op0 - step * diff
	if err = eval.Sub(op0, stepdiff, stepdiff); err != nil {
		return
	}

	return stepdiff, nil
}

func (eval ComparisonEvaluator) stepdiff(op0, op1 *rlwe.Ciphertext) (stepdiff *rlwe.Ciphertext, err error) {
	params := eval.Parameters

	// diff = op0 - op1
	var diff *rlwe.Ciphertext
	if diff, err = eval.SubNew(op0, op1); err != nil {
		return
	}

	// Required for the scale matching before the last multiplication.
	if diff.Level() < params.LevelsConsummedPerRescaling()*2 {
		if diff, err = eval.Bootstrap(diff); err != nil {
			return
		}
	}

	// step = 1 if diff > 0, 0 if diff < 0 else 0.5
	var step *rlwe.Ciphertext
	if step, err = eval.Step(diff); err != nil {
		return
	}

	// Required for the following multiplication
	if step.Level() < params.LevelsConsummedPerRescaling() {
		if step, err = eval.Bootstrap(step); err != nil {
			return
		}
	}

	// Extremum gate: op0 * step + op1 * (1 - step) = step * diff + op1
	level := utils.Min(diff.Level(), step.Level())

	ratio := rlwe.NewScale(1)
	for i := 0; i < params.LevelsConsummedPerRescaling(); i++ {
		ratio = ratio.Mul(rlwe.NewScale(params.Q()[level-i]))
	}

	ratio = ratio.Div(diff.Scale)
	if err = eval.Mul(diff, &ratio.Value, diff); err != nil {
		return
	}

	if err = eval.Rescale(diff, diff); err != nil {
		return
	}
	diff.Scale = diff.Scale.Mul(ratio)

	// max = step * diff
	if err = eval.MulRelin(diff, step, diff); err != nil {
		return
	}

	if err = eval.Rescale(diff, diff); err != nil {
		return
	}

	return diff, nil
}
