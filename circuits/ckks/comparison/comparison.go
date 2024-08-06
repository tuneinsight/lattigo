// Package comparison implements homomorphic comparisons for the CKKS scheme.
package comparison

import (
	"math/big"

	"github.com/tuneinsight/lattigo/v6/circuits/ckks/minimax"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
	"github.com/tuneinsight/lattigo/v6/utils"
	"github.com/tuneinsight/lattigo/v6/utils/bignum"
)

// Evaluator is an evaluator providing an API for homomorphic comparisons.
// All fields of this struct are public, enabling custom instantiations.
type Evaluator struct {
	Parameters ckks.Parameters
	*minimax.Evaluator
	MinimaxCompositeSignPolynomial minimax.Polynomial
}

// NewEvaluator instantiates a new ComparisonEvaluator.
// The default ckks.Evaluator is compliant with the EvaluatorForMinimaxCompositePolynomial interface.
//
// Giving a MinimaxCompositePolynomial is optional, but it is highly recommended to provide one that is optimized
// for the circuit requiring the comparisons as this polynomial will define the internal precision of all computations
// performed by this evaluator.
//
// The MinimaxCompositePolynomial must be a composite minimax approximation of the sign function:
// f(x) = 1 if x > 0, -1 if x < 0, else 0, in the interval [-1, 1].
// Such composite polynomial can be obtained with the function GenMinimaxCompositePolynomialForSign.
//
// If no MinimaxCompositePolynomial is given, then it will use by default the variable DefaultMinimaxCompositePolynomialForSign.
// See the doc of DefaultMinimaxCompositePolynomialForSign for additional information about the performance of this approximation.
//
// This method is allocation free if a MinimaxCompositePolynomial is given.
func NewEvaluator(params ckks.Parameters, eval *minimax.Evaluator, signPoly ...minimax.Polynomial) *Evaluator {
	if len(signPoly) == 1 {
		return &Evaluator{
			Parameters:                     params,
			Evaluator:                      eval,
			MinimaxCompositeSignPolynomial: signPoly[0],
		}
	} else {
		return &Evaluator{
			Parameters:                     params,
			Evaluator:                      eval,
			MinimaxCompositeSignPolynomial: minimax.NewPolynomial(DefaultCompositePolynomialForSign),
		}
	}
}

// DefaultCompositePolynomialForSign is an example of composite minimax polynomial
// for the sign function that is able to distinguish between value with a delta of up to
// 2^{-alpha=30}, tolerates a scheme error of 2^{-35} and outputs a binary value (-1, or 1)
// of up to 20x4 bits of precision.
//
// It was computed with GenMinimaxCompositePolynomialForSign(256, 30, 35, []int{15, 15, 15, 17, 31, 31, 31, 31})
// which outputs a minimax composite polynomial of precision 21.926741, which is further composed with
// CoeffsSignX4Cheby to bring it to ~80bits of precision.
var DefaultCompositePolynomialForSign = [][]string{
	{"0", "0.6371462957672043333", "0", "-0.2138032460610765328", "0", "0.1300439303835664499", "0", "-0.0948842756566191044", "0", "0.0760417811618939909", "0", "-0.0647714820920817557", "0", "0.0577904411211959048", "0", "-0.5275634328386103792"},
	{"0", "0.6371463830322414578", "0", "-0.2138032749880402509", "0", "0.1300439475440832118", "0", "-0.0948842877009570762", "0", "0.0760417903036533484", "0", "-0.0647714893343788749", "0", "0.0577904470018789283", "0", "-0.5275633669027163690"},
	{"0", "0.6371474873319408921", "0", "-0.2138036410457105809", "0", "0.1300441647026617059", "0", "-0.0948844401165889295", "0", "0.0760419059884502454", "0", "-0.0647715809823254389", "0", "0.0577905214191996406", "0", "-0.5275625325136631842"},
	{"0", "0.6370469776996076431", "0", "-0.2134526779726600620", "0", "0.1294300181775238920", "0", "-0.0939692999460324791", "0", "0.0747629355709698798", "0", "-0.0630298319949635571", "0", "0.0554299627688379896", "0", "-0.0504549111784642023", "0", "0.5242368268605847996"},
	{"0", "0.6371925153898374380", "0", "-0.2127272333844484291", "0", "0.1280350175397897124", "0", "-0.0918861831051024970", "0", "0.0719237384158242601", "0", "-0.0593247422790627989", "0", "0.0506973946536399213", "0", "-0.0444605229007162961", "0", "0.0397788020190944552", "0", "-0.0361705584687241925", "0", "0.0333397971860406254", "0", "-0.0310960060432036761", "0", "0.0293126335952747929", "0", "-0.0279042579223662982", "0", "0.0268135229627401517", "0", "-0.5128179323757194002"},
	{"0", "0.6484328404896112084", "0", "-0.2164688471885406655", "0", "0.1302737771018761402", "0", "-0.0934786176742356885", "0", "0.0731553324133884104", "0", "-0.0603252338481440981", "0", "0.0515366139595849853", "0", "-0.0451803385226980999", "0", "0.0404062758116036740", "0", "-0.0367241775307736352", "0", "0.0338327393147257876", "0", "-0.0315379870551266008", "0", "0.0297110181467332488", "0", "-0.0282647625290482803", "0", "0.0271406820054187399", "0", "-0.5041440308249296747"},
	{"0", "0.8988231150519633581", "0", "-0.2996064625122592138", "0", "0.1797645789317822353", "0", "-0.1284080039344265678", "0", "0.0998837306152582349", "0", "-0.0817422066647773587", "0", "0.0691963884439569899", "0", "-0.0600136111161848355", "0", "0.0530132660795356506", "0", "-0.0475133961913746909", "0", "0.0430936248086665091", "0", "-0.0394819050695222720", "0", "0.0364958013826412785", "0", "-0.0340100990129699835", "0", "0.0319381346687564699", "0", "-0.3095637759472512887"},
	{"0", "1.2654405107323937767", "0", "-0.4015427502443620045", "0", "0.2182109348265640036", "0", "-0.1341692540177466882", "0", "0.0852282854825304735", "0", "-0.0539043807248265057", "0", "0.0332611560159092728", "0", "-0.0197419082926337129", "0", "0.0111368708758574529", "0", "-0.0058990205011466309", "0", "0.0028925861201479251", "0", "-0.0012889673944941461", "0", "0.0005081425552893727", "0", "-0.0001696330470066833", "0", "0.0000440808328172753", "0", "-0.0000071549240608255"},
	minimax.CoeffsSignX4Cheby, // Quadruples the output precision (up to the scheme error)
}

// Sign evaluates f(x) = 1 if x > 0, -1 if x < 0, else 0.
// This will ensure that sign.Scale = params.DefaultScale().
func (eval Evaluator) Sign(op0 *rlwe.Ciphertext) (sign *rlwe.Ciphertext, err error) {
	return eval.Evaluate(op0, eval.MinimaxCompositeSignPolynomial)
}

// Step evaluates f(x) = 1 if x > 0, 0 if x < 0, else 0.5 (i.e. (sign+1)/2).
// This will ensure that step.Scale = params.DefaultScale().
func (eval Evaluator) Step(op0 *rlwe.Ciphertext) (step *rlwe.Ciphertext, err error) {

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
func (eval Evaluator) Max(op0, op1 *rlwe.Ciphertext) (max *rlwe.Ciphertext, err error) {

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
func (eval Evaluator) Min(op0, op1 *rlwe.Ciphertext) (min *rlwe.Ciphertext, err error) {

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

func (eval Evaluator) stepdiff(op0, op1 *rlwe.Ciphertext) (stepdiff *rlwe.Ciphertext, err error) {
	params := eval.Parameters

	// diff = op0 - op1
	var diff *rlwe.Ciphertext
	if diff, err = eval.SubNew(op0, op1); err != nil {
		return
	}

	// Required for the scale matching before the last multiplication.
	if diff.Level() < params.LevelsConsumedPerRescaling()*2 {
		if diff, err = eval.BtsEval.Bootstrap(diff); err != nil {
			return
		}
	}

	// step = 1 if diff > 0, 0 if diff < 0 else 0.5
	var step *rlwe.Ciphertext
	if step, err = eval.Step(diff); err != nil {
		return
	}

	// Required for the following multiplication
	if step.Level() < params.LevelsConsumedPerRescaling() {
		if step, err = eval.BtsEval.Bootstrap(step); err != nil {
			return
		}
	}

	// Extremum gate: op0 * step + op1 * (1 - step) = step * diff + op1
	level := utils.Min(diff.Level(), step.Level())

	ratio := rlwe.NewScale(1)
	for i := 0; i < params.LevelsConsumedPerRescaling(); i++ {
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
