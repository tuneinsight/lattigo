// Package inverse implements a homomorphic inversion circuit for the CKKS scheme.
package inverse

import (
	"fmt"
	"math"

	"github.com/tuneinsight/lattigo/v6/circuits/ckks/bootstrapping"
	"github.com/tuneinsight/lattigo/v6/circuits/ckks/comparison"
	"github.com/tuneinsight/lattigo/v6/circuits/ckks/minimax"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
	"github.com/tuneinsight/lattigo/v6/utils"
)

// Evaluator is an evaluator used to evaluate the inverses of ciphertexts.
// All fields of this struct are public, enabling custom instantiations.
type Evaluator struct {
	Parameters ckks.Parameters
	*minimax.Evaluator
}

// NewEvaluator instantiates a new InverseEvaluator.
// This method is allocation free.
func NewEvaluator(params ckks.Parameters, eval *minimax.Evaluator) Evaluator {
	return Evaluator{
		Parameters: params,
		Evaluator:  eval,
	}
}

// EvaluateFullDomainNew computes 1/x for x in [-2^{log2max}, -2^{log2min}] U [2^{log2min}, 2^{log2max}].
//  1. Reduce the interval from [-max, -min] U [min, max] to [-1, -min] U [min, 1] by computing an approximate
//     inverse c such that |c * x| <= 1. For |x| > 1, c tends to 1/x while for |x| < c tends to 1.
//     This is done by using the work Efficient Homomorphic Evaluation on Large Intervals (https://eprint.iacr.org/2022/280.pdf).
//  2. Compute |c * x| = sign(x * c) * (x * c), this is required for the next step, which can only accept positive values.
//  3. Compute y' = 1/(|c * x|) with the iterative Goldschmidt division algorithm.
//  4. Compute y = y' * c * sign(x * c)
//
// The user can provide a minimax composite polynomial (signMinimaxPoly) for the sign function in the interval
// [-1-e, -2^{log2min}] U [2^{log2min}, 1+e] (where e is an upperbound on the scheme error).
// If no such polynomial is provided, then the [DefaultMinimaxCompositePolynomialForSign] is used by default.
// Note that the precision of the output of sign(x * c) does not impact the circuit precision since this value ends up being both at
// the numerator and denominator, thus cancelling itself.
func (eval Evaluator) EvaluateFullDomainNew(ct *rlwe.Ciphertext, log2min, log2max float64, signMinimaxPoly ...minimax.Polynomial) (cInv *rlwe.Ciphertext, err error) {

	var poly minimax.Polynomial
	if len(signMinimaxPoly) == 1 {
		poly = signMinimaxPoly[0]
	} else {
		poly = minimax.NewPolynomial(comparison.DefaultCompositePolynomialForSign)
	}

	return eval.evaluateNew(ct, log2min, log2max, true, poly)
}

// EvaluatePositiveDomainNew computes 1/x for x in [2^{log2min}, 2^{log2max}].
//  1. Reduce the interval from [min, max] to [min, 1] by computing an approximate
//     inverse c such that |c * x| <= 1. For |x| > 1, c tends to 1/x while for |x| < c tends to 1.
//     This is done by using the work Efficient Homomorphic Evaluation on Large Intervals (https://eprint.iacr.org/2022/280.pdf).
//  2. Compute y' = 1/(c * x) with the iterative Goldschmidt division algorithm.
//  3. Compute y = y' * c
func (eval Evaluator) EvaluatePositiveDomainNew(ct *rlwe.Ciphertext, log2min, log2max float64) (cInv *rlwe.Ciphertext, err error) {
	return eval.evaluateNew(ct, log2min, log2max, false, nil)
}

// EvaluateNegativeDomainNew computes 1/x for x in [-2^{log2max}, -2^{log2min}].
//  1. Reduce the interval from [-max, -min] to [-1, -min] by computing an approximate
//     inverse c such that |c * x| <= 1. For |x| > 1, c tends to 1/x while for |x| < c tends to 1.
//     This is done by using the work Efficient Homomorphic Evaluation on Large Intervals (https://eprint.iacr.org/2022/280.pdf).
//  2. Compute y' = 1/(c * x) with the iterative Goldschmidt division algorithm.
//  3. Compute y = y' * c
func (eval Evaluator) EvaluateNegativeDomainNew(ct *rlwe.Ciphertext, log2min, log2max float64) (cInv *rlwe.Ciphertext, err error) {

	var ctNeg *rlwe.Ciphertext
	if ctNeg, err = eval.MulNew(ct, -1); err != nil {
		return
	}

	if cInv, err = eval.EvaluatePositiveDomainNew(ctNeg, log2min, log2max); err != nil {
		return
	}

	return cInv, eval.Mul(cInv, -1, cInv)
}

func (eval Evaluator) evaluateNew(ct *rlwe.Ciphertext, log2min, log2max float64, fulldomain bool, signMinimaxPoly minimax.Polynomial) (cInv *rlwe.Ciphertext, err error) {

	params := eval.Parameters

	levelsPerRescaling := params.LevelsConsumedPerRescaling()

	btp := eval.BtsEval

	var normalizationfactor *rlwe.Ciphertext

	// If max > 1, then normalizes the ciphertext interval from  [-max, -min] U [min, max]
	// to [-1, -min] U [min, 1], and returns the encrypted normalization factor.
	if log2max > 0 {

		if cInv, normalizationfactor, err = eval.IntervalNormalization(ct, log2max, btp); err != nil {
			return nil, fmt.Errorf("preprocessing: normalizationfactor: %w", err)
		}

	} else {
		cInv = ct.CopyNew()
	}

	var sign *rlwe.Ciphertext

	if fulldomain {

		if eval.Evaluator == nil {
			return nil, fmt.Errorf("preprocessing: cannot EvaluateNew: MinimaxCompositePolynomialEvaluator is nil but fulldomain is set to true")
		}

		// Computes the sign with precision [-1, -2^-a] U [2^-a, 1]
		if sign, err = eval.Evaluator.Evaluate(cInv, signMinimaxPoly); err != nil {
			return nil, fmt.Errorf("preprocessing: fulldomain: true -> sign: %w", err)
		}

		if sign.Level() < btp.MinimumInputLevel()+levelsPerRescaling {
			if sign, err = btp.Bootstrap(sign); err != nil {
				return nil, fmt.Errorf("preprocessing: fulldomain: true -> sign -> bootstrap(sign): %w", err)
			}
		}

		// Checks that cInv have at least one level remaining above the minimum
		// level required for the bootstrapping.
		if cInv.Level() < btp.MinimumInputLevel()+levelsPerRescaling {
			if cInv, err = btp.Bootstrap(cInv); err != nil {
				return nil, fmt.Errorf("preprocessing: fulldomain: true -> sign -> bootstrap(cInv): %w", err)
			}
		}

		// Gets |x| = x * sign(x)
		if err = eval.MulRelin(cInv, sign, cInv); err != nil {
			return nil, fmt.Errorf("preprocessing: fulldomain: true -> sign -> bootstrap -> mul(cInv, sign): %w", err)
		}

		if err = eval.Rescale(cInv, cInv); err != nil {
			return nil, fmt.Errorf("preprocessing: fulldomain: true -> sign -> bootstrap -> mul(cInv, sign) -> rescale: %w", err)
		}
	}

	// Computes the inverse of x in [min = 2^-a, 1]
	if cInv, err = eval.GoldschmidtDivisionNew(cInv, log2min); err != nil {
		return nil, fmt.Errorf("division: GoldschmidtDivisionNew: %w", err)
	}

	var postprocessdepth int

	if normalizationfactor != nil || fulldomain {
		postprocessdepth += levelsPerRescaling
	}

	if fulldomain {
		postprocessdepth += levelsPerRescaling
	}

	// If x > 1 then multiplies back with the encrypted normalization vector
	if normalizationfactor != nil {

		if cInv.Level() < btp.MinimumInputLevel()+postprocessdepth {
			if cInv, err = btp.Bootstrap(cInv); err != nil {
				return nil, fmt.Errorf("normalizationfactor: bootstrap(cInv): %w", err)
			}
		}

		if normalizationfactor.Level() < btp.MinimumInputLevel()+postprocessdepth {
			if normalizationfactor, err = btp.Bootstrap(normalizationfactor); err != nil {
				return nil, fmt.Errorf("normalizationfactor: bootstrap(normalizationfactor): %w", err)
			}
		}

		if err = eval.MulRelin(cInv, normalizationfactor, cInv); err != nil {
			return nil, fmt.Errorf("normalizationfactor: mul(cInv): %w", err)
		}

		if err = eval.Rescale(cInv, cInv); err != nil {
			return nil, fmt.Errorf("normalizationfactor: rescale(cInv): %w", err)
		}
	}

	if fulldomain {

		// Multiplies back with the encrypted sign
		if err = eval.MulRelin(cInv, sign, cInv); err != nil {
			return nil, fmt.Errorf("fulldomain: mul(cInv):  %w", err)
		}

		if err = eval.Rescale(cInv, cInv); err != nil {
			return nil, fmt.Errorf("fulldomain: rescale(cInv): %w", err)
		}
	}

	return cInv, nil
}

// GoldschmidtDivisionNew homomorphically computes 1/x in the domain [0, 2].
// input: ct: Enc(x) with values in the interval [0+2^{-log2min}, 2-2^{-log2min}].
// output: Enc(1/x - e), where |e| <= (1-x)^2^(#iterations+1) -> the bit-precision doubles after each iteration.
// This method automatically estimates how many iterations are needed to
// achieve the optimal precision, which is derived from the plaintext scale.
// This method will return an error if the input ciphertext does not have enough
// remaining level and if the InverseEvaluator was instantiated with no bootstrapper.
// This method will return an error if something goes wrong with the bootstrapping or the rescaling operations.
func (eval Evaluator) GoldschmidtDivisionNew(ct *rlwe.Ciphertext, log2min float64) (ctInv *rlwe.Ciphertext, err error) {

	btp := eval.BtsEval

	params := eval.Parameters

	// 2^{-(prec - LogN + 1)}
	prec := float64(params.N()/2) / ct.Scale.Float64()

	// Estimates the number of iterations required to achieve the desired precision, given the interval [min, 2-min]
	start := 1 - math.Exp2(log2min)
	var iters = 1
	for start >= prec {
		start *= start // Doubles the bit-precision at each iteration
		iters++
	}

	// Minimum of 3 iterations
	// This minimum is set in the case where log2min is close to 0.
	iters = utils.Max(iters, 3)

	levelsPerRescaling := params.LevelsConsumedPerRescaling()

	if depth := iters * levelsPerRescaling; btp == nil && depth > ct.Level() {
		return nil, fmt.Errorf("cannot GoldschmidtDivisionNew: ct.Level()=%d < depth=%d and rlwe.Bootstrapper is nil", ct.Level(), depth)
	}

	var a *rlwe.Ciphertext
	if a, err = eval.MulNew(ct, -1); err != nil {
		return nil, err
	}

	b := a.CopyNew()

	if err = eval.Add(a, 2, a); err != nil {
		return nil, err
	}

	if err = eval.Add(b, 1, b); err != nil {
		return nil, err
	}

	for i := 1; i < iters; i++ {

		if btp != nil && (b.Level() == btp.MinimumInputLevel() || b.Level() == levelsPerRescaling-1) {
			if b, err = btp.Bootstrap(b); err != nil {
				return
			}
		}

		if btp != nil && (a.Level() == btp.MinimumInputLevel() || a.Level() == levelsPerRescaling-1) {
			if a, err = btp.Bootstrap(a); err != nil {
				return
			}
		}

		if err = eval.MulRelin(b, b, b); err != nil {
			return
		}

		if err = eval.Rescale(b, b); err != nil {
			return
		}

		if btp != nil && (b.Level() == btp.MinimumInputLevel() || b.Level() == levelsPerRescaling-1) {
			if b, err = btp.Bootstrap(b); err != nil {
				return
			}
		}

		var tmp *rlwe.Ciphertext
		if tmp, err = eval.MulRelinNew(a, b); err != nil {
			return
		}

		if err = eval.Rescale(tmp, tmp); err != nil {
			return
		}

		// a is at a higher level than tmp but at the same scale magnitude
		// We consume a level to bring a to the same level as tmp
		if err = eval.SetScale(a, tmp.Scale); err != nil {
			return
		}

		if err = eval.Add(a, tmp, a); err != nil {
			return
		}
	}

	return a, nil
}

// IntervalNormalization applies a modified version of Algorithm 2 of Efficient Homomorphic Evaluation on Large Intervals (https://eprint.iacr.org/2022/280)
// to normalize the interval from [-max, max] to [-1, 1]. Also returns the encrypted normalization factor.
//
// The original algorithm of https://eprint.iacr.org/2022/280 works by successive evaluation of a function that compresses values greater than some threshold
// to this threshold and let values smaller than the threshold untouched (mostly). The process is iterated, each time reducing the threshold by a pre-defined
// factor L. We can modify the algorithm to keep track of the compression factor so that we can get back the original values (before the compression) afterward.
//
// Given ct with values [-max, max], the method will compute y such that ct * y has values in [-1, 1].
// The normalization factor is independant to each slot:
//   - values smaller than 1 will have a normalization factor that tends to 1
//   - values greater than 1 will have a normalization factor that tends to 1/x
func (eval Evaluator) IntervalNormalization(ct *rlwe.Ciphertext, log2Max float64, btp bootstrapping.Bootstrapper) (ctNorm, ctNormFac *rlwe.Ciphertext, err error) {

	ctNorm = ct.CopyNew()

	levelsPerRescaling := eval.Parameters.LevelsConsumedPerRescaling()

	L := 2.45 // Compression factor (experimental)

	// n = log_{L}(max)
	n := math.Ceil(log2Max / math.Log2(L))

	for i := 0; i < int(n); i++ {

		if ctNorm.Level() < btp.MinimumInputLevel()+4*levelsPerRescaling {
			if ctNorm, err = btp.Bootstrap(ctNorm); err != nil {
				return
			}
		}

		if ctNormFac != nil && (ctNormFac.Level() == btp.MinimumInputLevel() || ctNormFac.Level() == levelsPerRescaling-1) {
			if ctNormFac, err = btp.Bootstrap(ctNormFac); err != nil {
				return
			}
		}

		// c = 2/sqrt(27 * L^(2 * (n-1-i)))
		c := 2.0 / math.Sqrt(27*math.Pow(L, 2*(n-1-float64(i))))

		// Depth 2
		// Computes: z = 1 - (y * c)^2

		// 1 level
		var z *rlwe.Ciphertext
		// (c * y)
		if z, err = eval.MulNew(ctNorm, c); err != nil {
			return
		}

		if err = eval.Rescale(z, z); err != nil {
			return
		}

		// 1 level
		// (c * y)^2
		if err = eval.MulRelin(z, z, z); err != nil {
			return
		}

		if err = eval.Rescale(z, z); err != nil {
			return
		}

		// -(c * y)^2
		if err = eval.Mul(z, -1, z); err != nil {
			return
		}

		// 1-(c * y)^2
		if err = eval.Add(z, 1, z); err != nil {
			return
		}

		if z.Level() < btp.MinimumInputLevel()+levelsPerRescaling {
			if z, err = btp.Bootstrap(z); err != nil {
				return
			}
		}

		// Updates the normalization factor
		if ctNormFac == nil {
			ctNormFac = z
		} else {

			// 1 level
			if err = eval.MulRelin(ctNormFac, z, ctNormFac); err != nil {
				return
			}

			if err = eval.Rescale(ctNormFac, ctNormFac); err != nil {
				return
			}
		}

		// Updates the ciphertext
		// 1 level
		if err = eval.MulRelin(ctNorm, z, ctNorm); err != nil {
			return
		}

		if err = eval.Rescale(ctNorm, ctNorm); err != nil {
			return
		}
	}

	return ctNorm, ctNormFac, nil
}
