package float

import (
	"fmt"
	"math"

	"github.com/tuneinsight/lattigo/v4/circuits"
	"github.com/tuneinsight/lattigo/v4/ckks"
	"github.com/tuneinsight/lattigo/v4/rlwe"
)

// EvaluatorForInverse defines a set of common and scheme agnostic method that are necessary to instantiate an InverseEvaluator.
type EvaluatorForInverse interface {
	circuits.Evaluator
	SetScale(ct *rlwe.Ciphertext, scale rlwe.Scale) (err error)
}

// InverseEvaluator is an evaluator used to evaluate the inverses of ciphertexts.
type InverseEvaluator struct {
	EvaluatorForInverse
	*PieceWiseFunctionEvaluator
	Parameters ckks.Parameters
}

// NewInverseEvaluator instantiates a new InverseEvaluator from an EvaluatorForInverse.
// evalPWF can be nil and is not be required if 'canBeNegative' of EvaluateNew is set to false.
// This method is allocation free.
func NewInverseEvaluator(params ckks.Parameters, evalInv EvaluatorForInverse, evalPWF EvaluatorForPieceWiseFunction) InverseEvaluator {

	var PWFEval *PieceWiseFunctionEvaluator

	if evalPWF != nil {
		PWFEval = NewPieceWiseFunctionEvaluator(params, evalPWF)
	}

	return InverseEvaluator{
		EvaluatorForInverse:        evalInv,
		PieceWiseFunctionEvaluator: PWFEval,
		Parameters:                 params,
	}
}

// EvaluateNew computes 1/x for x in [-max, -min] U [min, max].
//  1. Reduce the interval from [-max, -min] U [min, max] to [-1, -min] U [min, 1] by computing an approximate
//     inverse c such that |c * x| <= 1. For |x| > 1, c tends to 1/x while for |x| < c tends to 1.
//     This is done by using the work Efficient Homomorphic Evaluation on Large Intervals (https://eprint.iacr.org/2022/280.pdf).
//  2. Compute |c * x| = sign(x * c) * (x * c), this is required for the next step, which can only accept positive values.
//  3. Compute y' = 1/(|c * x|) with the iterative Goldschmidt division algorithm.
//  4. Compute y = y' * c * sign(x * c)
//
// canBeNegative: if set to false, then step 2 is skipped.
// prec: the desired precision of the GoldschmidtDivisionNew given the interval [min, 1].
func (eval InverseEvaluator) EvaluateNew(ct *rlwe.Ciphertext, min, max float64, canBeNegative bool, btp rlwe.Bootstrapper) (cInv *rlwe.Ciphertext, err error) {

	params := eval.Parameters

	levelsPerRescaling := params.LevelsConsummedPerRescaling()

	var normalizationfactor *rlwe.Ciphertext

	// If max > 1, then normalizes the ciphertext interval from  [-max, -min] U [min, max]
	// to [-1, -min] U [min, 1], and returns the encrypted normalization factor.
	if max > 1.0 {

		if cInv, normalizationfactor, err = eval.IntervalNormalization(ct, max, btp); err != nil {
			return
		}

	} else {
		cInv = ct.CopyNew()
	}

	var sign *rlwe.Ciphertext

	if canBeNegative {

		if eval.PieceWiseFunctionEvaluator == nil {
			return nil, fmt.Errorf("cannot EvaluateNew: PieceWiseFunctionEvaluator is nil but canBeNegative is set to true")
		}

		// Computes the sign with precision [-1, -2^-a] U [2^-a, 1]
		if sign, err = eval.PieceWiseFunctionEvaluator.EvaluateSign(cInv, int(math.Ceil(math.Log2(1/min))), btp); err != nil {
			return nil, fmt.Errorf("canBeNegative: true -> sign: %w", err)
		}

		if sign, err = btp.Bootstrap(sign); err != nil {
			return
		}

		if cInv.Level() == btp.MinimumInputLevel() || cInv.Level() == levelsPerRescaling-1 {
			if cInv, err = btp.Bootstrap(cInv); err != nil {
				return nil, fmt.Errorf("canBeNegative: true -> sign -> bootstrap: %w", err)
			}
		}

		// Gets the absolute value
		if err = eval.MulRelin(cInv, sign, cInv); err != nil {
			return nil, fmt.Errorf("canBeNegative: true -> sign -> bootstrap -> mul(cInv, sign): %w", err)
		}

		if err = eval.Rescale(cInv, cInv); err != nil {
			return nil, fmt.Errorf("canBeNegative: true -> sign -> bootstrap -> mul(cInv, sign) -> rescale: %w", err)
		}
	}

	// Computes the inverse of x in [min = 2^-a, 1]
	if cInv, err = eval.GoldschmidtDivisionNew(cInv, min, btp); err != nil {
		return
	}

	if cInv, err = btp.Bootstrap(cInv); err != nil {
		return
	}

	// If x > 1 then multiplies back with the encrypted normalization vector
	if normalizationfactor != nil {

		if normalizationfactor, err = btp.Bootstrap(normalizationfactor); err != nil {
			return
		}

		if err = eval.MulRelin(cInv, normalizationfactor, cInv); err != nil {
			return
		}

		if err = eval.Rescale(cInv, cInv); err != nil {
			return
		}
	}

	if canBeNegative {
		// Multiplies back with the encrypted sign
		if err = eval.MulRelin(cInv, sign, cInv); err != nil {
			return
		}

		if err = eval.Rescale(cInv, cInv); err != nil {
			return
		}

		if cInv, err = btp.Bootstrap(cInv); err != nil {
			return
		}
	}

	return cInv, nil
}

// GoldschmidtDivisionNew homomorphically computes 1/x.
// input: ct: Enc(x) with values in the interval [0+minvalue, 2-minvalue].
// output: Enc(1/x - e), where |e| <= (1-x)^2^(#iterations+1) -> the bit-precision doubles after each iteration.
// The method automatically estimates how many iterations are needed to achieve the optimal precision, which is derived from the plaintext scale,
// and will returns an error if the input ciphertext does not have enough remaining level and if no bootstrapper was given.
// This method will return an error if something goes wrong with the bootstrapping or the rescaling operations.
func (eval InverseEvaluator) GoldschmidtDivisionNew(ct *rlwe.Ciphertext, minValue float64, btp rlwe.Bootstrapper) (ctInv *rlwe.Ciphertext, err error) {

	params := eval.Parameters

	// 2^{-(prec - LogN + 1)}
	prec := float64(params.N()/2) / ct.Scale.Float64()

	// Estimates the number of iterations required to achieve the desired precision, given the interval [min, 2-min]
	start := 1 - minValue
	var iters = 1
	for start >= prec {
		start *= start // Doubles the bit-precision at each iteration
		iters++
	}

	levelsPerRescaling := params.LevelsConsummedPerRescaling()

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
// Given ct with values [-max, max], the method will compute y such that ct * y has values in [-1, 1].
// The normalization factor is independant to each slot:
//   - values smaller than 1 will have a normalizes factor that tends to 1
//   - values greater than 1 will have a normalizes factor that tends to 1/x
func (eval InverseEvaluator) IntervalNormalization(ct *rlwe.Ciphertext, max float64, btp rlwe.Bootstrapper) (ctNorm, ctNormFac *rlwe.Ciphertext, err error) {

	ctNorm = ct.CopyNew()

	levelsPerRescaling := eval.Parameters.LevelsConsummedPerRescaling()

	L := 2.45 // Experimental

	n := math.Ceil(math.Log(max) / math.Log(L))

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
