package float

import (
	"fmt"

	"github.com/tuneinsight/lattigo/v4/circuits"
	"github.com/tuneinsight/lattigo/v4/ckks"
	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/rlwe"
	"github.com/tuneinsight/lattigo/v4/utils/bignum"
)

// EvaluatorForPieceWiseFunction defines a set of common and scheme agnostic method that are necessary to instantiate a PieceWiseFunctionEvaluator.
type EvaluatorForPieceWiseFunction interface {
	circuits.EvaluatorForPolynomialEvaluation
	circuits.Evaluator
	ConjugateNew(ct *rlwe.Ciphertext) (ctConj *rlwe.Ciphertext, err error)
}

// PieceWiseFunctionEvaluator is an evaluator used to evaluate piecewise functions on ciphertexts.
type PieceWiseFunctionEvaluator struct {
	EvaluatorForPieceWiseFunction
	*PolynomialEvaluator
	Parameters ckks.Parameters
}

// NewPieceWiseFunctionEvaluator instantiates a new PieceWiseFunctionEvaluator from an EvaluatorForPieceWiseFunction.
// This method is allocation free.
func NewPieceWiseFunctionEvaluator(params ckks.Parameters, eval EvaluatorForPieceWiseFunction) *PieceWiseFunctionEvaluator {
	return &PieceWiseFunctionEvaluator{eval, NewPolynomialEvaluator(params, eval), params}
}

// EvaluateSign takes a ciphertext with values in the interval [-1, -2^{-alpha}] U [2^{-alpha}, 1] and returns
// -  1 if x is in [2^{-alpha}, 1]
// -  a value between -1 and 1 if x is in [-2^{-alpha}, 2^{-alpha}]
// - -1 if x is in [-1, -2^{-alpha}]
func (eval PieceWiseFunctionEvaluator) EvaluateSign(ct *rlwe.Ciphertext, prec int, btp rlwe.Bootstrapper) (sign *rlwe.Ciphertext, err error) {

	params := eval.Parameters

	sign = ct.CopyNew()

	var polys [][]float64
	if polys, err = GetSignPoly30Coefficients(prec); err != nil {
		return
	}

	for _, coeffs := range polys {

		c128 := make([]complex128, len(coeffs))

		if params.RingType() == ring.Standard {
			for j := range c128 {
				c128[j] = complex(coeffs[j]/2, 0)
			}
		} else {
			for j := range c128 {
				c128[j] = complex(coeffs[j], 0)
			}
		}

		pol := bignum.NewPolynomial(bignum.Chebyshev, c128, nil)

		if sign.Level() < pol.Depth()+btp.MinimumInputLevel() {

			if params.MaxLevel() < pol.Depth()+btp.MinimumInputLevel() {
				return nil, fmt.Errorf("sign: parameters do not enable the evaluation of the circuit, missing %d levels", pol.Depth()+btp.MinimumInputLevel()-params.MaxLevel())
			}

			if sign, err = btp.Bootstrap(sign); err != nil {
				return
			}
		}

		if sign, err = eval.PolynomialEvaluator.Evaluate(sign, pol, ct.Scale); err != nil {
			return nil, fmt.Errorf("sign: polynomial: %w", err)
		}

		// Clean the imaginary part (else it tends to expload)
		if params.RingType() == ring.Standard {

			var signConj *rlwe.Ciphertext
			if signConj, err = eval.ConjugateNew(sign); err != nil {
				return
			}

			if err = eval.Add(sign, signConj, sign); err != nil {
				return
			}
		}
	}

	return
}

// EvaluateStep takes a ciphertext with values in the interval [0, 0.5-2^{-alpha}] U [0.5+2^{-alpha}, 1] and returns
// -  1 if x is in [0.5+2^{-alpha}, 1]
// -  a value between 0 and 1 if x is in [0.5-2^{-alpha}, 0.5+2^{-alpha}]
// -  0 if x is in [0, 0.5-2^{-alpha}]
func (eval PieceWiseFunctionEvaluator) EvaluateStep(ct *rlwe.Ciphertext, prec int, btp rlwe.Bootstrapper) (step *rlwe.Ciphertext, err error) {

	params := eval.Parameters

	step = ct.CopyNew()

	var polys [][]float64
	if polys, err = GetSignPoly30Coefficients(prec); err != nil {
		return
	}

	for i, coeffs := range polys {

		c128 := make([]complex128, len(coeffs))

		if params.RingType() == ring.Standard {
			for j := range c128 {
				c128[j] = complex(coeffs[j]/2, 0)
			}
		} else {
			for j := range c128 {
				c128[j] = complex(coeffs[j], 0)
			}
		}

		// Changes the last poly to scale the output by 0.5 and add 0.5
		if i == len(polys)-1 {
			for j := range c128 {
				c128[j] /= 2
			}
		}

		pol := bignum.NewPolynomial(bignum.Chebyshev, c128, nil)

		if step.Level() < pol.Depth()+btp.MinimumInputLevel() {

			if params.MaxLevel() < pol.Depth()+btp.MinimumInputLevel() {
				return nil, fmt.Errorf("step: parameters do not enable the evaluation of the circuit, missing %d levels", pol.Depth()+btp.MinimumInputLevel()-params.MaxLevel())
			}

			if step, err = btp.Bootstrap(step); err != nil {
				return
			}
		}

		if step, err = eval.PolynomialEvaluator.Evaluate(step, pol, ct.Scale); err != nil {
			return nil, fmt.Errorf("step: polynomial: %w", err)
		}

		// Clean the imaginary part (else it tends to expload)
		if params.RingType() == ring.Standard {
			var stepConj *rlwe.Ciphertext
			if stepConj, err = eval.ConjugateNew(step); err != nil {
				return
			}

			if err = eval.Add(step, stepConj, step); err != nil {
				return
			}
		}
	}

	if err = eval.Add(step, 0.5, step); err != nil {
		return
	}

	return step, nil
}
