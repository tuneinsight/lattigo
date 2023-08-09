package float

import (
	"fmt"
	"math/big"

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

	var polys []bignum.Polynomial
	if polys, err = GetSignPoly30Polynomials(prec); err != nil {
		return
	}

	two := new(big.Float).SetInt64(2)

	for _, poly := range polys {

		if params.RingType() == ring.Standard {
			for j := range poly.Coeffs {
				poly.Coeffs[j][0].Quo(poly.Coeffs[j][0], two)
			}
		}

		if sign.Level() < poly.Depth()+btp.MinimumInputLevel() {

			if params.MaxLevel() < poly.Depth()+btp.MinimumInputLevel() {
				return nil, fmt.Errorf("sign: parameters do not enable the evaluation of the circuit, missing %d levels", poly.Depth()+btp.MinimumInputLevel()-params.MaxLevel())
			}

			if sign, err = btp.Bootstrap(sign); err != nil {
				return
			}
		}

		if sign, err = eval.PolynomialEvaluator.Evaluate(sign, poly, ct.Scale); err != nil {
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

	var polys []bignum.Polynomial
	if polys, err = GetSignPoly30Polynomials(prec); err != nil {
		return
	}

	two := new(big.Float).SetInt64(2)

	for i, poly := range polys {

		if params.RingType() == ring.Standard {
			for j := range poly.Coeffs {
				poly.Coeffs[j][0].Quo(poly.Coeffs[j][0], two)
			}
		}

		// Changes the last poly to scale the output by 0.5 and add 0.5
		if i == len(polys)-1 {
			for j := range poly.Coeffs {
				poly.Coeffs[j][0].Quo(poly.Coeffs[j][0], two)
			}
		}

		if step.Level() < poly.Depth()+btp.MinimumInputLevel() {

			if params.MaxLevel() < poly.Depth()+btp.MinimumInputLevel() {
				return nil, fmt.Errorf("step: parameters do not enable the evaluation of the circuit, missing %d levels", poly.Depth()+btp.MinimumInputLevel()-params.MaxLevel())
			}

			if step, err = btp.Bootstrap(step); err != nil {
				return
			}
		}

		if step, err = eval.PolynomialEvaluator.Evaluate(step, poly, ct.Scale); err != nil {
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
