// Package mod1 implements a homomorphic mod1 circuit for the CKKS scheme.
package mod1

import (
	"fmt"
	"math/big"
	"math/cmplx"

	"github.com/tuneinsight/lattigo/v6/circuits/ckks/polynomial"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
	"github.com/tuneinsight/lattigo/v6/utils/bignum"
)

// Evaluator is an evaluator providing an API for homomorphic evaluations of scaled x mod 1.
// All fields of this struct are public, enabling custom instantiations.
type Evaluator struct {
	*ckks.Evaluator
	PolynomialEvaluator *polynomial.Evaluator
	Parameters          Parameters
}

// NewEvaluator instantiates a new [Evaluator] evaluator from [ckks.Evaluator].
// This method is allocation free.
func NewEvaluator(eval *ckks.Evaluator, evalPoly *polynomial.Evaluator, Mod1Parameters Parameters) *Evaluator {
	return &Evaluator{Evaluator: eval, PolynomialEvaluator: evalPoly, Parameters: Mod1Parameters}
}

// EvaluateAndScaleNew calls [EvaluateNew] and scales the output values by `scaling` (without consuming additional depth).
// If `scaling` set to 1, then this is equivalent to simply calling [EvaluateNew].
func (eval Evaluator) EvaluateAndScaleNew(ct *rlwe.Ciphertext, scaling complex128) (res *rlwe.Ciphertext, err error) {

	evm := eval.Parameters

	if ct.Level() < evm.LevelQ {
		return nil, fmt.Errorf("cannot Evaluate: ct.Level() < Mod1Parameters.LevelQ")
	}

	if ct.Level() > evm.LevelQ {
		eval.DropLevel(ct, ct.Level()-evm.LevelQ)
	}

	res = ct.CopyNew()

	// Normalize the modular reduction to mod by 1 (division by Q)
	res.Scale = evm.ScalingFactor()

	// Compute the scales that the ciphertext should have before the double angle
	// formula such that after it it has the scale it had before the polynomial
	// evaluation

	Qi := eval.GetParameters().Q()

	targetScale := res.Scale
	for i := 0; i < evm.DoubleAngle; i++ {
		targetScale = targetScale.Mul(rlwe.NewScale(Qi[ct.Level()-evm.Mod1Poly.Depth()-evm.DoubleAngle+i+1]))
		targetScale.Value.Sqrt(&targetScale.Value)
	}

	// Division by 1/2^r and change of variable for the Chebyshev evaluation
	if evm.Mod1Type == CosDiscrete || evm.Mod1Type == CosContinuous {
		offset := new(big.Float).Sub(&evm.Mod1Poly.B, &evm.Mod1Poly.A)
		offset.Mul(offset, new(big.Float).SetFloat64(evm.IntervalShrinkFactor()))
		offset.Quo(new(big.Float).SetFloat64(-0.5), offset)
		if err = eval.Add(res, offset, res); err != nil {
			return nil, fmt.Errorf("cannot Evaluate: %w", err)
		}
	}

	// Double angle
	sqrt2pi := complex(evm.Sqrt2Pi, 0)

	var mod1Poly bignum.Polynomial
	if evm.Mod1InvPoly == nil {

		scaling := cmplx.Pow(scaling, complex(1/evm.IntervalShrinkFactor(), 0))

		mul := bignum.NewComplexMultiplier().Mul

		mod1Poly = evm.Mod1Poly.Clone()

		scalingPowBig := bignum.NewComplex().SetComplex128(scaling)

		for i := range mod1Poly.Coeffs {
			if mod1Poly.Coeffs[i] != nil {
				mul(mod1Poly.Coeffs[i], scalingPowBig, mod1Poly.Coeffs[i])
			}
		}

		sqrt2pi *= scaling

	} else {
		mod1Poly = evm.Mod1Poly
	}

	// Chebyshev evaluation
	if res, err = eval.PolynomialEvaluator.Evaluate(res, mod1Poly, rlwe.NewScale(targetScale)); err != nil {
		return nil, fmt.Errorf("cannot Evaluate: %w", err)
	}

	for i := 0; i < evm.DoubleAngle; i++ {
		sqrt2pi *= sqrt2pi

		if err = eval.MulRelin(res, res, res); err != nil {
			return nil, fmt.Errorf("cannot Evaluate: %w", err)
		}

		if err = eval.Add(res, res, res); err != nil {
			return nil, fmt.Errorf("cannot Evaluate: %w", err)
		}

		if err = eval.Add(res, -sqrt2pi, res); err != nil {
			return nil, fmt.Errorf("cannot Evaluate: %w", err)
		}

		if err = eval.Rescale(res, res); err != nil {
			return nil, fmt.Errorf("cannot Evaluate: %w", err)
		}
	}

	// ArcSine
	if evm.Mod1InvPoly != nil {

		mul := bignum.NewComplexMultiplier().Mul

		mod1InvPoly := evm.Mod1InvPoly.Clone()

		scalingBig := bignum.NewComplex().SetComplex128(scaling)

		for i := range mod1InvPoly.Coeffs {
			if mod1InvPoly.Coeffs[i] != nil {
				mul(mod1InvPoly.Coeffs[i], scalingBig, mod1InvPoly.Coeffs[i])
			}
		}

		if res, err = eval.PolynomialEvaluator.Evaluate(res, mod1InvPoly, res.Scale); err != nil {
			return nil, fmt.Errorf("cannot Evaluate: %w", err)
		}
	}

	// Multiplies back by q
	res.Scale = ct.Scale
	return res, nil
}

// EvaluateNew applies an homomorphic mod Q on a vector scaled by Delta, scaled down to mod 1:
//
//  1. Delta * (Q/Delta * I(X) + m(X)) (Delta = scaling factor, I(X) integer poly, m(X) message)
//  2. Delta * (I(X) + Delta/Q * m(X)) (divide by Q/Delta)
//  3. Delta * (Delta/Q * m(X)) (x mod 1)
//  4. Delta * (m(X)) (multiply back by Q/Delta)
//
// Since Q is not a power of two, but Delta is, then does an approximate division by the closest
// power of two to Q instead. Hence, it assumes that the input plaintext is already scaled by
// the correcting factor Q/2^{round(log(Q))}.
//
// !! Assumes that the input is normalized by 1/K for K the range of the approximation.
//
// Scaling back error correction by 2^{round(log(Q))}/Q afterward is included in the polynomial
func (eval Evaluator) EvaluateNew(ct *rlwe.Ciphertext) (*rlwe.Ciphertext, error) {
	return eval.EvaluateAndScaleNew(ct, 1)
}
