package hefloat

import (
	"fmt"
	"math/big"

	"github.com/tuneinsight/lattigo/v5/core/rlwe"
	"github.com/tuneinsight/lattigo/v5/he"
)

// EvaluatorForMod1 defines a set of common and scheme agnostic
// methods that are necessary to instantiate a Mod1Evaluator.
// The default hefloat.Evaluator is compliant to this interface.
type EvaluatorForMod1 interface {
	he.Evaluator
	DropLevel(*rlwe.Ciphertext, int)
	GetParameters() *Parameters
}

// Mod1Evaluator is an evaluator providing an API for homomorphic evaluations of scaled x mod 1.
// All fields of this struct are public, enabling custom instantiations.
type Mod1Evaluator struct {
	EvaluatorForMod1
	PolynomialEvaluator *PolynomialEvaluator
	Mod1Parameters      Mod1Parameters
}

// NewMod1Evaluator instantiates a new Mod1Evaluator evaluator.
// The default hefloat.Evaluator is compliant to the EvaluatorForMod1 interface.
// This method is allocation free.
func NewMod1Evaluator(eval EvaluatorForMod1, evalPoly *PolynomialEvaluator, Mod1Parameters Mod1Parameters) *Mod1Evaluator {
	return &Mod1Evaluator{EvaluatorForMod1: eval, PolynomialEvaluator: evalPoly, Mod1Parameters: Mod1Parameters}
}

// EvaluateNew applies a homomorphic mod Q on a vector scaled by Delta, scaled down to mod 1 :
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
func (eval Mod1Evaluator) EvaluateNew(ct *rlwe.Ciphertext) (*rlwe.Ciphertext, error) {

	var err error

	evm := eval.Mod1Parameters

	if ct.Level() < evm.LevelStart() {
		return nil, fmt.Errorf("cannot Evaluate: ct.Level() < Mod1Parameters.LevelStart")
	}

	if ct.Level() > evm.LevelStart() {
		eval.DropLevel(ct, ct.Level()-evm.LevelStart())
	}

	// Stores default scales
	prevScaleCt := ct.Scale

	// Normalize the modular reduction to mod by 1 (division by Q)
	ct.Scale = evm.ScalingFactor()

	// Compute the scales that the ciphertext should have before the double angle
	// formula such that after it it has the scale it had before the polynomial
	// evaluation

	Qi := eval.GetParameters().Q()

	targetScale := ct.Scale
	for i := 0; i < evm.doubleAngle; i++ {
		targetScale = targetScale.Mul(rlwe.NewScale(Qi[evm.levelStart-evm.mod1Poly.Depth()-evm.doubleAngle+i+1]))
		targetScale.Value.Sqrt(&targetScale.Value)
	}

	// Division by 1/2^r and change of variable for the Chebyshev evaluation
	if evm.Mod1Type == CosDiscrete || evm.Mod1Type == CosContinuous {
		offset := new(big.Float).Sub(&evm.mod1Poly.B, &evm.mod1Poly.A)
		offset.Mul(offset, new(big.Float).SetFloat64(evm.scFac))
		offset.Quo(new(big.Float).SetFloat64(-0.5), offset)

		if err = eval.Add(ct, offset, ct); err != nil {
			return nil, fmt.Errorf("cannot Evaluate: %w", err)
		}
	}

	// Chebyshev evaluation
	if ct, err = eval.PolynomialEvaluator.Evaluate(ct, evm.mod1Poly, rlwe.NewScale(targetScale)); err != nil {
		return nil, fmt.Errorf("cannot Evaluate: %w", err)
	}

	// Double angle
	sqrt2pi := evm.sqrt2Pi
	for i := 0; i < evm.doubleAngle; i++ {
		sqrt2pi *= sqrt2pi

		if err = eval.MulRelin(ct, ct, ct); err != nil {
			return nil, fmt.Errorf("cannot Evaluate: %w", err)
		}

		if err = eval.Add(ct, ct, ct); err != nil {
			return nil, fmt.Errorf("cannot Evaluate: %w", err)
		}

		if err = eval.Add(ct, -sqrt2pi, ct); err != nil {
			return nil, fmt.Errorf("cannot Evaluate: %w", err)
		}

		if err = eval.Rescale(ct, ct); err != nil {
			return nil, fmt.Errorf("cannot Evaluate: %w", err)
		}
	}

	// ArcSine
	if evm.mod1InvPoly != nil {
		if ct, err = eval.PolynomialEvaluator.Evaluate(ct, *evm.mod1InvPoly, ct.Scale); err != nil {
			return nil, fmt.Errorf("cannot Evaluate: %w", err)
		}
	}

	// Multiplies back by q
	ct.Scale = prevScaleCt
	return ct, nil
}
