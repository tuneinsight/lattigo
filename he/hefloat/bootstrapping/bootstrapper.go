package bootstrapping

import (
	"fmt"

	"github.com/tuneinsight/lattigo/v5/core/rlwe"
	"github.com/tuneinsight/lattigo/v5/he"
	"github.com/tuneinsight/lattigo/v5/ring"
)

// Ensures that the Evaluator complies to the he.Bootstrapper interface
var _ he.Bootstrapper[rlwe.Ciphertext] = (*Evaluator)(nil)

// Bootstrap bootstraps a single ciphertext and returns the bootstrapped ciphertext.
func (eval Evaluator) Bootstrap(ct *rlwe.Ciphertext) (*rlwe.Ciphertext, error) {
	cts := []rlwe.Ciphertext{*ct}
	cts, err := eval.BootstrapMany(cts)
	if err != nil {
		return nil, err
	}
	return &cts[0], nil
}

// BootstrapMany bootstraps a list of ciphertext and returns the list of bootstrapped ciphertexts.
func (eval Evaluator) BootstrapMany(cts []rlwe.Ciphertext) ([]rlwe.Ciphertext, error) {

	var err error

	switch eval.ResidualParameters.RingType() {
	case ring.ConjugateInvariant:

		for i := 0; i < len(cts); i = i + 2 {

			even, odd := i, i+1

			ct0 := &cts[even]

			var ct1 *rlwe.Ciphertext
			if odd < len(cts) {
				ct1 = &cts[odd]
			}

			if ct0, ct1, err = eval.EvaluateConjugateInvariant(ct0, ct1); err != nil {
				return nil, fmt.Errorf("cannot BootstrapMany: %w", err)
			}

			cts[even] = *ct0

			if ct1 != nil {
				cts[odd] = *ct1
			}
		}

	default:

		LogSlots := cts[0].LogSlots()
		nbCiphertexts := len(cts)

		if cts, err = eval.PackAndSwitchN1ToN2(cts); err != nil {
			return nil, fmt.Errorf("cannot BootstrapMany: %w", err)
		}

		for i := range cts {
			var ct *rlwe.Ciphertext
			if ct, err = eval.Evaluate(&cts[i]); err != nil {
				return nil, fmt.Errorf("cannot BootstrapMany: %w", err)
			}
			cts[i] = *ct
		}

		if cts, err = eval.UnpackAndSwitchN2Tn1(cts, LogSlots, nbCiphertexts); err != nil {
			return nil, fmt.Errorf("cannot BootstrapMany: %w", err)
		}
	}

	for i := range cts {
		cts[i].Scale = eval.ResidualParameters.DefaultScale()
	}

	return cts, err
}

// Depth returns the multiplicative depth (number of levels consumed) of the bootstrapping circuit.
func (eval Evaluator) Depth() int {
	return eval.BootstrappingParameters.MaxLevel() - eval.ResidualParameters.MaxLevel()
}

// OutputLevel returns the output level after the evaluation of the bootstrapping circuit.
func (eval Evaluator) OutputLevel() int {
	return eval.ResidualParameters.MaxLevel()
}

// MinimumInputLevel returns the minimum level at which a ciphertext must be to be bootstrapped.
func (eval Evaluator) MinimumInputLevel() int {
	return eval.BootstrappingParameters.LevelsConsumedPerRescaling()
}
