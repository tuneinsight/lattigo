// Package bootstrapper implements the Bootstrapper struct which provides generic bootstrapping for the CKKS scheme (and RLWE ciphertexts by extension).
// It notably abstracts scheme switching and ring dimension switching, enabling efficient bootstrapping of ciphertexts in the Conjugate Invariant ring
// or multiple ciphertexts of a lower ring dimension.
package bootstrapper

import (
	"fmt"
	"math/big"
	"runtime"

	"github.com/tuneinsight/lattigo/v4/circuits/float/bootstrapper/bootstrapping"
	"github.com/tuneinsight/lattigo/v4/ckks"
	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/rlwe"
)

// Bootstrapper is a struct storing the bootstrapping
// parameters, the bootstrapping evaluation keys and
// pre-computed constant necessary to carry out the
// bootstrapping circuit.
type Bootstrapper struct {
	Parameters
	bridge       ckks.DomainSwitcher
	bootstrapper *bootstrapping.Bootstrapper

	xPow2N1    []ring.Poly
	xPow2InvN1 []ring.Poly
	xPow2N2    []ring.Poly
	xPow2InvN2 []ring.Poly

	evk *BootstrappingKeys
}

// NewBootstrapper instantiates a new bootstrapper.Bootstrapper from a set of bootstrapper.Parameters
// and a set of bootstrapper.BootstrappingKeys
func NewBootstrapper(btpParams Parameters, evk *BootstrappingKeys) (*Bootstrapper, error) {

	b := &Bootstrapper{}

	paramsN1 := btpParams.ResidualParameters
	paramsN2 := btpParams.Parameters.Parameters

	switch paramsN1.RingType() {
	case ring.Standard:
		if evk.EvkN1ToN2 == nil || evk.EvkN2ToN1 == nil {
			return nil, fmt.Errorf("cannot NewBootstrapper: evk.(BootstrappingKeys) is missing EvkN1ToN2 and EvkN2ToN1")
		}
	case ring.ConjugateInvariant:
		if evk.EvkCmplxToReal == nil || evk.EvkRealToCmplx == nil {
			return nil, fmt.Errorf("cannot NewBootstrapper: evk.(BootstrappingKeys) is missing EvkN1ToN2 and EvkN2ToN1")
		}

		var err error
		if b.bridge, err = ckks.NewDomainSwitcher(paramsN2, evk.EvkCmplxToReal, evk.EvkRealToCmplx); err != nil {
			return nil, fmt.Errorf("cannot NewBootstrapper: ckks.NewDomainSwitcher: %w", err)
		}

		// The switch to standard to conjugate invariant multiplies the scale by 2
		btpParams.SlotsToCoeffsParameters.Scaling = new(big.Float).SetFloat64(0.5)
	}

	b.Parameters = btpParams
	b.evk = evk

	b.xPow2N2 = rlwe.GenXPow2(paramsN2.RingQ().AtLevel(0), paramsN2.LogN(), false)
	b.xPow2InvN2 = rlwe.GenXPow2(paramsN2.RingQ(), paramsN2.LogN(), true)

	if paramsN1.N() != paramsN2.N() {
		b.xPow2N1 = b.xPow2N2
		b.xPow2InvN1 = b.xPow2InvN2
	} else {
		b.xPow2N1 = rlwe.GenXPow2(paramsN1.RingQ().AtLevel(0), paramsN2.LogN(), false)
		b.xPow2InvN1 = rlwe.GenXPow2(paramsN1.RingQ(), paramsN2.LogN(), true)
	}

	var err error
	if b.bootstrapper, err = bootstrapping.NewBootstrapper(btpParams.Parameters, evk.EvkBootstrapping); err != nil {
		return nil, err
	}

	return b, nil
}

// Depth returns the multiplicative depth (number of levels consumed) of the bootstrapping circuit.
func (b Bootstrapper) Depth() int {
	return b.Parameters.Parameters.MaxLevel() - b.ResidualParameters.MaxLevel()
}

// OutputLevel returns the output level after the evaluation of the bootstrapping circuit.
func (b Bootstrapper) OutputLevel() int {
	return b.ResidualParameters.MaxLevel()
}

// MinimumInputLevel returns the minimum level at which a ciphertext must be to be
// bootstrapped.
func (b Bootstrapper) MinimumInputLevel() int {
	return b.LevelsConsumedPerRescaling()
}

// Bootstrap bootstraps a single ciphertext and returns the bootstrapped ciphertext.
func (b Bootstrapper) Bootstrap(ct *rlwe.Ciphertext) (*rlwe.Ciphertext, error) {
	cts := []rlwe.Ciphertext{*ct}
	cts, err := b.BootstrapMany(cts)
	if err != nil {
		return nil, err
	}
	return &cts[0], nil
}

// BootstrapMany bootstraps a list of ciphertext and returns the list of bootstrapped ciphertexts.
func (b Bootstrapper) BootstrapMany(cts []rlwe.Ciphertext) ([]rlwe.Ciphertext, error) {

	var err error

	switch b.ResidualParameters.RingType() {
	case ring.ConjugateInvariant:

		for i := 0; i < len(cts); i = i + 2 {

			even, odd := i, i+1

			ct0 := &cts[even]

			var ct1 *rlwe.Ciphertext
			if odd < len(cts) {
				ct1 = &cts[odd]
			}

			if ct0, ct1, err = b.refreshConjugateInvariant(ct0, ct1); err != nil {
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

		if cts, err = b.PackAndSwitchN1ToN2(cts); err != nil {
			return nil, fmt.Errorf("cannot BootstrapMany: %w", err)
		}

		for i := range cts {
			var ct *rlwe.Ciphertext
			if ct, err = b.bootstrapper.Bootstrap(&cts[i]); err != nil {
				return nil, fmt.Errorf("cannot BootstrapMany: %w", err)
			}
			cts[i] = *ct
		}

		if cts, err = b.UnpackAndSwitchN2Tn1(cts, LogSlots, nbCiphertexts); err != nil {
			return nil, fmt.Errorf("cannot BootstrapMany: %w", err)
		}
	}

	runtime.GC()

	for i := range cts {
		cts[i].Scale = b.ResidualParameters.DefaultScale()
	}

	return cts, err
}

// refreshConjugateInvariant takes two ciphertext in the Conjugate Invariant ring, repacks them in a single ciphertext in the standard ring
// using the real and imaginary part, bootstrap both ciphertext, and then extract back the real and imaginary part before repacking them
// individually in two new ciphertexts in the Conjugate Invariant ring.
func (b Bootstrapper) refreshConjugateInvariant(ctLeftN1Q0, ctRightN1Q0 *rlwe.Ciphertext) (ctLeftN1QL, ctRightN1QL *rlwe.Ciphertext, err error) {

	if ctLeftN1Q0 == nil {
		return nil, nil, fmt.Errorf("ctLeftN1Q0 cannot be nil")
	}

	// Switches ring from ring.ConjugateInvariant to ring.Standard
	ctLeftN2Q0 := b.RealToComplexNew(ctLeftN1Q0)

	// Repacks ctRightN1Q0 into the imaginary part of ctLeftN1Q0
	// which is zero since it comes from the Conjugate Invariant ring)
	if ctRightN1Q0 != nil {
		ctRightN2Q0 := b.RealToComplexNew(ctRightN1Q0)

		if err = b.bootstrapper.Evaluator.Mul(ctRightN2Q0, 1i, ctRightN2Q0); err != nil {
			return nil, nil, fmt.Errorf("cannot BootstrapMany: %w", err)
		}

		if err = b.bootstrapper.Evaluator.Add(ctLeftN2Q0, ctRightN2Q0, ctLeftN2Q0); err != nil {
			return nil, nil, fmt.Errorf("cannot BootstrapMany: %w", err)
		}
	}

	// Refreshes in the ring.Sstandard
	var ctLeftAndRightN2QL *rlwe.Ciphertext
	if ctLeftAndRightN2QL, err = b.bootstrapper.Bootstrap(ctLeftN2Q0); err != nil {
		return nil, nil, fmt.Errorf("cannot BootstrapMany: %w", err)
	}

	// The SlotsToCoeffs transformation scales the ciphertext by 0.5
	// This is done to compensate for the 2x factor introduced by ringStandardToConjugate(*).
	ctLeftAndRightN2QL.Scale = ctLeftAndRightN2QL.Scale.Mul(rlwe.NewScale(1 / 2.0))

	// Switches ring from ring.Standard to ring.ConjugateInvariant
	ctLeftN1QL = b.ComplexToRealNew(ctLeftAndRightN2QL)

	// Extracts the imaginary part
	if ctRightN1Q0 != nil {
		if err = b.bootstrapper.Mul(ctLeftAndRightN2QL, -1i, ctLeftAndRightN2QL); err != nil {
			return nil, nil, fmt.Errorf("cannot BootstrapMany: %w", err)
		}
		ctRightN1QL = b.ComplexToRealNew(ctLeftAndRightN2QL)
	}

	return
}
