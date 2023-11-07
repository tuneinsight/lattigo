package bootstrapping

import (
	"fmt"
	"math/big"
	"math/bits"
	"runtime"

	"github.com/tuneinsight/lattigo/v4/core/rlwe"
	"github.com/tuneinsight/lattigo/v4/he"
	"github.com/tuneinsight/lattigo/v4/he/hefloat"
	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/schemes/ckks"
	"github.com/tuneinsight/lattigo/v4/utils"
)

// Bootstrapper is a high level wrapper of the bootstrapping circuit that
// stores the bootstrapping parameters, the bootstrapping evaluation keys and
// pre-computed constant necessary to carry out the bootstrapping circuit.
type Bootstrapper struct {
	*Parameters
	*CoreBootstrapper
	ckks.DomainSwitcher

	// [1, x, x^2, x^4, ..., x^N1/2] / (X^N1 +1)
	xPow2N1 []ring.Poly
	// [1, x, x^2, x^4, ..., x^N2/2] / (X^N2 +1)
	xPow2N2 []ring.Poly
	// [1, x^-1, x^-2, x^-4, ..., x^-N2/2] / (X^N2 +1)
	xPow2InvN2 []ring.Poly
}

// Ensures that the bootstrapper complies to the he.Bootstrapper interface
var _ he.Bootstrapper[rlwe.Ciphertext] = (*Bootstrapper)(nil)

// NewBootstrapper instantiates a new bootstrapper.Bootstrapper from a set
// of bootstrapping.Parameters and a set of bootstrapping.EvaluationKeys.
// It notably abstracts scheme switching and ring dimension switching,
// enabling efficient bootstrapping of ciphertexts in the Conjugate
// Invariant ring or multiple ciphertexts of a lower ring dimension.
func NewBootstrapper(btpParams Parameters, evk *EvaluationKeys) (*Bootstrapper, error) {

	b := &Bootstrapper{}

	paramsN1 := btpParams.ResidualParameters
	paramsN2 := btpParams.BootstrappingParameters

	switch paramsN1.RingType() {
	case ring.Standard:
		if paramsN1.N() != paramsN2.N() && (evk.EvkN1ToN2 == nil || evk.EvkN2ToN1 == nil) {
			return nil, fmt.Errorf("cannot NewBootstrapper: evk.(BootstrappingKeys) is missing EvkN1ToN2 and EvkN2ToN1")
		}
	case ring.ConjugateInvariant:
		if evk.EvkCmplxToReal == nil || evk.EvkRealToCmplx == nil {
			return nil, fmt.Errorf("cannot NewBootstrapper: evk.(BootstrappingKeys) is missing EvkN1ToN2 and EvkN2ToN1")
		}

		var err error
		if b.DomainSwitcher, err = ckks.NewDomainSwitcher(paramsN2.Parameters, evk.EvkCmplxToReal, evk.EvkRealToCmplx); err != nil {
			return nil, fmt.Errorf("cannot NewBootstrapper: ckks.NewDomainSwitcher: %w", err)
		}

		// The switch to standard to conjugate invariant multiplies the scale by 2
		btpParams.SlotsToCoeffsParameters.Scaling = new(big.Float).SetFloat64(0.5)
	}

	b.Parameters = &btpParams

	if paramsN1.N() != paramsN2.N() {
		b.xPow2N1 = rlwe.GenXPow2(paramsN1.RingQ().AtLevel(0), paramsN2.LogN(), false)
		b.xPow2N2 = rlwe.GenXPow2(paramsN2.RingQ().AtLevel(0), paramsN2.LogN(), false)
		b.xPow2InvN2 = rlwe.GenXPow2(paramsN2.RingQ(), paramsN2.LogN(), true)
	}

	var err error
	if b.CoreBootstrapper, err = NewCoreBootstrapper(btpParams, evk); err != nil {
		return nil, err
	}

	return b, nil
}

// Depth returns the multiplicative depth (number of levels consumed) of the bootstrapping circuit.
func (b Bootstrapper) Depth() int {
	return b.BootstrappingParameters.MaxLevel() - b.ResidualParameters.MaxLevel()
}

// OutputLevel returns the output level after the evaluation of the bootstrapping circuit.
func (b Bootstrapper) OutputLevel() int {
	return b.ResidualParameters.MaxLevel()
}

// MinimumInputLevel returns the minimum level at which a ciphertext must be to be
// bootstrapped.
func (b Bootstrapper) MinimumInputLevel() int {
	return b.BootstrappingParameters.LevelsConsumedPerRescaling()
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
			if ct, err = b.CoreBootstrapper.Bootstrap(&cts[i]); err != nil {
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

		if err = b.CoreBootstrapper.Evaluator.Mul(ctRightN2Q0, 1i, ctRightN2Q0); err != nil {
			return nil, nil, fmt.Errorf("cannot BootstrapMany: %w", err)
		}

		if err = b.CoreBootstrapper.Evaluator.Add(ctLeftN2Q0, ctRightN2Q0, ctLeftN2Q0); err != nil {
			return nil, nil, fmt.Errorf("cannot BootstrapMany: %w", err)
		}
	}

	// Refreshes in the ring.Sstandard
	var ctLeftAndRightN2QL *rlwe.Ciphertext
	if ctLeftAndRightN2QL, err = b.CoreBootstrapper.Bootstrap(ctLeftN2Q0); err != nil {
		return nil, nil, fmt.Errorf("cannot BootstrapMany: %w", err)
	}

	// The SlotsToCoeffs transformation scales the ciphertext by 0.5
	// This is done to compensate for the 2x factor introduced by ringStandardToConjugate(*).
	ctLeftAndRightN2QL.Scale = ctLeftAndRightN2QL.Scale.Mul(rlwe.NewScale(1 / 2.0))

	// Switches ring from ring.Standard to ring.ConjugateInvariant
	ctLeftN1QL = b.ComplexToRealNew(ctLeftAndRightN2QL)

	// Extracts the imaginary part
	if ctRightN1Q0 != nil {
		if err = b.CoreBootstrapper.Evaluator.Mul(ctLeftAndRightN2QL, -1i, ctLeftAndRightN2QL); err != nil {
			return nil, nil, fmt.Errorf("cannot BootstrapMany: %w", err)
		}
		ctRightN1QL = b.ComplexToRealNew(ctLeftAndRightN2QL)
	}

	return
}

func (b Bootstrapper) SwitchRingDegreeN1ToN2New(ctN1 *rlwe.Ciphertext) (ctN2 *rlwe.Ciphertext) {
	ctN2 = hefloat.NewCiphertext(b.BootstrappingParameters, 1, ctN1.Level())

	// Sanity check, this error should never happen unless this algorithm has been improperly
	// modified to pass invalid inputs.
	if err := b.CoreBootstrapper.ApplyEvaluationKey(ctN1, b.EvkN1ToN2, ctN2); err != nil {
		panic(err)
	}
	return
}

func (b Bootstrapper) SwitchRingDegreeN2ToN1New(ctN2 *rlwe.Ciphertext) (ctN1 *rlwe.Ciphertext) {
	ctN1 = hefloat.NewCiphertext(b.ResidualParameters, 1, ctN2.Level())

	// Sanity check, this error should never happen unless this algorithm has been improperly
	// modified to pass invalid inputs.
	if err := b.CoreBootstrapper.ApplyEvaluationKey(ctN2, b.EvkN2ToN1, ctN1); err != nil {
		panic(err)
	}
	return
}

func (b Bootstrapper) ComplexToRealNew(ctCmplx *rlwe.Ciphertext) (ctReal *rlwe.Ciphertext) {
	ctReal = hefloat.NewCiphertext(b.ResidualParameters, 1, ctCmplx.Level())

	// Sanity check, this error should never happen unless this algorithm has been improperly
	// modified to pass invalid inputs.
	if err := b.DomainSwitcher.ComplexToReal(&b.CoreBootstrapper.Evaluator.Evaluator, ctCmplx, ctReal); err != nil {
		panic(err)
	}
	return
}

func (b Bootstrapper) RealToComplexNew(ctReal *rlwe.Ciphertext) (ctCmplx *rlwe.Ciphertext) {
	ctCmplx = hefloat.NewCiphertext(b.BootstrappingParameters, 1, ctReal.Level())

	// Sanity check, this error should never happen unless this algorithm has been improperly
	// modified to pass invalid inputs.
	if err := b.DomainSwitcher.RealToComplex(&b.CoreBootstrapper.Evaluator.Evaluator, ctReal, ctCmplx); err != nil {
		panic(err)
	}
	return
}

func (b Bootstrapper) PackAndSwitchN1ToN2(cts []rlwe.Ciphertext) ([]rlwe.Ciphertext, error) {

	var err error

	if b.ResidualParameters.N() != b.BootstrappingParameters.N() {
		if cts, err = b.Pack(cts, b.ResidualParameters, b.xPow2N1); err != nil {
			return nil, fmt.Errorf("cannot PackAndSwitchN1ToN2: PackN1: %w", err)
		}

		for i := range cts {
			cts[i] = *b.SwitchRingDegreeN1ToN2New(&cts[i])
		}
	}

	if cts, err = b.Pack(cts, b.BootstrappingParameters, b.xPow2N2); err != nil {
		return nil, fmt.Errorf("cannot PackAndSwitchN1ToN2: PackN2: %w", err)
	}

	return cts, nil
}

func (b Bootstrapper) UnpackAndSwitchN2Tn1(cts []rlwe.Ciphertext, LogSlots, Nb int) ([]rlwe.Ciphertext, error) {

	var err error

	if b.ResidualParameters.N() != b.BootstrappingParameters.N() {
		if cts, err = b.UnPack(cts, b.BootstrappingParameters, LogSlots, Nb, b.xPow2InvN2); err != nil {
			return nil, fmt.Errorf("cannot UnpackAndSwitchN2Tn1: UnpackN2: %w", err)
		}

		for i := range cts {
			cts[i] = *b.SwitchRingDegreeN2ToN1New(&cts[i])
		}
	}

	for i := range cts {
		cts[i].LogDimensions.Cols = LogSlots
	}

	return cts, nil
}

func (b Bootstrapper) UnPack(cts []rlwe.Ciphertext, params hefloat.Parameters, LogSlots, Nb int, xPow2Inv []ring.Poly) ([]rlwe.Ciphertext, error) {
	LogGap := params.LogMaxSlots() - LogSlots

	if LogGap == 0 {
		return cts, nil
	}

	cts = append(cts, make([]rlwe.Ciphertext, Nb-1)...)

	for i := 1; i < len(cts); i++ {
		cts[i] = *cts[0].CopyNew()
	}

	r := params.RingQ().AtLevel(cts[0].Level())

	N := len(cts)

	for i := 0; i < utils.Min(bits.Len64(uint64(N-1)), LogGap); i++ {

		step := 1 << (i + 1)

		for j := 0; j < N; j += step {

			for k := step >> 1; k < step; k++ {

				if (j + k) >= N {
					break
				}

				r.MulCoeffsMontgomery(cts[j+k].Value[0], xPow2Inv[i], cts[j+k].Value[0])
				r.MulCoeffsMontgomery(cts[j+k].Value[1], xPow2Inv[i], cts[j+k].Value[1])
			}
		}
	}

	return cts, nil
}

func (b Bootstrapper) Pack(cts []rlwe.Ciphertext, params hefloat.Parameters, xPow2 []ring.Poly) ([]rlwe.Ciphertext, error) {

	var LogSlots = cts[0].LogSlots()
	RingDegree := params.N()

	for i, ct := range cts {
		if N := ct.LogSlots(); N != LogSlots {
			return nil, fmt.Errorf("cannot Pack: cts[%d].PlaintextLogSlots()=%d != cts[0].PlaintextLogSlots=%d", i, N, LogSlots)
		}

		if N := ct.Value[0].N(); N != RingDegree {
			return nil, fmt.Errorf("cannot Pack: cts[%d].Value[0].N()=%d != params.N()=%d", i, N, RingDegree)
		}
	}

	LogGap := params.LogMaxSlots() - LogSlots

	if LogGap == 0 {
		return cts, nil
	}

	for i := 0; i < LogGap; i++ {

		for j := 0; j < len(cts)>>1; j++ {

			eve := cts[j*2+0]
			odd := cts[j*2+1]

			level := utils.Min(eve.Level(), odd.Level())

			r := params.RingQ().AtLevel(level)

			r.MulCoeffsMontgomeryThenAdd(odd.Value[0], xPow2[i], eve.Value[0])
			r.MulCoeffsMontgomeryThenAdd(odd.Value[1], xPow2[i], eve.Value[1])

			cts[j] = eve
		}

		if len(cts)&1 == 1 {
			cts[len(cts)>>1] = cts[len(cts)-1]
			cts = cts[:len(cts)>>1+1]
		} else {
			cts = cts[:len(cts)>>1]
		}
	}

	LogMaxDimensions := params.LogMaxDimensions()
	for i := range cts {
		cts[i].LogDimensions = LogMaxDimensions
	}

	return cts, nil
}
