// Package bootstrapper implements the Bootstrapper struct which provides generic bootstrapping for the CKKS scheme (and RLWE ciphertexts by extension).
// It notably abstracts scheme switching and ring dimension switching, enabling efficient bootstrapping of ciphertexts in the Conjugate Invariant ring
// or multiple ciphertexts of a lower ring dimension.
package bootstrapper

import (
	"fmt"
	"math/big"
	"runtime"

	"github.com/tuneinsight/lattigo/v4/ckks"
	"github.com/tuneinsight/lattigo/v4/ckks/bootstrapper/bootstrapping"
	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/rlwe"
	"github.com/tuneinsight/lattigo/v4/utils"
)

type Bootstrapper struct {
	bridge       *ckks.DomainSwitcher
	bootstrapper *bootstrapping.Bootstrapper

	paramsN1    *ckks.Parameters
	paramsN2    ckks.Parameters
	btpParamsN2 bootstrapping.Parameters

	evk BootstrappingKeys
}

type BootstrappingKeys struct {
	SwkN1toN2      *rlwe.EvaluationKey
	SwkN2toN1      *rlwe.EvaluationKey
	SwkReals2Cmplx *rlwe.EvaluationKey
	SwkCmplx2Reals *rlwe.EvaluationKey
	BtpKeys        *bootstrapping.EvaluationKeySet
}

func (b BootstrappingKeys) BinarySize() (dLen int) {
	if b.SwkN1toN2 != nil {
		dLen += b.SwkN1toN2.BinarySize()
	}

	if b.SwkN2toN1 != nil {
		dLen += b.SwkN2toN1.BinarySize()
	}

	if b.SwkReals2Cmplx != nil {
		dLen += b.SwkReals2Cmplx.BinarySize()
	}

	if b.SwkCmplx2Reals != nil {
		dLen += b.SwkCmplx2Reals.BinarySize()
	}

	if b.BtpKeys != nil {
		dLen += b.BtpKeys.BinarySize()
	}

	return
}

func GenBootstrappingKeys(paramsN1 *ckks.Parameters, paramsN2 ckks.Parameters, btpParamsN2 bootstrapping.Parameters, skN1 *rlwe.SecretKey, skN2 rlwe.SecretKey) (BootstrappingKeys, error) {

	var swkN1toN2, swkN2toN1 *rlwe.EvaluationKey
	var SwkReals2Cmplx *rlwe.EvaluationKey
	var SwkCmplx2Reals *rlwe.EvaluationKey
	if paramsN1 != nil {

		// Checks that the maximum level of paramsN1 is equal to the remaining level after the bootstrapping of paramsN2
		if paramsN2.MaxLevel()-btpParamsN2.SlotsToCoeffsParameters.Depth(true)-btpParamsN2.EvalModParameters.Depth()-btpParamsN2.CoeffsToSlotsParameters.Depth(true) < paramsN1.MaxLevel() {
			return BootstrappingKeys{}, fmt.Errorf("GenBootstrappingKeys(*): bootstrapping depth is too large, level after bootstrapping is smaller than paramsN1.MaxLevel()")
		}

		// Checks that the overlapping primes between paramsN1 and paramsN2 are the same, i.e.
		// pN1: q0, q1, q2, ..., qL
		// pN2: q0, q1, q2, ..., qL, [bootstrapping primes]
		QN1 := paramsN1.Q()
		QN2 := paramsN2.Q()

		for i := range QN1 {
			if QN1[i] != QN2[i] {
				return BootstrappingKeys{}, fmt.Errorf("GenBootstrappingKeys(*): paramsN1.Q() is not a subset of paramsN2.Q()")
			}
		}

		kgen := ckks.NewKeyGenerator(paramsN2)

		switch paramsN1.RingType() {
		// In this case we need need generate the bridge switching keys between the two rings
		case ring.ConjugateInvariant:

			if paramsN1.LogN() != paramsN2.LogN()-1 {
				return BootstrappingKeys{}, fmt.Errorf("GenBootstrappingKeys(*): if paramsN1.RingType() == ring.ConjugateInvariant then must ensure that paramsN1.LogN()+1 == paramsN2.LogN()-1")
			}

			SwkCmplx2Reals, SwkReals2Cmplx = kgen.GenEvaluationKeysForRingSwapNew(&skN2, skN1)

		// Only regular key-switching is required in this case
		case ring.Standard:
			swkN1toN2 = kgen.GenEvaluationKeyNew(skN1, &skN2)
			swkN2toN1 = kgen.GenEvaluationKeyNew(&skN2, skN1)
		}
	}

	return BootstrappingKeys{
		SwkN1toN2:      swkN1toN2,
		SwkN2toN1:      swkN2toN1,
		SwkReals2Cmplx: SwkReals2Cmplx,
		SwkCmplx2Reals: SwkCmplx2Reals,
		BtpKeys:        bootstrapping.GenEvaluationKeySetNew(btpParamsN2, paramsN2, &skN2),
	}, nil
}

func NewBootstrapper(paramsN1 *ckks.Parameters, paramsN2 ckks.Parameters, btpParamsN2 bootstrapping.Parameters, evk BootstrappingKeys) (rlwe.Bootstrapper, error) {

	b := &Bootstrapper{}

	if paramsN1 != nil {
		switch paramsN1.RingType() {
		case ring.Standard:
			if evk.SwkN1toN2 == nil || evk.SwkN2toN1 == nil {
				return nil, fmt.Errorf("NewBootstrapper(*): evk.(BootstrappingKeys) is missing SwkN1toN2 and SwkN2toN1")
			}

		case ring.ConjugateInvariant:
			if evk.SwkCmplx2Reals == nil || evk.SwkReals2Cmplx == nil {
				return nil, fmt.Errorf("NewBootstrapper(*): evk.(BootstrappingKeys) is missing SwkN1toN2 and SwkN2toN1")
			}

			bridge, _ := ckks.NewDomainSwitcher(paramsN2, evk.SwkCmplx2Reals, evk.SwkReals2Cmplx)

			b.bridge = &bridge

			// The switch to standard to conjugate invariant multiplies the scale by 2
			btpParamsN2.SlotsToCoeffsParameters.Scaling = new(big.Float).SetFloat64(0.5)
		}
	}

	if paramsN1 == nil && (evk.SwkN1toN2 != nil || evk.SwkN2toN1 != nil) && (evk.SwkCmplx2Reals != nil || evk.SwkReals2Cmplx != nil) {
		return nil, fmt.Errorf("NewBootstrapper(*): missing argument paramsN1")
	}

	b.paramsN1 = paramsN1
	b.paramsN2 = paramsN2
	b.btpParamsN2 = btpParamsN2
	b.evk = evk

	var err error
	if b.bootstrapper, err = bootstrapping.NewBootstrapper(paramsN2, btpParamsN2, evk.BtpKeys); err != nil {
		return nil, err
	}

	return b, nil
}

func (b Bootstrapper) Depth() int {
	return b.btpParamsN2.SlotsToCoeffsParameters.Depth(true) + b.btpParamsN2.EvalModParameters.Depth() + b.btpParamsN2.CoeffsToSlotsParameters.Depth(true)
}

func (b Bootstrapper) OutputLevel() int {
	return b.paramsN2.MaxLevel() - b.Depth()
}

func (b Bootstrapper) MinimumInputLevel() int {
	return 0
}

// RingN1toRingN2 switches the ringdegree of ctN1 to the ring degree of ctN2.
func (b Bootstrapper) ringN1toRingN2(eval *rlwe.Evaluator, ctN1, ctN2 *rlwe.Ciphertext) {
	eval.ApplyEvaluationKey(ctN1, b.evk.SwkN1toN2, ctN2)
}

func (b Bootstrapper) ringStandardToConjugate(eval *rlwe.Evaluator, ctIn, ctOut *rlwe.Ciphertext) {
	level := utils.Min(ctIn.Level(), ctOut.Level())

	tmp := &rlwe.Ciphertext{}
	tmp.Value = []ring.Poly{eval.BuffQP[1].Q, eval.BuffQP[2].Q}
	tmp.IsNTT = true

	eval.GadgetProduct(level, ctIn.Value[1], &b.evk.SwkCmplx2Reals.GadgetCiphertext, tmp)
	b.paramsN2.RingQ().AtLevel(level).Add(eval.BuffQP[1].Q, ctIn.Value[0], eval.BuffQP[1].Q)

	b.paramsN1.RingQ().AtLevel(level).FoldStandardToConjugateInvariant(eval.BuffQP[1].Q, eval.AutomorphismIndex[uint64(2*b.paramsN2.N()-1)], ctOut.Value[0])
	b.paramsN1.RingQ().AtLevel(level).FoldStandardToConjugateInvariant(eval.BuffQP[2].Q, eval.AutomorphismIndex[uint64(2*b.paramsN2.N()-1)], ctOut.Value[1])
	ctOut.PlaintextScale = ctIn.PlaintextScale.Mul(rlwe.NewScale(2))
}

func (b Bootstrapper) ringConjugateToRingStandard(eval *rlwe.Evaluator, ctIn, ctOut *rlwe.Ciphertext) {

	level := utils.Min(ctIn.Level(), ctOut.Level())

	b.paramsN2.RingQ().AtLevel(level).UnfoldConjugateInvariantToStandard(ctIn.Value[0], ctOut.Value[0])
	b.paramsN2.RingQ().AtLevel(level).UnfoldConjugateInvariantToStandard(ctIn.Value[1], ctOut.Value[1])

	tmp := &rlwe.Ciphertext{}
	tmp.Value = []ring.Poly{eval.BuffQP[1].Q, eval.BuffQP[2].Q}
	tmp.IsNTT = true

	// Switches the RCKswitcher key [X+X^-1] to a CKswitcher key [X]
	eval.GadgetProduct(level, ctOut.Value[1], &b.evk.SwkReals2Cmplx.GadgetCiphertext, tmp)
	b.paramsN2.RingQ().AtLevel(level).Add(ctOut.Value[0], eval.BuffQP[1].Q, ctOut.Value[0])
	ring.Copy(eval.BuffQP[2].Q, ctOut.Value[1])
}

// RingN2toRingN1 switches the ringdegree of ctN2 to the ring degree of ctN1.
func (b Bootstrapper) ringN2toRingN1(eval *rlwe.Evaluator, ctN2, ctN1 *rlwe.Ciphertext) {
	eval.ApplyEvaluationKey(ctN2, b.evk.SwkN2toN1, ctN1)
}

func (b Bootstrapper) Bootstrap(ct *rlwe.Ciphertext) (*rlwe.Ciphertext, error) {
	cts := []*rlwe.Ciphertext{ct}
	cts, err := b.BootstrapMany(cts)
	return cts[0], err
}

func (b Bootstrapper) BootstrapMany(cts []*rlwe.Ciphertext) ([]*rlwe.Ciphertext, error) {

	var err error

	if b.evk.SwkN1toN2 != nil {
		for i := range cts {
			cts[i] = b.refreshRingDegreeSwitch(cts[i])
		}

	} else if b.evk.SwkReals2Cmplx != nil {

		for i := 0; i < len(cts); i = i + 2 {

			even, odd := i, i+1

			ct0 := cts[even]

			var ct1 *rlwe.Ciphertext
			if odd < len(cts) {
				ct1 = cts[odd]
			}

			ct0, ct1 = b.refreshConjugateInvariant(ct0, ct1)

			if ct0 != nil {
				cts[even] = ct0
			}

			if ct1 != nil {
				cts[odd] = ct1
			}
		}

	} else {

		for i := range cts {
			cts[i] = b.refreshStandard(cts[i])
		}
	}

	runtime.GC()

	for i := range cts {
		cts[i].PlaintextScale = b.paramsN1.PlaintextScale()
	}

	return cts, err
}

func (b Bootstrapper) refreshStandard(ctN2Q0 *rlwe.Ciphertext) (ctN2QL *rlwe.Ciphertext) {
	return b.bootstrapper.Bootstrap(ctN2Q0)
}

func (b Bootstrapper) refreshRingDegreeSwitch(ctN1Q0 *rlwe.Ciphertext) (ctN1QL *rlwe.Ciphertext) {

	if ctN1Q0 != nil {
		// Switches ring from N1 to N2
		ctLeftN2Q0 := ckks.NewCiphertext(b.paramsN2, 1, ctN1Q0.Level())
		ctLeftN2Q0.PlaintextScale = ctN1Q0.PlaintextScale
		b.ringN1toRingN2(b.bootstrapper.Evaluator.Evaluator, ctN1Q0, ctLeftN2Q0)

		// Refreshes in the ring N2
		ctN2QL := b.bootstrapper.Bootstrap(ctLeftN2Q0)

		// Switches ring from N2 to N1
		ctN1QL = ckks.NewCiphertext(*b.paramsN1, 1, b.paramsN1.MaxLevel())
		ctN1QL.PlaintextScale = ctN2QL.PlaintextScale
		b.ringN2toRingN1(b.bootstrapper.Evaluator.Evaluator, ctN2QL, ctN1QL)
	}

	return
}

// refreshConjugateInvariant takes two ciphertext in the Conjugate Invariant ring, repacks them in a single ciphertext in the standard ring
// using the real and imaginary part, bootstrap both ciphertext, and then extract back the real and imaginary part before repacking them
// individually in two new ciphertexts in the Conjugate Invariant ring.
func (b Bootstrapper) refreshConjugateInvariant(ctLeftN1Q0, ctRightN1Q0 *rlwe.Ciphertext) (ctLeftN1QL, ctRightN1QL *rlwe.Ciphertext) {

	var ctLeftN2Q0 *rlwe.Ciphertext
	if ctLeftN1Q0 != nil {
		// Switches ring from ring.ConjugateInvariant to ring.Standard
		ctLeftN2Q0 = ckks.NewCiphertext(b.paramsN2, 1, ctLeftN1Q0.Level())
		ctLeftN2Q0.PlaintextScale = ctLeftN1Q0.PlaintextScale
		b.ringConjugateToRingStandard(b.bootstrapper.Evaluator.Evaluator, ctLeftN1Q0, ctLeftN2Q0)
	}

	// Repacks ctRightN1Q0 into the imaginary part of ctLeftN1Q0
	// which is zero since it comes from the Conjugate Invariant ring)
	if ctRightN1Q0 != nil {
		ctRightN2Q0 := ckks.NewCiphertext(b.paramsN2, 1, ctRightN1Q0.Level())
		ctRightN2Q0.PlaintextScale = ctRightN1Q0.PlaintextScale
		b.ringConjugateToRingStandard(b.bootstrapper.Evaluator.Evaluator, ctRightN1Q0, ctRightN2Q0)

		if ctLeftN1Q0 != nil {
			b.bootstrapper.Evaluator.Mul(ctRightN2Q0, 1i, ctRightN2Q0)
			b.bootstrapper.Evaluator.Add(ctLeftN2Q0, ctRightN2Q0, ctLeftN2Q0)
		} else {
			ctLeftN2Q0 = ctRightN2Q0
		}
	}

	// Refreshes in the ring.Sstandard
	ctLeftAndRightN2QL := b.bootstrapper.Bootstrap(ctLeftN2Q0)

	// The SlotsToCoeffs transformation scales the ciphertext by 0.5
	// This is done to compensate for the 2x factor introduced by ringStandardToConjugate(*).
	ctLeftAndRightN2QL.PlaintextScale = ctLeftAndRightN2QL.PlaintextScale.Mul(rlwe.NewScale(1 / 2.0))

	if ctLeftN1Q0 != nil {
		// Switches ring from ring.Standard to ring.ConjugateInvariant
		ctLeftN1QL = ckks.NewCiphertext(*b.paramsN1, 1, b.paramsN1.MaxLevel())
		ctLeftN1QL.PlaintextScale = ctLeftAndRightN2QL.PlaintextScale
		b.ringStandardToConjugate(b.bootstrapper.Evaluator.Evaluator, ctLeftAndRightN2QL, ctLeftN1QL)
	}

	// Also extracts the imaginary part
	if ctRightN1Q0 != nil {

		ctRightN1QL = ckks.NewCiphertext(*b.paramsN1, 1, b.paramsN1.MaxLevel())
		ctRightN1QL.PlaintextScale = ctLeftAndRightN2QL.PlaintextScale

		// If the left ct is not nil, then extract as the imaginary part
		if ctLeftN1Q0 != nil {
			b.bootstrapper.Mul(ctLeftAndRightN2QL, -1i, ctLeftAndRightN2QL)
		}

		b.ringStandardToConjugate(b.bootstrapper.Evaluator.Evaluator, ctLeftAndRightN2QL, ctRightN1QL)
	}

	return
}
