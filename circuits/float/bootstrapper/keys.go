package bootstrapper

import (
	"fmt"

	"github.com/tuneinsight/lattigo/v4/circuits/float/bootstrapper/bootstrapping"
	"github.com/tuneinsight/lattigo/v4/ckks"
	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/rlwe"
)

// BootstrappingKeys is a struct storing the different
// evaluation keys required by the bootstrapper.
type BootstrappingKeys struct {
	// EvkN1ToN2 is an evaluation key to switch from the residual parameters'
	// ring degree (N1) to the bootstrapping parameters' ring degre (N2)
	EvkN1ToN2 *rlwe.EvaluationKey
	// EvkN2ToN1 is an evaluation key to switch from the bootstrapping parameters'
	// ring degre (N2) to the residual parameters' ring degree (N1)
	EvkN2ToN1 *rlwe.EvaluationKey
	// EvkRealToCmplx is an evaluation key to switch from the standard ring to the
	// conjugate invariant ring.
	EvkRealToCmplx *rlwe.EvaluationKey
	// EvkCmplxToReal is an evaluation key to switch from the conjugate invariant
	// ring to the standard ring.
	EvkCmplxToReal *rlwe.EvaluationKey
	// EvkBootstrapping is a set of evaluation keys for the bootstraping circuit.
	EvkBootstrapping *bootstrapping.EvaluationKeySet
}

// BinarySize returns the total binary size of the bootstrapper's keys.
func (b BootstrappingKeys) BinarySize() (dLen int) {
	if b.EvkN1ToN2 != nil {
		dLen += b.EvkN1ToN2.BinarySize()
	}

	if b.EvkN2ToN1 != nil {
		dLen += b.EvkN2ToN1.BinarySize()
	}

	if b.EvkRealToCmplx != nil {
		dLen += b.EvkRealToCmplx.BinarySize()
	}

	if b.EvkCmplxToReal != nil {
		dLen += b.EvkCmplxToReal.BinarySize()
	}

	if b.EvkBootstrapping != nil {
		dLen += b.EvkBootstrapping.BinarySize()
	}

	return
}

// GenBootstrappingKeys generates the bootstrapping keys, which include:
// - If the bootstrapping parameters' ring degree > residual parameters' ring degree:
//   - An evaluation key to switch from the residual parameters' ring to the bootstrapping parameters' ring
//   - An evaluation key to switch from the bootstrapping parameters' ring to the residual parameters' ring
//
// - If the residual parameters use the Conjugate Invariant ring:
//   - An evaluation key to switch from the conjugate invariant ring to the standard ring
//   - An evaluation key to switch from the standard ring to the conjugate invariant ring
//
// - The bootstrapping evaluation keys:
//   - Relinearization key
//   - Galois keys
//   - The encapsulation evaluation keys (https://eprint.iacr.org/2022/024)
//
// Note: These evaluation keys are generated under an ephemeral secret key using the distribution
//
//	specified in the  bootstrapping parameters.
func (p Parameters) GenBootstrappingKeys(skN1 *rlwe.SecretKey) (*BootstrappingKeys, error) {

	var EvkN1ToN2, EvkN2ToN1 *rlwe.EvaluationKey
	var EvkRealToCmplx *rlwe.EvaluationKey
	var EvkCmplxToReal *rlwe.EvaluationKey
	paramsN2 := p.Parameters.Parameters

	var skN2 *rlwe.SecretKey
	kgen := ckks.NewKeyGenerator(paramsN2)

	// Ephemeral secret-key used to generate the evaluation keys.

	ringQ := paramsN2.RingQ()
	ringP := paramsN2.RingP()

	switch p.ResidualParameters.RingType() {
	// In this case we need need generate the bridge switching keys between the two rings
	case ring.ConjugateInvariant:

		if skN1.Value.Q.N() != paramsN2.N()>>1 {
			return nil, fmt.Errorf("cannot GenBootstrappingKeys: if paramsN1.RingType() == ring.ConjugateInvariant then must ensure that paramsN1.LogN()+1 == paramsN2.LogN()-1")
		}

		skN2 = rlwe.NewSecretKey(paramsN2)
		buff := paramsN2.RingQ().NewPoly()

		// R[X+X^-1]/(X^N +1) -> R[X]/(X^2N + 1)
		ringQ.AtLevel(skN1.LevelQ()).UnfoldConjugateInvariantToStandard(skN1.Value.Q, skN2.Value.Q)

		// Extends basis Q0 -> QL
		rlwe.ExtendBasisSmallNormAndCenterNTTMontgomery(ringQ, ringQ, skN2.Value.Q, buff, skN2.Value.Q)

		// Extends basis Q0 -> P
		rlwe.ExtendBasisSmallNormAndCenterNTTMontgomery(ringQ, ringP, skN2.Value.Q, buff, skN2.Value.P)

		EvkCmplxToReal, EvkRealToCmplx = kgen.GenEvaluationKeysForRingSwapNew(skN2, skN1)

	// Only regular key-switching is required in this case
	case ring.Standard:

		if skN1.Value.Q.N() == paramsN2.N() {
			skN2 = skN1
		} else {
			skN2 = kgen.GenSecretKeyNew()
		}

		EvkN1ToN2 = kgen.GenEvaluationKeyNew(skN1, skN2)
		EvkN2ToN1 = kgen.GenEvaluationKeyNew(skN2, skN1)
	}

	return &BootstrappingKeys{
		EvkN1ToN2:        EvkN1ToN2,
		EvkN2ToN1:        EvkN2ToN1,
		EvkRealToCmplx:   EvkRealToCmplx,
		EvkCmplxToReal:   EvkCmplxToReal,
		EvkBootstrapping: p.Parameters.GenEvaluationKeySetNew(skN2),
	}, nil
}
