package bootstrapper

import (
	"fmt"

	"github.com/tuneinsight/lattigo/v4/circuits/float/bootstrapper/bootstrapping"
	"github.com/tuneinsight/lattigo/v4/ckks"
	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/rlwe"
)

type BootstrappingKeys struct {
	EvkN1ToN2        *rlwe.EvaluationKey
	EvkN2ToN1        *rlwe.EvaluationKey
	EvkRealToCmplx   *rlwe.EvaluationKey
	EvkCmplxToReal   *rlwe.EvaluationKey
	EvkBootstrapping *bootstrapping.EvaluationKeySet
}

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

func (p Parameters) GenBootstrappingKeys(paramsN1 ckks.Parameters, skN1 *rlwe.SecretKey) (*BootstrappingKeys, error) {

	var EvkN1ToN2, EvkN2ToN1 *rlwe.EvaluationKey
	var EvkRealToCmplx *rlwe.EvaluationKey
	var EvkCmplxToReal *rlwe.EvaluationKey
	paramsN2 := p.Parameters.Parameters

	// Checks that the maximum level of paramsN1 is equal to the remaining level after the bootstrapping of paramsN2
	if paramsN2.MaxLevel()-p.SlotsToCoeffsParameters.Depth(true)-p.Mod1ParametersLiteral.Depth()-p.CoeffsToSlotsParameters.Depth(true) < paramsN1.MaxLevel() {
		return nil, fmt.Errorf("cannot GenBootstrappingKeys: bootstrapping depth is too large, level after bootstrapping is smaller than paramsN1.MaxLevel()")
	}

	// Checks that the overlapping primes between paramsN1 and paramsN2 are the same, i.e.
	// pN1: q0, q1, q2, ..., qL
	// pN2: q0, q1, q2, ..., qL, [bootstrapping primes]
	QN1 := paramsN1.Q()
	QN2 := paramsN2.Q()

	for i := range QN1 {
		if QN1[i] != QN2[i] {
			return nil, fmt.Errorf("cannot GenBootstrappingKeys: paramsN1.Q() is not a subset of paramsN2.Q()")
		}
	}

	kgen := ckks.NewKeyGenerator(paramsN2)

	// Ephemeral secret-key used to generate the evaluation keys.
	skN2 := rlwe.NewSecretKey(paramsN2)
	buff := paramsN2.RingQ().NewPoly()
	ringQ := paramsN2.RingQ()
	ringP := paramsN2.RingP()

	switch paramsN1.RingType() {
	// In this case we need need generate the bridge switching keys between the two rings
	case ring.ConjugateInvariant:

		if paramsN1.LogN() != paramsN2.LogN()-1 {
			return nil, fmt.Errorf("cannot GenBootstrappingKeys: if paramsN1.RingType() == ring.ConjugateInvariant then must ensure that paramsN1.LogN()+1 == paramsN2.LogN()-1")
		}

		// R[X+X^-1]/(X^N +1) -> R[X]/(X^2N + 1)
		ringQ.AtLevel(skN1.LevelQ()).UnfoldConjugateInvariantToStandard(skN1.Value.Q, skN2.Value.Q)

		// Extends basis Q0 -> QL
		rlwe.ExtendBasisSmallNormAndCenterNTTMontgomery(ringQ, ringQ, skN2.Value.Q, buff, skN2.Value.Q)

		// Extends basis Q0 -> P
		rlwe.ExtendBasisSmallNormAndCenterNTTMontgomery(ringQ, ringP, skN2.Value.Q, buff, skN2.Value.P)

		EvkCmplxToReal, EvkRealToCmplx = kgen.GenEvaluationKeysForRingSwapNew(skN2, skN1)

	// Only regular key-switching is required in this case
	case ring.Standard:

		// Maps the smaller key to the largest with Y = X^{N/n}.
		ring.MapSmallDimensionToLargerDimensionNTT(skN1.Value.Q, skN2.Value.Q)

		// Extends basis Q0 -> QL
		rlwe.ExtendBasisSmallNormAndCenterNTTMontgomery(ringQ, ringQ, skN2.Value.Q, buff, skN2.Value.Q)

		// Extends basis Q0 -> P
		rlwe.ExtendBasisSmallNormAndCenterNTTMontgomery(ringQ, ringP, skN2.Value.Q, buff, skN2.Value.P)

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
