package bootstrapping

import (
	"github.com/tuneinsight/lattigo/v5/core/rlwe"
	"github.com/tuneinsight/lattigo/v5/ring"
)

// EvaluationKeys is a struct storing the different
// evaluation keys required by the bootstrapper.
type EvaluationKeys struct {

	// EvkN1ToN2 is the evaluation key to switch from the residual parameters'
	// ring degree (N1) to the bootstrapping parameters' ring degree (N2)
	EvkN1ToN2 *rlwe.EvaluationKey

	// EvkN2ToN1 is the evaluation key to switch from the bootstrapping parameters'
	// ring degree (N2) to the residual parameters' ring degree (N1)
	EvkN2ToN1 *rlwe.EvaluationKey

	// EvkRealToCmplx is the evaluation key to switch from the standard ring to the
	// conjugate invariant ring.
	EvkRealToCmplx *rlwe.EvaluationKey

	// EvkCmplxToReal is the evaluation key to switch from the conjugate invariant
	// ring to the standard ring.
	EvkCmplxToReal *rlwe.EvaluationKey

	// EvkDenseToSparse is the evaluation key to switch
	// from the dense secret to the sparse secret.
	// https://eprint.iacr.org/2022/024
	EvkDenseToSparse *rlwe.EvaluationKey

	// EvkSparseToDense is the evaluation key to switch
	// from the sparse secret to the dense secret.
	// https://eprint.iacr.org/2022/024
	EvkSparseToDense *rlwe.EvaluationKey

	// MemEvaluationKeySet is the evaluation key set storing the relinearization
	// key and the Galois keys necessary for the bootstrapping circuit.
	*rlwe.MemEvaluationKeySet
}

// BinarySize returns the total binary size of the bootstrapper's keys.
func (b EvaluationKeys) BinarySize() (dLen int) {
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

	if b.EvkDenseToSparse != nil {
		dLen += b.EvkDenseToSparse.BinarySize()
	}

	if b.EvkSparseToDense != nil {
		dLen += b.EvkSparseToDense.BinarySize()
	}

	if b.MemEvaluationKeySet != nil {
		dLen += b.MemEvaluationKeySet.BinarySize()
	}

	return
}

// GenEvaluationKeys generates the bootstrapping evaluation keys, which include:
//
// - If the bootstrapping parameters' ring degree > residual parameters' ring degree:
//   - An evaluation key to switch from the residual parameters' ring to the bootstrapping parameters' ring
//   - An evaluation key to switch from the bootstrapping parameters' ring to the residual parameters' ring
//
// - If the residual parameters use the Conjugate Invariant ring:
//   - An evaluation key to switch from the conjugate invariant ring to the standard ring
//   - An evaluation key to switch from the standard ring to the conjugate invariant ring
//
// - The core bootstrapping circuit evaluation keys:
//   - Relinearization key
//   - Galois keys
//   - The encapsulation evaluation keys (https://eprint.iacr.org/2022/024)
//
// Note:
//   - These evaluation keys are generated under an ephemeral secret key skN2 using the distribution
//     specified in the bootstrapping parameters.
//   - The ephemeral key used to generate the bootstrapping keys is returned by this method for debugging purposes.
//   - !WARNING! The bootstrapping parameters use their own and independent cryptographic parameters (i.e. float.Parameters)
//     and it is the user's responsibility to ensure that these parameters meet the target security and tweak them if necessary.
func (p Parameters) GenEvaluationKeys(skN1 *rlwe.SecretKey) (btpkeys *EvaluationKeys, skN2 *rlwe.SecretKey, err error) {

	var EvkN1ToN2, EvkN2ToN1 *rlwe.EvaluationKey
	var EvkRealToCmplx *rlwe.EvaluationKey
	var EvkCmplxToReal *rlwe.EvaluationKey
	paramsN2 := p.BootstrappingParameters

	kgen := rlwe.NewKeyGenerator(paramsN2)

	if p.ResidualParameters.N() != paramsN2.N() {
		// If the ring degree do not match
		// (if the residual parameters are Conjugate Invariant, N1 = N2/2)
		skN2 = kgen.GenSecretKeyNew()

		if p.ResidualParameters.RingType() == ring.ConjugateInvariant {
			EvkCmplxToReal, EvkRealToCmplx = kgen.GenEvaluationKeysForRingSwapNew(skN2, skN1)
		} else {
			EvkN1ToN2 = kgen.GenEvaluationKeyNew(skN1, skN2)
			EvkN2ToN1 = kgen.GenEvaluationKeyNew(skN2, skN1)
		}

	} else {

		ringQ := paramsN2.RingQ()
		ringP := paramsN2.RingP()

		// Else, keeps the same secret, but extends to the full modulus of the bootstrapping parameters.
		skN2 = rlwe.NewSecretKey(paramsN2)
		buff := ringQ.NewPoly()

		// Extends basis Q0 -> QL
		rlwe.ExtendBasisSmallNormAndCenterNTTMontgomery(ringQ, ringQ, skN1.Value.Q, buff, skN2.Value.Q)

		// Extends basis Q0 -> P
		rlwe.ExtendBasisSmallNormAndCenterNTTMontgomery(ringQ, ringP, skN1.Value.Q, buff, skN2.Value.P)
	}

	EvkDenseToSparse, EvkSparseToDense := p.genEncapsulationEvaluationKeysNew(skN2)

	rlk := kgen.GenRelinearizationKeyNew(skN2)
	gks := kgen.GenGaloisKeysNew(append(p.GaloisElements(paramsN2), paramsN2.GaloisElementForComplexConjugation()), skN2)

	return &EvaluationKeys{
		EvkN1ToN2:           EvkN1ToN2,
		EvkN2ToN1:           EvkN2ToN1,
		EvkRealToCmplx:      EvkRealToCmplx,
		EvkCmplxToReal:      EvkCmplxToReal,
		MemEvaluationKeySet: rlwe.NewMemEvaluationKeySet(rlk, gks...),
		EvkDenseToSparse:    EvkDenseToSparse,
		EvkSparseToDense:    EvkSparseToDense,
	}, skN2, nil
}

// GenEncapsulationEvaluationKeysNew generates the low level encapsulation EvaluationKeys for the bootstrapping.
func (p Parameters) genEncapsulationEvaluationKeysNew(skDense *rlwe.SecretKey) (EvkDenseToSparse, EvkSparseToDense *rlwe.EvaluationKey) {

	params := p.BootstrappingParameters

	if p.EphemeralSecretWeight == 0 {
		return
	}

	paramsSparse, _ := rlwe.NewParametersFromLiteral(rlwe.ParametersLiteral{
		LogN: params.LogN(),
		Q:    params.Q()[:1],
		P:    params.P()[:1],
	})

	kgenSparse := rlwe.NewKeyGenerator(paramsSparse)
	kgenDense := rlwe.NewKeyGenerator(params)
	skSparse := kgenSparse.GenSecretKeyWithHammingWeightNew(p.EphemeralSecretWeight)

	EvkDenseToSparse = kgenDense.GenEvaluationKeyNew(skDense, skSparse)
	EvkSparseToDense = kgenDense.GenEvaluationKeyNew(skSparse, skDense)
	return
}
