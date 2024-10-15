package rlwe

import (
	"fmt"

	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/ring/ringqp"
	"github.com/tuneinsight/lattigo/v6/utils"
	"github.com/tuneinsight/lattigo/v6/utils/sampling"
)

// KeyGenerator is a structure that stores the elements required to create new keys,
// as well as a memory buffer for intermediate values.
type KeyGenerator struct {
	*Encryptor
	bufSkIn  ring.Poly
	bufSkOut ringqp.Poly
}

// NewKeyGenerator creates a new KeyGenerator, from which the secret and public keys, as well as [EvaluationKey].
func NewKeyGenerator(params ParameterProvider) *KeyGenerator {
	return &KeyGenerator{
		Encryptor: NewEncryptor(params, nil),
		bufSkIn:   params.GetRLWEParameters().RingQ().NewPoly(),
		bufSkOut:  params.GetRLWEParameters().RingQP().NewPoly(),
	}
}

// GenSecretKeyNew generates a new [SecretKey].
// Distribution is set according to [rlwe.Parameters.HammingWeight].
func (kgen KeyGenerator) GenSecretKeyNew() (sk *SecretKey) {
	sk = NewSecretKey(kgen.params)
	kgen.GenSecretKey(sk)
	return
}

// GenSecretKey generates a [SecretKey].
// Distribution is set according to [rlwe.Parameters.HammingWeight].
func (kgen KeyGenerator) GenSecretKey(sk *SecretKey) {
	kgen.genSecretKeyFromSampler(kgen.xsSampler, sk)
}

// GenSecretKeyWithHammingWeightNew generates a new [SecretKey] with exactly hw non-zero coefficients.
func (kgen *KeyGenerator) GenSecretKeyWithHammingWeightNew(hw int) (sk *SecretKey) {
	sk = NewSecretKey(kgen.params)
	kgen.GenSecretKeyWithHammingWeight(hw, sk)
	return
}

// GenSecretKeyWithHammingWeight generates a [SecretKey] with exactly hw non-zero coefficients.
func (kgen KeyGenerator) GenSecretKeyWithHammingWeight(hw int, sk *SecretKey) {
	Xs, err := ring.NewSampler(kgen.prng, kgen.params.RingQ(), ring.Ternary{H: hw}, false)
	// Sanity check, this error should not happen.
	if err != nil {
		panic(err)
	}
	kgen.genSecretKeyFromSampler(Xs, sk)
}

func (kgen KeyGenerator) genSecretKeyFromSampler(sampler ring.Sampler, sk *SecretKey) {

	ringQP := kgen.params.RingQP().AtLevel(sk.LevelQ(), sk.LevelP())

	sampler.AtLevel(sk.LevelQ()).Read(sk.Value.Q)

	if levelP := sk.LevelP(); levelP > -1 {
		ringQP.ExtendBasisSmallNormAndCenter(sk.Value.Q, levelP, sk.Value.Q, sk.Value.P)
	}

	ringQP.NTT(sk.Value, sk.Value)
	ringQP.MForm(sk.Value, sk.Value)
}

// GenPublicKeyNew generates a new [PublicKey] from the provided [SecretKey].
func (kgen KeyGenerator) GenPublicKeyNew(sk *SecretKey) (pk *PublicKey) {
	pk = NewPublicKey(kgen.params)
	kgen.GenPublicKey(sk, pk)
	return
}

// GenPublicKey generates a [PublicKey] from the provided [SecretKey].
func (kgen KeyGenerator) GenPublicKey(sk *SecretKey, pk *PublicKey) {
	if err := kgen.WithKey(sk).EncryptZero(Element[ringqp.Poly]{
		MetaData: &MetaData{CiphertextMetaData: CiphertextMetaData{IsNTT: true, IsMontgomery: true}},
		Value:    []ringqp.Poly(pk.Value),
	}); err != nil {
		// Sanity check, this error should not happen.
		panic(err)
	}
}

// GenKeyPairNew generates a new [SecretKey] and a corresponding [PublicKey].
// Distribution of the [SecretKey] set according to [rlwe.Parameters.HammingWeight].
func (kgen KeyGenerator) GenKeyPairNew() (sk *SecretKey, pk *PublicKey) {
	sk = kgen.GenSecretKeyNew()
	pk = kgen.GenPublicKeyNew(sk)
	return
}

// GenRelinearizationKeyNew generates a new [EvaluationKey] that will be used to relinearize [Ciphertexts] during multiplication.
func (kgen KeyGenerator) GenRelinearizationKeyNew(sk *SecretKey, evkParams ...EvaluationKeyParameters) (rlk *RelinearizationKey) {
	levelQ, levelP, BaseTwoDecomposition, compressed := ResolveEvaluationKeyParameters(kgen.params, evkParams)

	degree := 1
	if compressed {
		degree = 0
	}

	rlk = &RelinearizationKey{EvaluationKey: EvaluationKey{GadgetCiphertext: *NewGadgetCiphertext(kgen.params, degree, levelQ, levelP, BaseTwoDecomposition)}}
	kgen.GenRelinearizationKey(sk, rlk)
	return
}

// GenRelinearizationKey generates an [EvaluationKey] that will be used to relinearize [Ciphertexts] during multiplication.
func (kgen KeyGenerator) GenRelinearizationKey(sk *SecretKey, rlk *RelinearizationKey) {
	sk2 := kgen.bufSkIn
	sk2.CopyLvl(rlk.LevelQ(), sk.Value.Q)
	kgen.params.RingQ().AtLevel(rlk.LevelQ()).MulCoeffsMontgomery(sk2, sk.Value.Q, sk2)
	kgen.genEvaluationKey(sk2, sk.Value, &rlk.EvaluationKey)
}

// GenGaloisKeyNew generates a new [GaloisKey], enabling the automorphism X^{i} -> X^{i * galEl}.
func (kgen KeyGenerator) GenGaloisKeyNew(galEl uint64, sk *SecretKey, evkParams ...EvaluationKeyParameters) (gk *GaloisKey) {
	levelQ, levelP, BaseTwoDecomposition, compressed := ResolveEvaluationKeyParameters(kgen.params, evkParams)

	degree := 1
	if compressed {
		degree = 0
	}

	gk = &GaloisKey{
		EvaluationKey: EvaluationKey{GadgetCiphertext: *NewGadgetCiphertext(kgen.params, degree, levelQ, levelP, BaseTwoDecomposition)},
		NthRoot:       kgen.params.GetRLWEParameters().RingQ().NthRoot(),
	}
	kgen.GenGaloisKey(galEl, sk, gk)
	return
}

// GenGaloisKey generates a GaloisKey, enabling the automorphism X^{i} -> X^{i * galEl}.
func (kgen KeyGenerator) GenGaloisKey(galEl uint64, sk *SecretKey, gk *GaloisKey) {

	skIn := sk.Value
	skOut := kgen.bufSkOut

	ringQP := kgen.params.RingQP().AtLevel(gk.LevelQ(), gk.LevelP())

	ringQ := ringQP.RingQ
	ringP := ringQP.RingP

	// We encrypt [-a * pi_{k^-1}(sk) + sk, a]
	// This enables to first apply the gadget product, re-encrypting
	// a ciphetext from sk to pi_{k^-1}(sk) and then we apply pi_{k}
	// on the ciphertext.
	galElInv := kgen.params.ModInvGaloisElement(galEl)

	index, err := ring.AutomorphismNTTIndex(ringQ.N(), ringQ.NthRoot(), galElInv)

	// Sanity check, this error should not happen unless the
	// evaluator's buffer thave been improperly tempered with.
	if err != nil {
		panic(err)
	}

	ringQ.AutomorphismNTTWithIndex(skIn.Q, index, skOut.Q)

	if ringP != nil {
		ringP.AutomorphismNTTWithIndex(skIn.P, index, skOut.P)
	}

	kgen.genEvaluationKey(skIn.Q, skOut, &gk.EvaluationKey)

	gk.GaloisElement = galEl
	gk.NthRoot = ringQ.NthRoot()
}

// GenGaloisKeys generates the [GaloisKey] objects for all galois elements in galEls, and stores
// the resulting key for galois element i in gks[i].
// The galEls and gks parameters must have the same length.
func (kgen KeyGenerator) GenGaloisKeys(galEls []uint64, sk *SecretKey, gks []*GaloisKey) {

	// Sanity check
	if len(galEls) != len(gks) {
		panic(fmt.Errorf("galEls and gks must have the same length"))
	}

	for i, galEl := range galEls {
		if gks[i] == nil {
			gks[i] = kgen.GenGaloisKeyNew(galEl, sk)
		} else {
			kgen.GenGaloisKey(galEl, sk, gks[i])
		}
	}
}

// GenGaloisKeysNew generates the [GaloisKey] objects for all galois elements in galEls, and
// returns the resulting keys in a newly allocated []*[GaloisKey].
func (kgen KeyGenerator) GenGaloisKeysNew(galEls []uint64, sk *SecretKey, evkParams ...EvaluationKeyParameters) (gks []*GaloisKey) {
	levelQ, levelP, BaseTwoDecomposition, compressed := ResolveEvaluationKeyParameters(kgen.params, evkParams)
	gks = make([]*GaloisKey, len(galEls))
	for i, galEl := range galEls {
		gks[i] = newGaloisKey(kgen.params, levelQ, levelP, BaseTwoDecomposition, compressed)
		kgen.GenGaloisKey(galEl, sk, gks[i])
	}
	return
}

// GenEvaluationKeysForRingSwapNew generates the necessary evaluation keys to switch from a standard ring to to a conjugate invariant ring and vice-versa.
func (kgen KeyGenerator) GenEvaluationKeysForRingSwapNew(skStd, skConjugateInvariant *SecretKey, evkParams ...EvaluationKeyParameters) (stdToci, ciToStd *EvaluationKey) {

	levelQ := utils.Min(skStd.Value.Q.Level(), skConjugateInvariant.Value.Q.Level())

	skCIMappedToStandard := &SecretKey{Value: kgen.params.RingQP().AtLevel(levelQ, kgen.params.MaxLevelP()).NewPoly()}
	kgen.params.RingQ().AtLevel(levelQ).UnfoldConjugateInvariantToStandard(skConjugateInvariant.Value.Q, skCIMappedToStandard.Value.Q)

	if kgen.params.PCount() != 0 {
		ExtendBasisSmallNormAndCenterNTTMontgomery(kgen.params.RingQ(), kgen.params.RingP(), skCIMappedToStandard.Value.Q, kgen.buffQP[0].Q, skCIMappedToStandard.Value.P)
	}

	levelQ, levelP, BaseTwoDecomposition, compressed := ResolveEvaluationKeyParameters(kgen.params, evkParams)

	stdToci = newEvaluationKey(kgen.params, levelQ, levelP, BaseTwoDecomposition, compressed)
	kgen.GenEvaluationKey(skStd, skCIMappedToStandard, stdToci)

	ciToStd = newEvaluationKey(kgen.params, levelQ, levelP, BaseTwoDecomposition, compressed)
	kgen.GenEvaluationKey(skCIMappedToStandard, skStd, ciToStd)

	return
}

// GenEvaluationKeyNew generates a new [EvaluationKey], that will re-encrypt a [Ciphertext] encrypted under the input key into the output key.
// If the ringDegree(skOutput) > ringDegree(skInput),  generates [-a*SkOut + w*P*skIn_{Y^{N/n}} + e, a] in X^{N}.
// If the ringDegree(skOutput) < ringDegree(skInput),  generates [-a*skOut_{Y^{N/n}} + w*P*skIn + e_{N}, a_{N}] in X^{N}.
// Else generates [-a*skOut + w*P*skIn + e, a] in X^{N}.
// The output [EvaluationKey] is always given in max(N, n) and in the moduli of the output [EvaluationKey].
// When re-encrypting a [Ciphertext] from Y^{N/n} to X^{N}, the Ciphertext must first be mapped to X^{N}
// using [SwitchCiphertextRingDegreeNTT](ctSmallDim, nil, ctLargeDim).
// When re-encrypting a [Ciphertext] from X^{N} to Y^{N/n}, the output of the re-encryption is in still X^{N} and
// must be mapped Y^{N/n} using [SwitchCiphertextRingDegreeNTT](ctLargeDim, ringQLargeDim, ctSmallDim).
func (kgen KeyGenerator) GenEvaluationKeyNew(skInput, skOutput *SecretKey, evkParams ...EvaluationKeyParameters) (evk *EvaluationKey) {
	levelQ, levelP, BaseTwoDecomposition, compressed := ResolveEvaluationKeyParameters(kgen.params, evkParams)
	evk = newEvaluationKey(kgen.params, levelQ, levelP, BaseTwoDecomposition, compressed)
	kgen.GenEvaluationKey(skInput, skOutput, evk)
	return
}

// GenEvaluationKey generates an [EvaluationKey], that will re-encrypt a [Ciphertext] encrypted under the input key into the output key.
// If the ringDegree(skOutput) > ringDegree(skInput),  generates [-a*SkOut + w*P*skIn_{Y^{N/n}} + e, a] in X^{N}.
// If the ringDegree(skOutput) < ringDegree(skInput),  generates [-a*skOut_{Y^{N/n}} + w*P*skIn + e_{N}, a_{N}] in X^{N}.
// Else generates [-a*skOut + w*P*skIn + e, a] in X^{N}.
// The output [EvaluationKey] is always given in max(N, n) and in the moduli of the output [EvaluationKey].
// When re-encrypting a [Ciphertext] from Y^{N/n} to X^{N}, the [Ciphertext] must first be mapped to X^{N}
// using [SwitchCiphertextRingDegreeNTT](ctSmallDim, nil, ctLargeDim).
// When re-encrypting a [Ciphertext] from X^{N} to Y^{N/n}, the output of the re-encryption is in still X^{N} and
// must be mapped Y^{N/n} using [SwitchCiphertextRingDegreeNTT](ctLargeDim, ringQLargeDim, ctSmallDim).
func (kgen KeyGenerator) GenEvaluationKey(skInput, skOutput *SecretKey, evk *EvaluationKey) {

	ringQ := kgen.params.RingQ()
	ringP := kgen.params.RingP()

	// Maps the smaller key to the largest with Y = X^{N/n}.
	ring.MapSmallDimensionToLargerDimensionNTT(skOutput.Value.Q, kgen.bufSkOut.Q)

	// Extends the modulus P of skOutput to the one of skInput
	if levelP := evk.LevelP(); levelP != -1 {
		ExtendBasisSmallNormAndCenterNTTMontgomery(ringQ, ringP.AtLevel(levelP), kgen.bufSkOut.Q, kgen.buffQP[0].Q, kgen.bufSkOut.P)
	}

	// Maps the smaller key to the largest dimension with Y = X^{N/n}.
	ring.MapSmallDimensionToLargerDimensionNTT(skInput.Value.Q, kgen.bufSkIn)
	ExtendBasisSmallNormAndCenterNTTMontgomery(ringQ, ringQ.AtLevel(skOutput.Value.Q.Level()), kgen.bufSkIn, kgen.buffQP[0].Q, kgen.bufSkIn)

	kgen.genEvaluationKey(kgen.bufSkIn, kgen.bufSkOut, evk)
}

func (kgen KeyGenerator) genEvaluationKey(skIn ring.Poly, skOut ringqp.Poly, evk *EvaluationKey) {

	enc := kgen.WithKey(&SecretKey{Value: skOut})

	// For a compressed evaluation key, a seed is created and stored in the EvaluationKey struct
	// struct while an uncompressed key uses an ephemeral seed.
	if evk.IsCompressed() {
		evk.Seed = make([]byte, 32)
		if n, err := kgen.prng.Read(evk.Seed); n != 32 || err != nil {
			panic(fmt.Errorf("unable to sample evaluation key seed"))
		}

		sampler, err := sampling.NewKeyedPRNG(evk.Seed)
		if err != nil {
			panic(fmt.Errorf("sampling.NewKeyedPRNG: %w", err))
		}

		enc = enc.WithPRNG(sampler)
	}

	// Samples an encryption of zero for each element of the EvaluationKey.
	for i := 0; i < len(evk.Value); i++ {
		for j := 0; j < len(evk.Value[i]); j++ {
			// evk[i][j] = (-a*sk + e, a) if the degree of degree of the GadgetCiphertext is 1
			// evk[i][j] = (-a*sk + e) if the degree is 0
			if err := enc.EncryptZero(Element[ringqp.Poly]{
				MetaData: &MetaData{CiphertextMetaData: CiphertextMetaData{IsNTT: true, IsMontgomery: true}},
				Value:    []ringqp.Poly(evk.Value[i][j]),
			}); err != nil {
				panic(err)
			}
		}
	}

	// Adds the plaintext (input-key) to the EvaluationKey.
	if err := AddPolyTimesGadgetVectorToGadgetCiphertext(skIn, []GadgetCiphertext{evk.GadgetCiphertext}, *kgen.params.RingQP(), kgen.buffQP[0].Q); err != nil {
		// Sanity check, this error should not happen.
		panic(err)
	}
}
