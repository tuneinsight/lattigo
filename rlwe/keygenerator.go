package rlwe

import (
	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/rlwe/ringqp"
	"github.com/tuneinsight/lattigo/v4/utils"
)

// KeyGenerator is a structure that stores the elements required to create new keys,
// as well as a memory buffer for intermediate values.
type KeyGenerator struct {
	*EncryptorSecretKey
}

// NewKeyGenerator creates a new KeyGenerator, from which the secret and public keys, as well as EvaluationKeys.
func NewKeyGenerator(params ParametersInterface) *KeyGenerator {
	return &KeyGenerator{
		EncryptorSecretKey: NewEncryptorSecretKey(params, NewSecretKey(params)),
	}
}

// GenSecretKeyNew generates a new SecretKey.
// Distribution is set according to `rlwe.Parameters.HammingWeight()`.
func (kgen KeyGenerator) GenSecretKeyNew() (sk *SecretKey) {
	sk = NewSecretKey(kgen.params)
	kgen.GenSecretKey(sk)
	return
}

// GenSecretKey generates a SecretKey.
// Distribution is set according to `rlwe.Parameters.HammingWeight()`.
func (kgen KeyGenerator) GenSecretKey(sk *SecretKey) {
	kgen.genSecretKeyFromSampler(kgen.xsSampler, sk)
}

// GenSecretKeyWithHammingWeightNew generates a new SecretKey with exactly hw non-zero coefficients.
func (kgen *KeyGenerator) GenSecretKeyWithHammingWeightNew(hw int) (sk *SecretKey) {
	sk = NewSecretKey(kgen.params)
	kgen.GenSecretKeyWithHammingWeight(hw, sk)
	return
}

// GenSecretKeyWithHammingWeight generates a SecretKey with exactly hw non-zero coefficients.
func (kgen KeyGenerator) GenSecretKeyWithHammingWeight(hw int, sk *SecretKey) {
	kgen.genSecretKeyFromSampler(ring.NewSampler(kgen.prng, kgen.params.RingQ(), ring.Ternary{H: hw}, false), sk)
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

// GenPublicKeyNew generates a new public key from the provided SecretKey.
func (kgen KeyGenerator) GenPublicKeyNew(sk *SecretKey) (pk *PublicKey) {
	pk = NewPublicKey(kgen.params)
	kgen.GenPublicKey(sk, pk)
	return
}

// GenPublicKey generates a public key from the provided SecretKey.
func (kgen KeyGenerator) GenPublicKey(sk *SecretKey, pk *PublicKey) {
	kgen.WithKey(sk).EncryptZero(OperandQP{
		MetaData: MetaData{IsNTT: true, IsMontgomery: true},
		Value:    []ringqp.Poly(pk.Value)})
}

// GenKeyPairNew generates a new SecretKey and a corresponding public key.
// Distribution is of the SecretKey set according to `rlwe.Parameters.HammingWeight()`.
func (kgen KeyGenerator) GenKeyPairNew() (sk *SecretKey, pk *PublicKey) {
	sk = kgen.GenSecretKeyNew()
	return sk, kgen.GenPublicKeyNew(sk)
}

// GenRelinearizationKeyNew generates a new EvaluationKey that will be used to relinearize Ciphertexts during multiplication.
func (kgen KeyGenerator) GenRelinearizationKeyNew(sk *SecretKey) (rlk *RelinearizationKey) {
	rlk = NewRelinearizationKey(kgen.params, kgen.params.MaxLevelQ(), kgen.params.MaxLevelP(), 0)
	kgen.GenRelinearizationKey(sk, rlk)
	return
}

// GenRelinearizationKey generates an EvaluationKey that will be used to relinearize Ciphertexts during multiplication.
func (kgen KeyGenerator) GenRelinearizationKey(sk *SecretKey, rlk *RelinearizationKey) {
	kgen.buffQP.Q.CopyValues(sk.Value.Q)
	kgen.params.RingQ().AtLevel(rlk.LevelQ()).MulCoeffsMontgomery(kgen.buffQP.Q, sk.Value.Q, kgen.buffQP.Q)
	kgen.genEvaluationKey(kgen.buffQP.Q, sk.Value, &rlk.EvaluationKey)
}

// GenGaloisKeyNew generates a new GaloisKey, enabling the automorphism X^{i} -> X^{i * galEl}.
func (kgen KeyGenerator) GenGaloisKeyNew(galEl uint64, sk *SecretKey) (gk *GaloisKey) {
	gk = &GaloisKey{EvaluationKey: *NewEvaluationKey(kgen.params, sk.LevelQ(), sk.LevelP(), 0)}
	kgen.GenGaloisKey(galEl, sk, gk)
	return
}

// GenGaloisKey generates a GaloisKey, enabling the automorphism X^{i} -> X^{i * galEl}.
func (kgen KeyGenerator) GenGaloisKey(galEl uint64, sk *SecretKey, gk *GaloisKey) {

	skIn := sk.Value
	skOut := kgen.buffQP

	ringQP := kgen.params.RingQP().AtLevel(gk.LevelQ(), gk.LevelP())

	ringQ := ringQP.RingQ
	ringP := ringQP.RingP

	// We encrypt [-a * pi_{k^-1}(sk) + sk, a]
	// This enables to first apply the gadget product, re-encrypting
	// a ciphetext from sk to pi_{k^-1}(sk) and then we apply pi_{k}
	// on the ciphertext.
	galElInv := kgen.params.ModInvGaloisElement(galEl)

	index := ring.AutomorphismNTTIndex(ringQ.N(), ringQ.NthRoot(), galElInv)

	ringQ.AutomorphismNTTWithIndex(skIn.Q, index, skOut.Q)

	if ringP != nil {
		ringP.AutomorphismNTTWithIndex(skIn.P, index, skOut.P)
	}

	kgen.genEvaluationKey(skIn.Q, skOut, &gk.EvaluationKey)

	gk.GaloisElement = galEl
	gk.NthRoot = ringQ.NthRoot()
}

// GenGaloisKeys generates the GaloisKey objects for all galois elements in galEls, and stores
// the resulting key for galois element i in gks[i].
// The galEls and gks parameters must have the same length.
func (kgen KeyGenerator) GenGaloisKeys(galEls []uint64, sk *SecretKey, gks []*GaloisKey) {
	if len(galEls) != len(gks) {
		panic("galEls and gks must have the same length")
	}
	for i, galEl := range galEls {
		if gks[i] == nil {
			gks[i] = kgen.GenGaloisKeyNew(galEl, sk)
		} else {
			kgen.GenGaloisKey(galEl, sk, gks[i])
		}
	}
}

// GenGaloisKeysNew generates the GaloisKey objects for all galois elements in galEls, and
// returns the resulting keys in a newly allocated []*GaloisKey.
func (kgen KeyGenerator) GenGaloisKeysNew(galEls []uint64, sk *SecretKey) []*GaloisKey {
	gks := make([]*GaloisKey, len(galEls))
	for i, galEl := range galEls {
		gks[i] = kgen.GenGaloisKeyNew(galEl, sk)
	}
	return gks
}

// GenEvaluationKeysForRingSwapNew generates the necessary EvaluationKeys to switch from a standard ring to to a conjugate invariant ring and vice-versa.
func (kgen KeyGenerator) GenEvaluationKeysForRingSwapNew(skStd, skConjugateInvariant *SecretKey) (stdToci, ciToStd *EvaluationKey) {

	levelQ := utils.Min(skStd.Value.Q.Level(), skConjugateInvariant.Value.Q.Level())

	skCIMappedToStandard := &SecretKey{Value: kgen.params.RingQP().AtLevel(levelQ, kgen.params.MaxLevelP()).NewPoly()}
	kgen.params.RingQ().AtLevel(levelQ).UnfoldConjugateInvariantToStandard(skConjugateInvariant.Value.Q, skCIMappedToStandard.Value.Q)

	if kgen.params.PCount() != 0 {
		kgen.extendQ2P2(kgen.params.MaxLevelP(), skCIMappedToStandard.Value.Q, kgen.buffQ[1], skCIMappedToStandard.Value.P)
	}

	return kgen.GenEvaluationKeyNew(skStd, skCIMappedToStandard), kgen.GenEvaluationKeyNew(skCIMappedToStandard, skStd)
}

// GenEvaluationKeyNew generates a new EvaluationKey, that will re-encrypt a Ciphertext encrypted under the input key into the output key.
// If the ringDegree(skOutput) > ringDegree(skInput),  generates [-a*SkOut + w*P*skIn_{Y^{N/n}} + e, a] in X^{N}.
// If the ringDegree(skOutput) < ringDegree(skInput),  generates [-a*skOut_{Y^{N/n}} + w*P*skIn + e_{N}, a_{N}] in X^{N}.
// Else generates [-a*skOut + w*P*skIn + e, a] in X^{N}.
// The output EvaluationKey is always given in max(N, n) and in the moduli of the output EvaluationKey.
// When re-encrypting a Ciphertext from Y^{N/n} to X^{N}, the Ciphertext must first be mapped to X^{N}
// using SwitchCiphertextRingDegreeNTT(ctSmallDim, nil, ctLargeDim).
// When re-encrypting a Ciphertext from X^{N} to Y^{N/n}, the output of the re-encryption is in still X^{N} and
// must be mapped Y^{N/n} using SwitchCiphertextRingDegreeNTT(ctLargeDim, ringQLargeDim, ctSmallDim).
func (kgen KeyGenerator) GenEvaluationKeyNew(skInput, skOutput *SecretKey) (evk *EvaluationKey) {
	levelQ := utils.Min(skOutput.LevelQ(), kgen.params.MaxLevelQ())
	levelP := utils.Min(skOutput.LevelP(), kgen.params.MaxLevelP())
	evk = NewEvaluationKey(kgen.params, levelQ, levelP, 0)
	kgen.GenEvaluationKey(skInput, skOutput, evk)
	return
}

// GenEvaluationKey generates an EvaluationKey, that will re-encrypt a Ciphertext encrypted under the input key into the output key.
// If the ringDegree(skOutput) > ringDegree(skInput),  generates [-a*SkOut + w*P*skIn_{Y^{N/n}} + e, a] in X^{N}.
// If the ringDegree(skOutput) < ringDegree(skInput),  generates [-a*skOut_{Y^{N/n}} + w*P*skIn + e_{N}, a_{N}] in X^{N}.
// Else generates [-a*skOut + w*P*skIn + e, a] in X^{N}.
// The output EvaluationKey is always given in max(N, n) and in the moduli of the output EvaluationKey.
// When re-encrypting a Ciphertext from Y^{N/n} to X^{N}, the Ciphertext must first be mapped to X^{N}
// using SwitchCiphertextRingDegreeNTT(ctSmallDim, nil, ctLargeDim).
// When re-encrypting a Ciphertext from X^{N} to Y^{N/n}, the output of the re-encryption is in still X^{N} and
// must be mapped Y^{N/n} using SwitchCiphertextRingDegreeNTT(ctLargeDim, ringQLargeDim, ctSmallDim).
func (kgen KeyGenerator) GenEvaluationKey(skInput, skOutput *SecretKey, evk *EvaluationKey) {

	ringQ := kgen.params.RingQ()
	ringP := kgen.params.RingP()

	// Maps the smaller key to the largest with Y = X^{N/n}.
	ring.MapSmallDimensionToLargerDimensionNTT(skOutput.Value.Q, kgen.buffQP.Q)

	// Extends the modulus P of skOutput to the one of skInput
	if levelP := evk.LevelP(); levelP != -1 {
		kgen.extendQ2P(ringQ, ringP.AtLevel(levelP), kgen.buffQP.Q, kgen.buffQ[0], kgen.buffQP.P)
	}

	// Maps the smaller key to the largest dimension with Y = X^{N/n}.
	ring.MapSmallDimensionToLargerDimensionNTT(skInput.Value.Q, kgen.buffQ[0])
	kgen.extendQ2P(ringQ, ringQ.AtLevel(skOutput.Value.Q.Level()), kgen.buffQ[0], kgen.buffQ[1], kgen.buffQ[0])

	kgen.genEvaluationKey(kgen.buffQ[0], kgen.buffQP, evk)
}

func (kgen KeyGenerator) extendQ2P2(levelP int, polQ, buff, polP ring.Poly) {
	ringQ := kgen.params.RingQ().AtLevel(0)
	ringP := kgen.params.RingP().AtLevel(levelP)

	// Switches Q[0] out of the NTT and Montgomery domain.
	ringQ.INTT(polQ, buff)
	ringQ.IMForm(buff, buff)

	// Reconstruct P from Q
	Q := ringQ.SubRings[0].Modulus
	QHalf := Q >> 1

	P := ringP.ModuliChain()
	N := ringQ.N()

	var sign uint64
	for j := 0; j < N; j++ {

		coeff := buff.Coeffs[0][j]

		sign = 1
		if coeff > QHalf {
			coeff = Q - coeff
			sign = 0
		}

		for i := 0; i < levelP+1; i++ {
			polP.Coeffs[i][j] = (coeff * sign) | (P[i]-coeff)*(sign^1)
		}
	}

	ringP.NTT(polP, polP)
	ringP.MForm(polP, polP)
}

func (kgen KeyGenerator) extendQ2P(rQ, rP *ring.Ring, polQ, buff, polP ring.Poly) {
	rQ = rQ.AtLevel(0)

	levelP := rP.Level()

	// Switches Q[0] out of the NTT and Montgomery domain.
	rQ.INTT(polQ, buff)
	rQ.IMForm(buff, buff)

	// Reconstruct P from Q
	Q := rQ.SubRings[0].Modulus
	QHalf := Q >> 1

	P := rP.ModuliChain()
	N := rQ.N()

	var sign uint64
	for j := 0; j < N; j++ {

		coeff := buff.Coeffs[0][j]

		sign = 1
		if coeff > QHalf {
			coeff = Q - coeff
			sign = 0
		}

		for i := 0; i < levelP+1; i++ {
			polP.Coeffs[i][j] = (coeff * sign) | (P[i]-coeff)*(sign^1)
		}
	}

	rP.NTT(polP, polP)
	rP.MForm(polP, polP)
}

func (kgen KeyGenerator) genEvaluationKey(skIn ring.Poly, skOut ringqp.Poly, evk *EvaluationKey) {

	enc := kgen.WithKey(&SecretKey{Value: skOut})
	// Samples an encryption of zero for each element of the EvaluationKey.
	for i := 0; i < len(evk.Value); i++ {
		for j := 0; j < len(evk.Value[0]); j++ {
			enc.EncryptZero(OperandQP{MetaData: MetaData{IsNTT: true, IsMontgomery: true}, Value: []ringqp.Poly(evk.Value[i][j])})
		}
	}

	// Adds the plaintext (input-key) to the EvaluationKey.
	AddPolyTimesGadgetVectorToGadgetCiphertext(skIn, []GadgetCiphertext{evk.GadgetCiphertext}, *kgen.params.RingQP(), kgen.buffQ[0])
}
