package rlwe

import (
	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/ring/distribution"
	"github.com/tuneinsight/lattigo/v4/utils"
)

// KeyGenerator is a structure that stores the elements required to create new keys,
// as well as a memory buffer for intermediate values.
type KeyGenerator struct {
	*skEncryptor
}

// NewKeyGenerator creates a new KeyGenerator, from which the secret and public keys, as well as EvaluationKeys.
func NewKeyGenerator(params ParametersInterface) *KeyGenerator {
	return &KeyGenerator{
		skEncryptor: newSkEncryptor(params, NewSecretKey(params)),
	}
}

// GenSecretKeyNew generates a new SecretKey.
// Distribution is set according to `rlwe.Parameters.HammingWeight()`.
func (kgen *KeyGenerator) GenSecretKeyNew() (sk *SecretKey) {
	sk = NewSecretKey(kgen.params)
	kgen.GenSecretKey(sk)
	return
}

// GenSecretKey generates a SecretKey.
// Distribution is set according to `rlwe.Parameters.HammingWeight()`.
func (kgen *KeyGenerator) GenSecretKey(sk *SecretKey) {
	kgen.genSecretKeyFromSampler(kgen.ternarySampler, sk)
}

// GenSecretKeyWithHammingWeightNew generates a new SecretKey with exactly hw non-zero coefficients.
func (kgen *KeyGenerator) GenSecretKeyWithHammingWeightNew(hw int) (sk *SecretKey) {
	sk = NewSecretKey(kgen.params)
	kgen.GenSecretKeyWithHammingWeight(hw, sk)
	return
}

// GenSecretKeyWithHammingWeight generates a SecretKey with exactly hw non-zero coefficients.
func (kgen *KeyGenerator) GenSecretKeyWithHammingWeight(hw int, sk *SecretKey) {
	kgen.genSecretKeyFromSampler(ring.NewSampler(kgen.prng, kgen.params.RingQ(), &distribution.Ternary{H: hw}, false), sk)
}

func (kgen *KeyGenerator) genSecretKeyFromSampler(sampler ring.Sampler, sk *SecretKey) {

	ringQP := kgen.params.RingQP().AtLevel(sk.LevelQ(), sk.LevelP())

	sampler.AtLevel(sk.LevelQ()).Read(sk.Value.Q)

	if levelP := sk.LevelP(); levelP > -1 {
		ringQP.ExtendBasisSmallNormAndCenter(sk.Value.Q, levelP, nil, sk.Value.P)
	}

	ringQP.NTT(&sk.Value, &sk.Value)
	ringQP.MForm(&sk.Value, &sk.Value)
}

// GenPublicKeyNew generates a new public key from the provided SecretKey.
func (kgen *KeyGenerator) GenPublicKeyNew(sk *SecretKey) (pk *PublicKey) {
	pk = NewPublicKey(kgen.params)
	kgen.GenPublicKey(sk, pk)
	return
}

// GenPublicKey generates a public key from the provided SecretKey.
func (kgen *KeyGenerator) GenPublicKey(sk *SecretKey, pk *PublicKey) {
	kgen.WithKey(sk).EncryptZero(&pk.OperandQP)
}

// GenKeyPairNew generates a new SecretKey and a corresponding public key.
// Distribution is of the SecretKey set according to `rlwe.Parameters.HammingWeight()`.
func (kgen *KeyGenerator) GenKeyPairNew() (sk *SecretKey, pk *PublicKey) {
	sk = kgen.GenSecretKeyNew()
	return sk, kgen.GenPublicKeyNew(sk)
}

// GenRelinearizationKeyNew generates a new EvaluationKey that will be used to relinearize Ciphertexts during multiplication.
func (kgen *KeyGenerator) GenRelinearizationKeyNew(sk *SecretKey) (rlk *RelinearizationKey) {
	rlk = NewRelinearizationKey(kgen.params)
	kgen.GenRelinearizationKey(sk, rlk)
	return
}

// GenRelinearizationKey generates an EvaluationKey that will be used to relinearize Ciphertexts during multiplication.
func (kgen *KeyGenerator) GenRelinearizationKey(sk *SecretKey, rlk *RelinearizationKey) {
	kgen.buffQP.Q.CopyValues(sk.Value.Q)
	kgen.params.RingQ().AtLevel(rlk.LevelQ()).MulCoeffsMontgomery(kgen.buffQP.Q, sk.Value.Q, kgen.buffQP.Q)
	kgen.genEvaluationKey(kgen.buffQP.Q, sk, &rlk.EvaluationKey)
}

// GenGaloisKeyNew generates a new GaloisKey, enabling the automorphism X^{i} -> X^{i * galEl}.
func (kgen *KeyGenerator) GenGaloisKeyNew(galEl uint64, sk *SecretKey) (gk *GaloisKey) {
	gk = &GaloisKey{EvaluationKey: *NewEvaluationKey(kgen.params, sk.LevelQ(), sk.LevelP())}
	kgen.GenGaloisKey(galEl, sk, gk)
	return
}

// GenGaloisKey generates a GaloisKey, enabling the automorphism X^{i} -> X^{i * galEl}.
func (kgen *KeyGenerator) GenGaloisKey(galEl uint64, sk *SecretKey, gk *GaloisKey) {

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

	kgen.genEvaluationKey(skIn.Q, &SecretKey{Value: skOut}, &gk.EvaluationKey)

	gk.GaloisElement = galEl
	gk.NthRoot = ringQ.NthRoot()
}

// GenEvaluationKeysForRingSwapNew generates the necessary EvaluationKeys to switch from a standard ring to to a conjugate invariant ring and vice-versa.
func (kgen *KeyGenerator) GenEvaluationKeysForRingSwapNew(skStd, skConjugateInvariant *SecretKey) (stdToci, ciToStd *EvaluationKey) {

	levelQ := utils.Min(skStd.Value.Q.Level(), skConjugateInvariant.Value.Q.Level())

	skCIMappedToStandard := &SecretKey{Value: kgen.buffQP}
	kgen.params.RingQ().AtLevel(levelQ).UnfoldConjugateInvariantToStandard(skConjugateInvariant.Value.Q, skCIMappedToStandard.Value.Q)

	if kgen.params.PCount() != 0 {
		kgen.extendQ2P(kgen.params.MaxLevelP(), skCIMappedToStandard.Value.Q, kgen.buffQ[0], skCIMappedToStandard.Value.P)
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
func (kgen *KeyGenerator) GenEvaluationKeyNew(skInput, skOutput *SecretKey) (evk *EvaluationKey) {
	levelQ := utils.Min(skOutput.LevelQ(), kgen.params.MaxLevelQ())
	levelP := utils.Min(skOutput.LevelP(), kgen.params.MaxLevelP())
	evk = NewEvaluationKey(kgen.params, levelQ, levelP)
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
func (kgen *KeyGenerator) GenEvaluationKey(skInput, skOutput *SecretKey, evk *EvaluationKey) {

	// N -> n (evk is to switch to a smaller dimension).
	if len(skInput.Value.Q.Coeffs[0]) > len(skOutput.Value.Q.Coeffs[0]) {

		// Maps the smaller key to the largest with Y = X^{N/n}.
		ring.MapSmallDimensionToLargerDimensionNTT(skOutput.Value.Q, kgen.buffQP.Q)

		// Extends the modulus P of skOutput to the one of skInput
		if levelP := evk.LevelP(); levelP != -1 {
			kgen.extendQ2P(levelP, kgen.buffQP.Q, kgen.buffQ[0], kgen.buffQP.P)
		}

		kgen.genEvaluationKey(skInput.Value.Q, &SecretKey{Value: kgen.buffQP}, evk)

	} else { // N -> N or n -> N (evk switch to the same or a larger dimension)

		// Maps the smaller key to the largest dimension with Y = X^{N/n}.
		ring.MapSmallDimensionToLargerDimensionNTT(skInput.Value.Q, kgen.buffQ[0])

		// Extends the modulus of the input key to the one of the output key
		// if the former is smaller.
		if skInput.Value.Q.Level() < skOutput.Value.Q.Level() {

			ringQ := kgen.params.RingQ().AtLevel(0)

			// Switches out of the NTT and Montgomery domain.
			ringQ.INTT(kgen.buffQ[0], kgen.buffQP.Q)
			ringQ.IMForm(kgen.buffQP.Q, kgen.buffQP.Q)

			// Extends the RNS basis of the small norm polynomial.
			Qi := ringQ.ModuliChain()
			Q := Qi[0]
			QHalf := Q >> 1

			polQ := kgen.buffQP.Q
			polP := kgen.buffQ[0]
			var sign uint64
			N := ringQ.N()
			for j := 0; j < N; j++ {

				coeff := polQ.Coeffs[0][j]

				sign = 1
				if coeff > QHalf {
					coeff = Q - coeff
					sign = 0
				}

				for i := skInput.LevelQ() + 1; i < skOutput.LevelQ()+1; i++ {
					polP.Coeffs[i][j] = (coeff * sign) | (Qi[i]-coeff)*(sign^1)
				}
			}

			// Switches back to the NTT and Montgomery domain.
			for i := skInput.Value.Q.Level() + 1; i < skOutput.Value.Q.Level()+1; i++ {
				ringQ.SubRings[i].NTT(polP.Coeffs[i], polP.Coeffs[i])
				ringQ.SubRings[i].MForm(polP.Coeffs[i], polP.Coeffs[i])
			}
		}

		kgen.genEvaluationKey(kgen.buffQ[0], skOutput, evk)
	}
}

func (kgen *KeyGenerator) extendQ2P(levelP int, polQ, buff, polP *ring.Poly) {
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

func (kgen *KeyGenerator) genEvaluationKey(skIn *ring.Poly, skOut *SecretKey, evk *EvaluationKey) {

	enc := kgen.WithKey(skOut)
	// Samples an encryption of zero for each element of the EvaluationKey.
	for i := 0; i < len(evk.Value); i++ {
		for j := 0; j < len(evk.Value[0]); j++ {
			enc.EncryptZero(evk.Value[i][j])
		}
	}

	// Adds the plaintext (input-key) to the EvaluationKey.
	AddPolyTimesGadgetVectorToGadgetCiphertext(skIn, []GadgetCiphertext{evk.GadgetCiphertext}, *kgen.params.RingQP(), kgen.params.Pow2Base(), kgen.buffQ[0])
}
