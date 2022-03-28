package rlwe

import (
<<<<<<< dev_bfv_poly
	"math"

=======
>>>>>>> First step for adding bit-decomp
	"github.com/tuneinsight/lattigo/v3/ring"
	"github.com/tuneinsight/lattigo/v3/rlwe/gadget"
	"github.com/tuneinsight/lattigo/v3/rlwe/ringqp"
)

// KeyGenerator is an interface implementing the methods of the KeyGenerator.
type KeyGenerator interface {
	GenSecretKey() (sk *SecretKey)
	GenSecretKeyGaussian() (sk *SecretKey)
	GenSecretKeyWithDistrib(p float64) (sk *SecretKey)
	GenSecretKeyWithHammingWeight(hw int) (sk *SecretKey)
	GenPublicKey(sk *SecretKey) (pk *PublicKey)
	GenKeyPair() (sk *SecretKey, pk *PublicKey)
	GenRelinearizationKey(sk *SecretKey, maxDegree int) (evk *RelinearizationKey)
	GenSwitchingKey(skInput, skOutput *SecretKey) (newevakey *SwitchingKey)
	GenSwitchingKeyForGalois(galEl uint64, sk *SecretKey) (swk *SwitchingKey)
	GenRotationKeys(galEls []uint64, sk *SecretKey) (rks *RotationKeySet)
	GenSwitchingKeyForRotationBy(k int, sk *SecretKey) (swk *SwitchingKey)
	GenRotationKeysForRotations(ks []int, inclueSwapRows bool, sk *SecretKey) (rks *RotationKeySet)
	GenSwitchingKeyForRowRotation(sk *SecretKey) (swk *SwitchingKey)
	GenRotationKeysForInnerSum(sk *SecretKey) (rks *RotationKeySet)
	GenSwitchingKeysForRingSwap(skCKKS, skCI *SecretKey) (swkStdToConjugateInvariant, swkConjugateInvariantToStd *SwitchingKey)
}

// KeyGenerator is a structure that stores the elements required to create new keys,
// as well as a memory buffer for intermediate values.
type keyGenerator struct {
<<<<<<< dev_bfv_poly
	params           Parameters
<<<<<<< dev_bfv_poly
	buffQ            *ring.Poly
	buffQP           PolyQP
=======
	poolQ            *ring.Poly
	poolQP           ringqp.Poly
>>>>>>> [rlwe]: complete refactoring
	ternarySampler   *ring.TernarySampler
	gaussianSamplerQ *ring.GaussianSampler
	uniformSamplerQ  *ring.UniformSampler
	uniformSamplerP  *ring.UniformSampler
=======
	encryptor
>>>>>>> [rlwe]: simplified KeyGenerator, which is now based on Encryptor.
}

// NewKeyGenerator creates a new KeyGenerator, from which the secret and public keys, as well as the evaluation,
// rotation and switching keys can be generated.
func NewKeyGenerator(params Parameters) KeyGenerator {
<<<<<<< dev_bfv_poly

	prng, err := utils.NewPRNG()
	if err != nil {
		panic(err)
	}

<<<<<<< dev_bfv_poly
	var buffQP PolyQP
=======
	var poolQP ringqp.Poly
>>>>>>> [rlwe]: complete refactoring
	var uniformSamplerP *ring.UniformSampler
	if params.PCount() > 0 {
		buffQP = params.RingQP().NewPoly()
		uniformSamplerP = ring.NewUniformSampler(prng, params.RingP())
	} else {
		poolQP = ringqp.Poly{Q: params.RingQ().NewPoly()}
	}

	return &keyGenerator{
		params:           params,
		buffQ:            params.RingQ().NewPoly(),
		buffQP:           buffQP,
		ternarySampler:   ring.NewTernarySamplerWithHammingWeight(prng, params.ringQ, params.h, false),
		gaussianSamplerQ: ring.NewGaussianSampler(prng, params.RingQ(), params.Sigma(), int(6*params.Sigma())),
		uniformSamplerQ:  ring.NewUniformSampler(prng, params.RingQ()),
		uniformSamplerP:  uniformSamplerP,
=======
	return &keyGenerator{
		encryptor: newEncryptor(params),
>>>>>>> [rlwe]: simplified KeyGenerator, which is now based on Encryptor.
	}
}

// GenSecretKey generates a new SecretKey with the distribution [1/3, 1/3, 1/3].
func (keygen *keyGenerator) GenSecretKey() (sk *SecretKey) {
	return keygen.genSecretKeyFromSampler(keygen.ternarySampler)
}

// GenSecretKey generates a new SecretKey with the error distribution.
func (keygen *keyGenerator) GenSecretKeyGaussian() (sk *SecretKey) {
	return keygen.genSecretKeyFromSampler(keygen.gaussianSampler)
}

// GenSecretKeyWithDistrib generates a new SecretKey with the distribution [(p-1)/2, p, (p-1)/2].
func (keygen *keyGenerator) GenSecretKeyWithDistrib(p float64) (sk *SecretKey) {
	return keygen.genSecretKeyFromSampler(ring.NewTernarySampler(keygen.prng, keygen.params.RingQ(), p, false))
}

// GenSecretKeyWithHammingWeight generates a new SecretKey with exactly hw non-zero coefficients.
func (keygen *keyGenerator) GenSecretKeyWithHammingWeight(hw int) (sk *SecretKey) {
	return keygen.genSecretKeyFromSampler(ring.NewTernarySamplerWithHammingWeight(keygen.prng, keygen.params.RingQ(), hw, false))
}

// genSecretKeyFromSampler generates a new SecretKey sampled from the provided Sampler.
func (keygen *keyGenerator) genSecretKeyFromSampler(sampler ring.Sampler) (sk *SecretKey) {
	sk = new(SecretKey)
	ringQP := keygen.params.RingQP()
	sk.Value = ringQP.NewPoly()
	levelQ, levelP := sk.LevelQ(), sk.LevelP()
	sampler.Read(sk.Value.Q)

	if levelP > -1 {
		ringQP.ExtendBasisSmallNormAndCenter(sk.Value.Q, levelP, nil, sk.Value.P)
	}

	ringQP.NTTLvl(levelQ, levelP, sk.Value, sk.Value)
	ringQP.MFormLvl(levelQ, levelP, sk.Value, sk.Value)

	return
}

// GenPublicKey generates a new public key from the provided SecretKey.
func (keygen *keyGenerator) GenPublicKey(sk *SecretKey) (pk *PublicKey) {
	pk = NewPublicKey(keygen.params)
	keygen.encryptZeroSymetricQPNTT(pk.LevelQ(), pk.LevelP(), sk.Value, true, false, pk.Value)
	return
}

// GenKeyPair generates a new SecretKey with distribution [1/3, 1/3, 1/3] and a corresponding public key.
func (keygen *keyGenerator) GenKeyPair() (sk *SecretKey, pk *PublicKey) {
	sk = keygen.GenSecretKey()
	return sk, keygen.GenPublicKey(sk)
}

// GenRelinKey generates a new EvaluationKey that will be used to relinearize Ciphertexts during multiplication.
<<<<<<< 83ae36f5f9908381fe0d957ce0daa4f037d38e6f
func (keygen *keyGenerator) GenRelinearizationKey(sk *SecretKey, maxDegree, bitDecomp int) (evk *RelinearizationKey) {

	if keygen.params.PCount() == 0 {
		panic("cannot GenRelinearizationKey: modulus P is empty")
	}
=======
func (keygen *keyGenerator) GenRelinearizationKey(sk *SecretKey, maxDegree int) (evk *RelinearizationKey) {
>>>>>>> wip

	levelQ := keygen.params.QCount() - 1
	levelP := keygen.params.PCount() - 1

	evk = new(RelinearizationKey)
	evk.Keys = make([]*SwitchingKey, maxDegree)
	for i := range evk.Keys {
		evk.Keys[i] = NewSwitchingKey(keygen.params, levelQ, levelP)
	}

	keygen.buffQP.Q.CopyValues(sk.Value.Q)
	ringQ := keygen.params.RingQ()
	for i := 0; i < maxDegree; i++ {
		ringQ.MulCoeffsMontgomery(keygen.buffQP.Q, sk.Value.Q, keygen.buffQP.Q)
		keygen.genSwitchingKey(keygen.buffQP.Q, sk.Value, evk.Keys[i])
	}

	return
}

// GenRotationKeys generates a RotationKeySet from a list of galois element corresponding to the desired rotations
// See also GenRotationKeysForRotations.
func (keygen *keyGenerator) GenRotationKeys(galEls []uint64, sk *SecretKey) (rks *RotationKeySet) {
	rks = NewRotationKeySet(keygen.params, galEls)
	for _, galEl := range galEls {
		keygen.genrotKey(sk.Value, keygen.params.InverseGaloisElement(galEl), rks.Keys[galEl])
	}
	return rks
}

func (keygen *keyGenerator) GenSwitchingKeyForRotationBy(k int, sk *SecretKey) (swk *SwitchingKey) {
	swk = NewSwitchingKey(keygen.params, keygen.params.QCount()-1, keygen.params.PCount()-1)
	galElInv := keygen.params.GaloisElementForColumnRotationBy(-int(k))
	keygen.genrotKey(sk.Value, galElInv, swk)
	return
}

// GenRotationKeysForRotations generates a RotationKeySet supporting left rotations by k positions for all k in ks.
// Negative k is equivalent to a right rotation by k positions
// If includeConjugate is true, the resulting set contains the conjugation key.
func (keygen *keyGenerator) GenRotationKeysForRotations(ks []int, includeConjugate bool, sk *SecretKey) (rks *RotationKeySet) {
	galEls := make([]uint64, len(ks), len(ks)+1)
	for i, k := range ks {
		galEls[i] = keygen.params.GaloisElementForColumnRotationBy(k)
	}
	if includeConjugate {
		galEls = append(galEls, keygen.params.GaloisElementForRowRotation())
	}
	return keygen.GenRotationKeys(galEls, sk)
}

func (keygen *keyGenerator) GenSwitchingKeyForRowRotation(sk *SecretKey) (swk *SwitchingKey) {
	swk = NewSwitchingKey(keygen.params, keygen.params.QCount()-1, keygen.params.PCount()-1)
	keygen.genrotKey(sk.Value, keygen.params.GaloisElementForRowRotation(), swk)
	return
}

func (keygen *keyGenerator) GenSwitchingKeyForGalois(galoisEl uint64, sk *SecretKey) (swk *SwitchingKey) {
	swk = NewSwitchingKey(keygen.params, keygen.params.QCount()-1, keygen.params.PCount()-1)
	keygen.genrotKey(sk.Value, keygen.params.InverseGaloisElement(galoisEl), swk)
	return
}

// GenRotationKeysForInnerSum generates a RotationKeySet supporting the InnerSum operation of the Evaluator
func (keygen *keyGenerator) GenRotationKeysForInnerSum(sk *SecretKey) (rks *RotationKeySet) {
	return keygen.GenRotationKeys(keygen.params.GaloisElementsForRowInnerSum(), sk)
}

func (keygen *keyGenerator) genrotKey(sk ringqp.Poly, galEl uint64, swk *SwitchingKey) {

	skIn := sk
	skOut := keygen.buffQP
	ringQ := keygen.params.RingQ()

	index := ringQ.PermuteNTTIndex(galEl)
	ringQ.PermuteNTTWithIndexLvl(keygen.params.QCount()-1, skIn.Q, index, skOut.Q)
	ringQ.PermuteNTTWithIndexLvl(keygen.params.PCount()-1, skIn.P, index, skOut.P)

	keygen.genSwitchingKey(skIn.Q, skOut, swk)
}

// GenSwitchingKeysForRingSwap generates the necessary switching keys to switch from a standard ring to to a conjugate invariant ring and vice-versa.
func (keygen *keyGenerator) GenSwitchingKeysForRingSwap(skStd, skConjugateInvariant *SecretKey) (swkStdToConjugateInvariant, swkConjugateInvariantToStd *SwitchingKey) {

	skCIMappedToStandard := &SecretKey{Value: keygen.buffQP}
	keygen.params.RingQ().UnfoldConjugateInvariantToStandard(skConjugateInvariant.Value.Q.Level(), skConjugateInvariant.Value.Q, skCIMappedToStandard.Value.Q)
	keygen.params.RingQ().UnfoldConjugateInvariantToStandard(skConjugateInvariant.Value.P.Level(), skConjugateInvariant.Value.P, skCIMappedToStandard.Value.P)

	swkConjugateInvariantToStd = keygen.GenSwitchingKey(skCIMappedToStandard, skStd)
	swkStdToConjugateInvariant = keygen.GenSwitchingKey(skStd, skCIMappedToStandard)
	return
}

// GenSwitchingKey generates a new key-switching key, that will re-encrypt a Ciphertext encrypted under the input key into the output key.
// If the ringDegree(skOutput) > ringDegree(skInput),  generates [-a*SkOut + w*P*skIn_{Y^{N/n}} + e, a] in X^{N}.
// If the ringDegree(skOutput) < ringDegree(skInput),  generates [-a*skOut_{Y^{N/n}} + w*P*skIn + e_{N}, a_{N}] in X^{N}.
// Else generates [-a*skOut + w*P*skIn + e, a] in X^{N}.
// The output switching key is always given in max(N, n) and in the moduli of the output switching key.
// When key-switching a ciphertext from Y^{N/n} to X^{N}, the ciphertext must first be mapped to X^{N}
// using SwitchCiphertextRingDegreeNTT(ctSmallDim, nil, ctLargeDim).
// When key-switching a ciphertext from X^{N} to Y^{N/n}, the output of the key-switch is in still X^{N} and
// must be mapped Y^{N/n} using SwitchCiphertextRingDegreeNTT(ctLargeDim, ringQLargeDim, ctSmallDim).
func (keygen *keyGenerator) GenSwitchingKey(skInput, skOutput *SecretKey) (swk *SwitchingKey) {

	levelP := skOutput.LevelP()

	// Allocates the switching-key.
	swk = NewSwitchingKey(keygen.params, skOutput.Value.Q.Level(), levelP)

<<<<<<< dev_bfv_poly:rlwe/keygen.go
	if len(skInput.Value.Q.Coeffs[0]) > len(skOutput.Value.Q.Coeffs[0]) { // N -> n
<<<<<<< dev_bfv_poly
		ring.MapSmallDimensionToLargerDimensionNTT(skOutput.Value.Q, keygen.buffQP.Q)
		ring.MapSmallDimensionToLargerDimensionNTT(skOutput.Value.P, keygen.buffQP.P)
		keygen.genSwitchingKey(skInput.Value.Q, keygen.buffQP, swk)
=======
=======
	// N -> n (swk is to switch to a smaller dimension).
	if len(skInput.Value.Q.Coeffs[0]) > len(skOutput.Value.Q.Coeffs[0]) {
>>>>>>> file name update & doc:rlwe/keygenerator.go

		// Maps the smaller key to the largest with Y = X^{N/n}.
		ring.MapSmallDimensionToLargerDimensionNTT(skOutput.Value.Q, keygen.poolQP.Q)
		if levelP != -1 {
			ring.MapSmallDimensionToLargerDimensionNTT(skOutput.Value.P, keygen.poolQP.P)
		}

		keygen.genSwitchingKey(skInput.Value.Q, keygen.poolQP, swk)
<<<<<<< dev_bfv_poly:rlwe/keygen.go
>>>>>>> wip
	} else { // N -> N or n -> N
<<<<<<< dev_bfv_poly
		ring.MapSmallDimensionToLargerDimensionNTT(skInput.Value.Q, keygen.buffQ)
=======
=======

	} else { // N -> N or n -> N (swk switch to the same or a larger dimension)

		// Maps the smaller key to the largest dimension with Y = X^{N/n}.
>>>>>>> file name update & doc:rlwe/keygenerator.go
		ring.MapSmallDimensionToLargerDimensionNTT(skInput.Value.Q, keygen.poolQ[0])
>>>>>>> [rlwe]: simplified KeyGenerator, which is now based on Encryptor.

		// Extends the modulus of the input key to the one of the output key
		// if the former is smaller.
		if skInput.Value.Q.Level() < skOutput.Value.Q.Level() {

			ringQ := keygen.params.RingQ()

<<<<<<< dev_bfv_poly:rlwe/keygen.go
<<<<<<< dev_bfv_poly
			ringQ.InvNTTLvl(0, keygen.buffQ, keygen.buffQP.Q)
			ringQ.InvMFormLvl(0, keygen.buffQP.Q, keygen.buffQP.Q)
=======
=======
			// Switches out of the NTT and Montgomery domain.
>>>>>>> file name update & doc:rlwe/keygenerator.go
			ringQ.InvNTTLvl(0, keygen.poolQ[0], keygen.poolQP.Q)
			ringQ.InvMFormLvl(0, keygen.poolQP.Q, keygen.poolQP.Q)
>>>>>>> [rlwe]: simplified KeyGenerator, which is now based on Encryptor.

			// Extends the RNS basis of the small norm polynomial.
			Q := ringQ.Modulus[0]
			QHalf := Q >> 1

<<<<<<< dev_bfv_poly
			polQ := keygen.buffQP.Q
			polP := keygen.buffQ
=======
			polQ := keygen.poolQP.Q
			polP := keygen.poolQ[0]
>>>>>>> [rlwe]: simplified KeyGenerator, which is now based on Encryptor.
			var sign uint64
			for j := 0; j < ringQ.N; j++ {

				coeff := polQ.Coeffs[0][j]

				sign = 1
				if coeff > QHalf {
					coeff = Q - coeff
					sign = 0
				}

				for i := skInput.LevelQ() + 1; i < skOutput.LevelQ()+1; i++ {
					polP.Coeffs[i][j] = (coeff * sign) | (ringQ.Modulus[i]-coeff)*(sign^1)
				}
			}

			// Switches back to the NTT and Montgomery domain.
			for i := skInput.Value.Q.Level() + 1; i < skOutput.Value.Q.Level()+1; i++ {
				ringQ.NTTSingle(i, polP.Coeffs[i], polP.Coeffs[i])
				ring.MFormVec(polP.Coeffs[i], polP.Coeffs[i], ringQ.Modulus[i], ringQ.BredParams[i])
			}
		}

<<<<<<< dev_bfv_poly
		keygen.genSwitchingKey(keygen.buffQ, skOutput.Value, swk)
=======
		keygen.genSwitchingKey(keygen.poolQ[0], skOutput.Value, swk)
>>>>>>> [rlwe]: simplified KeyGenerator, which is now based on Encryptor.
	}

	return
}

func (keygen *keyGenerator) genSwitchingKey(skIn *ring.Poly, skOut ringqp.Poly, swk *SwitchingKey) {

<<<<<<< dev_bfv_poly
	ringQP := keygen.params.RingQP()

<<<<<<< dev_bfv_poly
	levelQ := swk.Value[0][0][0].Q.Level()
	hasModulusP := swk.Value[0][0][0].P != nil

<<<<<<< dev_bfv_poly
	// Computes P * skIn
<<<<<<< 83ae36f5f9908381fe0d957ce0daa4f037d38e6f
	ringQ.MulScalarBigintLvl(levelQ, skIn, ringQP.RingP.ModulusAtLevel[levelP], keygen.buffQ)
=======
<<<<<<< dev_bfv_poly
	ringQ.MulScalarBigintLvl(levelQ, skIn, ringQP.RingP.ModulusBigint[levelP], keygen.buffQ)
=======
	ringQ.MulScalarBigintLvl(levelQ, skIn, ringQP.RingP.ModulusBigint[levelP], keygen.poolQ)
=======
	var levelP int
	if hasModulusP {

		levelP = swk.Value[0][0][0].P.Level()

		var pBigInt *big.Int
		if levelP == keygen.params.PCount()-1 {
			pBigInt = keygen.params.RingP().ModulusBigint
		} else {
			P := keygen.params.RingP().Modulus
			pBigInt = new(big.Int).SetUint64(P[0])
			for i := 1; i < levelP+1; i++ {
				pBigInt.Mul(pBigInt, ring.NewUint(P[i]))
			}
		}

		// Computes P * skIn
		ringQ.MulScalarBigintLvl(levelQ, skIn, pBigInt, keygen.poolQ)
	} else {
		levelP = 0
		ring.CopyLvl(levelQ, skIn, keygen.poolQ)
	}
>>>>>>> First step for adding bit-decomp
<<<<<<< dev_bfv_poly
>>>>>>> First step for adding bit-decomp
<<<<<<< 83ae36f5f9908381fe0d957ce0daa4f037d38e6f
>>>>>>> First step for adding bit-decomp
=======
=======
=======
=======
>>>>>>> [rlwe]: simplified KeyGenerator, which is now based on Encryptor.
	levelQ := swk.LevelQ()
	levelP := swk.LevelP()
>>>>>>> wip
>>>>>>> [rlwe]: complete refactoring
>>>>>>> [rlwe]: complete refactoring

<<<<<<< dev_bfv_poly:rlwe/keygen.go
<<<<<<< dev_bfv_poly
	RNSDecomp := len(swk.Value)
	BITDecomp := len(swk.Value[0])

	for j := 0; j < BITDecomp; j++ {
		for i := 0; i < RNSDecomp; i++ {

			// (e, 0)
			keygen.gaussianSamplerQ.ReadLvl(levelQ, swk.Value[i][j][0].Q)

			if levelP != -1 {
				ringQP.ExtendBasisSmallNormAndCenter(swk.Value[i][j][0].Q, levelP, nil, swk.Value[i][j][0].P)
				keygen.uniformSamplerP.ReadLvl(levelP, swk.Value[i][j][1].P)
			}

<<<<<<< dev_bfv_poly
			qi := ringQ.Modulus[index]
			p0tmp := keygen.buffQ.Coeffs[index]
			p1tmp := swk.Value[i][0].Q.Coeffs[index]
=======
			ringQP.NTTLazyLvl(levelQ, levelP, swk.Value[i][j][0], swk.Value[i][j][0])
			ringQP.MFormLvl(levelQ, levelP, swk.Value[i][j][0], swk.Value[i][j][0])

			// (e, a)
			keygen.uniformSamplerQ.ReadLvl(levelQ, swk.Value[i][j][1].Q)
>>>>>>> First step for adding bit-decomp

			// (- a * skOut + e, a)
			ringQP.MulCoeffsMontgomeryAndSubLvl(levelQ, levelP, swk.Value[i][j][1], skOut, swk.Value[i][j][0])
=======
=======
	// Samples an encryption of zero for each element of the switching-key.
>>>>>>> file name update & doc:rlwe/keygenerator.go
	for i := 0; i < len(swk.Value); i++ {
		for j := 0; j < len(swk.Value[0]); j++ {
			keygen.encryptZeroSymetricQPNTT(levelQ, levelP, skOut, true, true, swk.Value[i][j])
>>>>>>> [rlwe]: simplified KeyGenerator, which is now based on Encryptor.
		}
	}

	// Adds the plaintext (input-key) to the switching-key.
	gadget.AddPolyToGadgetMatrix(skIn, []gadget.Ciphertext{swk.Ciphertext}, *keygen.params.RingQP(), keygen.params.LogBase2(), keygen.poolQ[0])
}
