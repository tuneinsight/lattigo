package rlwe

import (
<<<<<<< dev_bfv_poly
	"math"

=======
>>>>>>> First step for adding bit-decomp
	"github.com/tuneinsight/lattigo/v3/ring"
	"github.com/tuneinsight/lattigo/v3/utils"
	"math/big"
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
	params           Parameters
	buffQ            *ring.Poly
	buffQP           PolyQP
	ternarySampler   *ring.TernarySampler
	gaussianSamplerQ *ring.GaussianSampler
	uniformSamplerQ  *ring.UniformSampler
	uniformSamplerP  *ring.UniformSampler
}

// NewKeyGenerator creates a new KeyGenerator, from which the secret and public keys, as well as the evaluation,
// rotation and switching keys can be generated.
func NewKeyGenerator(params Parameters) KeyGenerator {

	prng, err := utils.NewPRNG()
	if err != nil {
		panic(err)
	}

	var buffQP PolyQP
	var uniformSamplerP *ring.UniformSampler
	if params.PCount() > 0 {
		buffQP = params.RingQP().NewPoly()
		uniformSamplerP = ring.NewUniformSampler(prng, params.RingP())
	} else {
		poolQP = PolyQP{Q: params.RingQ().NewPoly()}
	}

	return &keyGenerator{
		params:           params,
		buffQ:            params.RingQ().NewPoly(),
		buffQP:           buffQP,
		ternarySampler:   ring.NewTernarySamplerWithHammingWeight(prng, params.ringQ, params.h, false),
		gaussianSamplerQ: ring.NewGaussianSampler(prng, params.RingQ(), params.Sigma(), int(6*params.Sigma())),
		uniformSamplerQ:  ring.NewUniformSampler(prng, params.RingQ()),
		uniformSamplerP:  uniformSamplerP,
	}
}

// GenSecretKey generates a new SecretKey with the distribution [1/3, 1/3, 1/3].
func (keygen *keyGenerator) GenSecretKey() (sk *SecretKey) {
	return keygen.genSecretKeyFromSampler(keygen.ternarySampler)
}

// GenSecretKey generates a new SecretKey with the error distribution.
func (keygen *keyGenerator) GenSecretKeyGaussian() (sk *SecretKey) {
	return keygen.genSecretKeyFromSampler(keygen.gaussianSamplerQ)
}

// GenSecretKeyWithDistrib generates a new SecretKey with the distribution [(p-1)/2, p, (p-1)/2].
func (keygen *keyGenerator) GenSecretKeyWithDistrib(p float64) (sk *SecretKey) {
	prng, err := utils.NewPRNG()
	if err != nil {
		panic(err)
	}
	ternarySamplerMontgomery := ring.NewTernarySampler(prng, keygen.params.RingQ(), p, false)
	return keygen.genSecretKeyFromSampler(ternarySamplerMontgomery)
}

// GenSecretKeyWithHammingWeight generates a new SecretKey with exactly hw non-zero coefficients.
func (keygen *keyGenerator) GenSecretKeyWithHammingWeight(hw int) (sk *SecretKey) {
	prng, err := utils.NewPRNG()
	if err != nil {
		panic(err)
	}
	ternarySamplerMontgomery := ring.NewTernarySamplerWithHammingWeight(prng, keygen.params.RingQ(), hw, false)
	return keygen.genSecretKeyFromSampler(ternarySamplerMontgomery)
}

// genSecretKeyFromSampler generates a new SecretKey sampled from the provided Sampler.
func (keygen *keyGenerator) genSecretKeyFromSampler(sampler ring.Sampler) (sk *SecretKey) {
	sk = new(SecretKey)
	if keygen.params.PCount() > 0 {
		ringQP := keygen.params.RingQP()
		sk.Value = ringQP.NewPoly()
		levelQ, levelP := keygen.params.QCount()-1, keygen.params.PCount()-1
		sampler.Read(sk.Value.Q)
		ringQP.ExtendBasisSmallNormAndCenter(sk.Value.Q, levelP, nil, sk.Value.P)
		ringQP.NTTLvl(levelQ, levelP, sk.Value, sk.Value)
		ringQP.MFormLvl(levelQ, levelP, sk.Value, sk.Value)
	} else {
		ringQ := keygen.params.RingQ()
		sk = new(SecretKey)
		sk.Value.Q = ringQ.NewPoly()
		sampler.Read(sk.Value.Q)
		ringQ.NTT(sk.Value.Q, sk.Value.Q)
		ringQ.MForm(sk.Value.Q, sk.Value.Q)
	}

	return
}

// GenPublicKey generates a new public key from the provided SecretKey.
func (keygen *keyGenerator) GenPublicKey(sk *SecretKey) (pk *PublicKey) {

	pk = new(PublicKey)

	if keygen.params.PCount() > 0 {

		ringQP := keygen.params.RingQP()
		levelQ, levelP := keygen.params.QCount()-1, keygen.params.PCount()-1

		//pk[0] = [-as + e]
		//pk[1] = [a]
		pk = NewPublicKey(keygen.params)
		keygen.gaussianSamplerQ.Read(pk.Value[0].Q)
		ringQP.ExtendBasisSmallNormAndCenter(pk.Value[0].Q, levelP, nil, pk.Value[0].P)
		ringQP.NTTLvl(levelQ, levelP, pk.Value[0], pk.Value[0])

		keygen.uniformSamplerQ.Read(pk.Value[1].Q)
		keygen.uniformSamplerP.Read(pk.Value[1].P)

		ringQP.MulCoeffsMontgomeryAndSubLvl(levelQ, levelP, sk.Value, pk.Value[1], pk.Value[0])
	} else {
		ringQ := keygen.params.RingQ()

		//pk[0] = [-as + e]
		//pk[1] = [a]
		pk = NewPublicKey(keygen.params)
		keygen.gaussianSamplerQ.Read(pk.Value[0].Q)

		ringQ.NTT(pk.Value[0].Q, pk.Value[0].Q)

		keygen.uniformSamplerQ.Read(pk.Value[1].Q)

		ringQ.MulCoeffsMontgomeryAndSub(sk.Value.Q, pk.Value[1].Q, pk.Value[0].Q)
	}
	return pk
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

func (keygen *keyGenerator) genrotKey(sk PolyQP, galEl uint64, swk *SwitchingKey) {

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

	var levelP int
	if skOutput.Value.P != nil {
		levelP = skOutput.Value.P.Level()
	} else {
		levelP = -1
	}

	swk = NewSwitchingKey(keygen.params, skOutput.Value.Q.Level(), levelP)

	if len(skInput.Value.Q.Coeffs[0]) > len(skOutput.Value.Q.Coeffs[0]) { // N -> n
<<<<<<< dev_bfv_poly
		ring.MapSmallDimensionToLargerDimensionNTT(skOutput.Value.Q, keygen.buffQP.Q)
		ring.MapSmallDimensionToLargerDimensionNTT(skOutput.Value.P, keygen.buffQP.P)
		keygen.genSwitchingKey(skInput.Value.Q, keygen.buffQP, swk)
=======

		ring.MapSmallDimensionToLargerDimensionNTT(skOutput.Value.Q, keygen.poolQP.Q)

		if levelP != -1 {
			ring.MapSmallDimensionToLargerDimensionNTT(skOutput.Value.P, keygen.poolQP.P)
		}

		keygen.genSwitchingKey(skInput.Value.Q, keygen.poolQP, swk)
>>>>>>> wip
	} else { // N -> N or n -> N
		ring.MapSmallDimensionToLargerDimensionNTT(skInput.Value.Q, keygen.buffQ)

		if skInput.Value.Q.Level() < skOutput.Value.Q.Level() {

			ringQ := keygen.params.RingQ()

			ringQ.InvNTTLvl(0, keygen.buffQ, keygen.buffQP.Q)
			ringQ.InvMFormLvl(0, keygen.buffQP.Q, keygen.buffQP.Q)

			Q := ringQ.Modulus[0]
			QHalf := Q >> 1

			polQ := keygen.buffQP.Q
			polP := keygen.buffQ
			var sign uint64
			for j := 0; j < ringQ.N; j++ {

				coeff := polQ.Coeffs[0][j]

				sign = 1
				if coeff > QHalf {
					coeff = Q - coeff
					sign = 0
				}

				for i := skInput.Value.Q.Level() + 1; i < skOutput.Value.Q.Level()+1; i++ {
					polP.Coeffs[i][j] = (coeff * sign) | (ringQ.Modulus[i]-coeff)*(sign^1)
				}
			}

			for i := skInput.Value.Q.Level() + 1; i < skOutput.Value.Q.Level()+1; i++ {
				ringQ.NTTSingle(i, polP.Coeffs[i], polP.Coeffs[i])
				ring.MFormVec(polP.Coeffs[i], polP.Coeffs[i], ringQ.Modulus[i], ringQ.BredParams[i])
			}
		}

		keygen.genSwitchingKey(keygen.buffQ, skOutput.Value, swk)
	}

	return
}

func (keygen *keyGenerator) genSwitchingKey(skIn *ring.Poly, skOut PolyQP, swk *SwitchingKey) {

	ringQ := keygen.params.RingQ()
	ringQP := keygen.params.RingQP()

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
>>>>>>> First step for adding bit-decomp
>>>>>>> First step for adding bit-decomp

	RNSDecomp := len(swk.Value)
	BITDecomp := len(swk.Value[0])

	var index int
	for j := 0; j < BITDecomp; j++ {
		for i := 0; i < RNSDecomp; i++ {

			// e
			keygen.gaussianSamplerQ.ReadLvl(levelQ, swk.Value[i][j][0].Q)

			if hasModulusP {
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

			// a (since a is uniform, we consider we already sample it in the NTT and Montgomery domain)
			keygen.uniformSamplerQ.ReadLvl(levelQ, swk.Value[i][j][1].Q)
>>>>>>> First step for adding bit-decomp

			// e + (skIn * P) * w^2 * (q_star * q_tild) mod QP
			//
			// q_prod = prod(q[i*alpha+j])
			// q_star = Q/qprod
			// q_tild = q_star^-1 mod q_prod
			//
			// Therefore : (skIn * P) * (q_star * q_tild) = sk*P mod q[i*alpha+j], else 0
			for k := 0; k < levelP+1; k++ {

				index = i*(levelP+1) + k

				// It handles the case where nb pj does not divide nb qi
				if index >= levelQ+1 {
					break
				}

				qi := ringQ.Modulus[index]
				p0tmp := keygen.poolQ.Coeffs[index]
				p1tmp := swk.Value[i][j][0].Q.Coeffs[index]

				for w := 0; w < ringQ.N; w++ {
					p1tmp[w] = ring.CRed(p1tmp[w]+p0tmp[w], qi)
				}
			}

			// (skIn * P) * (q_star * q_tild) - a * skOut + e mod QP
			ringQP.MulCoeffsMontgomeryAndSubLvl(levelQ, levelP, swk.Value[i][j][1], skOut, swk.Value[i][j][0])
		}

		ringQ.MulScalar(keygen.poolQ, 1<<keygen.params.logbase2, keygen.poolQ)
	}
}
