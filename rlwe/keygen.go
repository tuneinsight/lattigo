package rlwe

import (
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/utils"
	"math"
	"math/big"
)

// KeyGenerator is an interface implementing the methods of the KeyGenerator.
type KeyGenerator interface {
	GenSecretKey() (sk *SecretKey)
	GenSecretKeyGaussian() (sk *SecretKey)
	GenSecretKeyWithDistrib(p float64) (sk *SecretKey)
	GenSecretKeySparse(hw int) (sk *SecretKey)
	GenPublicKey(sk *SecretKey) (pk *PublicKey)
	GenKeyPair() (sk *SecretKey, pk *PublicKey)
	GenKeyPairSparse(hw int) (sk *SecretKey, pk *PublicKey)
	GenRelinearizationKey(sk *SecretKey, maxDegree int) (evk *RelinearizationKey)
	GenSwitchingKey(levelQ, levelP int, skInput, skOutput *SecretKey) (newevakey *SwitchingKey)
	GenSwitchingKeyForGalois(galEl uint64, sk *SecretKey) (swk *SwitchingKey)
	GenRotationKeys(galEls []uint64, sk *SecretKey) (rks *RotationKeySet)
	GenSwitchingKeyForRotationBy(k int, sk *SecretKey) (swk *SwitchingKey)
	GenRotationKeysForRotations(ks []int, inclueSwapRows bool, sk *SecretKey) (rks *RotationKeySet)
	GenSwitchingKeyForRowRotation(sk *SecretKey) (swk *SwitchingKey)
	GenRotationKeysForInnerSum(sk *SecretKey) (rks *RotationKeySet)
}

// KeyGenerator is a structure that stores the elements required to create new keys,
// as well as a small memory pool for intermediate values.
type keyGenerator struct {
	params           Parameters
	ringQ            *ring.Ring
	ringP            *ring.Ring
	poolQ            [2]*ring.Poly
	poolP            [2]*ring.Poly
	gaussianSamplerQ *ring.GaussianSampler
	uniformSamplerQ  *ring.UniformSampler
	uniformSamplerP  *ring.UniformSampler
}

// NewKeyGenerator creates a new KeyGenerator, from which the secret and public keys, as well as the evaluation,
// rotation and switching keys can be generated.
func NewKeyGenerator(params Parameters) KeyGenerator {

	ringQ := params.RingQ()
	ringP := params.RingP()

	prng, err := utils.NewPRNG()
	if err != nil {
		panic(err)
	}

	return &keyGenerator{
		params:           params,
		ringQ:            ringQ,
		ringP:            ringP,
		poolQ:            [2]*ring.Poly{ringQ.NewPoly(), ringQ.NewPoly()},
		poolP:            [2]*ring.Poly{ringP.NewPoly(), ringP.NewPoly()},
		gaussianSamplerQ: ring.NewGaussianSampler(prng, ringQ, params.Sigma(), int(6*params.Sigma())),
		uniformSamplerQ:  ring.NewUniformSampler(prng, ringQ),
		uniformSamplerP:  ring.NewUniformSampler(prng, ringP),
	}
}

// GenSecretKey generates a new SecretKey with the distribution [1/3, 1/3, 1/3].
func (keygen *keyGenerator) GenSecretKey() (sk *SecretKey) {
	return keygen.GenSecretKeyWithDistrib(1.0 / 3)
}

func (keygen *keyGenerator) GenSecretKeyGaussian() (sk *SecretKey) {
	sk = new(SecretKey)

	ringQ := keygen.ringQ
	ringP := keygen.ringP

	sk.Value[0] = keygen.gaussianSamplerQ.ReadNew()
	sk.Value[1] = ringP.NewPoly()
	extendBasisSmallNormAndCenter(ringQ, ringP, len(ringP.Modulus)-1, sk.Value[0], sk.Value[1])
	ringQ.NTT(sk.Value[0], sk.Value[0])
	ringP.NTT(sk.Value[1], sk.Value[1])
	ringQ.MForm(sk.Value[0], sk.Value[0])
	ringP.MForm(sk.Value[1], sk.Value[1])
	return sk
}

// GenSecretKeyWithDistrib generates a new SecretKey with the distribution [(p-1)/2, p, (p-1)/2].
func (keygen *keyGenerator) GenSecretKeyWithDistrib(p float64) (sk *SecretKey) {
	prng, err := utils.NewPRNG()
	if err != nil {
		panic(err)
	}

	ringQ := keygen.ringQ
	ringP := keygen.ringP

	ternarySamplerMontgomery := ring.NewTernarySampler(prng, ringQ, p, false)

	sk = new(SecretKey)
	sk.Value[0] = ternarySamplerMontgomery.ReadNew()
	sk.Value[1] = ringP.NewPoly()
	extendBasisSmallNormAndCenter(ringQ, ringP, len(ringP.Modulus)-1, sk.Value[0], sk.Value[1])
	ringQ.NTT(sk.Value[0], sk.Value[0])
	ringP.NTT(sk.Value[1], sk.Value[1])
	ringQ.MForm(sk.Value[0], sk.Value[0])
	ringP.MForm(sk.Value[1], sk.Value[1])
	return sk
}

// GenSecretKeySparse generates a new SecretKey with exactly hw non-zero coefficients.
func (keygen *keyGenerator) GenSecretKeySparse(hw int) (sk *SecretKey) {
	prng, err := utils.NewPRNG()
	if err != nil {
		panic(err)
	}

	ringQ := keygen.ringQ
	ringP := keygen.ringP

	ternarySamplerMontgomery := ring.NewTernarySamplerSparse(prng, ringQ, hw, false)

	sk = new(SecretKey)
	sk.Value[0] = ternarySamplerMontgomery.ReadNew()
	sk.Value[1] = ringP.NewPoly()
	extendBasisSmallNormAndCenter(ringQ, ringP, len(ringP.Modulus)-1, sk.Value[0], sk.Value[1])
	ringQ.NTT(sk.Value[0], sk.Value[0])
	ringP.NTT(sk.Value[1], sk.Value[1])
	ringQ.MForm(sk.Value[0], sk.Value[0])
	ringP.MForm(sk.Value[1], sk.Value[1])
	return sk
}

// GenPublicKey generates a new public key from the provided SecretKey.
func (keygen *keyGenerator) GenPublicKey(sk *SecretKey) (pk *PublicKey) {

	pk = new(PublicKey)

	ringQ := keygen.ringQ
	ringP := keygen.ringP

	//pk[0] = [-as + e]
	//pk[1] = [a]

	pk.Value[0][0] = keygen.gaussianSamplerQ.ReadNew()
	pk.Value[0][1] = ringP.NewPoly()
	extendBasisSmallNormAndCenter(ringQ, ringP, len(ringP.Modulus)-1, pk.Value[0][0], pk.Value[0][1])

	ringQ.NTT(pk.Value[0][0], pk.Value[0][0])
	ringP.NTT(pk.Value[0][1], pk.Value[0][1])

	pk.Value[1][0] = keygen.uniformSamplerQ.ReadNew()
	pk.Value[1][1] = keygen.uniformSamplerP.ReadNew()

	ringQ.MulCoeffsMontgomeryAndSub(sk.Value[0], pk.Value[1][0], pk.Value[0][0])
	ringP.MulCoeffsMontgomeryAndSub(sk.Value[1], pk.Value[1][1], pk.Value[0][1])

	return pk
}

// GenKeyPair generates a new SecretKey with distribution [1/3, 1/3, 1/3] and a corresponding public key.
func (keygen *keyGenerator) GenKeyPair() (sk *SecretKey, pk *PublicKey) {
	sk = keygen.GenSecretKey()
	return sk, keygen.GenPublicKey(sk)
}

// GenKeyPairSparse generates a new SecretKey with exactly hw non zero coefficients [1/2, 0, 1/2].
func (keygen *keyGenerator) GenKeyPairSparse(hw int) (sk *SecretKey, pk *PublicKey) {
	sk = keygen.GenSecretKeySparse(hw)
	return sk, keygen.GenPublicKey(sk)
}

// GenRelinKey generates a new EvaluationKey that will be used to relinearize Ciphertexts during multiplication.
func (keygen *keyGenerator) GenRelinearizationKey(sk *SecretKey, maxDegree int) (evk *RelinearizationKey) {

	if keygen.ringP == nil {
		panic("modulus P is empty")
	}

	levelQ := keygen.params.QCount() - 1
	levelP := keygen.params.PCount() - 1

	evk = new(RelinearizationKey)
	evk.Keys = make([]*SwitchingKey, maxDegree)
	for i := range evk.Keys {
		evk.Keys[i] = NewSwitchingKey(keygen.params, levelQ, levelP)
	}

	ringQ := keygen.ringQ
	ringP := keygen.ringP

	keygen.poolQ[1].CopyValues(sk.Value[0])
	keygen.poolP[1].CopyValues(sk.Value[1])
	for i := 0; i < maxDegree; i++ {
		ringQ.MulCoeffsMontgomery(keygen.poolQ[1], sk.Value[0], keygen.poolQ[1])
		ringP.MulCoeffsMontgomery(keygen.poolP[1], sk.Value[1], keygen.poolP[1])
		keygen.newSwitchingKey([2]*ring.Poly{keygen.poolQ[1], keygen.poolP[1]}, sk.Value, evk.Keys[i])
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

func (keygen *keyGenerator) genrotKey(sk [2]*ring.Poly, galEl uint64, swk *SwitchingKey) {

	skIn := sk
	skOut := [2]*ring.Poly{keygen.poolQ[1], keygen.poolP[1]}

	index := ring.PermuteNTTIndex(galEl, uint64(keygen.ringQ.N))
	ring.PermuteNTTWithIndexLvl(keygen.params.QCount()-1, skIn[0], index, skOut[0])
	ring.PermuteNTTWithIndexLvl(keygen.params.PCount()-1, skIn[1], index, skOut[1])

	keygen.newSwitchingKey(skIn, skOut, swk)
}

// GenSwitchingKey generates a new key-switching key, that will re-encrypt a Ciphertext encrypted under the input key into the output key.
// If the degree of the output key is larger than the input key, then generates a new key-switching key, that will re-encrypt a ciphertext encrypted
// under skIn of dimension n to a ciphertext encrypted under sKOut of dimension N > n.
// [-a*SkOut + w*P*skIn_{Y^{N/n}} + e, a] in X^{N}
// If the degree of the output key is smaller than the input key, then generates a new key-switching key, that will re-encrypt a ciphertext encrypted
// under skIn of dimension N to a ciphertext encrypted under sKOut of dimension n < N.
// [-a*skOut_{Y^{N/n}} + w*P*skIn + e_{N}, a_{N}] in X^{N}
// The output switching key is always given in max(N, n).
func (keygen *keyGenerator) GenSwitchingKey(levelQ, levelP int, skInput, skOutput *SecretKey) (swk *SwitchingKey) {

	if keygen.params.PCount() == 0 {
		panic("Cannot GenSwitchingKey: modulus P is empty")
	}

	swk = NewSwitchingKey(keygen.params, levelQ, levelP)

	// N -> N
	if len(skInput.Value[0].Coeffs[0]) == len(skOutput.Value[0].Coeffs[0]) {
		keygen.newSwitchingKey(skInput.Value, skOutput.Value, swk)
		// n -> N
	} else if len(skInput.Value[0].Coeffs[0]) < len(skOutput.Value[0].Coeffs[0]) {
		ring.MapSmallDimensionToLargerDimensionNTT(skInput.Value[0], keygen.poolQ[1])
		ring.MapSmallDimensionToLargerDimensionNTT(skInput.Value[1], keygen.poolP[1])
		keygen.newSwitchingKey([2]*ring.Poly{keygen.poolQ[1], keygen.poolP[1]}, skOutput.Value, swk)
		// N -> n
	} else {
		ring.MapSmallDimensionToLargerDimensionNTT(skOutput.Value[0], keygen.poolQ[1])
		ring.MapSmallDimensionToLargerDimensionNTT(skOutput.Value[1], keygen.poolP[1])
		keygen.newSwitchingKey(skInput.Value, [2]*ring.Poly{keygen.poolQ[1], keygen.poolP[1]}, swk)
	}

	return
}

func (keygen *keyGenerator) newSwitchingKey(skIn, skOut [2]*ring.Poly, swk *SwitchingKey) {

	ringQ := keygen.ringQ
	ringP := keygen.ringP

	levelQ := len(swk.Value[0][0][0].Coeffs) - 1
	levelP := len(swk.Value[0][0][1].Coeffs) - 1

	var pBigInt *big.Int
	if levelP == keygen.params.PCount()-1 {
		pBigInt = ringP.ModulusBigint
	} else {
		pBigInt = new(big.Int).SetUint64(ringP.Modulus[0])
		for i := 1; i < levelP+1; i++ {
			pBigInt.Mul(pBigInt, ring.NewUint(ringP.Modulus[i]))
		}
	}

	// Computes P * skIn
	ringQ.MulScalarBigintLvl(levelQ, skIn[0], pBigInt, keygen.poolQ[0])

	alpha := levelP + 1
	beta := int(math.Ceil(float64(levelQ+1) / float64(levelP+1)))

	var index int
	for i := 0; i < beta; i++ {

		// e
		keygen.gaussianSamplerQ.ReadLvl(levelQ, swk.Value[i][0][0])
		extendBasisSmallNormAndCenter(ringQ, ringP, levelP, swk.Value[i][0][0], swk.Value[i][0][1])
		ringQ.NTTLazyLvl(levelQ, swk.Value[i][0][0], swk.Value[i][0][0])
		ringP.NTTLazyLvl(levelP, swk.Value[i][0][1], swk.Value[i][0][1])
		ringQ.MFormLvl(levelQ, swk.Value[i][0][0], swk.Value[i][0][0])
		ringP.MFormLvl(levelP, swk.Value[i][0][1], swk.Value[i][0][1])

		// a (since a is uniform, we consider we already sample it in the NTT and Montgomery domain)
		keygen.uniformSamplerQ.ReadLvl(levelQ, swk.Value[i][1][0])
		keygen.uniformSamplerP.ReadLvl(levelP, swk.Value[i][1][1])

		// e + (skIn * P) * (q_star * q_tild) mod QP
		//
		// q_prod = prod(q[i*alpha+j])
		// q_star = Q/qprod
		// q_tild = q_star^-1 mod q_prod
		//
		// Therefore : (skIn * P) * (q_star * q_tild) = sk*P mod q[i*alpha+j], else 0
		for j := 0; j < alpha; j++ {

			index = i*alpha + j

			// It handles the case where nb pj does not divide nb qi
			if index >= levelQ+1 {
				break
			}

			qi := ringQ.Modulus[index]
			p0tmp := keygen.poolQ[0].Coeffs[index]
			p1tmp := swk.Value[i][0][0].Coeffs[index]

			for w := 0; w < ringQ.N; w++ {
				p1tmp[w] = ring.CRed(p1tmp[w]+p0tmp[w], qi)
			}
		}

		// (skIn * P) * (q_star * q_tild) - a * skOut + e mod QP
		ringQ.MulCoeffsMontgomeryAndSubLvl(levelQ, swk.Value[i][1][0], skOut[0], swk.Value[i][0][0])
		ringP.MulCoeffsMontgomeryAndSubLvl(levelP, swk.Value[i][1][1], skOut[1], swk.Value[i][0][1])
	}
}
