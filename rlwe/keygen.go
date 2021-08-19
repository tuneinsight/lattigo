package rlwe

import (
	"math"
	"math/big"

	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/utils"
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
	GenSwitchingKey(skInput, skOutput *SecretKey) (newevakey *SwitchingKey)
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
	poolQ            [2]*ring.Poly // TODO: as PolyQP ?
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

func (keygen *keyGenerator) genSecretKeyFromSampler(sampler ring.Sampler) *SecretKey {
	ringQP := keygen.params.ringQP
	sk := new(SecretKey)
	sk.Value = *ringQP.NewPoly()
	sampler.Read(sk.Value.Q)
	ringQP.ExtendBasisSmallNormAndCenter(sk.Value.Q, keygen.params.PCount()-1, &sk.Value)
	ringQP.NTT(&sk.Value, &sk.Value)
	ringQP.MForm(&sk.Value, &sk.Value)
	return sk
}

// GenSecretKey generates a new SecretKey with the distribution [1/3, 1/3, 1/3].
func (keygen *keyGenerator) GenSecretKey() (sk *SecretKey) {
	return keygen.GenSecretKeyWithDistrib(1.0 / 3)
}

func (keygen *keyGenerator) GenSecretKeyGaussian() (sk *SecretKey) {
	return keygen.genSecretKeyFromSampler(keygen.gaussianSamplerQ)
}

// GenSecretKeyWithDistrib generates a new SecretKey with the distribution [(p-1)/2, p, (p-1)/2].
func (keygen *keyGenerator) GenSecretKeyWithDistrib(p float64) (sk *SecretKey) {
	prng, err := utils.NewPRNG()
	if err != nil {
		panic(err)
	}
	ternarySamplerMontgomery := ring.NewTernarySampler(prng, keygen.ringQ, p, false)
	return keygen.genSecretKeyFromSampler(ternarySamplerMontgomery)
}

// GenSecretKeySparse generates a new SecretKey with exactly hw non-zero coefficients.
func (keygen *keyGenerator) GenSecretKeySparse(hw int) (sk *SecretKey) {
	prng, err := utils.NewPRNG()
	if err != nil {
		panic(err)
	}
	ternarySamplerMontgomery := ring.NewTernarySamplerSparse(prng, keygen.ringQ, hw, false)
	return keygen.genSecretKeyFromSampler(ternarySamplerMontgomery)
}

// GenPublicKey generates a new public key from the provided SecretKey.
func (keygen *keyGenerator) GenPublicKey(sk *SecretKey) (pk *PublicKey) {

	pk = new(PublicKey)
	ringQP := keygen.params.ringQP

	//pk[0] = [-as + e]
	//pk[1] = [a]
	pk = NewPublicKey(keygen.params)
	keygen.gaussianSamplerQ.Read(pk.Value[0].Q)
	ringQP.ExtendBasisSmallNormAndCenter(pk.Value[0].Q, keygen.params.PCount()-1, &pk.Value[0])
	ringQP.NTT(&pk.Value[0], &pk.Value[0])

	keygen.uniformSamplerQ.Read(pk.Value[1].Q)
	keygen.uniformSamplerP.Read(pk.Value[1].P)

	ringQP.MulCoeffsMontgomeryAndSub(&sk.Value, &pk.Value[1], &pk.Value[0])
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

	keygen.poolQ[1].CopyValues(sk.Value.Q)
	keygen.poolP[1].CopyValues(sk.Value.P)
	for i := 0; i < maxDegree; i++ {
		keygen.ringQ.MulCoeffsMontgomery(keygen.poolQ[1], sk.Value.Q, keygen.poolQ[1])
		keygen.genSwitchingKey(keygen.poolQ[1], &sk.Value, evk.Keys[i])
	}

	return
}

// GenRotationKeys generates a RotationKeySet from a list of galois element corresponding to the desired rotations
// See also GenRotationKeysForRotations.
func (keygen *keyGenerator) GenRotationKeys(galEls []uint64, sk *SecretKey) (rks *RotationKeySet) {
	rks = NewRotationKeySet(keygen.params, galEls)
	for _, galEl := range galEls {
		keygen.genrotKey(&sk.Value, keygen.params.InverseGaloisElement(galEl), rks.Keys[galEl])
	}
	return rks
}

func (keygen *keyGenerator) GenSwitchingKeyForRotationBy(k int, sk *SecretKey) (swk *SwitchingKey) {
	swk = NewSwitchingKey(keygen.params, keygen.params.QCount()-1, keygen.params.PCount()-1)
	galElInv := keygen.params.GaloisElementForColumnRotationBy(-int(k))
	keygen.genrotKey(&sk.Value, galElInv, swk)
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
	keygen.genrotKey(&sk.Value, keygen.params.GaloisElementForRowRotation(), swk)
	return
}

func (keygen *keyGenerator) GenSwitchingKeyForGalois(galoisEl uint64, sk *SecretKey) (swk *SwitchingKey) {
	swk = NewSwitchingKey(keygen.params, keygen.params.QCount()-1, keygen.params.PCount()-1)
	keygen.genrotKey(&sk.Value, keygen.params.InverseGaloisElement(galoisEl), swk)
	return
}

// GenRotationKeysForInnerSum generates a RotationKeySet supporting the InnerSum operation of the Evaluator
func (keygen *keyGenerator) GenRotationKeysForInnerSum(sk *SecretKey) (rks *RotationKeySet) {
	return keygen.GenRotationKeys(keygen.params.GaloisElementsForRowInnerSum(), sk)
}

func (keygen *keyGenerator) genrotKey(sk *ring.PolyQP, galEl uint64, swk *SwitchingKey) {

	skIn := sk
	skOut := &ring.PolyQP{Q: keygen.poolQ[1], P: keygen.poolP[1]}

	index := ring.PermuteNTTIndex(galEl, uint64(keygen.ringQ.N))
	ring.PermuteNTTWithIndexLvl(keygen.params.QCount()-1, skIn.Q, index, skOut.Q)
	ring.PermuteNTTWithIndexLvl(keygen.params.PCount()-1, skIn.P, index, skOut.P)

	keygen.genSwitchingKey(skIn.Q, skOut, swk)
}

// GenSwitchingKey generates a new key-switching key, that will re-encrypt a Ciphertext encrypted under the input key into the output key.
// If the degree of the output key is larger than the input key, then generates a new key-switching key, that will re-encrypt a ciphertext encrypted
// under skIn of dimension n to a ciphertext encrypted under sKOut of dimension N > n.
// [-a*SkOut + w*P*skIn_{Y^{N/n}} + e, a] in X^{N}
// If the degree of the output key is smaller than the input key, then generates a new key-switching key, that will re-encrypt a ciphertext encrypted
// under skIn of dimension N to a ciphertext encrypted under sKOut of dimension n < N.
// [-a*skOut_{Y^{N/n}} + w*P*skIn + e_{N}, a_{N}] in X^{N}
// The output switching key is always given in max(N, n) and in the moduli of the output switching key.
func (keygen *keyGenerator) GenSwitchingKey(skInput, skOutput *SecretKey) (swk *SwitchingKey) {

	if keygen.params.PCount() == 0 {
		panic("Cannot GenSwitchingKey: modulus P is empty")
	}

	swk = NewSwitchingKey(keygen.params, skOutput.Value.Q.Level(), skOutput.Value.P.Level())

	// n -> N
	if len(skInput.Value.Q.Coeffs[0]) > len(skOutput.Value.Q.Coeffs[0]) {

		ring.MapSmallDimensionToLargerDimensionNTT(skOutput.Value.Q, keygen.poolQ[1])
		ring.MapSmallDimensionToLargerDimensionNTT(skOutput.Value.P, keygen.poolP[1])
		keygen.genSwitchingKey(skInput.Value.Q, &ring.PolyQP{Q: keygen.poolQ[1], P: keygen.poolP[1]}, swk)

		// N -> N or N -> n
	} else {

		ring.MapSmallDimensionToLargerDimensionNTT(skInput.Value.Q, keygen.poolQ[1])

		if skInput.Value.Q.Level() < skOutput.Value.Q.Level() {

			ringQ := keygen.ringQ

			ringQ.InvNTTLvl(0, keygen.poolQ[1], keygen.poolQ[0])
			ringQ.InvMFormLvl(0, keygen.poolQ[0], keygen.poolQ[0])

			Q := keygen.ringQ.Modulus[0]
			QHalf := Q >> 1

			polQ := keygen.poolQ[0]
			polP := keygen.poolQ[1]
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
				ring.NTT(polP.Coeffs[i], polP.Coeffs[i], ringQ.N, ringQ.NttPsi[i], ringQ.Modulus[i], ringQ.MredParams[i], ringQ.BredParams[i])
				ring.MFormVec(polP.Coeffs[i], polP.Coeffs[i], ringQ.Modulus[i], ringQ.BredParams[i])
			}
		}

		keygen.genSwitchingKey(keygen.poolQ[1], &skOutput.Value, swk)
	}

	return
}

func (keygen *keyGenerator) genSwitchingKey(skIn *ring.Poly, skOut *ring.PolyQP, swk *SwitchingKey) {

	ringQ := keygen.ringQ
	ringP := keygen.ringP
	ringQP := keygen.params.ringQP

	levelQ := len(swk.Value[0][0].Q.Coeffs) - 1
	levelP := len(swk.Value[0][0].P.Coeffs) - 1

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
	ringQ.MulScalarBigintLvl(levelQ, skIn, pBigInt, keygen.poolQ[0])

	alpha := levelP + 1
	beta := int(math.Ceil(float64(levelQ+1) / float64(levelP+1)))

	var index int
	for i := 0; i < beta; i++ {

		// e
		keygen.gaussianSamplerQ.ReadLvl(levelQ, swk.Value[i][0].Q)
		ringQP.ExtendBasisSmallNormAndCenter(swk.Value[i][0].Q, levelP, &swk.Value[i][0])
		ringQ.NTTLazyLvl(levelQ, swk.Value[i][0].Q, swk.Value[i][0].Q)
		ringP.NTTLazyLvl(levelP, swk.Value[i][0].P, swk.Value[i][0].P)
		ringQ.MFormLvl(levelQ, swk.Value[i][0].Q, swk.Value[i][0].Q)
		ringP.MFormLvl(levelP, swk.Value[i][0].P, swk.Value[i][0].P)

		// a (since a is uniform, we consider we already sample it in the NTT and Montgomery domain)
		keygen.uniformSamplerQ.ReadLvl(levelQ, swk.Value[i][1].Q)
		keygen.uniformSamplerP.ReadLvl(levelP, swk.Value[i][1].P)

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
			p1tmp := swk.Value[i][0].Q.Coeffs[index]

			for w := 0; w < ringQ.N; w++ {
				p1tmp[w] = ring.CRed(p1tmp[w]+p0tmp[w], qi)
			}
		}

		// (skIn * P) * (q_star * q_tild) - a * skOut + e mod QP
		ringQ.MulCoeffsMontgomeryAndSubLvl(levelQ, swk.Value[i][1].Q, skOut.Q, swk.Value[i][0].Q)
		ringP.MulCoeffsMontgomeryAndSubLvl(levelP, swk.Value[i][1].P, skOut.P, swk.Value[i][0].P)
	}
}
