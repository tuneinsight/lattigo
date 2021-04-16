package bfv

import (
	"math/big"

	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/rlwe"
	"github.com/ldsec/lattigo/v2/utils"
)

// KeyGenerator is an interface implementing the methods of the keyGenerator.
type KeyGenerator interface {
	GenSecretKey() (sk *rlwe.SecretKey)
	GenSecretkeyWithDistrib(p float64) (sk *rlwe.SecretKey)
	GenPublicKey(sk *rlwe.SecretKey) (pk *rlwe.PublicKey)
	GenKeyPair() (sk *rlwe.SecretKey, pk *rlwe.PublicKey)
	GenSwitchingKey(skIn, skOut *rlwe.SecretKey) (evk *rlwe.SwitchingKey)
	GenRelinearizationKey(sk *rlwe.SecretKey, maxDegree int) (evk *rlwe.RelinearizationKey)
	GenSwitchingKeyForGalois(galEl uint64, sk *rlwe.SecretKey) (swk *rlwe.SwitchingKey)
	GenRotationKeys(galEls []uint64, sk *rlwe.SecretKey) (rks *rlwe.RotationKeySet)
	GenRotationKeysForRotations(ks []int, includeSwapRow bool, sk *rlwe.SecretKey) (rks *rlwe.RotationKeySet)
	GenRotationKeysForInnerSum(sk *rlwe.SecretKey) (rks *rlwe.RotationKeySet)
}

// keyGenerator is a structure that stores the elements required to create new keys,
// as well as a small memory pool for intermediate values.
type keyGenerator struct {
	params          *Parameters
	ringQP          *ring.Ring
	pBigInt         *big.Int
	polypool        [2]*ring.Poly
	gaussianSampler *ring.GaussianSampler
	uniformSampler  *ring.UniformSampler
}

// NewKeyGenerator creates a new KeyGenerator, from which the secret and public keys, as well as the evaluation,
// rotation and switching keys can be generated.
func NewKeyGenerator(params *Parameters) KeyGenerator {

	var ringQP *ring.Ring
	var err error
	if ringQP, err = ring.NewRing(params.N(), append(params.qi, params.pi...)); err != nil {
		panic(err)
	}

	var pBigInt *big.Int
	if len(params.pi) != 0 {
		pBigInt = ring.NewUint(1)
		for _, pi := range params.pi {
			pBigInt.Mul(pBigInt, ring.NewUint(pi))
		}
	}

	prng, err := utils.NewPRNG()
	if err != nil {
		panic(err)
	}

	return &keyGenerator{
		params:          params.Copy(),
		ringQP:          ringQP,
		pBigInt:         pBigInt,
		polypool:        [2]*ring.Poly{ringQP.NewPoly(), ringQP.NewPoly()},
		gaussianSampler: ring.NewGaussianSampler(prng, ringQP, params.Sigma(), int(6*params.Sigma())),
		uniformSampler:  ring.NewUniformSampler(prng, ringQP),
	}
}

// GenSecretKey creates a new SecretKey with the distribution [1/3, 1/3, 1/3].
func (keygen *keyGenerator) GenSecretKey() (sk *rlwe.SecretKey) {
	return keygen.GenSecretkeyWithDistrib(1.0 / 3)
}

// GenSecretkeyWithDistrib creates a new SecretKey with the distribution [(1-p)/2, p, (1-p)/2].
func (keygen *keyGenerator) GenSecretkeyWithDistrib(p float64) (sk *rlwe.SecretKey) {
	prng, err := utils.NewPRNG()
	if err != nil {
		panic(err)
	}
	ternarySamplerMontgomery := ring.NewTernarySampler(prng, keygen.ringQP, p, true)

	sk = new(rlwe.SecretKey)
	sk.Value = ternarySamplerMontgomery.ReadNew()
	keygen.ringQP.NTT(sk.Value, sk.Value)
	return sk
}

// GenPublicKey generates a new PublicKey from the provided SecretKey.
func (keygen *keyGenerator) GenPublicKey(sk *rlwe.SecretKey) (pk *rlwe.PublicKey) {

	pk = new(rlwe.PublicKey)

	ringQP := keygen.ringQP

	//pk[0] = [-(a*s + e)]
	//pk[1] = [a]

	pk.Value[0] = keygen.gaussianSampler.ReadNew()
	ringQP.NTT(pk.Value[0], pk.Value[0])
	pk.Value[1] = keygen.uniformSampler.ReadNew()

	ringQP.MulCoeffsMontgomeryAndSub(sk.Value, pk.Value[1], pk.Value[0])

	return pk
}

// NewKeyPair generates a new SecretKey with distribution [1/3, 1/3, 1/3] and a corresponding PublicKey.
func (keygen *keyGenerator) GenKeyPair() (sk *rlwe.SecretKey, pk *rlwe.PublicKey) {
	sk = keygen.GenSecretKey()
	return sk, keygen.GenPublicKey(sk)
}

// NewRelinKey generates a new evaluation key from the provided SecretKey. It will be used to relinearize a ciphertext (encrypted under a PublicKey generated from the provided SecretKey)
// of degree > 1 to a ciphertext of degree 1. Max degree is the maximum degree of the ciphertext allowed to relinearize.
func (keygen *keyGenerator) GenRelinearizationKey(sk *rlwe.SecretKey, maxDegree int) (evk *rlwe.RelinearizationKey) {

	if keygen.ringQP == nil {
		panic("modulus P is empty")
	}

	evk = new(rlwe.RelinearizationKey)
	evk.Keys = make([]*rlwe.SwitchingKey, maxDegree)
	for i := range evk.Keys {
		evk.Keys[i] = NewSwitchingKey(keygen.params)
	}

	keygen.polypool[0].Copy(sk.Value) // TODO Remove ?

	ringQP := keygen.ringQP

	keygen.polypool[1].Copy(sk.Value)
	for i := 0; i < maxDegree; i++ {
		ringQP.MulCoeffsMontgomery(keygen.polypool[1], sk.Value, keygen.polypool[1])
		keygen.newSwitchingKey(keygen.polypool[1], sk.Value, evk.Keys[i])
	}

	keygen.polypool[0].Zero()
	keygen.polypool[1].Zero()

	return
}

// GenSwitchingKey generates a new key-switching key, that will allow to re-encrypt under the output-key a ciphertext encrypted under the input-key.
func (keygen *keyGenerator) GenSwitchingKey(skInput, skOutput *rlwe.SecretKey) (swkOut *rlwe.SwitchingKey) {

	if keygen.ringQP == nil {
		panic("modulus P is empty")
	}

	swkOut = NewSwitchingKey(keygen.params)

	keygen.ringQP.Copy(skInput.Value, keygen.polypool[0]) // TODO: remove and pass skInput directly ?
	keygen.newSwitchingKey(keygen.polypool[0], skOutput.Value, swkOut)
	keygen.polypool[0].Zero()
	return
}

func (keygen *keyGenerator) GenSwitchingKeyForGalois(galoisEl uint64, sk *rlwe.SecretKey) (swk *rlwe.SwitchingKey) {
	swk = NewSwitchingKey(keygen.params)
	keygen.genrotKey(sk.Value, keygen.params.InverseGaloisElement(galoisEl), swk)
	return
}

// GenRotationKeys generates a RotationKeySet from a list of galois element corresponding to the desired rotations
// See also GenRotationKeysForRotations.
func (keygen *keyGenerator) GenRotationKeys(galEls []uint64, sk *rlwe.SecretKey) (rks *rlwe.RotationKeySet) {
	keys := make(map[uint64]*rlwe.SwitchingKey)
	for _, galEl := range galEls {
		keys[galEl] = keygen.GenSwitchingKeyForGalois(galEl, sk)
	}
	return &rlwe.RotationKeySet{Keys: keys}
}

// GenRotationKeysForRotations generates a RotationKeySet supporting left rotations by k positions for all k in ks.
// Negative k is equivalent to a right rotation by k positions
// If includeConjugate is true, the resulting set contains the conjugation key.
func (keygen *keyGenerator) GenRotationKeysForRotations(ks []int, includeConjugate bool, sk *rlwe.SecretKey) (rks *rlwe.RotationKeySet) {
	galEls := make([]uint64, len(ks), len(ks)+1)
	for i, k := range ks {
		galEls[i] = keygen.params.GaloisElementForColumnRotationBy(k)
	}
	if includeConjugate {
		galEls = append(galEls, keygen.params.GaloisElementForRowRotation())
	}
	return keygen.GenRotationKeys(galEls, sk)
}

// GenRotationKeysForInnerSum generates a RotationKeySet supporting the InnerSum operation of the Evaluator
func (keygen *keyGenerator) GenRotationKeysForInnerSum(sk *rlwe.SecretKey) (rks *rlwe.RotationKeySet) {
	return keygen.GenRotationKeys(keygen.params.GaloisElementsForRowInnerSum(), sk)
}

func (keygen *keyGenerator) genrotKey(sk *ring.Poly, gen uint64, swkOut *rlwe.SwitchingKey) {

	skIn := sk
	skOut := keygen.polypool[1]

	ring.PermuteNTT(skIn, gen, skOut)

	keygen.newSwitchingKey(skIn, skOut, swkOut)

	keygen.polypool[0].Zero()
	keygen.polypool[1].Zero()

	return
}

func (keygen *keyGenerator) newSwitchingKey(skIn, skOut *ring.Poly, swkOut *rlwe.SwitchingKey) {

	ringQP := keygen.ringQP

	alpha := keygen.params.Alpha()
	beta := keygen.params.Beta()

	var index int

	// delta_sk = skIn - skOut = GaloisEnd(skOut, rotation) - skOut

	ringQP.MulScalarBigint(skIn, keygen.pBigInt, keygen.polypool[0])

	for i := 0; i < beta; i++ {

		// e
		keygen.gaussianSampler.Read(swkOut.Value[i][0])
		ringQP.NTTLazy(swkOut.Value[i][0], swkOut.Value[i][0])
		ringQP.MForm(swkOut.Value[i][0], swkOut.Value[i][0])
		// a
		keygen.uniformSampler.Read(swkOut.Value[i][1])

		// e + skIn * (qiBarre*qiStar) * 2^w
		// (qiBarre*qiStar)%qi = 1, else 0

		for j := 0; j < alpha; j++ {

			index = i*alpha + j

			qi := ringQP.Modulus[index]
			p0tmp := keygen.polypool[0].Coeffs[index]
			p1tmp := swkOut.Value[i][0].Coeffs[index]

			for w := 0; w < ringQP.N; w++ {
				p1tmp[w] = ring.CRed(p1tmp[w]+p0tmp[w], qi)
			}

			// Handles the case where nb pj does not divide nb qi
			if index >= keygen.params.QCount() {
				break
			}

		}

		// skIn * (qiBarre*qiStar) * 2^w - a*sk + e
		ringQP.MulCoeffsMontgomeryAndSub(swkOut.Value[i][1], skOut, swkOut.Value[i][0])
	}

	return
}
