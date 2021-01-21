package bfv

import (
	"math/big"

	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/utils"
)

// KeyGenerator is an interface implementing the methods of the keyGenerator.
type KeyGenerator interface {
	GenSecretKey() (sk *SecretKey)
	GenSecretkeyWithDistrib(p float64) (sk *SecretKey)
	GenPublicKey(sk *SecretKey) (pk *PublicKey)
	GenKeyPair() (sk *SecretKey, pk *PublicKey)
	GenSwitchingKey(skIn, skOut *SecretKey) (evk *SwitchingKey)
	GenRelinKey(sk *SecretKey, maxDegree uint64) (evk *RelinearizationKey)
	GenSwitchingKeyForGalois(galEl uint64, sk *SecretKey) (swk *SwitchingKey)
	GenSwitchingKeyForRotationBy(k int, sk *SecretKey) (swk *SwitchingKey)
	GenSwitchingKeyForRowSwap(sk *SecretKey) (swk *SwitchingKey)
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
		gaussianSampler: ring.NewGaussianSampler(prng, ringQP, params.Sigma(), uint64(6*params.Sigma())),
		uniformSampler:  ring.NewUniformSampler(prng, ringQP),
	}
}

// GenSecretKey creates a new SecretKey with the distribution [1/3, 1/3, 1/3].
func (keygen *keyGenerator) GenSecretKey() (sk *SecretKey) {
	return keygen.GenSecretkeyWithDistrib(1.0 / 3)
}

// GenSecretkeyWithDistrib creates a new SecretKey with the distribution [(1-p)/2, p, (1-p)/2].
func (keygen *keyGenerator) GenSecretkeyWithDistrib(p float64) (sk *SecretKey) {
	prng, err := utils.NewPRNG()
	if err != nil {
		panic(err)
	}
	ternarySamplerMontgomery := ring.NewTernarySampler(prng, keygen.ringQP, p, true)

	sk = new(SecretKey)
	sk.sk = ternarySamplerMontgomery.ReadNew()
	keygen.ringQP.NTT(sk.sk, sk.sk)
	return sk
}

// NewSecretKey generates a new SecretKey with zero values.
func NewSecretKey(params *Parameters) *SecretKey {

	sk := new(SecretKey)
	sk.sk = ring.NewPoly(uint64(1<<params.logN), uint64(len(params.qi)+len(params.pi)))
	return sk
}

// Get returns the polynomial of the target SecretKey.
func (sk *SecretKey) Get() *ring.Poly {
	return sk.sk
}

// Set sets the polynomial of the target secret key as the input polynomial.
func (sk *SecretKey) Set(poly *ring.Poly) {
	sk.sk = poly.CopyNew()
}

// GenPublicKey generates a new PublicKey from the provided SecretKey.
func (keygen *keyGenerator) GenPublicKey(sk *SecretKey) (pk *PublicKey) {

	pk = new(PublicKey)

	ringQP := keygen.ringQP

	//pk[0] = [-(a*s + e)]
	//pk[1] = [a]

	pk.pk[0] = keygen.gaussianSampler.ReadNew()
	ringQP.NTT(pk.pk[0], pk.pk[0])
	pk.pk[1] = keygen.uniformSampler.ReadNew()

	ringQP.MulCoeffsMontgomeryAndSub(sk.sk, pk.pk[1], pk.pk[0])

	return pk
}

// NewPublicKey returns a new PublicKey with zero values.
func NewPublicKey(params *Parameters) (pk *PublicKey) {

	pk = new(PublicKey)

	pk.pk[0] = ring.NewPoly(uint64(1<<params.logN), uint64(len(params.qi)+len(params.pi)))
	pk.pk[1] = ring.NewPoly(uint64(1<<params.logN), uint64(len(params.qi)+len(params.pi)))

	return
}

// Get returns the polynomials of the PublicKey.
func (pk *PublicKey) Get() [2]*ring.Poly {
	return pk.pk
}

// Set sets the polynomial of the PublicKey as the input polynomials.
func (pk *PublicKey) Set(p [2]*ring.Poly) {
	pk.pk[0] = p[0].CopyNew()
	pk.pk[1] = p[1].CopyNew()
}

// NewKeyPair generates a new SecretKey with distribution [1/3, 1/3, 1/3] and a corresponding PublicKey.
func (keygen *keyGenerator) GenKeyPair() (sk *SecretKey, pk *PublicKey) {
	sk = keygen.GenSecretKey()
	return sk, keygen.GenPublicKey(sk)
}

// NewRelinKey generates a new evaluation key from the provided SecretKey. It will be used to relinearize a ciphertext (encrypted under a PublicKey generated from the provided SecretKey)
// of degree > 1 to a ciphertext of degree 1. Max degree is the maximum degree of the ciphertext allowed to relinearize.
func (keygen *keyGenerator) GenRelinKey(sk *SecretKey, maxDegree uint64) (evk *RelinearizationKey) {

	if keygen.ringQP == nil {
		panic("Cannot GenRelinKey: modulus P is empty")
	}

	evk = new(RelinearizationKey)
	evk.keys = make([]*SwitchingKey, maxDegree)
	for i := range evk.keys {
		evk.keys[i] = NewSwitchingKey(keygen.params)
	}

	keygen.polypool[0].Copy(sk.Get()) // TODO Remove ?

	ringQP := keygen.ringQP

	keygen.polypool[1].Copy(sk.Get())
	for i := uint64(0); i < maxDegree; i++ {
		ringQP.MulCoeffsMontgomery(keygen.polypool[1], sk.Get(), keygen.polypool[1])
		keygen.newSwitchingKey(keygen.polypool[1], sk.Get(), evk.keys[i])
	}

	keygen.polypool[0].Zero()
	keygen.polypool[1].Zero()

	return
}

// GenSwitchingKey generates a new key-switching key, that will allow to re-encrypt under the output-key a ciphertext encrypted under the input-key.
func (keygen *keyGenerator) GenSwitchingKey(skInput, skOutput *SecretKey) (swkOut *SwitchingKey) {

	if keygen.ringQP == nil {
		panic("Cannot GenRelinKey: modulus P is empty")
	}

	swkOut = NewSwitchingKey(keygen.params)

	keygen.ringQP.Copy(skInput.Get(), keygen.polypool[0]) // TODO: remove and pass skInput directly ?
	keygen.newSwitchingKey(keygen.polypool[0], skOutput.Get(), swkOut)
	keygen.polypool[0].Zero()
	return
}

// NewSwitchingKey returns a new SwitchingKey with zero values.
func NewSwitchingKey(params *Parameters) (evakey *SwitchingKey) {

	evakey = new(SwitchingKey)

	// delta_sk = skInput - skOutput = GaloisEnd(skOutput, rotation) - skOutput
	evakey.key = make([][2]*ring.Poly, params.Beta())

	for i := uint64(0); i < params.Beta(); i++ {
		evakey.key[i][0] = ring.NewPoly(uint64(1<<params.logN), uint64(len(params.qi)+len(params.pi)))
		evakey.key[i][1] = ring.NewPoly(uint64(1<<params.logN), uint64(len(params.qi)+len(params.pi)))
	}

	return
}

func (keygen *keyGenerator) GenSwitchingKeyForGalois(galoisEl uint64, sk *SecretKey) (swk *SwitchingKey) {
	swk = NewSwitchingKey(keygen.params)
	keygen.genrotKey(sk.sk, keygen.params.InverseGaloisElement(galoisEl), swk)
	return
}

func (keygen *keyGenerator) GenSwitchingKeyForRotationBy(k int, sk *SecretKey) (swk *SwitchingKey) {
	swk = NewSwitchingKey(keygen.params)
	galElInv := keygen.params.GaloisElementForColumnRotationBy(-int(k))
	keygen.genrotKey(sk.sk, galElInv, swk)
	return
}

func (keygen *keyGenerator) GenSwitchingKeyForRowSwap(sk *SecretKey) (swk *SwitchingKey) {
	swk = NewSwitchingKey(keygen.params)
	keygen.genrotKey(sk.sk, keygen.params.GaloisElementForRowRotation(), swk)
	return
}

// GenRotationKeysPow2 generates a new rotation key with all the power-of-two rotations to the left and right, as well as the conjugation.
// func (keygen *keyGenerator) GenRotationKeysPow2(skOutput *SecretKey) (rotKey *RotationKeySet) {

// 	if keygen.ringQP == nil {
// 		panic("Cannot GenRotationKeysPow2: modulus P is empty")
// 	}

// 	rotKey = NewRotationKeySet(keygen.params)

// 	for n := uint64(1); n < 1<<(keygen.params.LogN()-1); n <<= 1 {
// 		keygen.GenRot(RotationLeft, skOutput, n, rotKey)
// 		keygen.GenRot(RotationRight, skOutput, n, rotKey)
// 	}

// 	keygen.GenRot(RotationRow, skOutput, 0, rotKey)

// 	return
// }

func (keygen *keyGenerator) genrotKey(sk *ring.Poly, gen uint64, swkOut *SwitchingKey) {

	skIn := sk
	skOut := keygen.polypool[1]

	ring.PermuteNTT(skIn, gen, skOut)

	keygen.newSwitchingKey(skIn, skOut, swkOut)

	keygen.polypool[0].Zero()
	keygen.polypool[1].Zero()

	return
}

func (keygen *keyGenerator) newSwitchingKey(skIn, skOut *ring.Poly, swkOut *SwitchingKey) {

	ringQP := keygen.ringQP

	alpha := keygen.params.Alpha()
	beta := keygen.params.Beta()

	var index uint64

	// delta_sk = skIn - skOut = GaloisEnd(skOut, rotation) - skOut

	ringQP.MulScalarBigint(skIn, keygen.pBigInt, keygen.polypool[0])

	for i := uint64(0); i < beta; i++ {

		// e
		keygen.gaussianSampler.Read(swkOut.key[i][0])
		ringQP.NTTLazy(swkOut.key[i][0], swkOut.key[i][0])
		ringQP.MForm(swkOut.key[i][0], swkOut.key[i][0])
		// a
		keygen.uniformSampler.Read(swkOut.key[i][1])

		// e + skIn * (qiBarre*qiStar) * 2^w
		// (qiBarre*qiStar)%qi = 1, else 0

		for j := uint64(0); j < alpha; j++ {

			index = i*alpha + j

			qi := ringQP.Modulus[index]
			p0tmp := keygen.polypool[0].Coeffs[index]
			p1tmp := swkOut.key[i][0].Coeffs[index]

			for w := uint64(0); w < ringQP.N; w++ {
				p1tmp[w] = ring.CRed(p1tmp[w]+p0tmp[w], qi)
			}

			// Handles the case where nb pj does not divide nb qi
			if index >= keygen.params.QiCount() {
				break
			}

		}

		// skIn * (qiBarre*qiStar) * 2^w - a*sk + e
		ringQP.MulCoeffsMontgomeryAndSub(swkOut.key[i][1], skOut, swkOut.key[i][0])
	}

	return
}
