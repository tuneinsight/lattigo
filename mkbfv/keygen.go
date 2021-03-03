package mkbfv

import (
	"github.com/ldsec/lattigo/v2/bfv"
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/utils"
)

// KeyGen generated a secret key, a public key and a relinearization key
// given BFV paramters and the peer id.
func KeyGen(params *bfv.Parameters, peerID uint64) *MKKeys {

	generator := bfv.NewKeyGenerator(params)

	keyBag := new(MKKeys)

	// generate private and public BFV keys

	keyBag.secretKey.key = generator.GenSecretKey()
	keyBag.secretKey.peerID = peerID

	keyBag.publicKey.key = generator.GenPublicKey(keyBag.secretKey.key) //TODO: verify if format is in Rq or Rq^d
	keyBag.publicKey.peerID = peerID

	// generate evaluation key. The evaluation key is also used in the relinearization phase
	keyBag.evalKey = evaluationKeyGen(keyBag.secretKey, keyBag.publicKey, generator, params)

	return keyBag
}

// Symmetric encryption of a single ring element (mu) under the secret key (sk).
func uniEnc(mu *ring.Poly, sk MKSecretKey, pk MKPublicKey, generator bfv.KeyGenerator, params *bfv.Parameters) [3]*ring.Poly {

	random := generator.GenSecretKey() // random element as same distribution as the secret key

	// create an uniform sampler and a gaussian sampler in Rq
	var ringQP *ring.Ring
	var err error
	if ringQP, err = ring.NewRing(params.N(), append(params.Qi(), params.Pi()...)); err != nil {
		panic(err)
	}

	prng, err := utils.NewPRNG()
	if err != nil {
		panic(err)
	}

	uniformSampler := ring.NewUniformSampler(prng, ringQP)
	gaussianSampler := ring.NewGaussianSampler(prng, ringQP, params.Sigma(), uint64(6*params.Sigma()))

	// a  <- setup(1^\lambda)
	// e1 <- sample(\psi^d)
	// e2 <- sample(\psi^d)
	// r  <- sample(\chi)
	// d0 = -sk * d1 + e1 + r * g
	// d1 = U(Rq^d)
	// d2 = r * a + e2 + mu * g

	d1 := uniformSampler.ReadNew() // TODO: ask if format is correct for uniform sampling (Rq^d ?). No NTT on it ?

	e1 := gaussianSampler.ReadNew()
	e2 := gaussianSampler.ReadNew()
	ringQP.NTT(e1, e1)
	ringQP.NTT(e2, e2)

	a := pk.key.Value[1]

	d0 := e1
	d2 := e2

	multiplyByBaseAndAdd(random.Value, params, d0) //TODO: is it the correct way to perform these operations ?
	ringQP.MulCoeffsMontgomeryAndSub(sk.key.Value, d1, d0)
	ringQP.MulCoeffsMontgomeryAndAdd(random.Value, a, d2)
	multiplyByBaseAndAdd(mu, params, d2)

	return [3]*ring.Poly{d0, d1, d2}
}

// Function that multiply a ring element p1 by the decomposition basis and adds it to p2
func multiplyByBaseAndAdd(p1 *ring.Poly, params *bfv.Parameters, p2 *ring.Poly) {

	alpha := params.Alpha()
	beta := params.Beta()

	for j := alpha * beta; j < alpha*(beta+1)-1; j++ { //TODO: ask how to perform this operation

	}

}

// Function used to generate the evaluation key. The evaluation key is the encryption of the secret key under itself using uniEnc
func evaluationKeyGen(sk MKSecretKey, pk MKPublicKey, generator bfv.KeyGenerator, params *bfv.Parameters) MKEvaluationKey {

	return MKEvaluationKey{
		key:    uniEnc(sk.key.Value, sk, pk, generator, params),
		peerID: sk.peerID,
	}
}
