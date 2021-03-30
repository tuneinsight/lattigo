package mkbfv

import (
	"math/big"

	"github.com/ldsec/lattigo/v2/bfv"
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/utils"
)

// GetRingQP generates a RingQP from bfv parameters
func GetRingQP(params *bfv.Parameters) *ring.Ring {
	// create ring
	ringQP := new(ring.Ring)
	var err error
	if ringQP, err = ring.NewRing(params.N(), append(params.Qi(), params.Pi()...)); err != nil {
		panic(err)
	}
	return ringQP
}

// GetRingQ generates a RingQ from bfv parameters
func GetRingQ(params *bfv.Parameters) *ring.Ring {
	// create ring
	ringQ := new(ring.Ring)
	var err error
	if ringQ, err = ring.NewRing(params.N(), params.Qi()); err != nil {
		panic(err)
	}
	return ringQ
}

// GetRingP generates a RingP from bfv parameters
func GetRingP(params *bfv.Parameters) *ring.Ring {
	// create ring
	ringP := new(ring.Ring)
	var err error
	if ringP, err = ring.NewRing(params.N(), params.Pi()); err != nil {
		panic(err)
	}
	return ringP
}

// GetRingQMul generates a ringQMul from bfv parameters
func GetRingQMul(params *bfv.Parameters) *ring.Ring {

	qiMul := ring.GenerateNTTPrimesP(61, 2*params.N(), uint64(len(params.Qi())))

	ringQMul := new(ring.Ring)
	var err error
	if ringQMul, err = ring.NewRing(params.N(), qiMul); err != nil {
		panic(err)
	}

	return ringQMul
}

// KeyGen generated a secret key, a public key and a relinearization key
// given BFV paramters and the vector "a" common to all participants
func KeyGen(params *bfv.Parameters, a *MKDecomposedPoly) *MKKeys {

	// create ring
	ringQP := GetRingQP(params)

	generator := bfv.NewKeyGenerator(params)

	keyBag := new(MKKeys)

	// generate private and public BFV keys
	keyBag.secretKey = new(MKSecretKey)
	keyBag.secretKey.key = generator.GenSecretKey()

	//Public key = (b, a)
	keyBag.publicKey = new(MKPublicKey)
	keyBag.publicKey.key[0] = genPublicKey(keyBag.secretKey.key, params, generator, ringQP, a)
	keyBag.publicKey.key[1] = a

	// generate evaluation key. The evaluation key is also used in the relinearization phase.
	keyBag.evalKey = evaluationKeyGen(keyBag.secretKey, keyBag.publicKey, generator, params, ringQP)

	return keyBag
}

// Generate a public key in Rqp^d given an element a <- U(R_qp^d)
func genPublicKey(sk *bfv.SecretKey, params *bfv.Parameters, generator bfv.KeyGenerator, ringQP *ring.Ring, a *MKDecomposedPoly) *MKDecomposedPoly {

	//value in Rqp^d
	var res *MKDecomposedPoly

	// a <- U(Rqp^d)
	// e <- Gauss(Rqp^d)
	// b <- -s * a + e mod(qp) in Rqp^d

	beta := params.Beta()

	prng, err := utils.NewPRNG()
	if err != nil {
		panic(err)
	}

	res = GetGaussianDecomposed(GetGaussianSampler(params, ringQP, prng), beta) // e in Rqp^d

	for d := uint64(0); d < beta; d++ {
		current := res.poly[d]
		ringQP.NTT(current, current)                                   // Pass ei in NTT
		ringQP.MulCoeffsMontgomeryAndSub(sk.Value, a.poly[d], current) // bi = -s * ai + ei (mod qp)
	}

	return res
}

// Symmetric encryption of a single ring element (mu) under the secret key (sk).
func uniEnc(mu *ring.Poly, sk *MKSecretKey, pk *MKPublicKey, generator bfv.KeyGenerator, params *bfv.Parameters, ringQP *ring.Ring) []*MKDecomposedPoly {

	random := generator.GenSecretKey() // random element as same distribution as the secret key

	prng, err := utils.NewPRNG()
	if err != nil {
		panic(err)
	}
	randomValue := random.Value

	uniformSampler := GetUniformSampler(params, ringQP, prng)
	gaussianSampler := GetGaussianSampler(params, ringQP, prng)

	// a  <- setup(1^\lambda)
	// e1 <- sample(\psi^d)
	// e2 <- sample(\psi^d)
	// r  <- sample(\chi)
	// d0 = -sk * d1 + e1 + r * g
	// d1 = U(Rq^d)
	// d2 = r * a + e2 + mu * g

	// Size of decomposition (d)
	beta := params.Beta()

	d1 := GetUniformDecomposed(uniformSampler, beta)

	d0 := GetGaussianDecomposed(gaussianSampler, beta) // e1 <- Gauss(Rqp^d)
	d2 := GetGaussianDecomposed(gaussianSampler, beta) //e2 <- Gauss(Rqp^d)

	a := pk.key[1] // a <- U(Rqp^d) second component of the public key

	for d := uint64(0); d < beta; d++ {
		// Gaussian is not in NTT, so we convert it to NTT
		ringQP.NTT(d0.poly[d], d0.poly[d]) // pass e1_i in NTT
		ringQP.NTT(d2.poly[d], d2.poly[d]) // pass e2_i in NTT
		ringQP.MulCoeffsMontgomeryAndSub(sk.key.Value, d1.poly[d], d0.poly[d])
		ringQP.MulCoeffsMontgomeryAndAdd(randomValue, a.poly[d], d2.poly[d])
	}

	scaledMu := ringQP.NewPoly()
	scaledRandomValue := ringQP.NewPoly()

	var pBigInt *big.Int
	pis := params.Pi()
	if len(pis) != 0 {
		pBigInt = ring.NewUint(1)
		for _, pi := range pis {
			pBigInt.Mul(pBigInt, ring.NewUint(pi))
		}
	}

	ringQP.MulScalarBigint(randomValue, pBigInt, scaledRandomValue)
	ringQP.MulScalarBigint(mu, pBigInt, scaledMu)
	// the g_is mod q_i are either 0 or 1, so just need to compute sums of the correct random.Values
	MultiplyByBaseAndAdd(scaledRandomValue, params, d0)
	MultiplyByBaseAndAdd(scaledMu, params, d2)

	return []*MKDecomposedPoly{d0, d1, d2}
}

// MultiplyByBaseAndAdd multiplies a ring element p1 by the decomposition basis and adds it to p2
// Takes into account Modulus variant by scaling up by factor P
func MultiplyByBaseAndAdd(p1 *ring.Poly, params *bfv.Parameters, p2 *MKDecomposedPoly) {

	alpha := params.Alpha()
	// dimension of the vectors (d)
	beta := params.Beta()

	var index uint64
	ringQP := GetRingQP(params)

	for i := uint64(0); i < beta; i++ {

		for j := uint64(0); j < alpha; j++ {

			index = i*alpha + j

			qi := ringQP.Modulus[index]
			p0tmp := p1.Coeffs[index]
			p1tmp := p2.poly[i].Coeffs[index]

			for w := uint64(0); w < params.N(); w++ {
				p1tmp[w] = ring.CRed(p1tmp[w]+p0tmp[w], qi) // TODO: confirm code in next meeting (code review)
			}

			// Handles the case where nb pj does not divide nb qi
			if index >= params.QiCount() {
				break
			}
		}
	}
}

// Function used to generate the evaluation key. The evaluation key is the encryption of the secret key under itself using uniEnc
func evaluationKeyGen(sk *MKSecretKey, pk *MKPublicKey, generator bfv.KeyGenerator, params *bfv.Parameters, ringQP *ring.Ring) *MKEvaluationKey {

	return &MKEvaluationKey{
		key:    uniEnc(sk.key.Value, sk, pk, generator, params, ringQP),
		peerID: sk.peerID,
	}
}

// GetGaussianDecomposed samples from a gaussian distribution and build an element of Rq^d
func GetGaussianDecomposed(sampler *ring.GaussianSampler, dimension uint64) *MKDecomposedPoly {
	res := new(MKDecomposedPoly)

	for d := uint64(0); d < dimension; d++ {

		res.poly = append(res.poly, sampler.ReadNew())
	}

	return res
}

// GetUniformDecomposed samples from a uniform distribution and build an element of Rq^d
func GetUniformDecomposed(sampler *ring.UniformSampler, dimension uint64) *MKDecomposedPoly {
	res := new(MKDecomposedPoly)

	for d := uint64(0); d < dimension; d++ {

		res.poly = append(res.poly, sampler.ReadNew())
	}

	return res
}

// GetUniformSampler returns the uniform sampler associated to the given ring and prng
func GetUniformSampler(params *bfv.Parameters, r *ring.Ring, prng *utils.KeyedPRNG) *ring.UniformSampler {

	return ring.NewUniformSampler(prng, r)
}

// GetGaussianSampler returns the gaussian sampler associated to the given ring and prng
func GetGaussianSampler(params *bfv.Parameters, r *ring.Ring, prng *utils.KeyedPRNG) *ring.GaussianSampler {

	return ring.NewGaussianSampler(prng, r, params.Sigma(), uint64(6*params.Sigma()))
}

// GenCommonPublicParam generates the public parameter a <- U(R_qp^d) shared by all peers
func GenCommonPublicParam(params *bfv.Parameters) *MKDecomposedPoly {

	prng, err := utils.NewPRNG()
	if err != nil {
		panic(err)
	}
	ringQP := GetRingQP(params)

	uniformSampler := GetUniformSampler(params, ringQP, prng)

	return GetUniformDecomposed(uniformSampler, params.Beta())
}
