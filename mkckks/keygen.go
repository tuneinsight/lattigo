package mkckks

import (
	"math/big"

	"github.com/ldsec/lattigo/v2/ckks"
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/utils"
)

// GetRingQP generates a RingQP from ckks parameters
func GetRingQP(params *ckks.Parameters) *ring.Ring {
	// create ring
	ringQP := new(ring.Ring)
	var err error
	if ringQP, err = ring.NewRing(params.N(), append(params.Qi(), params.Pi()...)); err != nil {
		panic(err)
	}
	return ringQP
}

// GetRingQ generates a RingQ from ckks parameters
func GetRingQ(params *ckks.Parameters) *ring.Ring {
	// create ring
	ringQ := new(ring.Ring)
	var err error
	if ringQ, err = ring.NewRing(params.N(), params.Qi()); err != nil {
		panic(err)
	}
	return ringQ
}

// GetRingP generates a RingP from ckks parameters
func GetRingP(params *ckks.Parameters) *ring.Ring {
	// create ring
	ringP := new(ring.Ring)
	var err error
	if ringP, err = ring.NewRing(params.N(), params.Pi()); err != nil {
		panic(err)
	}
	return ringP
}

// GetRingQMul generates a ringQMul from ckks parameters
func GetRingQMul(params *ckks.Parameters) *ring.Ring {

	qiMul := ring.GenerateNTTPrimesP(61, 2*params.N(), uint64(len(params.Qi())))

	ringQMul := new(ring.Ring)
	var err error
	if ringQMul, err = ring.NewRing(params.N(), qiMul); err != nil {
		panic(err)
	}

	return ringQMul
}

// KeyGen generated a secret key, a public key and a relinearization key
// given ckks paramters, the peer id and the vector "a" common to all participants
func KeyGen(params *ckks.Parameters, a *MKDecomposedPoly) *MKKeys {

	// create ring
	ringQP := GetRingQP(params)

	generator := ckks.NewKeyGenerator(params)

	keyBag := new(MKKeys)

	// generate private and public ckks keys
	keyBag.secretKey = new(MKSecretKey)
	keyBag.secretKey.key = generator.GenSecretKey()

	//Public key = (b, a)
	keyBag.publicKey = new(MKPublicKey)
	keyBag.publicKey.key[0] = genPublicKey(keyBag.secretKey.key, params, generator, ringQP, a)
	keyBag.publicKey.key[1] = a

	// generate evaluation key. The evaluation key is used in the relinearization phase.
	keyBag.evalKey = evaluationKeyGen(keyBag.secretKey, keyBag.publicKey, generator, params, ringQP)

	return keyBag
}

// Generate a public key in Rqp^d given an element a <- U(R_qp^d)
func genPublicKey(sk *ckks.SecretKey, params *ckks.Parameters, generator ckks.KeyGenerator, ringQP *ring.Ring, a *MKDecomposedPoly) *MKDecomposedPoly {

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
func uniEnc(mu *ring.Poly, sk *MKSecretKey, pk *MKPublicKey, generator ckks.KeyGenerator, params *ckks.Parameters, ringQP *ring.Ring) []*MKDecomposedPoly {

	random := generator.GenSecretKey() // random element as same distribution as the secret key
	randomValue := random.Value

	prng, err := utils.NewPRNG()
	if err != nil {
		panic(err)
	}

	uniformSampler := GetUniformSampler(params, ringQP, prng)
	gaussianSampler := GetGaussianSampler(params, ringQP, prng)

	// a  <- setup(1^\lambda)
	// e1 <- sample(\psi^d)
	// e2 <- sample(\psi^d)
	// r  <- sample(\chi)
	// d0 = -sk * d1 + e1 + p * r * g
	// d1 = U(Rq^d)
	// d2 = r * a + e2 + p * mu * g

	// Size of decomposition (d)
	beta := params.Beta()

	d1 := GetUniformDecomposed(uniformSampler, beta)
	d0 := GetGaussianDecomposed(gaussianSampler, beta) // e1 <- Gauss(Rqp^d)
	d2 := GetGaussianDecomposed(gaussianSampler, beta) //e2 <- Gauss(Rqp^d)

	a := pk.key[1] // a <- U(Rqp^d) second component of the public key

	// multiply by P
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

	for i := uint64(0); i < beta; i++ {
		// Gaussian is not in NTT, so we convert it to NTT
		ringQP.NTTLazy(d0.poly[i], d0.poly[i]) // pass e1_i in NTT
		ringQP.NTTLazy(d2.poly[i], d2.poly[i]) // pass e2_i in NTT

		ringQP.MForm(d0.poly[i], d0.poly[i]) // pass e1_i in MForm
		ringQP.MForm(d2.poly[i], d2.poly[i]) // pass e2_i in Mform

		// the g_is mod q_i are either 0 or 1, so just need to compute sums
		MultiplyByBaseAndAdd(scaledRandomValue, params, d0.poly[i], i)
		MultiplyByBaseAndAdd(scaledMu, params, d2.poly[i], i)

		ringQP.MulCoeffsMontgomeryAndSub(sk.key.Value, d1.poly[i], d0.poly[i])
		ringQP.MulCoeffsMontgomeryAndAdd(randomValue, a.poly[i], d2.poly[i])

	}

	return []*MKDecomposedPoly{d0, d1, d2}
}

// MultiplyByBaseAndAdd multiplies a ring element p1 by the decomposition basis and adds it to p2
func MultiplyByBaseAndAdd(p1 *ring.Poly, params *ckks.Parameters, p2 *ring.Poly, beta uint64) {

	alpha := params.Alpha()

	var index uint64
	ringQP := GetRingQP(params)

	for j := uint64(0); j < alpha; j++ {

		index = beta*alpha + j

		qi := ringQP.Modulus[index]
		p0tmp := p1.Coeffs[index]
		p1tmp := p2.Coeffs[index]

		for w := uint64(0); w < ringQP.N; w++ {
			p1tmp[w] = ring.CRed(p1tmp[w]+p0tmp[w], qi)
		}

		// Handles the case where nb pj does not divide nb qi
		if index >= params.QiCount() {
			break
		}
	}

}

// Function used to generate the evaluation key. The evaluation key is the encryption of the secret key under itself using uniEnc
func evaluationKeyGen(sk *MKSecretKey, pk *MKPublicKey, generator ckks.KeyGenerator, params *ckks.Parameters, ringQ *ring.Ring) *MKEvaluationKey {

	return &MKEvaluationKey{
		key:    uniEnc(sk.key.Value, sk, pk, generator, params, ringQ),
		peerID: sk.peerID,
	}
}

// GaloisEvaluationKeyGen returns a galois evaluation key for a given automorphism
func GaloisEvaluationKeyGen(galEl uint64, sk *MKSecretKey, params *ckks.Parameters) *MKEvalGalKey {

	res := new(MKEvalGalKey)
	res.key = make([]*MKDecomposedPoly, 2)

	ringQP := GetRingQP(params)

	prng, err := utils.NewPRNG()
	if err != nil {
		panic(err)
	}

	h1 := GetUniformDecomposed(GetUniformSampler(params, ringQP, prng), params.Beta())
	h0 := GetGaussianDecomposed(GetGaussianSampler(params, ringQP, prng), params.Beta())

	permutedSecretKey := ringQP.NewPoly()

	index := ring.PermuteNTTIndex(galEl, ringQP.N)
	ring.PermuteNTTWithIndexLvl(params.QPiCount()-1, sk.key.Value, index, permutedSecretKey)

	for i := uint64(0); i < params.Beta(); i++ {
		ringQP.NTT(h0.poly[i], h0.poly[i])
		ringQP.MForm(h0.poly[i], h0.poly[i])

		MultiplyByBaseAndAdd(permutedSecretKey, params, h0.poly[i], i)

		ringQP.MulCoeffsMontgomeryAndSub(sk.key.Value, h1.poly[i], h0.poly[i])
	}

	res.key[0] = h0
	res.key[1] = h1

	return res
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
func GetUniformSampler(params *ckks.Parameters, r *ring.Ring, prng *utils.KeyedPRNG) *ring.UniformSampler {

	return ring.NewUniformSampler(prng, r)
}

// GetGaussianSampler returns the gaussian sampler associated to the given ring and prng
func GetGaussianSampler(params *ckks.Parameters, r *ring.Ring, prng *utils.KeyedPRNG) *ring.GaussianSampler {

	return ring.NewGaussianSampler(prng, r, params.Sigma(), uint64(6*params.Sigma()))
}

// GenCommonPublicParam generates the public parameter a <- U(R_qp^d) shared by all peers
func GenCommonPublicParam(params *ckks.Parameters, prng *utils.KeyedPRNG) *MKDecomposedPoly {

	if prng == nil {
		panic("Uninitialized prng")
	}

	ringQP := GetRingQP(params)

	uniformSampler := GetUniformSampler(params, ringQP, prng)

	return GetUniformDecomposed(uniformSampler, params.Beta())
}

// GenSwitchingKey generates a new key-switching key, that will re-encrypt a Ciphertext encrypted under the input key into the output key.
func GenSwitchingKey(skInput, skOutput *MKSecretKey, params *ckks.Parameters) (switchingKey *MKSwitchingKey) {

	ringQP := GetRingQP(params)
	keygenPool := ringQP.NewPoly()
	ringQP.Copy(skInput.key.Value, keygenPool)
	switchingKey = NewMKSwitchingKey(ringQP, params, 2, skInput.peerID)
	NewSwitchingKey(keygenPool, skOutput.key.Value, switchingKey, params)
	keygenPool.Zero()
	return switchingKey
}

// NewSwitchingKey generates a new switching key based on the secret key input and output, and stores it in swk
func NewSwitchingKey(skIn, skOut *ring.Poly, swk *MKSwitchingKey, params *ckks.Parameters) {

	ringQP := GetRingQP(params)
	var pBigInt *big.Int
	pis := params.Pi()
	if len(pis) != 0 {
		pBigInt = ring.NewUint(1)
		for _, pi := range pis {
			pBigInt.Mul(pBigInt, ring.NewUint(pi))
		}
	}

	// Computes P * skIn
	ringQP.MulScalarBigint(skIn, pBigInt, skIn)

	beta := params.Beta()

	prng, err := utils.NewPRNG()
	if err != nil {
		panic(err)
	}

	uniformSampler := GetUniformSampler(params, ringQP, prng)
	gaussianSampler := GetGaussianSampler(params, ringQP, prng)

	for i := uint64(0); i < beta; i++ {

		// e
		gaussianSampler.Read(swk.key[0].poly[i])
		ringQP.NTTLazy(swk.key[0].poly[i], swk.key[0].poly[i])
		ringQP.MForm(swk.key[0].poly[i], swk.key[0].poly[i])

		// a (since a is uniform, we consider we already sample it in the NTT and Montgomery domain)
		uniformSampler.Read(swk.key[1].poly[i])

		// e + (skIn * P) * (q_star * q_tild) mod QP
		//
		// q_prod = prod(q[i*alpha+j])
		// q_star = Q/qprod
		// q_tild = q_star^-1 mod q_prod
		//
		// Therefore : (skIn * P) * (q_star * q_tild) = sk*P mod q[i*alpha+j], else 0
		MultiplyByBaseAndAdd(skIn, params, swk.key[0].poly[i], i)

		// (skIn * P) * (q_star * q_tild) - a * skOut + e mod QP
		ringQP.MulCoeffsMontgomeryAndSub(swk.key[1].poly[i], skOut, swk.key[0].poly[i])
	}

	return
}
