package mkrlwe

import (
	"math/big"

	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/rlwe"
	"github.com/ldsec/lattigo/v2/utils"
)

// GetRingQP generates a RingQP from rlwe parameters
func GetRingQP(params *rlwe.Parameters) *ring.Ring {
	// create ring
	ringQP := new(ring.Ring)
	var err error
	if ringQP, err = ring.NewRing(params.N(), append(params.Q(), params.P()...)); err != nil {
		panic(err)
	}
	return ringQP
}

// GetRingQ generates a RingQ from rlwe parameters
func GetRingQ(params *rlwe.Parameters) *ring.Ring {
	// create ring
	ringQ := new(ring.Ring)
	var err error
	if ringQ, err = ring.NewRing(params.N(), params.Q()); err != nil {
		panic(err)
	}
	return ringQ
}

// GetRingP generates a RingP from rlwe parameters
func GetRingP(params *rlwe.Parameters) *ring.Ring {
	// create ring
	ringP := new(ring.Ring)
	var err error
	if ringP, err = ring.NewRing(params.N(), params.P()); err != nil {
		panic(err)
	}
	return ringP
}

// GetRingQMul generates a ringQMul from rlwe parameters
func GetRingQMul(params *rlwe.Parameters) *ring.Ring {

	qiMul := ring.GenerateNTTPrimesP(61, 2*params.N(), uint64(len(params.Q())))

	ringQMul := new(ring.Ring)
	var err error
	if ringQMul, err = ring.NewRing(params.N(), qiMul); err != nil {
		panic(err)
	}

	return ringQMul
}

// KeyGen generated a secret key, a public key and a relinearization key
// given rlwe paramters, the peer id and the vector "a" common to all participants
func KeyGen(params *rlwe.Parameters, a *MKDecomposedPoly) *MKKeys {

	// create ring
	ringQP := GetRingQP(params)

	keyBag := new(MKKeys)

	// generate private and public mkrlwe keys
	keyBag.SecretKey = new(MKSecretKey)
	keyBag.SecretKey.Key = GenSecretKey(ringQP)

	//Public key = (b, a)
	keyBag.PublicKey = new(MKPublicKey)
	keyBag.PublicKey.Key[0] = genPublicKey(keyBag.SecretKey.Key, params, ringQP, a)
	keyBag.PublicKey.Key[1] = a

	// generate evaluation key. The evaluation key is used in the relinearization phase.
	keyBag.EvalKey = evaluationKeyGen(keyBag.SecretKey, keyBag.PublicKey, params, ringQP)

	return keyBag
}

// GenSecretKey generates a new SecretKey with the distribution [1/3, 1/3, 1/3].
func GenSecretKey(ringQP *ring.Ring) *rlwe.SecretKey {
	return GenSecretKeyWithDistrib(1.0/3, ringQP)
}

// GenSecretKeyWithDistrib generates a new SecretKey with the distribution [(p-1)/2, p, (p-1)/2].
func GenSecretKeyWithDistrib(p float64, ringQP *ring.Ring) (sk *rlwe.SecretKey) {
	prng, err := utils.NewPRNG()
	if err != nil {
		panic(err)
	}
	ternarySamplerMontgomery := ring.NewTernarySampler(prng, ringQP, p, true)

	sk = new(rlwe.SecretKey)
	sk.Value = ternarySamplerMontgomery.ReadNew()
	ringQP.NTT(sk.Value, sk.Value)
	return sk
}

// Generate a public key in Rqp^d given an element a <- U(R_qp^d)
func genPublicKey(sk *rlwe.SecretKey, params *rlwe.Parameters, ringQP *ring.Ring, a *MKDecomposedPoly) *MKDecomposedPoly {

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
		current := res.Poly[d]
		ringQP.NTT(current, current)                                   // Pass ei in NTT
		ringQP.MulCoeffsMontgomeryAndSub(sk.Value, a.Poly[d], current) // bi = -s * ai + ei (mod qp)
	}

	return res
}

// Symmetric encryption of a single ring element (mu) under the secret key (sk).
// the output is not in MForm
func uniEnc(mu *ring.Poly, sk *MKSecretKey, pk *MKPublicKey, params *rlwe.Parameters, ringQP *ring.Ring) []*MKDecomposedPoly {

	random := GenSecretKey(ringQP) // random element as same distribution as the secret key
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

	a := pk.Key[1] // a <- U(Rqp^d) second component of the public key

	// multiply by P
	scaledMu := ringQP.NewPoly()
	scaledRandomValue := ringQP.NewPoly()

	var pBigInt *big.Int
	pis := params.P()
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
		ringQP.NTTLazy(d0.Poly[i], d0.Poly[i]) // pass e1_i in NTT
		ringQP.NTTLazy(d2.Poly[i], d2.Poly[i]) // pass e2_i in NTT

		ringQP.MForm(d0.Poly[i], d0.Poly[i]) // pass e1_i in MForm
		ringQP.MForm(d2.Poly[i], d2.Poly[i]) // pass e2_i in MForm

		// the g_is mod q_i are either 0 or 1, so just need to compute sums
		MultiplyByBaseAndAdd(scaledRandomValue, params, d0.Poly[i], i)
		MultiplyByBaseAndAdd(scaledMu, params, d2.Poly[i], i)

		ringQP.InvMForm(d0.Poly[i], d0.Poly[i])
		ringQP.InvMForm(d2.Poly[i], d2.Poly[i])

		ringQP.MulCoeffsMontgomeryAndSub(sk.Key.Value, d1.Poly[i], d0.Poly[i])
		ringQP.MulCoeffsMontgomeryAndAdd(randomValue, a.Poly[i], d2.Poly[i])

	}

	return []*MKDecomposedPoly{d0, d1, d2}
}

// MultiplyByBaseAndAdd multiplies a ring element p1 by the decomposition basis and adds it to p2
func MultiplyByBaseAndAdd(p1 *ring.Poly, params *rlwe.Parameters, p2 *ring.Poly, beta uint64) {

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
		if index >= params.QCount() {
			break
		}
	}

}

// Function used to generate the evaluation key. The evaluation key is the encryption of the secret key under itself using uniEnc
func evaluationKeyGen(sk *MKSecretKey, pk *MKPublicKey, params *rlwe.Parameters, ringQ *ring.Ring) *MKEvaluationKey {

	return &MKEvaluationKey{
		Key:    uniEnc(sk.Key.Value, sk, pk, params, ringQ),
		PeerID: sk.PeerID,
	}
}

// GaloisEvaluationKeyGen returns a galois evaluation key for a given automorphism
// the output is not in MForm
func GaloisEvaluationKeyGen(galEl uint64, sk *MKSecretKey, params *rlwe.Parameters) *MKEvalGalKey {

	res := new(MKEvalGalKey)
	res.Key = make([]*MKDecomposedPoly, 2)

	ringQP := GetRingQP(params)

	prng, err := utils.NewPRNG()
	if err != nil {
		panic(err)
	}

	var pBigInt *big.Int
	pis := params.P()
	if len(pis) != 0 {
		pBigInt = ring.NewUint(1)
		for _, pi := range pis {
			pBigInt.Mul(pBigInt, ring.NewUint(pi))
		}
	}

	h1 := GetUniformDecomposed(GetUniformSampler(params, ringQP, prng), params.Beta())
	h0 := GetGaussianDecomposed(GetGaussianSampler(params, ringQP, prng), params.Beta())

	permutedSecretKey := ringQP.NewPoly()
	index := ring.PermuteNTTIndex(galEl, ringQP.N)

	ring.PermuteNTTWithIndexLvl(uint64(len(ringQP.Modulus)-1), sk.Key.Value, index, permutedSecretKey)

	ringQP.MulScalarBigint(permutedSecretKey, pBigInt, permutedSecretKey)

	for i := uint64(0); i < params.Beta(); i++ {
		ringQP.NTTLazy(h0.Poly[i], h0.Poly[i])
		ringQP.MForm(h0.Poly[i], h0.Poly[i])
		MultiplyByBaseAndAdd(permutedSecretKey, params, h0.Poly[i], i)

		ringQP.InvMForm(h0.Poly[i], h0.Poly[i])

		ringQP.MulCoeffsMontgomeryAndSub(sk.Key.Value, h1.Poly[i], h0.Poly[i])

	}

	res.Key[0] = h0
	res.Key[1] = h1

	return res
}

// GetGaussianDecomposed samples from a gaussian distribution and build an element of Rq^d
func GetGaussianDecomposed(sampler *ring.GaussianSampler, dimension uint64) *MKDecomposedPoly {
	res := new(MKDecomposedPoly)

	for d := uint64(0); d < dimension; d++ {

		res.Poly = append(res.Poly, sampler.ReadNew())
	}

	return res
}

// GetUniformDecomposed samples from a uniform distribution and build an element of Rq^d
func GetUniformDecomposed(sampler *ring.UniformSampler, dimension uint64) *MKDecomposedPoly {
	res := new(MKDecomposedPoly)

	for d := uint64(0); d < dimension; d++ {

		res.Poly = append(res.Poly, sampler.ReadNew())
	}

	return res
}

// GetUniformSampler returns the uniform sampler associated to the given ring and prng
func GetUniformSampler(params *rlwe.Parameters, r *ring.Ring, prng *utils.KeyedPRNG) *ring.UniformSampler {

	return ring.NewUniformSampler(prng, r)
}

// GetGaussianSampler returns the gaussian sampler associated to the given ring and prng
func GetGaussianSampler(params *rlwe.Parameters, r *ring.Ring, prng *utils.KeyedPRNG) *ring.GaussianSampler {

	return ring.NewGaussianSampler(prng, r, params.Sigma(), uint64(6*params.Sigma()))
}

// GenCommonPublicParam generates the public parameter a <- U(R_qp^d) shared by all peers
func GenCommonPublicParam(params *rlwe.Parameters, prng *utils.KeyedPRNG) *MKDecomposedPoly {

	if prng == nil {
		panic("Uninitialized prng")
	}

	ringQP := GetRingQP(params)

	uniformSampler := GetUniformSampler(params, ringQP, prng)

	return GetUniformDecomposed(uniformSampler, params.Beta())
}

// GenSwitchingKey generates a new key-switching key, that will re-encrypt a Ciphertext encrypted under the input key into the output key.
func GenSwitchingKey(skInput, skOutput *MKSecretKey, params *rlwe.Parameters) (switchingKey *MKSwitchingKey) {

	ringQP := GetRingQP(params)
	keygenPool := ringQP.NewPoly()
	ringQP.Copy(skInput.Key.Value, keygenPool)
	switchingKey = NewMKSwitchingKey(ringQP, params, 2, skInput.PeerID)
	NewSwitchingKey(keygenPool, skOutput.Key.Value, switchingKey, params)
	keygenPool.Zero()
	return switchingKey
}

// NewSwitchingKey generates a new switching key based on the secret key input and output, and stores it in swk
func NewSwitchingKey(skIn, skOut *ring.Poly, swk *MKSwitchingKey, params *rlwe.Parameters) {

	ringQP := GetRingQP(params)
	var pBigInt *big.Int
	pis := params.P()
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
		gaussianSampler.Read(swk.Key[0].Poly[i])
		ringQP.NTTLazy(swk.Key[0].Poly[i], swk.Key[0].Poly[i])
		ringQP.MForm(swk.Key[0].Poly[i], swk.Key[0].Poly[i])

		// a (since a is uniform, we consider we already sample it in the NTT and Montgomery domain)
		uniformSampler.Read(swk.Key[1].Poly[i])

		// e + (skIn * P) * (q_star * q_tild) mod QP
		//
		// q_prod = prod(q[i*alpha+j])
		// q_star = Q/qprod
		// q_tild = q_star^-1 mod q_prod
		//
		// Therefore : (skIn * P) * (q_star * q_tild) = sk*P mod q[i*alpha+j], else 0
		MultiplyByBaseAndAdd(skIn, params, swk.Key[0].Poly[i], i)

		// (skIn * P) * (q_star * q_tild) - a * skOut + e mod QP
		ringQP.MulCoeffsMontgomeryAndSub(swk.Key[1].Poly[i], skOut, swk.Key[0].Poly[i])
	}

	return
}
