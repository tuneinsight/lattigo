package ckks

import (
	"math"
	"math/big"

	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/rlwe"
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
	GenSwitchingKey(skInput, skOutput *SecretKey) (newevakey *SwitchingKey)
	GenRelinearizationKey(sk *SecretKey) (evakey *RelinearizationKey)
	GenSwitchingKeyForGalois(galEl uint64, sk *SecretKey) (swk *SwitchingKey)

	GenRotationKeys(galEls []uint64, sk *SecretKey) (rks *RotationKeySet)

	GenRotationKeysForRotations(ks []int, includeConjugate bool, sk *SecretKey) (rks *RotationKeySet)

	GenRotationIndexesForBootstrapping(logSlots int, btpParams *BootstrappingParameters) []int

	GenRotationIndexesForInnerSum(batch, n int) []int

	GenRotationIndexesForInnerSumNaive(batch, n int) []int

	GenRotationIndexesForDiagMatrix(matrix *PtDiagMatrix) []int
}

// KeyGenerator is a structure that stores the elements required to create new keys,
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

	var qp *ring.Ring
	var err error
	if qp, err = ring.NewRing(params.N(), append(params.qi, params.pi...)); err != nil {
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
		ringQP:          qp,
		pBigInt:         pBigInt,
		polypool:        [2]*ring.Poly{qp.NewPoly(), qp.NewPoly()},
		gaussianSampler: ring.NewGaussianSampler(prng),
		uniformSampler:  ring.NewUniformSampler(prng, qp),
	}
}

// GenSecretKey generates a new SecretKey with the distribution [1/3, 1/3, 1/3].
func (keygen *keyGenerator) GenSecretKey() (sk *SecretKey) {
	return keygen.GenSecretKeyWithDistrib(1.0 / 3)
}

func (keygen *keyGenerator) GenSecretKeyGaussian() (sk *SecretKey) {
	sk = new(SecretKey)

	sk.Value = keygen.gaussianSampler.ReadNew(keygen.ringQP, keygen.params.sigma, int(6*keygen.params.sigma))
	keygen.ringQP.NTT(sk.Value, sk.Value)
	return sk
}

// GenSecretKeyWithDistrib generates a new SecretKey with the distribution [(p-1)/2, p, (p-1)/2].
func (keygen *keyGenerator) GenSecretKeyWithDistrib(p float64) (sk *SecretKey) {
	prng, err := utils.NewPRNG()
	if err != nil {
		panic(err)
	}
	ternarySamplerMontgomery := ring.NewTernarySampler(prng, keygen.ringQP, p, true)

	sk = new(SecretKey)
	sk.Value = ternarySamplerMontgomery.ReadNew()
	keygen.ringQP.NTT(sk.Value, sk.Value)
	return sk
}

// GenSecretKeySparse generates a new SecretKey with exactly hw non-zero coefficients.
func (keygen *keyGenerator) GenSecretKeySparse(hw int) (sk *SecretKey) {
	prng, err := utils.NewPRNG()
	if err != nil {
		panic(err)
	}
	ternarySamplerMontgomery := ring.NewTernarySamplerSparse(prng, keygen.ringQP, hw, true)

	sk = new(SecretKey)
	sk.Value = ternarySamplerMontgomery.ReadNew()
	keygen.ringQP.NTT(sk.Value, sk.Value)
	return sk
}

// GenPublicKey generates a new public key from the provided SecretKey.
func (keygen *keyGenerator) GenPublicKey(sk *SecretKey) (pk *PublicKey) {

	pk = new(PublicKey)

	ringQP := keygen.ringQP

	//pk[0] = [-(a*s + e)]
	//pk[1] = [a]

	pk.Value[0] = keygen.gaussianSampler.ReadNew(keygen.ringQP, keygen.params.sigma, int(6*keygen.params.sigma))
	ringQP.NTT(pk.Value[0], pk.Value[0])
	pk.Value[1] = keygen.uniformSampler.ReadNew()

	ringQP.MulCoeffsMontgomeryAndSub(sk.Value, pk.Value[1], pk.Value[0])

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
func (keygen *keyGenerator) GenRelinearizationKey(sk *SecretKey) (rlk *RelinearizationKey) {

	if len(keygen.params.pi) == 0 {
		panic("Cannot GenRelinKey: modulus P is empty")
	}

	rlk = NewRelinearizationKey(keygen.params)
	keygen.ringQP.MulCoeffsMontgomery(sk.Value, sk.Value, keygen.polypool[0])
	rlk.Keys[0] = &NewSwitchingKey(keygen.params).SwitchingKey
	keygen.newSwitchingKey(keygen.polypool[0], sk.Value, rlk.Keys[0])
	keygen.polypool[0].Zero()

	return
}

// GenSwitchingKey generates a new key-switching key, that will re-encrypt a Ciphertext encrypted under the input key into the output key.
func (keygen *keyGenerator) GenSwitchingKey(skInput, skOutput *SecretKey) (newevakey *SwitchingKey) {

	if len(keygen.params.pi) == 0 {
		panic("Cannot GenSwitchingKey: modulus P is empty")
	}

	keygen.ringQP.Copy(skInput.Value, keygen.polypool[0])
	newevakey = NewSwitchingKey(keygen.params)
	keygen.newSwitchingKey(keygen.polypool[0], skOutput.Value, &newevakey.SwitchingKey)
	keygen.polypool[0].Zero()
	return
}

func (keygen *keyGenerator) GenSwitchingKeyForGalois(galoisEl uint64, sk *SecretKey) (swk *SwitchingKey) {
	swk = NewSwitchingKey(keygen.params)
	keygen.genrotKey(sk.Value, keygen.params.InverseGaloisElement(galoisEl), &swk.SwitchingKey)
	return
}

func (keygen *keyGenerator) GenSwitchingKeyForRotationBy(k int, sk *SecretKey) (swk *SwitchingKey) {
	swk = NewSwitchingKey(keygen.params)
	galElInv := keygen.params.GaloisElementForColumnRotationBy(-int(k))
	keygen.genrotKey(sk.Value, galElInv, &swk.SwitchingKey)
	return
}

func (keygen *keyGenerator) GenSwitchingKeyForConjugate(sk *SecretKey) (swk *SwitchingKey) {
	swk = NewSwitchingKey(keygen.params)
	keygen.genrotKey(sk.Value, keygen.params.GaloisElementForRowRotation(), &swk.SwitchingKey)
	return
}

func (keygen *keyGenerator) genrotKey(sk *ring.Poly, galEl uint64, swk *rlwe.SwitchingKey) {

	skIn := sk
	skOut := keygen.polypool[1]

	index := ring.PermuteNTTIndex(galEl, uint64(keygen.ringQP.N))
	ring.PermuteNTTWithIndexLvl(keygen.params.QPiCount()-1, skIn, index, skOut)

	keygen.newSwitchingKey(skIn, skOut, swk)

	keygen.polypool[0].Zero()
	keygen.polypool[1].Zero()

	return
}

func (keygen *keyGenerator) newSwitchingKey(skIn, skOut *ring.Poly, swk *rlwe.SwitchingKey) {

	ringQP := keygen.ringQP

	// Computes P * skIn
	ringQP.MulScalarBigint(skIn, keygen.pBigInt, keygen.polypool[0])

	alpha := keygen.params.Alpha()
	beta := keygen.params.Beta()

	var index int
	for i := 0; i < beta; i++ {

		// e

		keygen.gaussianSampler.Read(swk.Value[i][0], keygen.ringQP, keygen.params.sigma, int(6*keygen.params.sigma))
		ringQP.NTTLazy(swk.Value[i][0], swk.Value[i][0])
		ringQP.MForm(swk.Value[i][0], swk.Value[i][0])

		// a (since a is uniform, we consider we already sample it in the NTT and Montgomery domain)
		keygen.uniformSampler.Read(swk.Value[i][1])

		// e + (skIn * P) * (q_star * q_tild) mod QP
		//
		// q_prod = prod(q[i*alpha+j])
		// q_star = Q/qprod
		// q_tild = q_star^-1 mod q_prod
		//
		// Therefore : (skIn * P) * (q_star * q_tild) = sk*P mod q[i*alpha+j], else 0
		for j := 0; j < alpha; j++ {

			index = i*alpha + j

			qi := ringQP.Modulus[index]
			p0tmp := keygen.polypool[0].Coeffs[index]
			p1tmp := swk.Value[i][0].Coeffs[index]

			for w := 0; w < ringQP.N; w++ {
				p1tmp[w] = ring.CRed(p1tmp[w]+p0tmp[w], qi)
			}

			// It handles the case where nb pj does not divide nb qi
			if index >= keygen.params.QiCount() {
				break
			}
		}

		// (skIn * P) * (q_star * q_tild) - a * skOut + e mod QP
		ringQP.MulCoeffsMontgomeryAndSub(swk.Value[i][1], skOut, swk.Value[i][0])
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

// GenRotationIndexesForInnerSumNaive generates the rotation indexes for the
// InnerSumNaive. To be then used with GenRotationKeysForRotations to generate
// the RotationKeySet.
func (keygen *keyGenerator) GenRotationIndexesForInnerSumNaive(batch, n int) (rotations []int) {
	rotations = []int{}
	for i := 1; i < n; i++ {
		rotations = append(rotations, i*batch)
	}
	return
}

// GenRotationIndexesForInnerSum generates the rotation indexes for the
// InnerSum. To be then used with GenRotationKeysForRotations to generate
// the RotationKeySet.
func (keygen *keyGenerator) GenRotationIndexesForInnerSum(batch, n int) (rotations []int) {

	rotations = []int{}
	var k int
	for i := 1; i < n; i <<= 1 {

		k = i
		k *= batch

		if !utils.IsInSliceInt(k, rotations) && k != 0 {
			rotations = append(rotations, k)
		}

		k = n - (n & ((i << 1) - 1))
		k *= batch

		if !utils.IsInSliceInt(k, rotations) && k != 0 {
			rotations = append(rotations, k)
		}
	}

	return
}

// GetRotationIndexForDiagMatrix generates of all the rotations needed for a the multiplication
// with the diagonal plaintext matrix.
func (keygen *keyGenerator) GenRotationIndexesForDiagMatrix(matrix *PtDiagMatrix) []int {
	slots := 1 << matrix.LogSlots

	rotKeyIndex := []int{}

	var index int

	N1 := matrix.N1

	if len(matrix.Vec) < 3 {

		for j := range matrix.Vec {

			if !utils.IsInSliceInt(j, rotKeyIndex) {
				rotKeyIndex = append(rotKeyIndex, j)
			}
		}

	} else {

		for j := range matrix.Vec {

			index = ((j / N1) * N1) & (slots - 1)

			if index != 0 && !utils.IsInSliceInt(index, rotKeyIndex) {
				rotKeyIndex = append(rotKeyIndex, index)
			}

			index = j & (N1 - 1)

			if index != 0 && !utils.IsInSliceInt(index, rotKeyIndex) {
				rotKeyIndex = append(rotKeyIndex, index)
			}
		}
	}

	return rotKeyIndex
}

func addMatrixRotToList(pVec map[int]bool, rotations []int, N1, slots int, repack bool) []int {

	if len(pVec) < 3 {
		for j := range pVec {
			if !utils.IsInSliceInt(j, rotations) {
				rotations = append(rotations, j)
			}
		}
	} else {
		var index int
		for j := range pVec {

			index = (j / N1) * N1

			if repack {
				// Sparse repacking, occurring during the first DFT matrix of the CoeffsToSlots.
				index &= (2*slots - 1)
			} else {
				// Other cases
				index &= (slots - 1)
			}

			if index != 0 && !utils.IsInSliceInt(index, rotations) {
				rotations = append(rotations, index)
			}

			index = j & (N1 - 1)

			if index != 0 && !utils.IsInSliceInt(index, rotations) {
				rotations = append(rotations, index)
			}
		}
	}

	return rotations
}

func (keygen *keyGenerator) GenRotationIndexesForBootstrapping(logSlots int, btpParams *BootstrappingParameters) (rotations []int) {

	// List of the rotation key values to needed for the bootstrapp
	rotations = []int{}

	logN := int(keygen.params.logN)

	slots := 1 << logSlots
	dslots := slots
	if logSlots < logN-1 {
		dslots <<= 1
	}

	//SubSum rotation needed X -> Y^slots rotations
	for i := logSlots; i < logN-1; i++ {
		if !utils.IsInSliceInt(1<<i, rotations) {
			rotations = append(rotations, 1<<i)
		}
	}

	indexCtS := computeBootstrappingDFTIndexMap(logN, logSlots, btpParams.CtSDepth(false), true)

	// Coeffs to Slots rotations
	for _, pVec := range indexCtS {
		N1 := findbestbabygiantstepsplit(pVec, dslots, btpParams.MaxN1N2Ratio)
		rotations = addMatrixRotToList(pVec, rotations, N1, slots, false)
	}

	indexStC := computeBootstrappingDFTIndexMap(logN, logSlots, btpParams.StCDepth(false), false)

	// Slots to Coeffs rotations
	for i, pVec := range indexStC {
		N1 := findbestbabygiantstepsplit(pVec, dslots, btpParams.MaxN1N2Ratio)
		rotations = addMatrixRotToList(pVec, rotations, N1, slots, logSlots < logN-1 && i == 0)
	}

	return
}

func computeBootstrappingDFTIndexMap(logN, logSlots, maxDepth int, forward bool) (rotationMap []map[int]bool) {

	bitreversed := false

	var level, depth, nextLevel int

	level = logSlots

	rotationMap = make([]map[int]bool, maxDepth)

	// We compute the chain of merge in order or reverse order depending if its DFT or InvDFT because
	// the way the levels are collapsed has an impact on the total number of rotations and keys to be
	// stored. Ex. instead of using 255 + 64 plaintext vectors, we can use 127 + 128 plaintext vectors
	// by reversing the order of the merging.
	merge := make([]int, maxDepth)
	for i := 0; i < maxDepth; i++ {

		depth = int(math.Ceil(float64(level) / float64(maxDepth-i)))

		if forward {
			merge[i] = depth
		} else {
			merge[len(merge)-i-1] = depth

		}

		level -= depth
	}

	level = logSlots
	for i := 0; i < maxDepth; i++ {

		if logSlots < logN-1 && !forward && i == 0 {

			// Special initial matrix for the repacking before SlotsToCoeffs
			rotationMap[i] = genWfftRepackIndexMap(logSlots, level)

			// Merges this special initial matrix with the first layer of SlotsToCoeffs DFT
			rotationMap[i] = nextLevelfftIndexMap(rotationMap[i], logSlots, 2<<logSlots, level, forward, bitreversed)

			// Continues the merging with the next layers if the total depth requires it.
			nextLevel = level - 1
			for j := 0; j < merge[i]-1; j++ {
				rotationMap[i] = nextLevelfftIndexMap(rotationMap[i], logSlots, 2<<logSlots, nextLevel, forward, bitreversed)
				nextLevel--
			}

		} else {
			// First layer of the i-th level of the DFT
			rotationMap[i] = genWfftIndexMap(logSlots, level, forward, bitreversed)

			// Merges the layer with the next levels of the DFT if the total depth requires it.
			nextLevel = level - 1
			for j := 0; j < merge[i]-1; j++ {
				rotationMap[i] = nextLevelfftIndexMap(rotationMap[i], logSlots, 1<<logSlots, nextLevel, forward, bitreversed)
				nextLevel--
			}
		}

		level -= merge[i]
	}

	return
}

func genWfftIndexMap(logL, level int, forward, bitreversed bool) (vectors map[int]bool) {

	var rot int

	if forward && !bitreversed || !forward && bitreversed {
		rot = 1 << (level - 1)
	} else {
		rot = 1 << (logL - level)
	}

	vectors = make(map[int]bool)
	vectors[0] = true
	vectors[rot] = true
	vectors[(1<<logL)-rot] = true
	return
}

func genWfftRepackIndexMap(logL, level int) (vectors map[int]bool) {
	vectors = make(map[int]bool)
	vectors[0] = true
	vectors[(1 << logL)] = true
	return
}

func nextLevelfftIndexMap(vec map[int]bool, logL, N, nextLevel int, forward, bitreversed bool) (newVec map[int]bool) {

	var rot int

	newVec = make(map[int]bool)

	if forward && !bitreversed || !forward && bitreversed {
		rot = (1 << (nextLevel - 1)) & (N - 1)
	} else {
		rot = (1 << (logL - nextLevel)) & (N - 1)
	}

	for i := range vec {
		newVec[i] = true
		newVec[(i+rot)&(N-1)] = true
		newVec[(i-rot)&(N-1)] = true
	}

	return
}
