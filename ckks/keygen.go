package ckks

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
	GenSecretKeySparse(hw uint64) (sk *SecretKey)
	GenPublicKey(sk *SecretKey) (pk *PublicKey)
	GenKeyPair() (sk *SecretKey, pk *PublicKey)
	GenKeyPairSparse(hw uint64) (sk *SecretKey, pk *PublicKey)
	GenSwitchingKey(skInput, skOutput *SecretKey) (newevakey *SwitchingKey)
	GenRelinearizationKey(sk *SecretKey) (evakey *RelinearizationKey)
	GenSwitchingKeyForGalois(galEl uint64, sk *SecretKey) (swk *SwitchingKey)
	GenRotationKeys(galEls []uint64, sk *SecretKey) (rks *RotationKeySet)
	GenRotationKeysForRotations(ks []int, includeConjugate bool, sk *SecretKey) (rks *RotationKeySet)
	GenRotationKeysForInnerSum(sk *SecretKey) (rks *RotationKeySet)
	GenBootstrappingKey(logSlots uint64, btpParams *BootstrappingParameters, sk *SecretKey) (btpKey *BootstrappingKey)
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
		gaussianSampler: ring.NewGaussianSampler(prng, qp, params.sigma, uint64(6*params.sigma)),
		uniformSampler:  ring.NewUniformSampler(prng, qp),
	}
}

// GenSecretKey generates a new SecretKey with the distribution [1/3, 1/3, 1/3].
func (keygen *keyGenerator) GenSecretKey() (sk *SecretKey) {
	return keygen.GenSecretKeyWithDistrib(1.0 / 3)
}

func (keygen *keyGenerator) GenSecretKeyGaussian() (sk *SecretKey) {
	sk = new(SecretKey)

	sk.sk = keygen.gaussianSampler.ReadNew()
	keygen.ringQP.NTT(sk.sk, sk.sk)
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
	sk.sk = ternarySamplerMontgomery.ReadNew()
	keygen.ringQP.NTT(sk.sk, sk.sk)
	return sk
}

// GenSecretKeySparse generates a new SecretKey with exactly hw non-zero coefficients.
func (keygen *keyGenerator) GenSecretKeySparse(hw uint64) (sk *SecretKey) {
	prng, err := utils.NewPRNG()
	if err != nil {
		panic(err)
	}
	ternarySamplerMontgomery := ring.NewTernarySamplerSparse(prng, keygen.ringQP, hw, true)

	sk = new(SecretKey)
	sk.sk = ternarySamplerMontgomery.ReadNew()
	keygen.ringQP.NTT(sk.sk, sk.sk)
	return sk
}

// GenPublicKey generates a new public key from the provided SecretKey.
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

// GenKeyPair generates a new SecretKey with distribution [1/3, 1/3, 1/3] and a corresponding public key.
func (keygen *keyGenerator) GenKeyPair() (sk *SecretKey, pk *PublicKey) {
	sk = keygen.GenSecretKey()
	return sk, keygen.GenPublicKey(sk)
}

// GenKeyPairSparse generates a new SecretKey with exactly hw non zero coefficients [1/2, 0, 1/2].
func (keygen *keyGenerator) GenKeyPairSparse(hw uint64) (sk *SecretKey, pk *PublicKey) {
	sk = keygen.GenSecretKeySparse(hw)
	return sk, keygen.GenPublicKey(sk)
}

// GenRelinKey generates a new EvaluationKey that will be used to relinearize Ciphertexts during multiplication.
func (keygen *keyGenerator) GenRelinearizationKey(sk *SecretKey) (evakey *RelinearizationKey) {

	if len(keygen.params.pi) == 0 {
		panic("Cannot GenRelinKey: modulus P is empty")
	}

	evakey = new(RelinearizationKey)
	keygen.ringQP.MulCoeffsMontgomery(sk.Get(), sk.Get(), keygen.polypool[0])
	evakey.evakey = NewSwitchingKey(keygen.params)
	keygen.newSwitchingKey(keygen.polypool[0], sk.Get(), evakey.evakey)
	keygen.polypool[0].Zero()

	return
}

// GenSwitchingKey generates a new key-switching key, that will re-encrypt a Ciphertext encrypted under the input key into the output key.
func (keygen *keyGenerator) GenSwitchingKey(skInput, skOutput *SecretKey) (newevakey *SwitchingKey) {

	if len(keygen.params.pi) == 0 {
		panic("Cannot GenSwitchingKey: modulus P is empty")
	}

	keygen.ringQP.Copy(skInput.Get(), keygen.polypool[0])
	newevakey = NewSwitchingKey(keygen.params)
	keygen.newSwitchingKey(keygen.polypool[0], skOutput.Get(), newevakey)
	keygen.polypool[0].Zero()
	return
}

// NewSwitchingKey returns a new SwitchingKey with zero values.
func NewSwitchingKey(params *Parameters) (evakey *SwitchingKey) {

	evakey = new(SwitchingKey)

	// delta_sk = skInput - skOutput = GaloisEnd(skOutput, rotation) - skOutput

	evakey.key = make([][2]*ring.Poly, params.Beta())

	for i := uint64(0); i < params.Beta(); i++ {
		evakey.key[i][0] = params.NewPolyQP()
		evakey.key[i][1] = params.NewPolyQP()
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

func (keygen *keyGenerator) GenSwitchingKeyForConjugate(sk *SecretKey) (swk *SwitchingKey) {
	swk = NewSwitchingKey(keygen.params)
	keygen.genrotKey(sk.sk, keygen.params.GaloisElementForRowRotation(), swk)
	return
}

// GenRot populates the input RotationKeys with a SwitchingKey for the given rotation type and amount.
// func (keygen *keyGenerator) GenRotationKey(rotType Rotation, sk *SecretKey, k uint64, rotKey *RotationKeys) {

// 	if len(keygen.params.pi) == 0 {
// 		panic("Cannot GenRot: modulus P is empty")
// 	}

// 	ringQP := keygen.ringQP

// 	if rotType != Conjugate {

// 		if rotKey.permuteNTTLeftIndex == nil {
// 			rotKey.permuteNTTLeftIndex = make(map[uint64][]uint64)
// 		}

// 		if rotKey.permuteNTTRightIndex == nil {
// 			rotKey.permuteNTTRightIndex = make(map[uint64][]uint64)
// 		}

// 		if _, inMap := rotKey.permuteNTTLeftIndex[k]; !inMap {
// 			rotKey.permuteNTTLeftIndex[k] = ring.PermuteNTTIndex(GaloisGen, k, ringQP.N)
// 		}

// 		if _, inMap := rotKey.permuteNTTRightIndex[k]; !inMap {
// 			rotKey.permuteNTTRightIndex[k] = ring.PermuteNTTIndex(GaloisGen, 2*ringQP.N-k, ringQP.N)
// 		}
// 	}

// 	switch rotType {
// 	case RotationLeft:

// 		if rotKey.evakeyRotColLeft == nil {
// 			rotKey.evakeyRotColLeft = make(map[uint64]*SwitchingKey)
// 		}

// 		if rotKey.evakeyRotColLeft[k] == nil && k != 0 {
// 			rotKey.evakeyRotColLeft[k] = keygen.genrotKey(sk.Get(), rotKey.permuteNTTRightIndex[k])
// 		}

// 	case RotationRight:

// 		if rotKey.evakeyRotColRight == nil {
// 			rotKey.evakeyRotColRight = make(map[uint64]*SwitchingKey)
// 		}

// 		if rotKey.evakeyRotColRight[k] == nil && k != 0 {
// 			rotKey.evakeyRotColRight[k] = keygen.genrotKey(sk.Get(), rotKey.permuteNTTLeftIndex[k])
// 		}

// 	case Conjugate:
// 		rotKey.permuteNTTConjugateIndex = ring.PermuteNTTIndex(2*ringQP.N-1, 1, ringQP.N)
// 		rotKey.evakeyConjugate = keygen.genrotKey(sk.Get(), rotKey.permuteNTTConjugateIndex)
// 	}
// }

func (keygen *keyGenerator) genrotKey(sk *ring.Poly, galEl uint64, swk *SwitchingKey) {

	skIn := sk
	skOut := keygen.polypool[1]

	index := ring.PermuteNTTIndex(galEl, keygen.ringQP.N)
	ring.PermuteNTTWithIndexLvl(keygen.params.QPiCount()-1, skIn, index, skOut)

	keygen.newSwitchingKey(skIn, skOut, swk)

	keygen.polypool[0].Zero()
	keygen.polypool[1].Zero()

	return
}

func (keygen *keyGenerator) newSwitchingKey(skIn, skOut *ring.Poly, swk *SwitchingKey) {

	ringQP := keygen.ringQP

	// Computes P * skIn
	ringQP.MulScalarBigint(skIn, keygen.pBigInt, keygen.polypool[0])

	alpha := keygen.params.Alpha()
	beta := keygen.params.Beta()

	var index uint64
	for i := uint64(0); i < beta; i++ {

		// e
		keygen.gaussianSampler.Read(swk.key[i][0])
		ringQP.NTTLazy(swk.key[i][0], swk.key[i][0])
		ringQP.MForm(swk.key[i][0], swk.key[i][0])

		// a (since a is uniform, we consider we already sample it in the NTT and Montgomery domain)
		keygen.uniformSampler.Read(swk.key[i][1])

		// e + (skIn * P) * (q_star * q_tild) mod QP
		//
		// q_prod = prod(q[i*alpha+j])
		// q_star = Q/qprod
		// q_tild = q_star^-1 mod q_prod
		//
		// Therefore : (skIn * P) * (q_star * q_tild) = sk*P mod q[i*alpha+j], else 0
		for j := uint64(0); j < alpha; j++ {

			index = i*alpha + j

			qi := ringQP.Modulus[index]
			p0tmp := keygen.polypool[0].Coeffs[index]
			p1tmp := swk.key[i][0].Coeffs[index]

			for w := uint64(0); w < ringQP.N; w++ {
				p1tmp[w] = ring.CRed(p1tmp[w]+p0tmp[w], qi)
			}

			// It handles the case where nb pj does not divide nb qi
			if index >= keygen.params.QiCount() {
				break
			}
		}

		// (skIn * P) * (q_star * q_tild) - a * skOut + e mod QP
		ringQP.MulCoeffsMontgomeryAndSub(swk.key[i][1], skOut, swk.key[i][0])
	}

	return
}

// GenRotationKeys generates a RotationKeySet from a list of galois element corresponding to the desired rotations
// See also GenRotationKeysForRotations.
func (keygen *keyGenerator) GenRotationKeys(galEls []uint64, sk *SecretKey) (rks *RotationKeySet) {
	rks = NewRotationKeySet(keygen.params)
	for _, galEl := range galEls {
		rks.keys[galEl] = keygen.GenSwitchingKeyForGalois(galEl, sk)
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

// GenRotationKeysForInnerSum generates a RotationKeySet supporting the InnerSum operation of the Evaluator
func (keygen *keyGenerator) GenRotationKeysForInnerSum(sk *SecretKey) (rks *RotationKeySet) {
	galEls := make([]uint64, keygen.params.logN+1, keygen.params.logN+1)
	galEls[0] = keygen.params.GaloisElementForRowRotation()
	for i := 0; i < int(keygen.params.logN)-1; i++ {
		galEls[i+1] = keygen.params.GaloisElementForColumnRotationBy(1 << i)
	}
	return keygen.GenRotationKeys(galEls, sk)
}

// GenKeys generates the bootstrapping keys
func (keygen *keyGenerator) GenBootstrappingKey(logSlots uint64, btpParams *BootstrappingParameters, sk *SecretKey) (btpKey *BootstrappingKey) {

	rotUint := computeBootstrappingDFTRotationList(keygen.params.logN, logSlots, btpParams)
	rotInt := make([]int, len(rotUint), len(rotUint))
	for i, r := range rotUint {
		rotInt[i] = int(r)
	}

	btpKey = &BootstrappingKey{
		Rlk:  keygen.GenRelinearizationKey(sk),
		Rtks: keygen.GenRotationKeysForRotations(rotInt, true, sk),
	}

	/*
		nbKeys := uint64(len(rotKeyIndex)) + 2 //rot keys + conj key + relin key
		nbPoly := keygen.params.Beta()
		nbCoefficients := 2 * keygen.params.N() * keygen.params.QPiCount()
		bytesPerCoeff := uint64(8)
		log.Println("Switching-Keys size (GB) :", float64(nbKeys*nbPoly*nbCoefficients*bytesPerCoeff)/float64(1000000000), "(", nbKeys, "keys)")
	*/

	return
}

func addMatrixRotToList(pVec map[uint64]bool, rotations []uint64, N1, slots uint64, repack bool) []uint64 {

	var index uint64
	for j := range pVec {

		index = (j / N1) * N1

		if repack {
			// Sparse repacking, occurring during the first DFT matrix of the CoeffsToSlots.
			index &= (2*slots - 1)
		} else {
			// Other cases
			index &= (slots - 1)
		}

		if index != 0 && !utils.IsInSliceUint64(index, rotations) {
			rotations = append(rotations, index)
		}

		index = j & (N1 - 1)

		if index != 0 && !utils.IsInSliceUint64(index, rotations) {
			rotations = append(rotations, index)
		}
	}

	return rotations
}

func computeBootstrappingDFTRotationList(logN, logSlots uint64, btpParams *BootstrappingParameters) (rotKeyIndex []uint64) {

	// List of the rotation key values to needed for the bootstrapp
	rotKeyIndex = []uint64{}

	var slots uint64 = 1 << logSlots
	var dslots uint64 = slots
	if logSlots < logN-1 {
		dslots <<= 1
	}

	//SubSum rotation needed X -> Y^slots rotations
	for i := logSlots; i < logN-1; i++ {
		if !utils.IsInSliceUint64(1<<i, rotKeyIndex) {
			rotKeyIndex = append(rotKeyIndex, 1<<i)
		}
	}

	indexCtS := computeBootstrappingDFTIndexMap(logN, logSlots, btpParams.CtSDepth(), true)

	// Coeffs to Slots rotations
	for _, pVec := range indexCtS {

		N1 := findbestbabygiantstepsplitIndexMap(pVec, dslots, btpParams.MaxN1N2Ratio)

		rotKeyIndex = addMatrixRotToList(pVec, rotKeyIndex, N1, slots, false)
	}

	indexStC := computeBootstrappingDFTIndexMap(logN, logSlots, btpParams.StCDepth(), false)

	// Slots to Coeffs rotations
	for i, pVec := range indexStC {

		N1 := findbestbabygiantstepsplitIndexMap(pVec, dslots, btpParams.MaxN1N2Ratio)

		if logSlots < logN-1 && i == 0 {
			rotKeyIndex = addMatrixRotToList(pVec, rotKeyIndex, N1, slots, true)
		} else {
			rotKeyIndex = addMatrixRotToList(pVec, rotKeyIndex, N1, slots, false)
		}
	}

	return rotKeyIndex
}

func computeBootstrappingDFTIndexMap(logN, logSlots, maxDepth uint64, forward bool) (rotationMap []map[uint64]bool) {

	var level, depth, nextLevel uint64

	level = logSlots

	rotationMap = make([]map[uint64]bool, maxDepth)

	// We compute the chain of merge in order or reverse order depending if its DFT or InvDFT because
	// the way the levels are collapsed has an inpact on the total number of rotations and keys to be
	// stored. Ex. instead of using 255 + 64 plaintext vectors, we can use 127 + 128 plaintext vectors
	// by reversing the order of the merging.
	merge := make([]uint64, maxDepth)
	for i := uint64(0); i < maxDepth; i++ {

		depth = uint64(math.Ceil(float64(level) / float64(maxDepth-i)))

		if forward {
			merge[i] = depth
		} else {
			merge[uint64(len(merge))-i-1] = depth

		}

		level -= depth
	}

	level = logSlots
	for i := uint64(0); i < maxDepth; i++ {

		if logSlots < logN-1 && !forward && i == 0 {

			// Special initial matrix for the repacking before SlotsToCoeffs
			rotationMap[i] = genWfftRepackIndexMap(logSlots, level)

			// Merges this special initial matrix with the first layer of SlotsToCoeffs DFT
			rotationMap[i] = nextLevelfftIndexMap(rotationMap[i], logSlots, 2<<logSlots, level, forward)

			// Continues the merging with the next layers if the total depth requires it.
			nextLevel = level - 1
			for j := uint64(0); j < merge[i]-1; j++ {
				rotationMap[i] = nextLevelfftIndexMap(rotationMap[i], logSlots, 2<<logSlots, nextLevel, forward)
				nextLevel--
			}

		} else {
			// First layer of the i-th level of the DFT
			rotationMap[i] = genWfftIndexMap(logSlots, level, forward)

			// Merges the layer with the next levels of the DFT if the total depth requires it.
			nextLevel = level - 1
			for j := uint64(0); j < merge[i]-1; j++ {
				rotationMap[i] = nextLevelfftIndexMap(rotationMap[i], logSlots, 1<<logSlots, nextLevel, forward)
				nextLevel--
			}
		}

		level -= merge[i]
	}

	return
}

// Finds the best N1*N2 = N for the baby-step giant-step algorithm for matrix multiplication.
func findbestbabygiantstepsplitIndexMap(vector map[uint64]bool, maxN uint64, maxRatio float64) (minN uint64) {

	for N1 := uint64(1); N1 < maxN; N1 <<= 1 {

		index := make(map[uint64][]uint64)

		for key := range vector {

			idx1 := key / N1
			idx2 := key & (N1 - 1)

			if index[idx1] == nil {
				index[idx1] = []uint64{idx2}
			} else {
				index[idx1] = append(index[idx1], idx2)
			}
		}

		if len(index[0]) > 0 {

			hoisted := len(index[0]) - 1
			normal := len(index) - 1

			// The matrice is very sparse already
			if normal == 0 {
				return N1 / 2
			}

			if hoisted > normal {
				// Finds the next split that has a ratio hoisted/normal greater or equal to maxRatio
				for float64(hoisted)/float64(normal) < maxRatio {

					if normal/2 == 0 {
						break
					}
					N1 *= 2
					hoisted = hoisted*2 + 1
					normal = normal / 2
				}
				return N1
			}
		}
	}

	return 1
}

func genWfftIndexMap(logL, level uint64, forward bool) (vectors map[uint64]bool) {

	var rot uint64

	if forward {
		rot = 1 << (level - 1)
	} else {
		rot = 1 << (logL - level)
	}

	vectors = make(map[uint64]bool)
	vectors[0] = true
	vectors[rot] = true
	vectors[(1<<logL)-rot] = true
	return
}

func genWfftRepackIndexMap(logL, level uint64) (vectors map[uint64]bool) {
	vectors = make(map[uint64]bool)
	vectors[0] = true
	vectors[(1 << logL)] = true
	return
}

func nextLevelfftIndexMap(vec map[uint64]bool, logL, N, nextLevel uint64, forward bool) (newVec map[uint64]bool) {

	var rot uint64

	newVec = make(map[uint64]bool)

	if forward {
		rot = (1 << (nextLevel - 1)) & (N - 1)
	} else {
		rot = (1 << (logL - nextLevel)) & (N - 1)
	}

	for i := range vec {
		newVec[i] = true
		newVec[(i+rot)&(N-1)] = true
		newVec[(i+N-rot)&(N-1)] = true
	}

	return
}
