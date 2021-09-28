package ckks

import (
	"github.com/ldsec/lattigo/ring"
	"math"
)

// KeyGenerator is an interface implementing the methods of the KeyGenerator.
type KeyGenerator interface {
	GenSecretKey() (sk *SecretKey)
	GenSecretKeyWithDistrib(p float64) (sk *SecretKey)
	GenSecretKeySparse(hw uint64) (sk *SecretKey)
	GenPublicKey(sk *SecretKey) (pk *PublicKey)
	GenKeyPair() (sk *SecretKey, pk *PublicKey)
	GenKeyPairSparse(hw uint64) (sk *SecretKey, pk *PublicKey)
	GenRelinKey(sk *SecretKey) (evakey *EvaluationKey)
	GenSwitchingKey(skInput, skOutput *SecretKey) (newevakey *SwitchingKey)
	GenRotationKeysPow2(skOutput *SecretKey) (rotKey *RotationKeys)
	GenRot(rotType Rotation, sk *SecretKey, k uint64, rotKey *RotationKeys)
}

// KeyGenerator is a structure that stores the elements required to create new keys,
// as well as a small memory pool for intermediate values.
type keyGenerator struct {
	params      *Parameters
	ckksContext *Context
	ringContext *ring.Context
	polypool    *ring.Poly
}

// SecretKey is a structure that stores the SecretKey
type SecretKey struct {
	sk *ring.Poly
}

// PublicKey is a structure that stores the PublicKey
type PublicKey struct {
	pk [2]*ring.Poly
}

// Rotation is a type used to represent the rotations types.
type Rotation int

// Constants for rotation types
const (
	RotationRight = iota + 1
	RotationLeft
	Conjugate
)

// RotationKeys is a structure that stores the switching-keys required during the homomorphic rotations.
type RotationKeys struct {
	permuteNTTRightIndex     map[uint64][]uint64
	permuteNTTLeftIndex      map[uint64][]uint64
	permuteNTTConjugateIndex []uint64

	evakeyRotColLeft  map[uint64]*SwitchingKey
	evakeyRotColRight map[uint64]*SwitchingKey
	evakeyConjugate   *SwitchingKey
}

// EvaluationKey is a structure that stores the switching-keys required during the relinearization.
type EvaluationKey struct {
	evakey *SwitchingKey
}

// SwitchingKey is a structure that stores the switching-keys required during the key-switching.
type SwitchingKey struct {
	evakey [][2]*ring.Poly
}

// Get returns the switching key backing slice
func (swk *SwitchingKey) Get() [][2]*ring.Poly {
	return swk.evakey
}

// NewKeyGenerator creates a new KeyGenerator, from which the secret and public keys, as well as the evaluation,
// rotation and switching keys can be generated.
func NewKeyGenerator(params *Parameters) KeyGenerator {

	if !params.isValid {
		panic("cannot NewKeyGenerator: parameters are invalid (check if the generation was done properly)")
	}

	ckksContext := newContext(params)
	ringContext := ckksContext.contextQP

	return &keyGenerator{
		params:      params.Copy(),
		ckksContext: ckksContext,
		ringContext: ringContext,
		polypool:    ringContext.NewPoly(),
	}
}

// GenSecretKey generates a new SecretKey with the distribution [1/3, 1/3, 1/3].
func (keygen *keyGenerator) GenSecretKey() (sk *SecretKey) {
	return keygen.GenSecretKeyWithDistrib(1.0 / 3)
}

// GenSecretKeyWithDistrib generates a new SecretKey with the distribution [(p-1)/2, p, (p-1)/2].
func (keygen *keyGenerator) GenSecretKeyWithDistrib(p float64) (sk *SecretKey) {
	sk = new(SecretKey)
	sk.sk = keygen.ckksContext.contextQP.SampleTernaryMontgomeryNTTNew(p)
	return sk
}

// GenSecretKeySparse generates a new SecretKey with exactly hw non-zero coefficients.
func (keygen *keyGenerator) GenSecretKeySparse(hw uint64) (sk *SecretKey) {
	sk = new(SecretKey)
	sk.sk = keygen.ckksContext.contextQP.SampleTernarySparseMontgomeryNTTNew(hw)
	return sk
}

// NewSecretKey generates a new SecretKey with zero values.
func NewSecretKey(params *Parameters) *SecretKey {

	if !params.isValid {
		panic("cannot NewSecretkey: parameters are invalid (check if the generation was done properly)")
	}

	sk := new(SecretKey)
	sk.sk = ring.NewPoly(1<<params.LogN, uint64(len(params.Qi)+len(params.Pi)))
	return sk
}

// Get returns the value of the SecretKey.
func (sk *SecretKey) Get() *ring.Poly {
	return sk.sk
}

// Set sets the value of the SecretKey to the provided value.
func (sk *SecretKey) Set(poly *ring.Poly) {
	sk.sk = poly.CopyNew()
}

// GenPublicKey generates a new public key from the provided SecretKey.
func (keygen *keyGenerator) GenPublicKey(sk *SecretKey) (pk *PublicKey) {

	pk = new(PublicKey)

	//pk[0] = [-(a*s + e)]
	//pk[1] = [a]
	pk.pk[0] = keygen.ckksContext.gaussianSampler.SampleNTTNew()
	pk.pk[1] = keygen.ringContext.NewUniformPoly()

	keygen.ringContext.MulCoeffsMontgomeryAndAdd(sk.sk, pk.pk[1], pk.pk[0])
	keygen.ringContext.Neg(pk.pk[0], pk.pk[0])

	return pk
}

// NewPublicKey returns a new PublicKey with zero values.
func NewPublicKey(params *Parameters) (pk *PublicKey) {

	if !params.isValid {
		panic("cannot NewPublicKey: parameters are invalid (check if the generation was done properly)")
	}

	pk = new(PublicKey)

	pk.pk[0] = ring.NewPoly(1<<params.LogN, uint64(len(params.Qi)+len(params.Pi)))
	pk.pk[1] = ring.NewPoly(1<<params.LogN, uint64(len(params.Qi)+len(params.Pi)))

	return
}

// Get returns the value of the the public key.
func (pk *PublicKey) Get() [2]*ring.Poly {
	return pk.pk
}

// Set sets the value of the public key to the provided value.
func (pk *PublicKey) Set(poly [2]*ring.Poly) {
	pk.pk[0] = poly[0].CopyNew()
	pk.pk[1] = poly[1].CopyNew()
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
func (keygen *keyGenerator) GenRelinKey(sk *SecretKey) (evakey *EvaluationKey) {

	if keygen.ckksContext.contextP == nil {
		panic("Cannot GenRelinKey: modulus P is empty")
	}

	evakey = new(EvaluationKey)
	keygen.polypool.Copy(sk.Get())
	keygen.ringContext.MulCoeffsMontgomery(keygen.polypool, sk.Get(), keygen.polypool)
	evakey.evakey = keygen.newSwitchingKey(keygen.polypool, sk.Get())
	keygen.polypool.Zero()

	return
}

// NewRelinKey returns a new EvaluationKey with zero values.
func NewRelinKey(params *Parameters) (evakey *EvaluationKey) {

	if !params.isValid {
		panic("cannot NewRelinKey: parameters are invalid (check if the generation was done properly)")
	}

	evakey = new(EvaluationKey)
	evakey.evakey = new(SwitchingKey)

	beta := uint64(math.Ceil(float64(len(params.Qi)) / float64(len(params.Pi))))

	// delta_sk = skInput - skOutput = GaloisEnd(skOutput, rotation) - skOutput
	evakey.evakey.evakey = make([][2]*ring.Poly, beta)
	for i := uint64(0); i < beta; i++ {

		evakey.evakey.evakey[i][0] = ring.NewPoly(1<<params.LogN, uint64(len(params.Qi)+len(params.Pi)))
		evakey.evakey.evakey[i][1] = ring.NewPoly(1<<params.LogN, uint64(len(params.Qi)+len(params.Pi)))
	}

	return
}

// Get returns the slice of switching keys of the evaluation-key.
func (evk *EvaluationKey) Get() *SwitchingKey {
	return evk.evakey
}

// Set sets the target Evaluation key with the input polynomials.
func (evk *EvaluationKey) Set(rlk [][2]*ring.Poly) {

	evk.evakey = new(SwitchingKey)
	evk.evakey.evakey = make([][2]*ring.Poly, len(rlk))
	for j := range rlk {
		evk.evakey.evakey[j][0] = rlk[j][0].CopyNew()
		evk.evakey.evakey[j][1] = rlk[j][1].CopyNew()
	}
}

// GenSwitchingKey generates a new key-switching key, that will re-encrypt a Ciphertext encrypted under the input key into the output key.
func (keygen *keyGenerator) GenSwitchingKey(skInput, skOutput *SecretKey) (newevakey *SwitchingKey) {

	if keygen.ckksContext.contextP == nil {
		panic("Cannot GenSwitchingKey: modulus P is empty")
	}

	//keygen.ringContext.Sub(skInput.Get(), skOutput.Get(), keygen.polypool)
	keygen.ringContext.Copy(skInput.Get(), keygen.polypool)
	newevakey = keygen.newSwitchingKey(keygen.polypool, skOutput.Get())
	keygen.polypool.Zero()
	return
}

// NewSwitchingKey returns a new SwitchingKey with zero values.
func NewSwitchingKey(params *Parameters) (evakey *SwitchingKey) {

	if !params.isValid {
		panic("cannot NewSwitchingKey: parameters are invalid (check if the generation was done properly)")
	}

	evakey = new(SwitchingKey)

	beta := uint64(math.Ceil(float64(len(params.Qi)) / float64(len(params.Pi))))

	// delta_sk = skInput - skOutput = GaloisEnd(skOutput, rotation) - skOutput
	evakey.evakey = make([][2]*ring.Poly, beta)

	for i := uint64(0); i < beta; i++ {
		evakey.evakey[i][0] = ring.NewPoly(1<<params.LogN, uint64(len(params.Qi)+len(params.Pi)))
		evakey.evakey[i][1] = ring.NewPoly(1<<params.LogN, uint64(len(params.Qi)+len(params.Pi)))
	}

	return
}

func (keygen *keyGenerator) newSwitchingKey(skIn, skOut *ring.Poly) (switchingkey *SwitchingKey) {

	switchingkey = new(SwitchingKey)

	context := keygen.ckksContext.contextQP

	// Computes P * skIn

	context.MulScalarBigint(skIn, keygen.ckksContext.contextP.ModulusBigint, skIn)

	alpha := keygen.params.Alpha()
	beta := keygen.params.Beta()

	var index uint64

	switchingkey.evakey = make([][2]*ring.Poly, beta)

	for i := uint64(0); i < beta; i++ {

		// e
		switchingkey.evakey[i][0] = keygen.ckksContext.gaussianSampler.SampleNTTNew()
		context.MForm(switchingkey.evakey[i][0], switchingkey.evakey[i][0])

		// a (since a is uniform, we consider we already sample it in the NTT and Montgomery domain)
		switchingkey.evakey[i][1] = keygen.ringContext.NewUniformPoly()

		// e + (skIn * P) * (q_star * q_tild) mod QP
		//
		// q_prod = prod(q[i*alpha+j])
		// q_star = Q/qprod
		// q_tild = q_star^-1 mod q_prod
		//
		// Therefore : (skIn * P) * (q_star * q_tild) = sk*P mod q[i*alpha+j], else 0
		for j := uint64(0); j < alpha; j++ {

			index = i*alpha + j

			qi := context.Modulus[index]
			p0tmp := skIn.Coeffs[index]
			p1tmp := switchingkey.evakey[i][0].Coeffs[index]

			for w := uint64(0); w < context.N; w++ {
				p1tmp[w] = ring.CRed(p1tmp[w]+p0tmp[w], qi)
			}

			// It handles the case where nb pj does not divide nb qi
			if index >= keygen.ckksContext.levels-1 {
				break
			}
		}

		// (skIn * P) * (q_star * q_tild) - a * skOut + e mod QP
		context.MulCoeffsMontgomeryAndSub(switchingkey.evakey[i][1], skOut, switchingkey.evakey[i][0])
	}

	return
}

// NewRotationKeys generates a new instance of RotationKeys, with the provided rotation to the left, right and conjugation if requested.
func NewRotationKeys() (rotKey *RotationKeys) {
	rotKey = new(RotationKeys)
	return
}

// GenRot populates the input RotationKeys with a SwitchingKey for the given rotation type and amount.
func (keygen *keyGenerator) GenRot(rotType Rotation, sk *SecretKey, k uint64, rotKey *RotationKeys) {

	if keygen.ckksContext.contextP == nil {
		panic("Cannot GenRot: modulus P is empty")
	}

	switch rotType {
	case RotationLeft:

		if rotKey.evakeyRotColLeft == nil {
			rotKey.evakeyRotColLeft = make(map[uint64]*SwitchingKey)
		}

		if rotKey.permuteNTTLeftIndex == nil {
			rotKey.permuteNTTLeftIndex = make(map[uint64][]uint64)
		}

		if rotKey.evakeyRotColLeft[k] == nil && k != 0 {
			rotKey.permuteNTTLeftIndex[k] = ring.PermuteNTTIndex(GaloisGen, k, keygen.ringContext.N)
			rotKey.evakeyRotColLeft[k] = keygen.genrotKey(sk.Get(), keygen.ckksContext.galElRotColLeft[k])
		}

	case RotationRight:

		if rotKey.evakeyRotColRight == nil {
			rotKey.evakeyRotColRight = make(map[uint64]*SwitchingKey)
		}

		if rotKey.permuteNTTRightIndex == nil {
			rotKey.permuteNTTRightIndex = make(map[uint64][]uint64)
		}

		if rotKey.evakeyRotColRight[k] == nil && k != 0 {
			rotKey.permuteNTTRightIndex[k] = ring.PermuteNTTIndex(GaloisGen, 2*keygen.ringContext.N-k, keygen.ringContext.N)
			rotKey.evakeyRotColRight[k] = keygen.genrotKey(sk.Get(), keygen.ckksContext.galElRotColRight[k])
		}

	case Conjugate:
		rotKey.permuteNTTConjugateIndex = ring.PermuteNTTIndex(2*keygen.ringContext.N-1, 1, keygen.ringContext.N)
		rotKey.evakeyConjugate = keygen.genrotKey(sk.Get(), keygen.ckksContext.galElConjugate)
	}
}

// GenRotationKeysPow2 generates a new rotation key with all the power-of-two rotations to the left and right, as well as the conjugation.
func (keygen *keyGenerator) GenRotationKeysPow2(skOutput *SecretKey) (rotKey *RotationKeys) {

	if keygen.ckksContext.contextP == nil {
		panic("Cannot GenRotationKeysPow2: modulus P is empty")
	}

	rotKey = new(RotationKeys)

	rotKey.evakeyRotColLeft = make(map[uint64]*SwitchingKey)
	rotKey.evakeyRotColRight = make(map[uint64]*SwitchingKey)

	rotKey.permuteNTTLeftIndex = make(map[uint64][]uint64)
	rotKey.permuteNTTRightIndex = make(map[uint64][]uint64)

	for n := uint64(1); n < 1<<(keygen.params.LogN-1); n <<= 1 {

		rotKey.permuteNTTLeftIndex[n] = ring.PermuteNTTIndex(GaloisGen, n, keygen.ringContext.N)
		rotKey.permuteNTTRightIndex[n] = ring.PermuteNTTIndex(GaloisGen, 2*keygen.ringContext.N-n, keygen.ringContext.N)

		rotKey.evakeyRotColLeft[n] = keygen.genrotKey(skOutput.Get(), keygen.ckksContext.galElRotColLeft[n])
		rotKey.evakeyRotColRight[n] = keygen.genrotKey(skOutput.Get(), keygen.ckksContext.galElRotColRight[n])
	}

	rotKey.permuteNTTConjugateIndex = ring.PermuteNTTIndex(2*keygen.ringContext.N-1, 1, keygen.ringContext.N)
	rotKey.evakeyConjugate = keygen.genrotKey(skOutput.Get(), keygen.ckksContext.galElConjugate)
	return
}

// SetRotKey sets the target RotationKeys' SwitchingKey for the specified rotation type and amount with the input polynomials.
func (rotKey *RotationKeys) SetRotKey(params *Parameters, evakey [][2]*ring.Poly, rotType Rotation, k uint64) {

	if !params.isValid {
		panic("cannot SetRotKey: parameters are invalid (check if the generation was done properly)")
	}

	switch rotType {
	case RotationLeft:

		if rotKey.evakeyRotColLeft == nil {
			rotKey.evakeyRotColLeft = make(map[uint64]*SwitchingKey)
		}

		if rotKey.permuteNTTLeftIndex == nil {
			rotKey.permuteNTTLeftIndex = make(map[uint64][]uint64)
		}

		if rotKey.evakeyRotColLeft[k] == nil && k != 0 {

			rotKey.permuteNTTLeftIndex[k] = ring.PermuteNTTIndex(GaloisGen, k, 1<<params.LogN)

			rotKey.evakeyRotColLeft[k] = new(SwitchingKey)
			rotKey.evakeyRotColLeft[k].evakey = make([][2]*ring.Poly, len(evakey))
			for j := range evakey {
				rotKey.evakeyRotColLeft[k].evakey[j][0] = evakey[j][0].CopyNew()
				rotKey.evakeyRotColLeft[k].evakey[j][1] = evakey[j][1].CopyNew()
			}
		}

	case RotationRight:

		if rotKey.evakeyRotColRight == nil {
			rotKey.evakeyRotColRight = make(map[uint64]*SwitchingKey)
		}

		if rotKey.permuteNTTRightIndex == nil {
			rotKey.permuteNTTRightIndex = make(map[uint64][]uint64)
		}

		if rotKey.evakeyRotColRight[k] == nil && k != 0 {

			rotKey.permuteNTTRightIndex[k] = ring.PermuteNTTIndex(GaloisGen, (2<<params.LogN)-1-k, 1<<params.LogN)

			rotKey.evakeyRotColRight[k] = new(SwitchingKey)
			rotKey.evakeyRotColRight[k].evakey = make([][2]*ring.Poly, len(evakey))
			for j := range evakey {
				rotKey.evakeyRotColRight[k].evakey[j][0] = evakey[j][0].CopyNew()
				rotKey.evakeyRotColRight[k].evakey[j][1] = evakey[j][1].CopyNew()
			}
		}

	case Conjugate:

		if rotKey.evakeyConjugate == nil {

			rotKey.permuteNTTConjugateIndex = ring.PermuteNTTIndex((2<<params.LogN)-1, 1, 1<<params.LogN)

			rotKey.evakeyConjugate = new(SwitchingKey)
			rotKey.evakeyConjugate.evakey = make([][2]*ring.Poly, len(evakey))
			for j := range evakey {
				rotKey.evakeyConjugate.evakey[j][0] = evakey[j][0].CopyNew()
				rotKey.evakeyConjugate.evakey[j][1] = evakey[j][1].CopyNew()
			}
		}
	}
}

func (keygen *keyGenerator) genrotKey(skOutput *ring.Poly, gen uint64) (switchingkey *SwitchingKey) {

	ring.PermuteNTT(skOutput, gen, keygen.polypool)
	switchingkey = keygen.newSwitchingKey(keygen.polypool, skOutput)
	keygen.polypool.Zero()

	return
}
