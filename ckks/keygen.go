package ckks

import (
	"github.com/ldsec/lattigo/ring"
)

// KeyGenerator is a structure that stores the elements required to create new keys,
// as well as a small memory pool for intermediate values.
type KeyGenerator struct {
	ckkscontext *Context
	context     *ring.Context
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

// NewKeyGenerator creates a new keygenerator, from which the secret and public keys, as well as the evaluation,
// rotation and switching keys can be generated.
func (ckkscontext *Context) NewKeyGenerator() (keygen *KeyGenerator) {
	keygen = new(KeyGenerator)
	keygen.ckkscontext = ckkscontext
	keygen.context = ckkscontext.contextKeys
	keygen.polypool = keygen.context.NewPoly()
	return
}

// NewSecretKey generates a new SecretKey with the distribution [1/3, 1/3, 1/3].
func (keygen *KeyGenerator) NewSecretKey() (sk *SecretKey) {
	return keygen.NewSecretKeyWithDistrib(1.0 / 3)
}

// NewSecretKeyWithDistrib generates a new SecretKey with the distribution [(p-1)/2, p, (p-1)/2].
func (keygen *KeyGenerator) NewSecretKeyWithDistrib(p float64) (sk *SecretKey) {
	sk = new(SecretKey)
	sk.sk = keygen.ckkscontext.contextKeys.SampleTernaryMontgomeryNTTNew(p)
	return sk
}

// NewSecretKeySparse generates a new SecretKey with exactly hw non zero coefficients.
func (keygen *KeyGenerator) NewSecretKeySparse(hw uint64) (sk *SecretKey) {
	sk = new(SecretKey)
	sk.sk = keygen.ckkscontext.contextKeys.SampleTernarySparseMontgomeryNTTNew(hw)
	return sk
}

// NewSecretKey generates a new SecretKey with zero values.
func (ckkscontext *Context) NewSecretKey() *SecretKey {
	sk := new(SecretKey)
	sk.sk = ckkscontext.contextKeys.NewPoly()
	return sk
}

// Get returns the SecretKey value of the SecretKey.
func (sk *SecretKey) Get() *ring.Poly {
	return sk.sk
}

// Set sets the value of the SecretKey to the provided value.
func (sk *SecretKey) Set(poly *ring.Poly) {
	sk.sk = poly.CopyNew()
}

// NewPublicKey generates a new public key from the provided SecretKey.
func (keygen *KeyGenerator) NewPublicKey(sk *SecretKey) (pk *PublicKey) {

	pk = new(PublicKey)

	//pk[0] = [-(a*s + e)]
	//pk[1] = [a]
	pk.pk[0] = keygen.ckkscontext.gaussianSampler.SampleNTTNew()
	pk.pk[1] = keygen.context.NewUniformPoly()

	keygen.context.MulCoeffsMontgomeryAndAdd(sk.sk, pk.pk[1], pk.pk[0])
	keygen.context.Neg(pk.pk[0], pk.pk[0])

	return pk
}

// NewPublicKey returns a new PublicKey with zero values.
func (ckkscontext *Context) NewPublicKey() (pk *PublicKey) {
	pk = new(PublicKey)

	pk.pk[0] = ckkscontext.contextKeys.NewPoly()
	pk.pk[1] = ckkscontext.contextKeys.NewPoly()

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

// NewKeyPair generates a new secretkey with distribution [1/3, 1/3, 1/3] and a corresponding public key.
func (keygen *KeyGenerator) NewKeyPair() (sk *SecretKey, pk *PublicKey) {
	sk = keygen.NewSecretKey()
	return sk, keygen.NewPublicKey(sk)
}

// NewKeyPairSparse generates a new secretkey with exactly hw non zero coefficients [1/2, 0, 1/2].
func (keygen *KeyGenerator) NewKeyPairSparse(hw uint64) (sk *SecretKey, pk *PublicKey) {
	sk = keygen.NewSecretKeySparse(hw)
	return sk, keygen.NewPublicKey(sk)
}

// NewRelinKey generates a new evaluation key that will be used to relinearize the ciphertexts during multiplication.
func (keygen *KeyGenerator) NewRelinKey(sk *SecretKey) (evakey *EvaluationKey) {

	evakey = new(EvaluationKey)
	keygen.polypool.Copy(sk.Get())
	keygen.context.MulCoeffsMontgomery(keygen.polypool, sk.Get(), keygen.polypool)
	evakey.evakey = keygen.newSwitchingKey(keygen.polypool, sk.Get())
	keygen.polypool.Zero()

	return
}

// NewRelinKey returns  new EvaluationKey with zero values.
func (ckkscontext *Context) NewRelinKey() (evakey *EvaluationKey) {
	evakey = new(EvaluationKey)
	evakey.evakey = new(SwitchingKey)

	// delta_sk = skInput - skOutput = GaloisEnd(skOutput, rotation) - skOutput
	evakey.evakey.evakey = make([][2]*ring.Poly, ckkscontext.beta)
	for i := uint64(0); i < ckkscontext.beta; i++ {

		evakey.evakey.evakey[i][0] = ckkscontext.contextKeys.NewPoly()
		evakey.evakey.evakey[i][1] = ckkscontext.contextKeys.NewPoly()
	}

	return
}

// Get returns the slice of switchintkeys of the evaluation-key.
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

// NewSwitchingKey generated a new keyswitching key, that will re-encrypt a ciphertext encrypted under the input key to the output key.
func (keygen *KeyGenerator) NewSwitchingKey(skInput, skOutput *SecretKey) (newevakey *SwitchingKey) {
	keygen.context.Sub(skInput.Get(), skOutput.Get(), keygen.polypool)
	newevakey = keygen.newSwitchingKey(keygen.polypool, skOutput.Get())
	keygen.polypool.Zero()
	return
}

// NewSwitchingKey returns a new SwitchingKey with zero values.
func (ckkscontext *Context) NewSwitchingKey() (evakey *SwitchingKey) {
	evakey = new(SwitchingKey)

	// delta_sk = skInput - skOutput = GaloisEnd(skOutput, rotation) - skOutput
	evakey.evakey = make([][2]*ring.Poly, ckkscontext.beta)

	for i := uint64(0); i < ckkscontext.beta; i++ {
		evakey.evakey[i][0] = ckkscontext.contextKeys.NewPoly()
		evakey.evakey[i][1] = ckkscontext.contextKeys.NewPoly()
	}

	return
}

func (keygen *KeyGenerator) newSwitchingKey(skIn, skOut *ring.Poly) (switchingkey *SwitchingKey) {

	switchingkey = new(SwitchingKey)

	context := keygen.ckkscontext.contextKeys

	// Computes P * skIn
	for _, pj := range keygen.ckkscontext.specialprimes {
		context.MulScalar(skIn, pj, skIn)
	}

	alpha := keygen.ckkscontext.alpha
	beta := keygen.ckkscontext.beta

	var index uint64

	switchingkey.evakey = make([][2]*ring.Poly, beta)

	for i := uint64(0); i < beta; i++ {

		// e
		switchingkey.evakey[i][0] = keygen.ckkscontext.gaussianSampler.SampleNTTNew()
		context.MForm(switchingkey.evakey[i][0], switchingkey.evakey[i][0])

		// a (since a is uniform, we consider we already sample it in the NTT and montgomery domain)
		switchingkey.evakey[i][1] = keygen.context.NewUniformPoly()

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

			// Handles the case where nb pj does not divides nb qi
			if index >= keygen.ckkscontext.levels-1 {
				break
			}
		}

		// (skIn * P) * (q_star * q_tild) - a * skOut + e mod QP
		context.MulCoeffsMontgomeryAndSub(switchingkey.evakey[i][1], skOut, switchingkey.evakey[i][0])
	}

	return
}

// NewRotationKeys generates a new instance of rotationkeys, with the provided rotation to the left, right and conjugation if asked.
func (ckkscontext *Context) NewRotationKeys() (rotKey *RotationKeys) {
	rotKey = new(RotationKeys)
	return
}

// GenRot populates input RotationKeys with a SwitchingKey for the given rotation type and amount.
func (keygen *KeyGenerator) GenRot(rotType Rotation, sk *SecretKey, k uint64, rotKey *RotationKeys) {
	switch rotType {
	case RotationLeft:

		if rotKey.evakeyRotColLeft == nil {
			rotKey.evakeyRotColLeft = make(map[uint64]*SwitchingKey)
		}

		if rotKey.permuteNTTLeftIndex == nil {
			rotKey.permuteNTTLeftIndex = make(map[uint64][]uint64)
		}

		if rotKey.evakeyRotColLeft[k] == nil && k != 0 {
			rotKey.permuteNTTLeftIndex[k] = ring.PermuteNTTIndex(keygen.ckkscontext.galElRotColLeft[k], 1<<keygen.ckkscontext.logN)
			rotKey.evakeyRotColLeft[k] = keygen.genrotKey(sk.Get(), keygen.ckkscontext.galElRotColLeft[k])
		}

	case RotationRight:

		if rotKey.evakeyRotColRight == nil {
			rotKey.evakeyRotColRight = make(map[uint64]*SwitchingKey)
		}

		if rotKey.permuteNTTRightIndex == nil {
			rotKey.permuteNTTRightIndex = make(map[uint64][]uint64)
		}

		if rotKey.evakeyRotColRight[k] == nil && k != 0 {
			rotKey.permuteNTTRightIndex[k] = ring.PermuteNTTIndex(keygen.ckkscontext.galElRotColRight[k], 1<<keygen.ckkscontext.logN)
			rotKey.evakeyRotColRight[k] = keygen.genrotKey(sk.Get(), keygen.ckkscontext.galElRotColRight[k])
		}

	case Conjugate:
		rotKey.permuteNTTConjugateIndex = ring.PermuteNTTIndex(keygen.ckkscontext.galElRotRow, 1<<keygen.ckkscontext.logN)
		rotKey.evakeyConjugate = keygen.genrotKey(sk.Get(), keygen.ckkscontext.galElRotRow)
	}
}

// NewRotationKeysPow2 generates a new rotation key with all the power of two rotation to the left and right, as well as the conjugation.
func (keygen *KeyGenerator) NewRotationKeysPow2(skOutput *SecretKey) (rotKey *RotationKeys) {

	rotKey = new(RotationKeys)

	rotKey.evakeyRotColLeft = make(map[uint64]*SwitchingKey)
	rotKey.evakeyRotColRight = make(map[uint64]*SwitchingKey)

	rotKey.permuteNTTLeftIndex = make(map[uint64][]uint64)
	rotKey.permuteNTTRightIndex = make(map[uint64][]uint64)

	for n := uint64(1); n < keygen.ckkscontext.n>>1; n <<= 1 {

		rotKey.permuteNTTLeftIndex[n] = ring.PermuteNTTIndex(keygen.ckkscontext.galElRotColLeft[n], 1<<keygen.ckkscontext.logN)
		rotKey.permuteNTTRightIndex[n] = ring.PermuteNTTIndex(keygen.ckkscontext.galElRotColRight[n], 1<<keygen.ckkscontext.logN)

		rotKey.evakeyRotColLeft[n] = keygen.genrotKey(skOutput.Get(), keygen.ckkscontext.galElRotColLeft[n])
		rotKey.evakeyRotColRight[n] = keygen.genrotKey(skOutput.Get(), keygen.ckkscontext.galElRotColRight[n])
	}

	rotKey.permuteNTTConjugateIndex = ring.PermuteNTTIndex(keygen.ckkscontext.galElRotRow, 1<<keygen.ckkscontext.logN)
	rotKey.evakeyConjugate = keygen.genrotKey(skOutput.Get(), keygen.ckkscontext.galElRotRow)
	return
}

// SetRotKey sets the target RotationKeys' SwitchingKey for the specified rotation type and amount with the input polynomials.
func (ckkscontext *Context) SetRotKey(evakey [][2]*ring.Poly, rotType Rotation, k uint64, rotKey *RotationKeys) {
	switch rotType {
	case RotationLeft:

		if rotKey.evakeyRotColLeft == nil {
			rotKey.evakeyRotColLeft = make(map[uint64]*SwitchingKey)
		}

		if rotKey.permuteNTTLeftIndex == nil {
			rotKey.permuteNTTLeftIndex = make(map[uint64][]uint64)
		}

		if rotKey.evakeyRotColLeft[k] == nil && k != 0 {

			rotKey.permuteNTTLeftIndex[k] = ring.PermuteNTTIndex(ckkscontext.galElRotColLeft[k], 1<<ckkscontext.logN)

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

		if rotKey.permuteNTTLeftIndex == nil {
			rotKey.permuteNTTRightIndex = make(map[uint64][]uint64)
		}

		if rotKey.evakeyRotColRight[k] == nil && k != 0 {

			rotKey.permuteNTTRightIndex[k] = ring.PermuteNTTIndex(ckkscontext.galElRotColRight[k], 1<<ckkscontext.logN)

			rotKey.evakeyRotColRight[k] = new(SwitchingKey)
			rotKey.evakeyRotColRight[k].evakey = make([][2]*ring.Poly, len(evakey))
			for j := range evakey {
				rotKey.evakeyRotColRight[k].evakey[j][0] = evakey[j][0].CopyNew()
				rotKey.evakeyRotColRight[k].evakey[j][1] = evakey[j][1].CopyNew()
			}
		}

	case Conjugate:

		if rotKey.evakeyConjugate == nil {

			rotKey.permuteNTTConjugateIndex = ring.PermuteNTTIndex(ckkscontext.galElRotRow, 1<<ckkscontext.logN)

			rotKey.evakeyConjugate = new(SwitchingKey)
			rotKey.evakeyConjugate.evakey = make([][2]*ring.Poly, len(evakey))
			for j := range evakey {
				rotKey.evakeyConjugate.evakey[j][0] = evakey[j][0].CopyNew()
				rotKey.evakeyConjugate.evakey[j][1] = evakey[j][1].CopyNew()
			}
		}
	}
}

func (keygen *KeyGenerator) genrotKey(skOutput *ring.Poly, gen uint64) (switchingkey *SwitchingKey) {

	ring.PermuteNTT(skOutput, gen, keygen.polypool)
	keygen.context.Sub(keygen.polypool, skOutput, keygen.polypool)
	switchingkey = keygen.newSwitchingKey(keygen.polypool, skOutput)
	keygen.polypool.Zero()

	return
}
