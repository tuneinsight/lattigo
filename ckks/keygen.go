package ckks

import (
	"github.com/ldsec/lattigo/ring"
)

// Keygenerator is a structure that stores the elements required to create new keys,
// as well as a small memory pool for intermediate values.
type KeyGenerator struct {
	ckksContext *CkksContext
	context     *ring.Context
	polypool    *ring.Poly
}

// Secretkey is a structure that stores the secret-key
type SecretKey struct {
	sk *ring.Poly
}

// Publickey is a structure that stores the public-key
type PublicKey struct {
	pk [2]*ring.Poly
}

type Rotation int

const (
	RotationRight = iota + 1
	RotationLeft
	Conjugate
)

// Rotationkeys is a structure that stores the switching-keys required during the homomorphic rotations.
type RotationKeys struct {
	ckksContext      *CkksContext
	evakey_rot_col_L map[uint64]*SwitchingKey
	evakey_rot_col_R map[uint64]*SwitchingKey
	evakey_rot_row   *SwitchingKey
}

// Evaluationkey is a structure that stores the switching-keys required during the relinearization.
type EvaluationKey struct {
	evakey *SwitchingKey
}

// Switchingkey is a structure that stores the switching-keys required during the key-switching.
type SwitchingKey struct {
	evakey [][2]*ring.Poly
}

// Get returns the switching key backing slice
func (swk *SwitchingKey) Get() [][2]*ring.Poly {
	return swk.evakey
}

// NewKeyGenerator creates a new keygenerator, from which the secret and public keys, as well as the evaluation,
// rotation and switching keys can be generated.
func (ckksContext *CkksContext) NewKeyGenerator() (keygen *KeyGenerator) {
	keygen = new(KeyGenerator)
	keygen.ckksContext = ckksContext
	keygen.context = ckksContext.contextKeys
	keygen.polypool = keygen.context.NewPoly()
	return
}

// NewSecretKey generates a new secret key with the distribution [1/3, 1/3, 1/3].
func (keygen *KeyGenerator) NewSecretKey() (sk *SecretKey) {
	return keygen.NewSecretKeyWithDistrib(1.0 / 3)
}

// NewSecretKey generates a new secret key with the distribution [(p-1)/2, p, (p-1)/2].
func (keygen *KeyGenerator) NewSecretKeyWithDistrib(p float64) (sk *SecretKey) {
	sk = new(SecretKey)
	sk.sk = keygen.ckksContext.contextKeys.SampleTernaryMontgomeryNTTNew(p)
	return sk
}

func (keygen *KeyGenerator) NewSecretKeySparse(hw uint64) (sk *SecretKey) {
	sk = new(SecretKey)
	sk.sk = keygen.ckksContext.contextKeys.SampleTernarySparseMontgomeryNTTNew(hw)
	return sk
}

func (keygen *KeyGenerator) NewSecretKeyEmpty() *SecretKey {
	sk := new(SecretKey)
	sk.sk = keygen.context.NewPoly()
	return sk
}

// Get returns the secret key value of the secret key.
func (sk *SecretKey) Get() *ring.Poly {
	return sk.sk
}

// Set sets the value of the secret key to the provided value.
func (sk *SecretKey) Set(poly *ring.Poly) {
	sk.sk = poly.CopyNew()
}

// NewPublicKey generates a new public key from the provided secret key.
func (keygen *KeyGenerator) NewPublicKey(sk *SecretKey) (pk *PublicKey) {

	pk = new(PublicKey)

	//pk[0] = [-(a*s + e)]
	//pk[1] = [a]
	pk.pk[0] = keygen.ckksContext.gaussianSampler.SampleNTTNew()
	pk.pk[1] = keygen.context.NewUniformPoly()

	keygen.context.MulCoeffsMontgomeryAndAdd(sk.sk, pk.pk[1], pk.pk[0])
	keygen.context.Neg(pk.pk[0], pk.pk[0])

	keygen.context.MForm(pk.pk[0], pk.pk[0])
	keygen.context.MForm(pk.pk[1], pk.pk[1])

	return pk
}

func (keygen *KeyGenerator) NewPublicKeyEmpty() (pk *PublicKey) {
	pk = new(PublicKey)

	pk.pk[0] = keygen.context.NewPoly()
	pk.pk[1] = keygen.context.NewPoly()

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

func (keygen *KeyGenerator) NewKeyPairSparse(hw uint64) (sk *SecretKey, pk *PublicKey) {
	sk = keygen.NewSecretKeySparse(hw)
	return sk, keygen.NewPublicKey(sk)
}

// NewRelinkey generates a new evaluation key that will be used to relinearize the ciphertexts during multiplication.
// Bitdecomposition aims at reducing the added noise at the expense of more storage needed for the keys and more computation
// during the relinearization. However for relinearization this bitdecomp value can be set to maximum as the encrypted value
// are also scaled up during the multiplication.
func (keygen *KeyGenerator) NewRelinKey(sk *SecretKey) (evakey *EvaluationKey) {

	evakey = new(EvaluationKey)
	keygen.polypool.Copy(sk.Get())
	keygen.context.MulCoeffsMontgomery(keygen.polypool, sk.Get(), keygen.polypool)
	evakey.evakey = keygen.newSwitchingKey(keygen.polypool, sk.Get())
	keygen.polypool.Zero()

	return
}

func (ckksContext *CkksContext) NewRelinKeyEmpty() (evakey *EvaluationKey) {
	evakey = new(EvaluationKey)
	evakey.evakey = new(SwitchingKey)

	// delta_sk = sk_input - sk_output = GaloisEnd(sk_output, rotation) - sk_output
	evakey.evakey.evakey = make([][2]*ring.Poly, ckksContext.beta)
	for i := uint64(0); i < ckksContext.beta; i++ {

		evakey.evakey.evakey[i][0] = ckksContext.contextKeys.NewPoly()
		evakey.evakey.evakey[i][1] = ckksContext.contextKeys.NewPoly()
	}

	return
}

// Get returns the slice of switchintkeys of the evaluation-key.
func (evk *EvaluationKey) Get() *SwitchingKey {
	return evk.evakey
}

func (evk *EvaluationKey) Set(rlk [][2]*ring.Poly) {

	evk.evakey = new(SwitchingKey)
	evk.evakey.evakey = make([][2]*ring.Poly, len(rlk))
	for j := range rlk {
		evk.evakey.evakey[j][0] = rlk[j][0].CopyNew()
		evk.evakey.evakey[j][1] = rlk[j][1].CopyNew()
	}
}

// NewSwitchingKey generated a new keyswitching key, that will re-encrypt a ciphertext encrypted under the input key to the output key.
// Here bitdecomp plays a role in the added noise if the scale of the input is smaller than the maximum size between the modulies.
func (keygen *KeyGenerator) NewSwitchingKey(sk_input, sk_output *SecretKey) (newevakey *SwitchingKey) {
	keygen.context.Sub(sk_input.Get(), sk_output.Get(), keygen.polypool)
	newevakey = keygen.newSwitchingKey(keygen.polypool, sk_output.Get())
	keygen.polypool.Zero()
	return
}

func (keygen *KeyGenerator) NewSwitchingKeyEmpty() (evakey *SwitchingKey) {
	evakey = new(SwitchingKey)

	// delta_sk = sk_input - sk_output = GaloisEnd(sk_output, rotation) - sk_output
	evakey.evakey = make([][2]*ring.Poly, keygen.ckksContext.beta)

	for i := uint64(0); i < keygen.ckksContext.beta; i++ {
		evakey.evakey[i][0] = keygen.context.NewPoly()
		evakey.evakey[i][1] = keygen.context.NewPoly()
	}

	return
}

// NewRotationKeys generates a new instance of rotationkeys, with the provided rotation to the left, right and conjugation if asked.
// Here bitdecomp plays a role in the added noise if the scale of the input is smaller than the maximum size between the modulies.
func (ckksContext *CkksContext) NewRotationKeys() (rotKey *RotationKeys) {
	rotKey = new(RotationKeys)
	rotKey.ckksContext = ckksContext
	return
}

// Newrotationkeys generates a new struct of rotationkeys storing the keys for the specified rotations. The provided secret-key must be the secret-key used to generate the public-key under
// which the ciphertexts to rotate are encrypted under. Bitdecomp is the power of two binary decomposition of the key. A higher bigdecomp will induce smaller keys, faster key-switching,
// but at the cost of more noise. rotLeft and rotRight must be a slice of uint64 rotations, row is a boolean value indicating if the key for the row rotation must be generated.
func (keygen *KeyGenerator) GenRot(rotType Rotation, sk *SecretKey, k uint64, rotKey *RotationKeys) {
	switch rotType {
	case RotationLeft:
		if rotKey.evakey_rot_col_L == nil {
			rotKey.evakey_rot_col_L = make(map[uint64]*SwitchingKey)
		}
		if rotKey.evakey_rot_col_L[k] == nil && k != 0 {
			rotKey.evakey_rot_col_L[k] = keygen.genrotKey(sk.Get(), keygen.ckksContext.galElRotColLeft[k])
		}
	case RotationRight:
		if rotKey.evakey_rot_col_R == nil {
			rotKey.evakey_rot_col_R = make(map[uint64]*SwitchingKey)
		}
		if rotKey.evakey_rot_col_R[k] == nil && k != 0 {
			rotKey.evakey_rot_col_R[k] = keygen.genrotKey(sk.Get(), keygen.ckksContext.galElRotColRight[k])
		}
	case Conjugate:
		rotKey.evakey_rot_row = keygen.genrotKey(sk.Get(), keygen.ckksContext.galElRotRow)
	}
}

// NewRotationkeysPow2 generates a new rotation key with all the power of two rotation to the left and right, as well as the conjugation
// key if asked. Here bitdecomp plays a role in the added noise if the scale of the input is smaller than the maximum size between the modulies.
func (keygen *KeyGenerator) NewRotationKeysPow2(sk_output *SecretKey) (rotKey *RotationKeys) {

	rotKey = new(RotationKeys)
	rotKey.ckksContext = keygen.ckksContext

	rotKey.evakey_rot_col_L = make(map[uint64]*SwitchingKey)
	rotKey.evakey_rot_col_R = make(map[uint64]*SwitchingKey)

	for n := uint64(1); n < rotKey.ckksContext.n>>1; n <<= 1 {

		rotKey.evakey_rot_col_L[n] = keygen.genrotKey(sk_output.Get(), keygen.ckksContext.galElRotColLeft[n])
		rotKey.evakey_rot_col_R[n] = keygen.genrotKey(sk_output.Get(), keygen.ckksContext.galElRotColRight[n])
	}

	rotKey.evakey_rot_row = keygen.genrotKey(sk_output.Get(), keygen.ckksContext.galElRotRow)
	return
}

func (rotKey *RotationKeys) SetRotKey(rotType Rotation, k uint64, evakey [][2]*ring.Poly) {
	switch rotType {
	case RotationLeft:
		if rotKey.evakey_rot_col_L == nil {
			rotKey.evakey_rot_col_L = make(map[uint64]*SwitchingKey)
		}
		if rotKey.evakey_rot_col_L[k] == nil && k != 0 {
			rotKey.evakey_rot_col_L[k] = new(SwitchingKey)
			rotKey.evakey_rot_col_L[k].evakey = make([][2]*ring.Poly, len(evakey))
			for j := range evakey {
				rotKey.evakey_rot_col_L[k].evakey[j][0] = evakey[j][0].CopyNew()
				rotKey.evakey_rot_col_L[k].evakey[j][1] = evakey[j][1].CopyNew()
			}
		}
	case RotationRight:
		if rotKey.evakey_rot_col_R == nil {
			rotKey.evakey_rot_col_R = make(map[uint64]*SwitchingKey)
		}
		if rotKey.evakey_rot_col_R[k] == nil && k != 0 {
			rotKey.evakey_rot_col_R[k] = new(SwitchingKey)
			rotKey.evakey_rot_col_R[k].evakey = make([][2]*ring.Poly, len(evakey))
			for j := range evakey {
				rotKey.evakey_rot_col_R[k].evakey[j][0] = evakey[j][0].CopyNew()
				rotKey.evakey_rot_col_R[k].evakey[j][1] = evakey[j][1].CopyNew()
			}
		}
	case Conjugate:
		if rotKey.evakey_rot_row == nil {
			rotKey.evakey_rot_row = new(SwitchingKey)
			rotKey.evakey_rot_row.evakey = make([][2]*ring.Poly, len(evakey))
			for j := range evakey {
				rotKey.evakey_rot_row.evakey[j][0] = evakey[j][0].CopyNew()
				rotKey.evakey_rot_row.evakey[j][1] = evakey[j][1].CopyNew()
			}
		}
	}
}

func (keygen *KeyGenerator) genrotKey(sk_output *ring.Poly, gen uint64) (switchingkey *SwitchingKey) {

	ring.PermuteNTT(sk_output, gen, keygen.polypool)
	keygen.context.Sub(keygen.polypool, sk_output, keygen.polypool)
	switchingkey = keygen.newSwitchingKey(keygen.polypool, sk_output)
	keygen.polypool.Zero()

	return
}

func (keygen *KeyGenerator) newSwitchingKey(sk_in, sk_out *ring.Poly) (switchingkey *SwitchingKey) {

	switchingkey = new(SwitchingKey)

	context := keygen.ckksContext.contextKeys

	// Computes P * sk_in
	for _, pj := range keygen.ckksContext.specialprimes {
		context.MulScalar(sk_in, pj, sk_in)
	}

	alpha := keygen.ckksContext.alpha
	beta := keygen.ckksContext.beta

	var index uint64

	switchingkey.evakey = make([][2]*ring.Poly, beta)

	for i := uint64(0); i < beta; i++ {

		// e
		switchingkey.evakey[i][0] = keygen.ckksContext.gaussianSampler.SampleNTTNew()
		context.MForm(switchingkey.evakey[i][0], switchingkey.evakey[i][0])

		// a (since a is uniform, we consider we already sample it in the NTT and montgomery domain)
		switchingkey.evakey[i][1] = keygen.context.NewUniformPoly()

		// e + (sk_in * P) * (q_star * q_tild) mod QP
		//
		// q_prod = prod(q[i*alpha+j])
		// q_star = Q/qprod
		// q_tild = q_star^-1 mod q_prod
		//
		// Therefore : (sk_in * P) * (q_star * q_tild) = sk*P mod q[i*alpha+j], else 0
		for j := uint64(0); j < alpha; j++ {

			index = i*alpha + j

			qi := context.Modulus[index]
			p0tmp := sk_in.Coeffs[index]
			p1tmp := switchingkey.evakey[i][0].Coeffs[index]

			for w := uint64(0); w < context.N; w++ {
				p1tmp[w] = ring.CRed(p1tmp[w]+p0tmp[w], qi)
			}

			// Handles the case where nb pj does not divides nb qi
			if index >= keygen.ckksContext.levels-1 {
				break
			}
		}

		// (sk_in * P) * (q_star * q_tild) - a * sk_out + e mod QP
		context.MulCoeffsMontgomeryAndSub(switchingkey.evakey[i][1], sk_out, switchingkey.evakey[i][0])
	}

	return
}
