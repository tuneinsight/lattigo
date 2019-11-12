package bfv

import (
	"github.com/ldsec/lattigo/ring"
)

// KeyGenerator is a structure that stores the elements required to create new keys,
// as well as a small memory pool for intermediate values.
type KeyGenerator struct {
	bfvcontext *BfvContext
	context    *ring.Context
	polypool   *ring.Poly
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
	RotationRow
)

// Rotationkeys is a structure that stores the switching-keys required during the homomorphic rotations.
type RotationKeys struct {
	bfvcontext       *BfvContext
	evakey_rot_col_L map[uint64]*SwitchingKey
	evakey_rot_col_R map[uint64]*SwitchingKey
	evakey_rot_row   *SwitchingKey
}

// Evaluationkey is a structure that stores the switching-keys required during the relinearization.
type EvaluationKey struct {
	evakey []*SwitchingKey
}

// Switchingkey is a structure that stores the switching-keys required during the key-switching.
type SwitchingKey struct {
	evakey [][2]*ring.Poly
}

// Get returns the switching key backing slice
func (swk *SwitchingKey) Get() [][2]*ring.Poly {
	return swk.evakey
}

// NewKeyGenerator creates a new KeyGenerator, from which the secret and public keys, as well as the evaluation,
// rotation and switching keys can be generated.
func (bfvcontext *BfvContext) NewKeyGenerator() (keygen *KeyGenerator) {
	keygen = new(KeyGenerator)
	keygen.bfvcontext = bfvcontext
	keygen.context = bfvcontext.contextKeys
	keygen.polypool = keygen.context.NewPoly()
	return
}

// Newsecretkey creates a new SecretKey with the distribution [1/3, 1/3, 1/3]
func (keygen *KeyGenerator) NewSecretKey() (sk *SecretKey) {
	sk, _ = keygen.NewSecretkeyWithDistrib(1.0 / 3)
	return sk
}

// Newsecretkey creates a new SecretKey with the distribution [(p-1)/2, p, (p-1)/2]
func (keygen *KeyGenerator) NewSecretkeyWithDistrib(p float64) (sk *SecretKey, err error) {

	sk = new(SecretKey)
	if sk.sk, err = keygen.bfvcontext.ternarySampler.SampleMontgomeryNTTNew(p); err != nil {
		return nil, err
	}

	return sk, nil
}

// NewSecretKeyEmpty creates a new SecretKey with all coeffcients set to zero, ready to received a marshaled SecretKey.
func (keygen *KeyGenerator) NewSecretKeyEmpty() *SecretKey {
	sk := new(SecretKey)
	sk.sk = keygen.context.NewPoly()
	return sk
}

// Get returns the polynomial of the target secret-key.
func (sk *SecretKey) Get() *ring.Poly {
	return sk.sk
}

// Set sets the polynomial of the target secret key as the input polynomial.
func (sk *SecretKey) Set(poly *ring.Poly) {
	sk.sk = poly.CopyNew()
}

// Newpublickey generates a new publickkey from the provided secret-key
func (keygen *KeyGenerator) NewPublicKey(sk *SecretKey) (pk *PublicKey) {

	pk = new(PublicKey)

	//pk[0] = [-(a*s + e)]
	//pk[1] = [a]
	pk.pk[0] = keygen.bfvcontext.gaussianSampler.SampleNTTNew()
	pk.pk[1] = keygen.context.NewUniformPoly()

	keygen.context.MulCoeffsMontgomeryAndAdd(sk.sk, pk.pk[1], pk.pk[0])
	keygen.context.Neg(pk.pk[0], pk.pk[0])

	return pk
}

func (bfvContext *BfvContext) NewPublicKey() (pk *PublicKey) {
	pk = new(PublicKey)
	pk.pk[0] = bfvContext.contextQ.NewPoly()
	pk.pk[1] = bfvContext.contextQ.NewPoly()
	return
}

// Get returns the polynomials of the public-key.
func (pk *PublicKey) Get() [2]*ring.Poly {
	return pk.pk
}

// Set sets the polynomial of the public-key as the input polynomials.
func (pk *PublicKey) Set(p [2]*ring.Poly) {
	pk.pk[0] = p[0].CopyNew()
	pk.pk[1] = p[1].CopyNew()
}

// NewKeyPair generates a new secret-key with distribution [1/3, 1/3, 1/3] and a corresponding public-key.
func (keygen *KeyGenerator) NewKeyPair() (sk *SecretKey, pk *PublicKey) {
	sk = keygen.NewSecretKey()
	return sk, keygen.NewPublicKey(sk)
}

// NewRelinKey generates a new evaluation key from the provided secret-key. It will be used to relinearize a ciphertext (encrypted under a public-key generated from the provided secret-key)
// of degree > 1 to a ciphertext of degree 1. Max degree is the maximum degree of the ciphertext allowed to relinearize. Bitdecomp is the power of two binary decomposition of the key.
// A higher bigdecomp will induce smaller keys, faster key-switching, but at the cost of more noise.
func (keygen *KeyGenerator) NewRelinKey(sk *SecretKey, maxDegree uint64) (newEvakey *EvaluationKey) {

	newEvakey = new(EvaluationKey)

	newEvakey.evakey = make([]*SwitchingKey, maxDegree)

	keygen.polypool.Copy(sk.Get())

	for _, pj := range keygen.bfvcontext.specialprimes {
		keygen.bfvcontext.contextKeys.MulScalar(keygen.polypool, pj, keygen.polypool)
	}

	for i := uint64(0); i < maxDegree; i++ {
		keygen.context.MulCoeffsMontgomery(keygen.polypool, sk.Get(), keygen.polypool)
		newEvakey.evakey[i] = newswitchintkey(keygen.bfvcontext, keygen.polypool, sk.Get())
	}

	keygen.polypool.Zero()

	return newEvakey
}

func (bfvcontext *BfvContext) NewRelinKeyEmpty(maxDegree uint64) (evakey *EvaluationKey) {

	evakey = new(EvaluationKey)

	context := bfvcontext.contextKeys

	// delta_sk = sk_input - sk_output = GaloisEnd(sk_output, rotation) - sk_output

	evakey.evakey = make([]*SwitchingKey, maxDegree)

	for w := uint64(0); w < maxDegree; w++ {

		evakey.evakey[w] = new(SwitchingKey)
		evakey.evakey[w].evakey = make([][2]*ring.Poly, len(context.Modulus))

		for i := range context.Modulus {
			evakey.evakey[w].evakey[i][0] = context.NewPoly()
			evakey.evakey[w].evakey[i][1] = context.NewPoly()

		}
	}

	return
}

// Get returns the slice of switchintkeys of the evaluation-key.
func (evk *EvaluationKey) Get() []*SwitchingKey {
	return evk.evakey
}

// SetRelinKeys sets the polynomial of the target evaluation-key as the input polynomials.
func (newevakey *EvaluationKey) SetRelinKeys(rlk [][][2]*ring.Poly) {

	newevakey.evakey = make([]*SwitchingKey, len(rlk))
	for i := range rlk {
		newevakey.evakey[i] = new(SwitchingKey)
		newevakey.evakey[i].evakey = make([][2]*ring.Poly, len(rlk[i]))
		for j := range rlk[i] {
			newevakey.evakey[i].evakey[j][0] = rlk[i][j][0].CopyNew()
			newevakey.evakey[i].evakey[j][1] = rlk[i][j][1].CopyNew()
		}
	}
}

// Newswitchintkey generates a new key-switching key, that will allow to re-encrypt under the output-key a ciphertext encrypted under the input-key. Bitdecomp
// is the power of two binary decomposition of the key. A higher bigdecomp will induce smaller keys, faster key-switching, but at the cost of more noise.
func (keygen *KeyGenerator) NewSwitchingKey(sk_input, sk_output *SecretKey) (newevakey *SwitchingKey) {

	keygen.context.Sub(sk_input.Get(), sk_output.Get(), keygen.polypool)

	for _, pj := range keygen.bfvcontext.specialprimes {
		keygen.bfvcontext.contextKeys.MulScalar(keygen.polypool, pj, keygen.polypool)
	}

	newevakey = newswitchintkey(keygen.bfvcontext, keygen.polypool, sk_output.Get())
	keygen.polypool.Zero()

	return
}

func (keygen *KeyGenerator) NewSwitchingKeyEmpty() (evakey *SwitchingKey) {
	evakey = new(SwitchingKey)

	context := keygen.context

	// delta_sk = sk_input - sk_output = GaloisEnd(sk_output, rotation) - sk_output
	evakey.evakey = make([][2]*ring.Poly, len(context.Modulus))

	for i := range context.Modulus {
		evakey.evakey[i][0] = context.NewPoly()
		evakey.evakey[i][1] = context.NewPoly()
	}

	return
}

func (bfvContext *BfvContext) NewRotationKeys() (rotKey *RotationKeys) {
	rotKey = new(RotationKeys)
	rotKey.bfvcontext = bfvContext
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
			rotKey.evakey_rot_col_L[k] = genrotkey(keygen, sk.Get(), keygen.bfvcontext.galElRotColLeft[k])
		}
	case RotationRight:
		if rotKey.evakey_rot_col_R == nil {
			rotKey.evakey_rot_col_R = make(map[uint64]*SwitchingKey)
		}
		if rotKey.evakey_rot_col_R[k] == nil && k != 0 {
			rotKey.evakey_rot_col_R[k] = genrotkey(keygen, sk.Get(), keygen.bfvcontext.galElRotColRight[k])
		}
	case RotationRow:
		rotKey.evakey_rot_row = genrotkey(keygen, sk.Get(), keygen.bfvcontext.galElRotRow)
	}
}

// Newrotationkeys generates a new struct of rotationkeys storing the keys of all the left and right powers of two rotations. The provided secret-key must be the secret-key used to generate the public-key under
// which the ciphertexts to rotate are encrypted under. rows is a boolean value indicatig if the keys for the row rotation have to be generated. Bitdecomp is the power of two binary decomposition of the key.
// A higher bigdecomp will induce smaller keys, faster key-switching, but at the cost of more noise.
func (keygen *KeyGenerator) NewRotationKeysPow2(sk *SecretKey) (rotKey *RotationKeys) {

	rotKey = new(RotationKeys)
	rotKey.bfvcontext = keygen.bfvcontext

	rotKey.evakey_rot_col_L = make(map[uint64]*SwitchingKey)
	rotKey.evakey_rot_col_R = make(map[uint64]*SwitchingKey)

	for n := uint64(1); n < rotKey.bfvcontext.n>>1; n <<= 1 {

		rotKey.evakey_rot_col_L[n] = genrotkey(keygen, sk.Get(), keygen.bfvcontext.galElRotColLeft[n])
		rotKey.evakey_rot_col_R[n] = genrotkey(keygen, sk.Get(), keygen.bfvcontext.galElRotColRight[n])
	}

	rotKey.evakey_rot_row = genrotkey(keygen, sk.Get(), keygen.bfvcontext.galElRotRow)

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
	case RotationRow:
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

// genrotkey is a methode used in the rotation-keys generation.
func genrotkey(keygen *KeyGenerator, sk *ring.Poly, gen uint64) (switchkey *SwitchingKey) {

	ring.PermuteNTT(sk, gen, keygen.polypool)
	keygen.context.Sub(keygen.polypool, sk, keygen.polypool)

	for _, pj := range keygen.bfvcontext.specialprimes {
		keygen.bfvcontext.contextKeys.MulScalar(keygen.polypool, pj, keygen.polypool)
	}

	switchkey = newswitchintkey(keygen.bfvcontext, keygen.polypool, sk)
	keygen.polypool.Zero()

	return
}

// newswitchintkey is a generic methode to generate key-switching keys used in the evaluation, key-switching and rotation-keys generation.
func newswitchintkey(bfvcontext *BfvContext, sk_in, sk_out *ring.Poly) (switchkey *SwitchingKey) {

	switchkey = new(SwitchingKey)

	context := bfvcontext.contextKeys

	var index uint64

	// delta_sk = sk_input - sk_output = GaloisEnd(sk_output, rotation) - sk_output

	switchkey.evakey = make([][2]*ring.Poly, bfvcontext.beta)

	for i := uint64(0); i < bfvcontext.beta; i++ {

		// e
		switchkey.evakey[i][0] = bfvcontext.gaussianSampler.SampleNTTNew()
		context.MForm(switchkey.evakey[i][0], switchkey.evakey[i][0])
		// a
		switchkey.evakey[i][1] = context.NewUniformPoly()

		// e + sk_in * (qiBarre*qiStar) * 2^w
		// (qiBarre*qiStar)%qi = 1, else 0

		for j := uint64(0); j < bfvcontext.alpha; j++ {

			index = i*bfvcontext.alpha + j

			for w := uint64(0); w < context.N; w++ {
				switchkey.evakey[i][0].Coeffs[index][w] = ring.CRed(switchkey.evakey[i][0].Coeffs[index][w]+sk_in.Coeffs[index][w], context.Modulus[index])
			}

			// Handles the case where nb pj does not divides nb qi
			if index >= uint64(len(context.Modulus)-1) {
				break
			}

		}

		// sk_in * (qiBarre*qiStar) * 2^w - a*sk + e
		context.MulCoeffsMontgomeryAndSub(switchkey.evakey[i][1], sk_out, switchkey.evakey[i][0])
	}

	return
}
