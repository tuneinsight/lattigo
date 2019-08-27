package bfv

import (
	"errors"
	"github.com/lca1/lattigo/ring"
	"math"
	"math/bits"
)

// SecretKey is a structure that stores the secret-key
type SecretKey struct {
	sk *ring.Poly
}

// PublicKey is a structure that stores the public-key
type PublicKey struct {
	pk [2]*ring.Poly
}

// KeyGenerator is a structure that stores the elements required to create new keys,
// as well as a small memory pool for intermediate values.
type KeyGenerator struct {
	bfvcontext *BfvContext
	context    *ring.Context
	polypool   *ring.Poly
}

// RotationKeys is a structure that stores the switching-keys required during the homomorphic rotations.
type RotationKeys struct {
	bfvcontext       *BfvContext
	bitDecomp        uint64
	evakey_rot_col_L map[uint64]*SwitchingKey
	evakey_rot_col_R map[uint64]*SwitchingKey
	evakey_rot_row   *SwitchingKey
}

// EvaluationKey is a structure that stores the switching-keys required during the relinearization.
type EvaluationKey struct {
	evakey []*SwitchingKey
}

// SwitchingKey is a structure that stores the switching-keys required during the key-switching.
type SwitchingKey struct {
	bitDecomp uint64
	evakey    [][][2]*ring.Poly
}

// NewKeyGenerator creates a new KeyGenerator from the target bfvcontext.
func (bfvcontext *BfvContext) NewKeyGenerator() (keygen *KeyGenerator) {
	keygen = new(KeyGenerator)
	keygen.bfvcontext = bfvcontext
	keygen.context = bfvcontext.contextQ
	keygen.polypool = keygen.context.NewPoly()
	return
}

// NewSecretKey creates a new secretkey with uniform distribution in [-1, 0, 1].
func (keygen *KeyGenerator) NewSecretKey() *SecretKey {

	sk := new(SecretKey)
	sk.sk = keygen.bfvcontext.ternarySampler.SampleMontgomeryNTTNew()

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

// check_sk checks if the input secret-key complies with the keygenerator context.
func (keygen *KeyGenerator) check_sk(sk_output *SecretKey) error {

	if sk_output.Get().GetDegree() != int(keygen.context.N) {
		return errors.New("error : pol degree sk != bfvcontext.n")
	}

	if len(sk_output.Get().Coeffs) != len(keygen.context.Modulus) {
		return errors.New("error : nb modulus sk != nb modulus bfvcontext")
	}

	return nil
}

// NewPublicKey generates a new publickkey from the provided secret-key
func (keygen *KeyGenerator) NewPublicKey(sk *SecretKey) (pk *PublicKey, err error) {

	if err = keygen.check_sk(sk); err != nil {
		return nil, err
	}

	pk = new(PublicKey)

	//pk[0] = [-(a*s + e)]
	//pk[1] = [a]
	pk.pk[0] = keygen.bfvcontext.gaussianSampler.SampleNTTNew()
	pk.pk[1] = keygen.context.NewUniformPoly()

	keygen.context.MulCoeffsMontgomeryAndAdd(sk.sk, pk.pk[1], pk.pk[0])
	keygen.context.Neg(pk.pk[0], pk.pk[0])

	return pk, nil
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

// NewKeyPair generates a new (secret-key, public-key) pair.
func (keygen *KeyGenerator) NewKeyPair() (sk *SecretKey, pk *PublicKey, err error) {
	sk = keygen.NewSecretKey()
	pk, err = keygen.NewPublicKey(sk)
	return
}

// NewRelinKey generates a new evaluation key from the provided secret-key. It will be used to relinearize a ciphertext (encrypted under a public-key generated from the provided secret-key)
// of degree > 1 to a ciphertext of degree 1. Max degree is the maximum degree of the ciphertext allowed to relinearize. Bitdecomp is the power of two binary decomposition of the key.
// A higher bigdecomp will induce smaller keys, faster key-switching, but at the cost of more noise.
func (keygen *KeyGenerator) NewRelinKey(sk *SecretKey, maxDegree, bitDecomp uint64) (newEvakey *EvaluationKey, err error) {

	newEvakey = new(EvaluationKey)
	newEvakey.evakey = make([]*SwitchingKey, maxDegree)
	sk.Get().Copy(keygen.polypool)
	for i := uint64(0); i < maxDegree; i++ {
		keygen.context.MulCoeffsMontgomery(keygen.polypool, sk.Get(), keygen.polypool)
		newEvakey.evakey[i] = newswitchingkey(keygen.bfvcontext, keygen.polypool, sk.Get(), bitDecomp)
	}
	keygen.polypool.Zero()

	return newEvakey, nil
}

// Get returns the slice of switchingkeys of the evaluation-key.
func (evk *EvaluationKey) Get() []*SwitchingKey {
	return evk.evakey
}

// SetRelinKeys sets the polynomial of the target evaluation-key as the input polynomials.
func (newevakey *EvaluationKey) SetRelinKeys(rlk [][][][2]*ring.Poly, bitDecomp uint64) {

	newevakey.evakey = make([]*SwitchingKey, len(rlk))
	for i := range rlk {
		newevakey.evakey[i] = new(SwitchingKey)
		newevakey.evakey[i].bitDecomp = bitDecomp
		newevakey.evakey[i].evakey = make([][][2]*ring.Poly, len(rlk[i]))
		for j := range rlk[i] {
			newevakey.evakey[i].evakey[j] = make([][2]*ring.Poly, len(rlk[i][j]))
			for u := range rlk[i][j] {
				newevakey.evakey[i].evakey[j][u][0] = rlk[i][j][u][0].CopyNew()
				newevakey.evakey[i].evakey[j][u][1] = rlk[i][j][u][1].CopyNew()
			}
		}
	}
}

// NewSwitchingKey generates a new key-switching key, that will allow to re-encrypt under the output-key a ciphertext encrypted under the input-key. Bitdecomp
// is the power of two binary decomposition of the key. A higher bigdecomp will induce smaller keys, faster key-switching, but at the cost of more noise.
func (keygen *KeyGenerator) NewSwitchingKey(sk_input, sk_output *SecretKey, bitDecomp uint64) (newevakey *SwitchingKey, err error) {

	if err = keygen.check_sk(sk_input); err != nil {
		return nil, err
	}

	if err = keygen.check_sk(sk_output); err != nil {
		return nil, err
	}

	keygen.context.Sub(sk_input.Get(), sk_output.Get(), keygen.polypool)
	newevakey = newswitchingkey(keygen.bfvcontext, keygen.polypool, sk_output.Get(), bitDecomp)
	keygen.polypool.Zero()

	return
}

// NewRotationKeys generates a new struct of rotationkeys storing the keys for the specified rotations. The provided secret-key must be the secret-key used to generate the public-key under
// which the ciphertexts to rotate are encrypted under. Bitdecomp is the power of two binary decomposition of the key. A higher bigdecomp will induce smaller keys, faster key-switching,
// but at the cost of more noise. rotLeft and rotRight must be a slice of uint64 rotations, row is a boolean value indicating if the key for the row rotation must be generated.
func (keygen *KeyGenerator) NewRotationKeys(sk *SecretKey, bitDecomp uint64, rotLeft []uint64, rotRight []uint64, row bool) (rotKey *RotationKeys, err error) {

	if err = keygen.check_sk(sk); err != nil {
		return nil, err
	}

	rotKey = new(RotationKeys)
	rotKey.bfvcontext = keygen.bfvcontext
	rotKey.bitDecomp = bitDecomp

	if rotLeft != nil {
		rotKey.evakey_rot_col_L = make(map[uint64]*SwitchingKey)
		for _, n := range rotLeft {
			if rotKey.evakey_rot_col_L[n] == nil && n != 0 {
				rotKey.evakey_rot_col_L[n] = genrotkey(keygen, sk.Get(), keygen.bfvcontext.galElRotColLeft[n], bitDecomp)
			}
		}
	}

	if rotRight != nil {
		rotKey.evakey_rot_col_R = make(map[uint64]*SwitchingKey)
		for _, n := range rotRight {
			if rotKey.evakey_rot_col_R[n] == nil && n != 0 {
				rotKey.evakey_rot_col_R[n] = genrotkey(keygen, sk.Get(), keygen.bfvcontext.galElRotColRight[n], bitDecomp)
			}
		}
	}

	if row {
		rotKey.evakey_rot_row = genrotkey(keygen, sk.Get(), keygen.bfvcontext.galElRotRow, bitDecomp)
	}

	return rotKey, nil

}

// NewRotationKeys generates a new struct of rotationkeys storing the keys of all the left and right powers of two rotations. The provided secret-key must be the secret-key used to generate the public-key under
// which the ciphertexts to rotate are encrypted under. rows is a boolean value indicatig if the keys for the row rotation have to be generated. Bitdecomp is the power of two binary decomposition of the key.
// A higher bigdecomp will induce smaller keys, faster key-switching, but at the cost of more noise.
func (keygen *KeyGenerator) NewRotationKeysPow2(sk *SecretKey, bitDecomp uint64, row bool) (rotKey *RotationKeys, err error) {

	if err = keygen.check_sk(sk); err != nil {
		return nil, err
	}

	rotKey = new(RotationKeys)
	rotKey.bfvcontext = keygen.bfvcontext
	rotKey.bitDecomp = bitDecomp

	rotKey.evakey_rot_col_L = make(map[uint64]*SwitchingKey)
	rotKey.evakey_rot_col_R = make(map[uint64]*SwitchingKey)

	for n := uint64(1); n < rotKey.bfvcontext.n>>1; n <<= 1 {

		rotKey.evakey_rot_col_L[n] = genrotkey(keygen, sk.Get(), keygen.bfvcontext.galElRotColLeft[n], bitDecomp)
		rotKey.evakey_rot_col_R[n] = genrotkey(keygen, sk.Get(), keygen.bfvcontext.galElRotColRight[n], bitDecomp)
	}

	if row {
		rotKey.evakey_rot_row = genrotkey(keygen, sk.Get(), keygen.bfvcontext.galElRotRow, bitDecomp)
	}

	return
}

// genrotkey is a methode used in the rotation-keys generation.
func genrotkey(keygen *KeyGenerator, sk *ring.Poly, gen, bitDecomp uint64) (switchingkey *SwitchingKey) {

	ring.PermuteNTT(sk, gen, keygen.polypool)
	keygen.context.Sub(keygen.polypool, sk, keygen.polypool)
	switchingkey = newswitchingkey(keygen.bfvcontext, keygen.polypool, sk, bitDecomp)
	keygen.polypool.Zero()

	return
}

// newswitchingkey is a generic methode to generate key-switching keys used in the evaluation, key-switching and rotation-keys generation.
func newswitchingkey(bfvcontext *BfvContext, sk_in, sk_out *ring.Poly, bitDecomp uint64) (switchingkey *SwitchingKey) {

	if bitDecomp > bfvcontext.maxBit || bitDecomp == 0 {
		bitDecomp = bfvcontext.maxBit
	}

	switchingkey = new(SwitchingKey)

	context := bfvcontext.contextQ

	switchingkey.bitDecomp = uint64(bitDecomp)

	mredParams := context.GetMredParams()

	// delta_sk = sk_input - sk_output = GaloisEnd(sk_output, rotation) - sk_output

	var bitLog uint64

	switchingkey.evakey = make([][][2]*ring.Poly, len(context.Modulus))

	for i, qi := range context.Modulus {

		bitLog = uint64(math.Ceil(float64(bits.Len64(qi)) / float64(bitDecomp)))

		switchingkey.evakey[i] = make([][2]*ring.Poly, bitLog)

		for j := uint64(0); j < bitLog; j++ {

			// e
			switchingkey.evakey[i][j][0] = bfvcontext.gaussianSampler.SampleNTTNew()
			// a
			switchingkey.evakey[i][j][1] = context.NewUniformPoly()

			// e + sk_in * (qiBarre*qiStar) * 2^w
			// (qiBarre*qiStar)%qi = 1, else 0
			for w := uint64(0); w < context.N; w++ {
				switchingkey.evakey[i][j][0].Coeffs[i][w] += ring.PowerOf2(sk_in.Coeffs[i][w], bitDecomp*j, qi, mredParams[i])
			}

			// sk_in * (qiBarre*qiStar) * 2^w - a*sk + e
			context.MulCoeffsMontgomeryAndSub(switchingkey.evakey[i][j][1], sk_out, switchingkey.evakey[i][j][0])

			context.MForm(switchingkey.evakey[i][j][0], switchingkey.evakey[i][j][0])
			context.MForm(switchingkey.evakey[i][j][1], switchingkey.evakey[i][j][1])
		}
	}

	return
}
