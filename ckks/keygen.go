package ckks

import (
	"errors"
	"github.com/lca1/lattigo/ring"
	"math"
	"math/bits"
)

// Keygenerator is a structure that stores the elements required to create new keys,
// as well as a small memory pool for intermediate values.
type keygenerator struct {
	ckkscontext *CkksContext
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

// Rotationkeys is a structure that stores the switching-keys required during the homomorphic rotations.
type RotationKey struct {
	ckkscontext      *CkksContext
	bitDecomp        uint64
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
	bitDecomp uint64
	evakey    [][][2]*ring.Poly
}

// NewKeyGenerator creates a new keygenerator, from which the secret and public keys, as well as the evaluation,
// rotation and switching keys can be generated.
func (ckkscontext *CkksContext) NewKeyGenerator() (keygen *keygenerator) {
	keygen = new(keygenerator)
	keygen.ckkscontext = ckkscontext
	keygen.context = ckkscontext.keyscontext
	keygen.polypool = ckkscontext.keyscontext.NewPoly()
	return
}

// check_sk checks if the input secret-key complies with the keygenerator context.
func (keygen *keygenerator) check_sk(sk_output *SecretKey) error {

	if sk_output.Get().GetDegree() != int(keygen.context.N) {
		return errors.New("error : pol degree sk != ckkscontext.n")
	}

	if len(sk_output.Get().Coeffs) != len(keygen.context.Modulus) {
		return errors.New("error : nb modulus sk != nb modulus ckkscontext")
	}

	return nil
}

// NewSecretKey generates a new secret key.
func (keygen *keygenerator) NewSecretKey() *SecretKey {
	sk := new(SecretKey)
	sk.sk = keygen.ckkscontext.ternarySampler.SampleMontgomeryNTTNew()
	return sk
}

func (keygen *keygenerator) NewSecretKeyEmpty() *SecretKey {
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
func (keygen *keygenerator) NewPublicKey(sk *SecretKey) (pk *PublicKey, err error) {

	if err = keygen.check_sk(sk); err != nil {
		return nil, err
	}

	pk = new(PublicKey)

	//pk[0] = [-(a*s + e)]
	//pk[1] = [a]
	pk.pk[0] = keygen.ckkscontext.gaussianSampler.SampleNTTNew()
	pk.pk[1] = keygen.context.NewUniformPoly()

	keygen.context.MulCoeffsMontgomeryAndAdd(sk.sk, pk.pk[1], pk.pk[0])
	keygen.context.Neg(pk.pk[0], pk.pk[0])

	return pk, nil
}

func (keygen *keygenerator) NewPublicKeyEmpty() (pk *PublicKey) {
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

// NewKeyPair generates a new secretkey and a corresponding public key.
func (keygen *keygenerator) NewKeyPair() (sk *SecretKey, pk *PublicKey, err error) {
	sk = keygen.NewSecretKey()
	pk, err = keygen.NewPublicKey(sk)
	return
}

// NewRelinkey generates a new evaluation key that will be used to relinearize the ciphertexts during multiplication.
// Bitdecomposition aims at reducing the added noise at the expense of more storage needed for the keys and more computation
// during the relinearization. However for relinearization this bitdecomp value can be set to maximum as the encrypted value
// are also scaled up during the multiplication.
func (keygen *keygenerator) NewRelinKey(sk *SecretKey, bitDecomp uint64) (evakey *EvaluationKey, err error) {

	if err = keygen.check_sk(sk); err != nil {
		return nil, err
	}
	evakey = new(EvaluationKey)
	sk.Get().Copy(keygen.polypool)
	keygen.context.MulCoeffsMontgomery(keygen.polypool, sk.Get(), keygen.polypool)
	evakey.evakey = newswitchingkey(keygen.ckkscontext, keygen.polypool, sk.Get(), bitDecomp)
	keygen.polypool.Zero()

	return
}

func (keygen *keygenerator) NewRelinKeyEmpty(bitDecomp uint64) (evakey *EvaluationKey) {
	evakey = new(EvaluationKey)
	evakey.evakey = new(SwitchingKey)

	if bitDecomp > keygen.ckkscontext.maxBit || bitDecomp == 0 {
		bitDecomp = keygen.ckkscontext.maxBit
	}

	context := keygen.ckkscontext.keyscontext

	evakey.evakey.bitDecomp = bitDecomp

	// delta_sk = sk_input - sk_output = GaloisEnd(sk_output, rotation) - sk_output
	var bitLog uint64

	evakey.evakey.evakey = make([][][2]*ring.Poly, len(context.Modulus))

	for i, qi := range context.Modulus {

		bitLog = uint64(math.Ceil(float64(bits.Len64(qi)) / float64(bitDecomp)))

		evakey.evakey.evakey[i] = make([][2]*ring.Poly, bitLog)

		for j := uint64(0); j < bitLog; j++ {
			evakey.evakey.evakey[i][j][0] = context.NewPoly()
			evakey.evakey.evakey[i][j][1] = context.NewPoly()
		}
	}

	return
}

func (keygen *keygenerator) SetRelinKeys(rlk [][][2]*ring.Poly, bitDecomp uint64) (*EvaluationKey, error) {

	newevakey := new(EvaluationKey)

	newevakey.evakey = new(SwitchingKey)
	newevakey.evakey.bitDecomp = bitDecomp
	newevakey.evakey.evakey = make([][][2]*ring.Poly, len(rlk))
	for j := range rlk {
		newevakey.evakey.evakey[j] = make([][2]*ring.Poly, len(rlk[j]))
		for u := range rlk[j] {
			newevakey.evakey.evakey[j][u][0] = rlk[j][u][0].CopyNew()
			newevakey.evakey.evakey[j][u][1] = rlk[j][u][1].CopyNew()
		}
	}

	return newevakey, nil
}

// NewSwitchingKey generated a new keyswitching key, that will re-encrypt a ciphertext encrypted under the input key to the output key.
// Here bitdecomp plays a role in the added noise if the scale of the input is smaller than the maximum size between the modulies.
func (keygen *keygenerator) NewSwitchingKey(sk_input, sk_output *SecretKey, bitDecomp uint64) (newevakey *SwitchingKey, err error) {

	if err = keygen.check_sk(sk_input); err != nil {
		return nil, err
	}

	if err = keygen.check_sk(sk_output); err != nil {
		return nil, err
	}

	keygen.context.Sub(sk_input.Get(), sk_output.Get(), keygen.polypool)
	newevakey = newswitchingkey(keygen.ckkscontext, keygen.polypool, sk_output.Get(), bitDecomp)
	keygen.polypool.Zero()

	return
}

func (keygen *keygenerator) NewSwitchingKeyEmpty(bitDecomp uint64) (evakey *SwitchingKey) {
	evakey = new(SwitchingKey)

	if bitDecomp > keygen.ckkscontext.maxBit || bitDecomp == 0 {
		bitDecomp = keygen.ckkscontext.maxBit
	}

	context := keygen.ckkscontext.keyscontext

	evakey.bitDecomp = bitDecomp

	// delta_sk = sk_input - sk_output = GaloisEnd(sk_output, rotation) - sk_output
	var bitLog uint64

	evakey.evakey = make([][][2]*ring.Poly, len(context.Modulus))

	for i, qi := range context.Modulus {

		bitLog = uint64(math.Ceil(float64(bits.Len64(qi)) / float64(bitDecomp)))

		evakey.evakey[i] = make([][2]*ring.Poly, bitLog)

		for j := uint64(0); j < bitLog; j++ {
			evakey.evakey[i][j][0] = context.NewPoly()
			evakey.evakey[i][j][1] = context.NewPoly()
		}
	}

	return
}

// NewRotationKeys generates a new instance of rotationkeys, with the provided rotation to the left, right and conjugation if asked.
// Here bitdecomp plays a role in the added noise if the scale of the input is smaller than the maximum size between the modulies.
func (keygen *keygenerator) NewRotationKeys(sk_output *SecretKey, bitDecomp uint64, rotLeft []uint64, rotRight []uint64, conjugate bool) (rotKey *RotationKey, err error) {

	if err = keygen.check_sk(sk_output); err != nil {
		return nil, err
	}

	if bitDecomp > keygen.ckkscontext.maxBit || bitDecomp == 0 {
		bitDecomp = keygen.ckkscontext.maxBit
	}

	rotKey = new(RotationKey)
	rotKey.ckkscontext = keygen.ckkscontext

	if rotLeft != nil {
		rotKey.evakey_rot_col_L = make(map[uint64]*SwitchingKey)
		for _, n := range rotLeft {
			if rotKey.evakey_rot_col_L[n] == nil && n != 0 {
				rotKey.evakey_rot_col_L[n] = genrotkey(keygen, sk_output.Get(), keygen.ckkscontext.galElRotColLeft[n], bitDecomp)
			}
		}
	}

	if rotRight != nil {
		rotKey.evakey_rot_col_R = make(map[uint64]*SwitchingKey)
		for _, n := range rotRight {
			if rotKey.evakey_rot_col_R[n] == nil && n != 0 {
				rotKey.evakey_rot_col_R[n] = genrotkey(keygen, sk_output.Get(), keygen.ckkscontext.galElRotColRight[n], bitDecomp)
			}
		}
	}

	if conjugate {
		rotKey.evakey_rot_row = genrotkey(keygen, sk_output.Get(), keygen.ckkscontext.galElRotRow, bitDecomp)
	}

	return rotKey, nil

}

// NewRotationKeys generates a new instance of rotationkeys, with the provided rotation to the left, right and conjugation if asked.
// Here bitdecomp plays a role in the added noise if the scale of the input is smaller than the maximum size between the modulies.
func (keygen *keygenerator) NewRotationKeysEmpty() (rotKey *RotationKey) {

	rotKey = new(RotationKey)
	rotKey.ckkscontext = keygen.ckkscontext

	return rotKey

}

// NewRotationkeysPow2 generates a new rotation key with all the power of two rotation to the left and right, as well as the conjugation
// key if asked. Here bitdecomp plays a role in the added noise if the scale of the input is smaller than the maximum size between the modulies.
func (keygen *keygenerator) NewRotationKeysPow2(sk_output *SecretKey, bitDecomp uint64, conjugate bool) (rotKey *RotationKey, err error) {

	if err = keygen.check_sk(sk_output); err != nil {
		return nil, err
	}

	if bitDecomp > keygen.ckkscontext.maxBit || bitDecomp == 0 {
		bitDecomp = keygen.ckkscontext.maxBit
	}

	rotKey = new(RotationKey)
	rotKey.ckkscontext = keygen.ckkscontext

	rotKey.evakey_rot_col_L = make(map[uint64]*SwitchingKey)
	rotKey.evakey_rot_col_R = make(map[uint64]*SwitchingKey)

	for n := uint64(1); n < rotKey.ckkscontext.n>>1; n <<= 1 {

		rotKey.evakey_rot_col_L[n] = genrotkey(keygen, sk_output.Get(), keygen.ckkscontext.galElRotColLeft[n], bitDecomp)
		rotKey.evakey_rot_col_R[n] = genrotkey(keygen, sk_output.Get(), keygen.ckkscontext.galElRotColRight[n], bitDecomp)
	}

	if conjugate {
		rotKey.evakey_rot_row = genrotkey(keygen, sk_output.Get(), keygen.ckkscontext.galElRotRow, bitDecomp)
	}

	return
}

func genrotkey(keygen *keygenerator, sk_output *ring.Poly, gen, bitDecomp uint64) (switchingkey *SwitchingKey) {

	ring.PermuteNTT(sk_output, gen, keygen.polypool)
	keygen.context.Sub(keygen.polypool, sk_output, keygen.polypool)
	switchingkey = newswitchingkey(keygen.ckkscontext, keygen.polypool, sk_output, bitDecomp)
	keygen.polypool.Zero()

	return
}

func newswitchingkey(ckkscontext *CkksContext, sk_in, sk_out *ring.Poly, bitDecomp uint64) (switchingkey *SwitchingKey) {

	if bitDecomp > ckkscontext.maxBit || bitDecomp == 0 {
		bitDecomp = ckkscontext.maxBit
	}

	switchingkey = new(SwitchingKey)

	context := ckkscontext.keyscontext

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
			switchingkey.evakey[i][j][0] = ckkscontext.gaussianSampler.SampleNTTNew()
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
