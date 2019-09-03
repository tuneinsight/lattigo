package bfv

import (
	"errors"
	"github.com/lca1/lattigo/ring"
	"math"
	"math/bits"
)

// Keygenerator is a structure that stores the elements required to create new keys,
// as well as a small memory pool for intermediate values.
type keygenerator struct {
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

// Rotationkeys is a structure that stores the switching-keys required during the homomorphic rotations.
type RotationKeys struct {
	bfvcontext       *BfvContext
	bitDecomp        uint64
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
	bitDecomp uint64
	evakey    [][][2]*ring.Poly
}

// Newkeygenerator creates a new KeyGenerator from the target bfvcontext.
func (bfvcontext *BfvContext) NewKeyGenerator() (keygen *keygenerator) {
	keygen = new(keygenerator)
	keygen.bfvcontext = bfvcontext
	keygen.context = bfvcontext.contextQ
	keygen.polypool = keygen.context.NewPoly()
	return
}

// Newsecretkey creates a new SecretKey with uniform distribution in [-1, 0, 1].
func (keygen *keygenerator) NewSecretKey() *SecretKey {

	sk := new(SecretKey)
	sk.sk = keygen.bfvcontext.ternarySampler.SampleMontgomeryNTTNew()

	return sk
}

// NewSecretKeyEmpty creates a new SecretKey with all coeffcients set to zero, ready to received a marshaled SecretKey.
func (keygen *keygenerator) NewSecretKeyEmpty() *SecretKey {
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

// check_sk checks if the input secret-key complies with the keygenerator context.
func (keygen *keygenerator) check_sk(sk_output *SecretKey) (err error) {

	if sk_output.Get().GetDegree() != int(keygen.context.N) {
		return errors.New("error : pol degree sk != bfvcontext.n")
	}

	if len(sk_output.Get().Coeffs) != len(keygen.context.Modulus) {
		return errors.New("error : nb modulus sk != nb modulus bfvcontext")
	}

	return nil
}

// Newpublickey generates a new publickkey from the provided secret-key
func (keygen *keygenerator) NewPublicKey(sk *SecretKey) (pk *PublicKey, err error) {

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

func (keygen *keygenerator) NewPublicKeyEmpty() (pk *PublicKey) {
	pk = new(PublicKey)

	pk.pk[0] = keygen.context.NewPoly()
	pk.pk[1] = keygen.context.NewPoly()

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

// NewKeyPair generates a new (secret-key, public-key) pair.
func (keygen *keygenerator) NewKeyPair() (sk *SecretKey, pk *PublicKey, err error) {
	sk = keygen.NewSecretKey()
	pk, err = keygen.NewPublicKey(sk)
	return
}

// NewRelinKey generates a new evaluation key from the provided secret-key. It will be used to relinearize a ciphertext (encrypted under a public-key generated from the provided secret-key)
// of degree > 1 to a ciphertext of degree 1. Max degree is the maximum degree of the ciphertext allowed to relinearize. Bitdecomp is the power of two binary decomposition of the key.
// A higher bigdecomp will induce smaller keys, faster key-switching, but at the cost of more noise.
func (keygen *keygenerator) NewRelinKey(sk *SecretKey, maxDegree, bitDecomp uint64) (newEvakey *EvaluationKey, err error) {

	newEvakey = new(EvaluationKey)
	newEvakey.evakey = make([]*SwitchingKey, maxDegree)
	sk.Get().Copy(keygen.polypool)
	for i := uint64(0); i < maxDegree; i++ {
		keygen.context.MulCoeffsMontgomery(keygen.polypool, sk.Get(), keygen.polypool)
		newEvakey.evakey[i] = newswitchintkey(keygen.bfvcontext, keygen.polypool, sk.Get(), bitDecomp)
	}
	keygen.polypool.Zero()

	return newEvakey, nil
}

func (keygen *keygenerator) NewRelinKeyEmpty(maxDegree, bitDecomp uint64) (evakey *EvaluationKey) {
	evakey = new(EvaluationKey)

	if bitDecomp > keygen.bfvcontext.maxBit || bitDecomp == 0 {
		bitDecomp = keygen.bfvcontext.maxBit
	}

	context := keygen.bfvcontext.contextQ

	// delta_sk = sk_input - sk_output = GaloisEnd(sk_output, rotation) - sk_output
	var bitLog uint64

	evakey.evakey = make([]*SwitchingKey, maxDegree)

	for w := uint64(0); w < maxDegree; w++ {

		evakey.evakey[w] = new(SwitchingKey)
		evakey.evakey[w].bitDecomp = bitDecomp
		evakey.evakey[w].evakey = make([][][2]*ring.Poly, len(context.Modulus))

		for i, qi := range context.Modulus {

			bitLog = uint64(math.Ceil(float64(bits.Len64(qi)) / float64(bitDecomp)))

			evakey.evakey[w].evakey[i] = make([][2]*ring.Poly, bitLog)

			for j := uint64(0); j < bitLog; j++ {
				evakey.evakey[w].evakey[i][j][0] = context.NewPoly()
				evakey.evakey[w].evakey[i][j][1] = context.NewPoly()
			}
		}
	}

	return
}

// Get returns the slice of switchintkeys of the evaluation-key.
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

// Newswitchintkey generates a new key-switching key, that will allow to re-encrypt under the output-key a ciphertext encrypted under the input-key. Bitdecomp
// is the power of two binary decomposition of the key. A higher bigdecomp will induce smaller keys, faster key-switching, but at the cost of more noise.
func (keygen *keygenerator) NewSwitchingKey(sk_input, sk_output *SecretKey, bitDecomp uint64) (newevakey *SwitchingKey, err error) {

	if err = keygen.check_sk(sk_input); err != nil {
		return nil, err
	}

	if err = keygen.check_sk(sk_output); err != nil {
		return nil, err
	}

	keygen.context.Sub(sk_input.Get(), sk_output.Get(), keygen.polypool)
	newevakey = newswitchintkey(keygen.bfvcontext, keygen.polypool, sk_output.Get(), bitDecomp)
	keygen.polypool.Zero()

	return
}

func (keygen *keygenerator) NewSwitchingKeyEmpty(bitDecomp uint64) (evakey *SwitchingKey) {
	evakey = new(SwitchingKey)

	if bitDecomp > keygen.bfvcontext.maxBit || bitDecomp == 0 {
		bitDecomp = keygen.bfvcontext.maxBit
	}

	context := keygen.bfvcontext.contextQ

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

// Newrotationkeys generates a new struct of rotationkeys storing the keys for the specified rotations. The provided secret-key must be the secret-key used to generate the public-key under
// which the ciphertexts to rotate are encrypted under. Bitdecomp is the power of two binary decomposition of the key. A higher bigdecomp will induce smaller keys, faster key-switching,
// but at the cost of more noise. rotLeft and rotRight must be a slice of uint64 rotations, row is a boolean value indicating if the key for the row rotation must be generated.
func (keygen *keygenerator) NewRotationKeys(sk *SecretKey, bitDecomp uint64, rotLeft []uint64, rotRight []uint64, row bool) (rotKey *RotationKeys, err error) {

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

func (keygen *keygenerator) NewRotationKeysEmpty() (rotKey *RotationKeys) {

	rotKey = new(RotationKeys)
	rotKey.bfvcontext = keygen.bfvcontext

	return rotKey

}

// Newrotationkeys generates a new struct of rotationkeys storing the keys of all the left and right powers of two rotations. The provided secret-key must be the secret-key used to generate the public-key under
// which the ciphertexts to rotate are encrypted under. rows is a boolean value indicatig if the keys for the row rotation have to be generated. Bitdecomp is the power of two binary decomposition of the key.
// A higher bigdecomp will induce smaller keys, faster key-switching, but at the cost of more noise.
func (keygen *keygenerator) NewRotationKeysPow2(sk *SecretKey, bitDecomp uint64, row bool) (rotKey *RotationKeys, err error) {

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
func genrotkey(keygen *keygenerator, sk *ring.Poly, gen, bitDecomp uint64) (switchkey *SwitchingKey) {

	ring.PermuteNTT(sk, gen, keygen.polypool)
	keygen.context.Sub(keygen.polypool, sk, keygen.polypool)
	switchkey = newswitchintkey(keygen.bfvcontext, keygen.polypool, sk, bitDecomp)
	keygen.polypool.Zero()

	return
}

// newswitchintkey is a generic methode to generate key-switching keys used in the evaluation, key-switching and rotation-keys generation.
func newswitchintkey(bfvcontext *BfvContext, sk_in, sk_out *ring.Poly, bitDecomp uint64) (switchkey *SwitchingKey) {

	if bitDecomp > bfvcontext.maxBit || bitDecomp == 0 {
		bitDecomp = bfvcontext.maxBit
	}

	switchkey = new(SwitchingKey)

	context := bfvcontext.contextQ

	switchkey.bitDecomp = uint64(bitDecomp)

	mredParams := context.GetMredParams()

	// delta_sk = sk_input - sk_output = GaloisEnd(sk_output, rotation) - sk_output

	var bitLog uint64

	switchkey.evakey = make([][][2]*ring.Poly, len(context.Modulus))

	for i, qi := range context.Modulus {

		bitLog = uint64(math.Ceil(float64(bits.Len64(qi)) / float64(bitDecomp)))

		switchkey.evakey[i] = make([][2]*ring.Poly, bitLog)

		for j := uint64(0); j < bitLog; j++ {

			// e
			switchkey.evakey[i][j][0] = bfvcontext.gaussianSampler.SampleNTTNew()
			// a
			switchkey.evakey[i][j][1] = context.NewUniformPoly()

			// e + sk_in * (qiBarre*qiStar) * 2^w
			// (qiBarre*qiStar)%qi = 1, else 0
			for w := uint64(0); w < context.N; w++ {
				switchkey.evakey[i][j][0].Coeffs[i][w] += ring.PowerOf2(sk_in.Coeffs[i][w], bitDecomp*j, qi, mredParams[i])
			}

			// sk_in * (qiBarre*qiStar) * 2^w - a*sk + e
			context.MulCoeffsMontgomeryAndSub(switchkey.evakey[i][j][1], sk_out, switchkey.evakey[i][j][0])

			context.MForm(switchkey.evakey[i][j][0], switchkey.evakey[i][j][0])
			context.MForm(switchkey.evakey[i][j][1], switchkey.evakey[i][j][1])
		}
	}

	return
}
