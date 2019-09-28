package ckks

import (
	"encoding/binary"
	"errors"
	"github.com/ldsec/lattigo/ring"
	"math"
	"math/bits"
)

// Keygenerator is a structure that stores the elements required to create new keys,
// as well as a small memory pool for intermediate values.
type KeyGenerator struct {
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
func (ckkscontext *CkksContext) NewKeyGenerator() (keygen *KeyGenerator) {
	keygen = new(KeyGenerator)
	keygen.ckkscontext = ckkscontext
	keygen.context = ckkscontext.keyscontext
	keygen.polypool = ckkscontext.keyscontext.NewPoly()
	return
}

// check_sk checks if the input secret-key complies with the keygenerator context.
func (keygen *KeyGenerator) check_sk(sk_output *SecretKey) error {

	if sk_output.Get().GetDegree() != int(keygen.context.N) {
		return errors.New("error : pol degree sk != ckkscontext.n")
	}

	if len(sk_output.Get().Coeffs) != len(keygen.context.Modulus) {
		return errors.New("error : nb modulus sk != nb modulus ckkscontext")
	}

	return nil
}

// NewSecretKey generates a new secret key.
func (keygen *KeyGenerator) NewSecretKey(p float64) (sk *SecretKey, err error) {
	sk = new(SecretKey)
	if sk.sk, err = keygen.ckkscontext.ternarySampler.SampleMontgomeryNTTNew(p); err != nil {
		return nil, err
	}
	return sk, nil
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
func (keygen *KeyGenerator) NewPublicKey(sk *SecretKey) (pk *PublicKey, err error) {

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

// NewKeyPair generates a new secretkey and a corresponding public key.
func (keygen *KeyGenerator) NewKeyPair(p float64) (sk *SecretKey, pk *PublicKey, err error) {
	if sk, err = keygen.NewSecretKey(p); err != nil {
		return nil, nil, err
	}
	if pk, err = keygen.NewPublicKey(sk); err != nil {
		return nil, nil, err
	}
	return
}

// NewRelinkey generates a new evaluation key that will be used to relinearize the ciphertexts during multiplication.
// Bitdecomposition aims at reducing the added noise at the expense of more storage needed for the keys and more computation
// during the relinearization. However for relinearization this bitdecomp value can be set to maximum as the encrypted value
// are also scaled up during the multiplication.
func (keygen *KeyGenerator) NewRelinKey(sk *SecretKey, bitDecomp uint64) (evakey *EvaluationKey, err error) {

	if err = keygen.check_sk(sk); err != nil {
		return nil, err
	}
	evakey = new(EvaluationKey)
	sk.Get().Copy(keygen.polypool)
	keygen.context.MulCoeffsMontgomery(keygen.polypool, sk.Get(), keygen.polypool)
	evakey.evakey = keygen.newSwitchingKey(keygen.polypool, sk.Get(), bitDecomp)
	keygen.polypool.Zero()

	return
}

func (keygen *KeyGenerator) NewRelinKeyEmpty(bitDecomp uint64) (evakey *EvaluationKey) {
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

func (keygen *KeyGenerator) SetRelinKeys(rlk [][][2]*ring.Poly, bitDecomp uint64) (*EvaluationKey, error) {

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
func (keygen *KeyGenerator) NewSwitchingKey(sk_input, sk_output *SecretKey, bitDecomp uint64) (newevakey *SwitchingKey, err error) {

	if err = keygen.check_sk(sk_input); err != nil {
		return nil, err
	}

	if err = keygen.check_sk(sk_output); err != nil {
		return nil, err
	}

	keygen.context.Sub(sk_input.Get(), sk_output.Get(), keygen.polypool)
	newevakey = keygen.newSwitchingKey(keygen.polypool, sk_output.Get(), bitDecomp)
	keygen.polypool.Zero()

	return
}

func (keygen *KeyGenerator) NewSwitchingKeyEmpty(bitDecomp uint64) (evakey *SwitchingKey) {
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
func (keygen *KeyGenerator) NewRotationKeys(sk_output *SecretKey, bitDecomp uint64, rotLeft []uint64, rotRight []uint64, conjugate bool) (rotKey *RotationKey, err error) {

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
				rotKey.evakey_rot_col_L[n] = keygen.genrotKey(sk_output.Get(), keygen.ckkscontext.galElRotColLeft[n], bitDecomp)
			}
		}
	}

	if rotRight != nil {
		rotKey.evakey_rot_col_R = make(map[uint64]*SwitchingKey)
		for _, n := range rotRight {
			if rotKey.evakey_rot_col_R[n] == nil && n != 0 {
				rotKey.evakey_rot_col_R[n] = keygen.genrotKey(sk_output.Get(), keygen.ckkscontext.galElRotColRight[n], bitDecomp)
			}
		}
	}

	if conjugate {
		rotKey.evakey_rot_row = keygen.genrotKey(sk_output.Get(), keygen.ckkscontext.galElRotRow, bitDecomp)
	}

	return rotKey, nil

}

// NewRotationKeys generates a new instance of rotationkeys, with the provided rotation to the left, right and conjugation if asked.
// Here bitdecomp plays a role in the added noise if the scale of the input is smaller than the maximum size between the modulies.
func (keygen *KeyGenerator) NewRotationKeysEmpty() (rotKey *RotationKey) {

	rotKey = new(RotationKey)
	rotKey.ckkscontext = keygen.ckkscontext
	return
}

// NewRotationkeysPow2 generates a new rotation key with all the power of two rotation to the left and right, as well as the conjugation
// key if asked. Here bitdecomp plays a role in the added noise if the scale of the input is smaller than the maximum size between the modulies.
func (keygen *KeyGenerator) NewRotationKeysPow2(sk_output *SecretKey, bitDecomp uint64, conjugate bool) (rotKey *RotationKey, err error) {

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
		rotKey.evakey_rot_col_L[n] = keygen.genrotKey(sk_output.Get(), keygen.ckkscontext.galElRotColLeft[n], bitDecomp)
		rotKey.evakey_rot_col_R[n] = keygen.genrotKey(sk_output.Get(), keygen.ckkscontext.galElRotColRight[n], bitDecomp)
	}

	if conjugate {
		rotKey.evakey_rot_row = keygen.genrotKey(sk_output.Get(), keygen.ckkscontext.galElRotRow, bitDecomp)
	}

	return
}

func (keygen *KeyGenerator) genrotKey(sk_output *ring.Poly, gen, bitDecomp uint64) (switchingkey *SwitchingKey) {

	ring.PermuteNTT(sk_output, gen, keygen.polypool)
	keygen.context.Sub(keygen.polypool, sk_output, keygen.polypool)
	switchingkey = keygen.newSwitchingKey(keygen.polypool, sk_output, bitDecomp)
	keygen.polypool.Zero()

	return
}

func (keygen *KeyGenerator) newSwitchingKey(sk_in, sk_out *ring.Poly, bitDecomp uint64) (switchingkey *SwitchingKey) {

	if bitDecomp > keygen.ckkscontext.maxBit || bitDecomp == 0 {
		bitDecomp = keygen.ckkscontext.maxBit
	}

	switchingkey = new(SwitchingKey)

	context := keygen.ckkscontext.keyscontext

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
			switchingkey.evakey[i][j][0] = keygen.ckkscontext.gaussianSampler.SampleNTTNew()
			// a
			switchingkey.evakey[i][j][1] = context.NewUniformPoly()

			// e + sk_in * (qiBarre*qiStar) * 2^w
			// (qiBarre*qiStar)%qi = 1 mod qi, else 0
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

// MarshalBinary encodes a secret-key on a byte slice. The total size in byte is 1 + N/4.
func (sk *SecretKey) MarshalBinary() (data []byte, err error) {

	N := uint64(len(sk.sk.Coeffs[0]))
	levels := uint64(len(sk.sk.Coeffs))

	if levels > 0xFF {
		return nil, errors.New("error : max degree uint8 overflow")
	}

	data = make([]byte, 2+((N*levels)<<3))

	data[0] = uint8((bits.Len64(uint64(N)) - 1))
	data[1] = uint8(levels)

	pointer := uint64(2)

	if _, err = ring.WriteCoeffsTo(pointer, N, levels, sk.sk.Coeffs, data); err != nil {
		return nil, err
	}

	return
}

// UnMarshalBinary decode a previously marshaled secret-key on the target secret-key.
// The target secret-key must be of the appropriate format, it can be created with the methode NewSecretKeyEmpty().
func (sk *SecretKey) UnMarshalBinary(data []byte) (err error) {

	N := uint64(1 << data[0])
	levels := uint64(data[1])

	pointer := uint64(2)

	if uint64(len(sk.sk.Coeffs[0])) != N {
		return errors.New("error : invalid publickey[0] receiver (logN do not match)")
	}

	if uint64(len(sk.sk.Coeffs)) != levels {
		return errors.New("error : invalid SecretKey receiver (level do not match data)")
	}

	if ((uint64(len(data)) - pointer) >> 3) != (N * levels) {
		return errors.New("error : invalid SecretKey encoding")
	}

	ring.DecodeCoeffs(pointer, N, levels, sk.sk.Coeffs, data)

	return
}

// MarshalBinary encodes a public-key on a byte slice. The total size is 2 + 16 * N * (level + 1).
func (pk *PublicKey) MarshalBinary() (data []byte, err error) {

	N := uint64(len(pk.pk[0].Coeffs[0]))
	levels := uint64(len(pk.pk[0].Coeffs))

	if levels > 0xFF {
		return nil, errors.New("error : max degree uint8 overflow")
	}

	data = make([]byte, 2+((N*levels)<<4))

	data[0] = uint8((bits.Len64(uint64(N)) - 1))
	data[1] = uint8(levels)

	pointer := uint64(2)

	if pointer, err = ring.WriteCoeffsTo(pointer, N, levels, pk.pk[0].Coeffs, data); err != nil {
		return nil, err
	}

	if pointer, err = ring.WriteCoeffsTo(pointer, N, levels, pk.pk[1].Coeffs, data); err != nil {
		return nil, err
	}

	return
}

// UnMarshalBinary decodes a previously marshaled public-key on the target public-key.
// The target public-key must have the appropriate format and size, it can be created with
// the methode NewPublicKeyEmpty().
func (pk *PublicKey) UnMarshalBinary(data []byte) (err error) {

	N := uint64(1 << data[0])
	levels := uint64(data[1])

	pointer := uint64(2)

	if uint64(len(pk.pk[0].Coeffs[0])) != N {
		return errors.New("error : invalid publickey[0] receiver (logN do not match)")
	}

	if uint64(len(pk.pk[0].Coeffs[1])) != N {
		return errors.New("error : invalid publickey[1] receiver (logN do not match)")
	}

	if uint64(len(pk.pk[0].Coeffs)) != levels {
		return errors.New("error : invalid publickey[0] receiver (level do not match data)")
	}

	if uint64(len(pk.pk[1].Coeffs)) != levels {
		return errors.New("error : invalid publickey[1] receiver (level do not match data)")
	}

	if ((uint64(len(data)) - pointer) >> 4) != (N * levels) {
		return errors.New("error : invalid PublicKey encoding")
	}

	pointer, _ = ring.DecodeCoeffs(pointer, N, levels, pk.pk[0].Coeffs, data)
	pointer, _ = ring.DecodeCoeffs(pointer, N, levels, pk.pk[1].Coeffs, data)

	return
}

// MarshalBinary encodes an evaluation key on a byte slice. The total size depends on each modulus size and the bit decomp.
// It will approximately be 5 + (level + 1) * ( 1 + 2 * 8 * N * (level + 1) * logQi/bitDecomp).
func (evaluationkey *EvaluationKey) MarshalBinary() (data []byte, err error) {

	N := uint64(len(evaluationkey.evakey.evakey[0][0][0].Coeffs[0]))
	levels := uint64(len(evaluationkey.evakey.evakey[0][0][0].Coeffs))
	decomposition := levels
	bitDecomp := evaluationkey.evakey.bitDecomp

	if levels > 0xFF {
		return nil, errors.New("error : max number modulis uint8 overflow")
	}

	if decomposition > 0xFF {
		return nil, errors.New("error : max decomposition uint8 overflow")
	}

	if bitDecomp > 0xFF {
		return nil, errors.New("error : max bitDecomp uint8 overflow")
	}

	var dataLen uint64
	dataLen = 4

	for j := uint64(0); j < decomposition; j++ {
		dataLen += 1                                                                                //Information about the size of the bitdecomposition
		dataLen += 2 * 8 * N * levels * decomposition * uint64(len(evaluationkey.evakey.evakey[j])) // nb coefficients * 8
	}

	data = make([]byte, dataLen)

	data[0] = uint8(bits.Len64(uint64(N)) - 1)
	data[1] = uint8(levels)
	data[2] = uint8(decomposition)
	data[3] = uint8(bitDecomp)

	pointer := uint64(4)

	var bitLog uint8

	for j := uint64(0); j < decomposition; j++ {
		bitLog = uint8(len(evaluationkey.evakey.evakey[j]))
		data[pointer] = bitLog
		pointer += 1
		for x := uint8(0); x < bitLog; x++ {
			if pointer, err = ring.WriteCoeffsTo(pointer, N, levels, evaluationkey.evakey.evakey[j][x][0].Coeffs, data); err != nil {
				return nil, err
			}

			if pointer, err = ring.WriteCoeffsTo(pointer, N, levels, evaluationkey.evakey.evakey[j][x][1].Coeffs, data); err != nil {
				return nil, err
			}
		}
	}

	return
}

// UnMarshalBinary decodes a previously marshaled evaluation-key on the target evaluation-key. The target evaluation-key
// must have the appropriate format and size, it can be created with the methode NewRelinKeyEmpty(uint64, uint64).
func (evaluationkey *EvaluationKey) UnMarshalBinary(data []byte) (err error) {

	N := uint64(1 << data[0])
	levels := uint64(data[1])
	decomposition := uint64(data[2])
	bitDecomp := uint64(data[3])

	pointer := uint64(4)
	var bitLog uint64

	evaluationkey.evakey.bitDecomp = bitDecomp

	for j := uint64(0); j < decomposition; j++ {

		bitLog = uint64(data[pointer])
		pointer += 1

		for x := uint64(0); x < bitLog; x++ {

			if uint64(len(evaluationkey.evakey.evakey[j][x][0].Coeffs)) != levels {
				return errors.New("error : evaluationkey receiver (level do not match data)")
			}

			if uint64(len(evaluationkey.evakey.evakey[j][x][1].Coeffs)) != levels {
				return errors.New("error : evaluationkey receiver (level do not match data)")
			}

			pointer, _ = ring.DecodeCoeffs(pointer, N, levels, evaluationkey.evakey.evakey[j][x][0].Coeffs, data)
			pointer, _ = ring.DecodeCoeffs(pointer, N, levels, evaluationkey.evakey.evakey[j][x][1].Coeffs, data)
		}
	}

	return
}

// MarshalBinary encodes an switching-key on a byte slice. The total size in byte will be approximately 5 + (level + 1) * ( 1 + 2 * 8 * N * (level + 1) * logQi/bitDecomp).
func (switchingkey *SwitchingKey) MarshalBinary() (data []byte, err error) {

	N := uint64(len(switchingkey.evakey[0][0][0].Coeffs[0]))
	level := uint64(len(switchingkey.evakey[0][0][0].Coeffs))
	decomposition := level
	bitDecomp := switchingkey.bitDecomp

	if level > 0xFF {
		return nil, errors.New("error : max number modulis uint8 overflow")
	}

	if decomposition > 0xFF {
		return nil, errors.New("error : max decomposition uint8 overflow")
	}

	if bitDecomp > 0xFF {
		return nil, errors.New("error : max bitDecomp uint8 overflow")
	}

	var dataLen uint64
	dataLen = 4

	for j := uint64(0); j < decomposition; j++ {
		dataLen += 1                                                                       //Information about the size of the bitdecomposition
		dataLen += 2 * 8 * N * level * decomposition * uint64(len(switchingkey.evakey[j])) // nb coefficients * 8
	}

	data = make([]byte, dataLen)

	data[0] = uint8(bits.Len64(uint64(N)) - 1)
	data[1] = uint8(level)
	data[2] = uint8(decomposition)
	data[3] = uint8(bitDecomp)

	pointer := uint64(4)

	var bitLog uint8

	for j := uint64(0); j < decomposition; j++ {
		bitLog = uint8(len(switchingkey.evakey[j]))
		data[pointer] = bitLog
		pointer += 1
		for x := uint8(0); x < bitLog; x++ {
			if pointer, err = ring.WriteCoeffsTo(pointer, N, level, switchingkey.evakey[j][x][0].Coeffs, data); err != nil {
				return nil, err
			}

			if pointer, err = ring.WriteCoeffsTo(pointer, N, level, switchingkey.evakey[j][x][1].Coeffs, data); err != nil {
				return nil, err
			}
		}
	}

	return
}

// UnMarshalBinary decode a previously marshaled switching-key on the target switching-key.
// The target switching-key must have the appropriate format and size, it can be created with the methode NewSwitchingKeyEmpty(uint64).
func (switchingkey *SwitchingKey) UnMarshalBinary(data []byte) (err error) {

	N := uint64(1 << data[0])
	level := uint64(data[1])
	decomposition := uint64(data[2])
	bitDecomp := uint64(data[3])

	pointer := uint64(4)
	var bitLog uint64

	switchingkey.bitDecomp = bitDecomp

	for j := uint64(0); j < decomposition; j++ {

		bitLog = uint64(data[pointer])
		pointer += 1

		for x := uint64(0); x < bitLog; x++ {
			pointer, _ = ring.DecodeCoeffs(pointer, N, level, switchingkey.evakey[j][x][0].Coeffs, data)
			pointer, _ = ring.DecodeCoeffs(pointer, N, level, switchingkey.evakey[j][x][1].Coeffs, data)
		}
	}

	return
}

// MarshalBinary encodes a rotationkeys structure on a byte slice. The total size in byte is approximately
// 5 + 4*(nb left rot + num right rot) + (nb left rot + num right rot + 1 (if rotate row)) * (level + 1) * ( 1 + 2 * 8 * N * (level + 1) * logQi/bitDecomp).
func (rotationkey *RotationKey) MarshalBinary() (data []byte, err error) {

	N := uint64(rotationkey.ckkscontext.n)
	level := uint64(len(rotationkey.ckkscontext.keyscontext.Modulus))
	decomposition := level
	bitDecomp := rotationkey.bitDecomp
	mappingRow := 0
	mappingColL := []uint64{}
	mappingColR := []uint64{}

	if level > 0xFF {
		return nil, errors.New("error : max number modulis uint8 overflow")
	}

	if decomposition > 0xFF {
		return nil, errors.New("error : max decomposition uint8 overflow")
	}

	if bitDecomp > 0xFF {
		return nil, errors.New("error : max bitDecomp uint8 overflow")
	}

	var dataLen uint64
	dataLen = 13

	for i := uint64(1); i < N>>1; i++ {
		if rotationkey.evakey_rot_col_L[i] != nil {

			mappingColL = append(mappingColL, i)

			for j := uint64(0); j < decomposition; j++ {
				dataLen += 1                                                                                          //Information about the size of the bitdecomposition
				dataLen += 2 * 8 * N * level * decomposition * uint64(len(rotationkey.evakey_rot_col_L[i].evakey[j])) // nb coefficients * 8
			}
		}

		if rotationkey.evakey_rot_col_L[i] != nil {

			mappingColR = append(mappingColR, i)

			for j := uint64(0); j < decomposition; j++ {
				dataLen += 1                                                                                          //Information about the size of the bitdecomposition
				dataLen += 2 * 8 * N * level * decomposition * uint64(len(rotationkey.evakey_rot_col_L[i].evakey[j])) // nb coefficients * 8
			}
		}
	}

	if rotationkey.evakey_rot_row != nil {
		mappingRow = 1
		for j := uint64(0); j < decomposition; j++ {
			dataLen += 1                                                                                     //Information about the size of the bitdecomposition
			dataLen += 2 * 8 * N * level * decomposition * uint64(len(rotationkey.evakey_rot_row.evakey[j])) // nb coefficients * 8
		}
	}

	dataLen += uint64(len(mappingColL)+len(mappingColR)) << 2 // size needed to encode what rotation are present

	data = make([]byte, dataLen)

	data[0] = uint8(bits.Len64(uint64(N)) - 1)
	data[1] = uint8(level)
	data[2] = uint8(decomposition)
	data[3] = uint8(bitDecomp)
	data[4] = uint8(mappingRow)

	pointer := uint64(5)

	binary.BigEndian.PutUint32(data[pointer:pointer+4], uint32(len(mappingColL)))
	pointer += 4

	binary.BigEndian.PutUint32(data[pointer:pointer+4], uint32(len(mappingColR)))
	pointer += 4

	for _, i := range mappingColL {

		binary.BigEndian.PutUint32(data[pointer:pointer+4], uint32(i))

		pointer += 4
	}

	for _, i := range mappingColR {

		binary.BigEndian.PutUint32(data[pointer:pointer+4], uint32(i))

		pointer += 4
	}

	// Encodes the different rotation key indexes
	var bitLog uint8
	if mappingRow == 1 {
		for j := uint64(0); j < decomposition; j++ {
			bitLog = uint8(len(rotationkey.evakey_rot_row.evakey[j]))
			data[pointer] = bitLog
			pointer += 1
			for x := uint8(0); x < bitLog; x++ {
				if pointer, err = ring.WriteCoeffsTo(pointer, N, level, rotationkey.evakey_rot_row.evakey[j][x][0].Coeffs, data); err != nil {
					return nil, err
				}

				if pointer, err = ring.WriteCoeffsTo(pointer, N, level, rotationkey.evakey_rot_row.evakey[j][x][1].Coeffs, data); err != nil {
					return nil, err
				}
			}
		}
	}

	for _, i := range mappingColL {
		for j := uint64(0); j < decomposition; j++ {
			bitLog = uint8(len(rotationkey.evakey_rot_col_L[i].evakey[j]))
			data[pointer] = bitLog
			pointer += 1
			for x := uint8(0); x < bitLog; x++ {
				if pointer, err = ring.WriteCoeffsTo(pointer, N, level, rotationkey.evakey_rot_col_L[i].evakey[j][x][0].Coeffs, data); err != nil {
					return nil, err
				}

				if pointer, err = ring.WriteCoeffsTo(pointer, N, level, rotationkey.evakey_rot_col_L[i].evakey[j][x][1].Coeffs, data); err != nil {
					return nil, err
				}
			}
		}
	}

	for _, i := range mappingColR {
		for j := uint64(0); j < decomposition; j++ {
			bitLog = uint8(len(rotationkey.evakey_rot_col_R[i].evakey[j]))
			data[pointer] = bitLog
			pointer += 1
			for x := uint8(0); x < bitLog; x++ {
				if pointer, err = ring.WriteCoeffsTo(pointer, N, level, rotationkey.evakey_rot_col_R[i].evakey[j][x][0].Coeffs, data); err != nil {
					return nil, err
				}

				if pointer, err = ring.WriteCoeffsTo(pointer, N, level, rotationkey.evakey_rot_col_R[i].evakey[j][x][1].Coeffs, data); err != nil {
					return nil, err
				}
			}
		}
	}

	return
}

// UnMarshalBinary decodes a previously marshaled rotation-keys on the target rotation-keys. In contrary to all
// the other structures, the unmarshaling for rotationkeys only need an empty receiver, as it is not possible to
// create receiver of the correct format and size without knowing all the content of the marshaled rotationkeys. The memory
// will be allocated on the fly.
func (rotationkey *RotationKey) UnMarshalBinary(data []byte) (err error) {

	N := uint64(1 << data[0])
	level := uint64(data[1])
	decomposition := uint64(data[2])
	bitDecomp := uint64(data[3])
	mappingRow := uint64(data[4])
	mappingColL := make([]uint64, binary.BigEndian.Uint32(data[5:9]))
	mappingColR := make([]uint64, binary.BigEndian.Uint32(data[9:13]))

	rotationkey.bitDecomp = uint64(bitDecomp)

	rotationkey.evakey_rot_col_L = make(map[uint64]*SwitchingKey)
	//rotationkey.evakey_rot_col_R = make(map[uint64][][][2]*ring.Poly)

	pointer := uint64(13)

	for i := 0; i < len(mappingColL); i++ {
		mappingColL[i] = uint64(binary.BigEndian.Uint32(data[pointer : pointer+4]))
		pointer += 4
	}

	for i := 0; i < len(mappingColR); i++ {
		mappingColR[i] = uint64(binary.BigEndian.Uint32(data[pointer : pointer+4]))
		pointer += 4
	}

	var bitLog uint64
	if mappingRow == 1 {

		rotationkey.evakey_rot_row = new(SwitchingKey)
		rotationkey.evakey_rot_row.bitDecomp = bitDecomp
		rotationkey.evakey_rot_row.evakey = make([][][2]*ring.Poly, decomposition)

		for j := uint64(0); j < decomposition; j++ {

			bitLog = uint64(data[pointer])
			pointer += 1

			rotationkey.evakey_rot_row.evakey[j] = make([][2]*ring.Poly, bitLog)

			for x := uint64(0); x < bitLog; x++ {

				rotationkey.evakey_rot_row.evakey[j][x][0] = new(ring.Poly)
				rotationkey.evakey_rot_row.evakey[j][x][0].Coeffs = make([][]uint64, level)
				pointer, _ = ring.DecodeCoeffsNew(pointer, N, level, rotationkey.evakey_rot_row.evakey[j][x][0].Coeffs, data)

				rotationkey.evakey_rot_row.evakey[j][x][1] = new(ring.Poly)
				rotationkey.evakey_rot_row.evakey[j][x][1].Coeffs = make([][]uint64, level)
				pointer, _ = ring.DecodeCoeffsNew(pointer, N, level, rotationkey.evakey_rot_row.evakey[j][x][1].Coeffs, data)
			}
		}
	}

	if len(mappingColL) > 0 {

		rotationkey.evakey_rot_col_L = make(map[uint64]*SwitchingKey)

		for _, i := range mappingColL {

			rotationkey.evakey_rot_col_L[i] = new(SwitchingKey)
			rotationkey.evakey_rot_col_L[i].bitDecomp = bitDecomp
			rotationkey.evakey_rot_col_L[i].evakey = make([][][2]*ring.Poly, decomposition)

			for j := uint64(0); j < decomposition; j++ {

				bitLog = uint64(data[pointer])
				pointer += 1

				rotationkey.evakey_rot_col_L[i].evakey[j] = make([][2]*ring.Poly, bitLog)

				for x := uint64(0); x < bitLog; x++ {

					rotationkey.evakey_rot_col_L[i].evakey[j][x][0] = new(ring.Poly)
					rotationkey.evakey_rot_col_L[i].evakey[j][x][0].Coeffs = make([][]uint64, level)
					pointer, _ = ring.DecodeCoeffsNew(pointer, N, level, rotationkey.evakey_rot_col_L[i].evakey[j][x][0].Coeffs, data)

					rotationkey.evakey_rot_col_L[i].evakey[j][x][1] = new(ring.Poly)
					rotationkey.evakey_rot_col_L[i].evakey[j][x][1].Coeffs = make([][]uint64, level)
					pointer, _ = ring.DecodeCoeffsNew(pointer, N, level, rotationkey.evakey_rot_col_L[i].evakey[j][x][1].Coeffs, data)
				}
			}
		}
	}

	if len(mappingColR) > 0 {

		rotationkey.evakey_rot_col_R = make(map[uint64]*SwitchingKey)

		for _, i := range mappingColR {

			rotationkey.evakey_rot_col_R[i] = new(SwitchingKey)
			rotationkey.evakey_rot_col_R[i].bitDecomp = bitDecomp
			rotationkey.evakey_rot_col_R[i].evakey = make([][][2]*ring.Poly, decomposition)

			for j := uint64(0); j < decomposition; j++ {

				bitLog = uint64(data[pointer])
				pointer += 1

				rotationkey.evakey_rot_col_R[i].evakey[j] = make([][2]*ring.Poly, bitLog)

				for x := uint64(0); x < bitLog; x++ {

					rotationkey.evakey_rot_col_R[i].evakey[j][x][0] = new(ring.Poly)
					rotationkey.evakey_rot_col_R[i].evakey[j][x][0].Coeffs = make([][]uint64, level)
					pointer, _ = ring.DecodeCoeffsNew(pointer, N, level, rotationkey.evakey_rot_col_R[i].evakey[j][x][0].Coeffs, data)

					rotationkey.evakey_rot_col_R[i].evakey[j][x][1] = new(ring.Poly)
					rotationkey.evakey_rot_col_R[i].evakey[j][x][1].Coeffs = make([][]uint64, level)
					pointer, _ = ring.DecodeCoeffsNew(pointer, N, level, rotationkey.evakey_rot_col_R[i].evakey[j][x][1].Coeffs, data)
				}
			}
		}
	}

	return
}
