package ckks

import (
	"encoding/binary"
	"errors"
	"github.com/ldsec/lattigo/ring"
	"math"
	"math/bits"
)

// KeyGenerator is a structure that stores the elements required to create new keys,
// as well as a small memory pool for intermediate values.
type KeyGenerator struct {
	ckkscontext *Context
	context     *ring.Context
	polypool    *ring.Poly
}

// SecretKey is a structure that stores the secret-key
type SecretKey struct {
	sk *ring.Poly
}

// PublicKey is a structure that stores the public-key
type PublicKey struct {
	pk [2]*ring.Poly
}

// RotationKeys is a structure that stores the switching-keys required during the homomorphic rotations.
type RotationKeys struct {
	ckkscontext       *Context
	bitDecomp         uint64
	evakeyRotColLeft  map[uint64]*SwitchingKey
	evakeyRotColRight map[uint64]*SwitchingKey
	evakeyRotRows     *SwitchingKey
}

// EvaluationKey is a structure that stores the switching-keys required during the relinearization.
type EvaluationKey struct {
	evakey *SwitchingKey
}

// SwitchingKey is a structure that stores the switching-keys required during the key-switching.
type SwitchingKey struct {
	bitDecomp uint64
	evakey    [][][2]*ring.Poly
}

// NewKeyGenerator creates a new keygenerator, from which the secret and public keys, as well as the evaluation,
// rotation and switching keys can be generated.
func (ckkscontext *Context) NewKeyGenerator() (keygen *KeyGenerator) {
	keygen = new(KeyGenerator)
	keygen.ckkscontext = ckkscontext
	keygen.context = ckkscontext.keyscontext
	keygen.polypool = ckkscontext.keyscontext.NewPoly()
	return
}

// NewSecretKey generates a new secret key with the distribution [1/3, 1/3, 1/3].
func (keygen *KeyGenerator) NewSecretKey() (sk *SecretKey) {
	sk, _ = keygen.NewSecretKeyWithDistrib(1.0 / 3)
	return sk
}

// NewSecretKeyWithDistrib generates a new secret key with the distribution [(p-1)/2, p, (p-1)/2].
func (keygen *KeyGenerator) NewSecretKeyWithDistrib(p float64) (sk *SecretKey, err error) {
	sk = new(SecretKey)
	if sk.sk, err = keygen.ckkscontext.ternarySampler.SampleMontgomeryNTTNew(p); err != nil {
		return nil, err
	}
	return sk, nil
}

// NewSecretKeyEmpty creates a new SecretKey struct initialized to zero.
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
	pk.pk[0] = keygen.ckkscontext.gaussianSampler.SampleNTTNew()
	pk.pk[1] = keygen.context.NewUniformPoly()

	keygen.context.MulCoeffsMontgomeryAndAdd(sk.sk, pk.pk[1], pk.pk[0])
	keygen.context.Neg(pk.pk[0], pk.pk[0])

	return pk
}

// NewPublicKeyEmpty creates a new PublicKey struct initialized to zero.
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
	pk = keygen.NewPublicKey(sk)
	return
}

// NewRelinKey generates a new evaluation key that will be used to relinearize the ciphertexts during multiplication.
// Bitdecomposition aims at reducing the added noise at the expense of more storage needed for the keys and more computation
// during the relinearization. However for relinearization this bitdecomp value can be set to maximum as the encrypted value
// are also scaled up during the multiplication.
func (keygen *KeyGenerator) NewRelinKey(sk *SecretKey, bitDecomp uint64) (evakey *EvaluationKey) {

	evakey = new(EvaluationKey)
	keygen.polypool.Copy(sk.Get())
	keygen.context.MulCoeffsMontgomery(keygen.polypool, sk.Get(), keygen.polypool)
	evakey.evakey = keygen.newSwitchingKey(keygen.polypool, sk.Get(), bitDecomp)
	keygen.polypool.Zero()

	return
}

// NewRelinKeyEmpty creates a new EvaluationKey struct initialized to zero.
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

// SetRelinKeys creates a new EvaluationKey struct with the input polynomials as value.
func (keygen *KeyGenerator) SetRelinKeys(rlk [][][2]*ring.Poly, bitDecomp uint64) *EvaluationKey {

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

	return newevakey
}

// NewSwitchingKey generated a new keyswitching key, that will re-encrypt a ciphertext encrypted under the input key to the output key.
// Here bitdecomp plays a role in the added noise if the scale of the input is smaller than the maximum size between the modulies.
func (keygen *KeyGenerator) NewSwitchingKey(skInput, skOutput *SecretKey, bitDecomp uint64) (newevakey *SwitchingKey) {

	keygen.context.Sub(skInput.Get(), skOutput.Get(), keygen.polypool)
	newevakey = keygen.newSwitchingKey(keygen.polypool, skOutput.Get(), bitDecomp)
	keygen.polypool.Zero()

	return
}

// NewSwitchingKeyEmpty creates a new SwitchingKey struct initialized to zero.
func (keygen *KeyGenerator) NewSwitchingKeyEmpty(bitDecomp uint64) (evakey *SwitchingKey) {
	evakey = new(SwitchingKey)

	if bitDecomp > keygen.ckkscontext.maxBit || bitDecomp == 0 {
		bitDecomp = keygen.ckkscontext.maxBit
	}

	context := keygen.ckkscontext.keyscontext

	evakey.bitDecomp = bitDecomp

	// delta_sk = skInput - skOutput = GaloisEnd(skOutput, rotation) - skOutput
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
func (keygen *KeyGenerator) NewRotationKeys(skOutput *SecretKey, bitDecomp uint64, rotLeft []uint64, rotRight []uint64, conjugate bool) (rotKey *RotationKeys) {

	if bitDecomp > keygen.ckkscontext.maxBit || bitDecomp == 0 {
		bitDecomp = keygen.ckkscontext.maxBit
	}

	rotKey = new(RotationKeys)
	rotKey.ckkscontext = keygen.ckkscontext

	if rotLeft != nil {
		rotKey.evakeyRotColLeft = make(map[uint64]*SwitchingKey)
		for _, n := range rotLeft {
			if rotKey.evakeyRotColLeft[n] == nil && n != 0 {
				rotKey.evakeyRotColLeft[n] = keygen.genrotKey(skOutput.Get(), keygen.ckkscontext.galElRotColLeft[n], bitDecomp)
			}
		}
	}

	if rotRight != nil {
		rotKey.evakeyRotColRight = make(map[uint64]*SwitchingKey)
		for _, n := range rotRight {
			if rotKey.evakeyRotColRight[n] == nil && n != 0 {
				rotKey.evakeyRotColRight[n] = keygen.genrotKey(skOutput.Get(), keygen.ckkscontext.galElRotColRight[n], bitDecomp)
			}
		}
	}

	if conjugate {
		rotKey.evakeyRotRows = keygen.genrotKey(skOutput.Get(), keygen.ckkscontext.galElRotRow, bitDecomp)
	}

	return rotKey

}

// NewRotationKeysEmpty creates a new empty RotationKeys struct.
func (keygen *KeyGenerator) NewRotationKeysEmpty() (rotKey *RotationKeys) {

	rotKey = new(RotationKeys)
	rotKey.ckkscontext = keygen.ckkscontext
	return
}

// NewRotationKeysPow2 generates a new rotation key with all the power of two rotation to the left and right, as well as the conjugation
// key if asked. Here bitdecomp plays a role in the added noise if the scale of the input is smaller than the maximum size between the modulies.
func (keygen *KeyGenerator) NewRotationKeysPow2(skOutput *SecretKey, bitDecomp uint64, conjugate bool) (rotKey *RotationKeys) {

	if bitDecomp > keygen.ckkscontext.maxBit || bitDecomp == 0 {
		bitDecomp = keygen.ckkscontext.maxBit
	}

	rotKey = new(RotationKeys)
	rotKey.ckkscontext = keygen.ckkscontext

	rotKey.evakeyRotColLeft = make(map[uint64]*SwitchingKey)
	rotKey.evakeyRotColRight = make(map[uint64]*SwitchingKey)

	for n := uint64(1); n < rotKey.ckkscontext.n>>1; n <<= 1 {
		rotKey.evakeyRotColLeft[n] = keygen.genrotKey(skOutput.Get(), keygen.ckkscontext.galElRotColLeft[n], bitDecomp)
		rotKey.evakeyRotColRight[n] = keygen.genrotKey(skOutput.Get(), keygen.ckkscontext.galElRotColRight[n], bitDecomp)
	}

	if conjugate {
		rotKey.evakeyRotRows = keygen.genrotKey(skOutput.Get(), keygen.ckkscontext.galElRotRow, bitDecomp)
	}

	return
}

func (keygen *KeyGenerator) genrotKey(skOutput *ring.Poly, gen, bitDecomp uint64) (switchingkey *SwitchingKey) {

	ring.PermuteNTT(skOutput, gen, keygen.polypool)
	keygen.context.Sub(keygen.polypool, skOutput, keygen.polypool)
	switchingkey = keygen.newSwitchingKey(keygen.polypool, skOutput, bitDecomp)
	keygen.polypool.Zero()

	return
}

func (keygen *KeyGenerator) newSwitchingKey(skIn, skOut *ring.Poly, bitDecomp uint64) (switchingkey *SwitchingKey) {

	if bitDecomp > keygen.ckkscontext.maxBit || bitDecomp == 0 {
		bitDecomp = keygen.ckkscontext.maxBit
	}

	switchingkey = new(SwitchingKey)

	context := keygen.ckkscontext.keyscontext

	switchingkey.bitDecomp = uint64(bitDecomp)

	mredParams := context.GetMredParams()

	// delta_sk = sk_input - skOutput = GaloisEnd(skOutput, rotation) - skOutput
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

			// e + skIn * (qiBarre*qiStar) * 2^w
			// (qiBarre*qiStar)%qi = 1 mod qi, else 0
			for w := uint64(0); w < context.N; w++ {
				switchingkey.evakey[i][j][0].Coeffs[i][w] += ring.PowerOf2(skIn.Coeffs[i][w], bitDecomp*j, qi, mredParams[i])
			}

			// skIn * (qiBarre*qiStar) * 2^w - a*sk + e
			context.MulCoeffsMontgomeryAndSub(switchingkey.evakey[i][j][1], skOut, switchingkey.evakey[i][j][0])

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
// The target secret-key must be of the appropriate format, it can be created with the method NewSecretKeyEmpty().
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
// the method NewPublicKeyEmpty().
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
		dataLen++                                                                                   //Information about the size of the bitdecomposition
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
		pointer++
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
// must have the appropriate format and size, it can be created with the method NewRelinKeyEmpty(uint64, uint64).
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
		pointer++

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
		dataLen++                                                                          //Information about the size of the bitdecomposition
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
		pointer++
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
// The target switching-key must have the appropriate format and size, it can be created with the method NewSwitchingKeyEmpty(uint64).
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
		pointer++

		for x := uint64(0); x < bitLog; x++ {
			pointer, _ = ring.DecodeCoeffs(pointer, N, level, switchingkey.evakey[j][x][0].Coeffs, data)
			pointer, _ = ring.DecodeCoeffs(pointer, N, level, switchingkey.evakey[j][x][1].Coeffs, data)
		}
	}

	return
}

// MarshalBinary encodes a rotationkeys structure on a byte slice. The total size in byte is approximately
// 5 + 4*(nb left rot + num right rot) + (nb left rot + num right rot + 1 (if rotate row)) * (level + 1) * ( 1 + 2 * 8 * N * (level + 1) * logQi/bitDecomp).
func (rotationkey *RotationKeys) MarshalBinary() (data []byte, err error) {

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
		if rotationkey.evakeyRotColLeft[i] != nil {

			mappingColL = append(mappingColL, i)

			for j := uint64(0); j < decomposition; j++ {
				dataLen++                                                                                             //Information about the size of the bitdecomposition
				dataLen += 2 * 8 * N * level * decomposition * uint64(len(rotationkey.evakeyRotColLeft[i].evakey[j])) // nb coefficients * 8
			}
		}

		if rotationkey.evakeyRotColLeft[i] != nil {

			mappingColR = append(mappingColR, i)

			for j := uint64(0); j < decomposition; j++ {
				dataLen++                                                                                             //Information about the size of the bitdecomposition
				dataLen += 2 * 8 * N * level * decomposition * uint64(len(rotationkey.evakeyRotColLeft[i].evakey[j])) // nb coefficients * 8
			}
		}
	}

	if rotationkey.evakeyRotRows != nil {
		mappingRow = 1
		for j := uint64(0); j < decomposition; j++ {
			dataLen++                                                                                       //Information about the size of the bitdecomposition
			dataLen += 2 * 8 * N * level * decomposition * uint64(len(rotationkey.evakeyRotRows.evakey[j])) // nb coefficients * 8
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
			bitLog = uint8(len(rotationkey.evakeyRotRows.evakey[j]))
			data[pointer] = bitLog
			pointer++
			for x := uint8(0); x < bitLog; x++ {
				if pointer, err = ring.WriteCoeffsTo(pointer, N, level, rotationkey.evakeyRotRows.evakey[j][x][0].Coeffs, data); err != nil {
					return nil, err
				}

				if pointer, err = ring.WriteCoeffsTo(pointer, N, level, rotationkey.evakeyRotRows.evakey[j][x][1].Coeffs, data); err != nil {
					return nil, err
				}
			}
		}
	}

	for _, i := range mappingColL {
		for j := uint64(0); j < decomposition; j++ {
			bitLog = uint8(len(rotationkey.evakeyRotColLeft[i].evakey[j]))
			data[pointer] = bitLog
			pointer++
			for x := uint8(0); x < bitLog; x++ {
				if pointer, err = ring.WriteCoeffsTo(pointer, N, level, rotationkey.evakeyRotColLeft[i].evakey[j][x][0].Coeffs, data); err != nil {
					return nil, err
				}

				if pointer, err = ring.WriteCoeffsTo(pointer, N, level, rotationkey.evakeyRotColLeft[i].evakey[j][x][1].Coeffs, data); err != nil {
					return nil, err
				}
			}
		}
	}

	for _, i := range mappingColR {
		for j := uint64(0); j < decomposition; j++ {
			bitLog = uint8(len(rotationkey.evakeyRotColRight[i].evakey[j]))
			data[pointer] = bitLog
			pointer++
			for x := uint8(0); x < bitLog; x++ {
				if pointer, err = ring.WriteCoeffsTo(pointer, N, level, rotationkey.evakeyRotColRight[i].evakey[j][x][0].Coeffs, data); err != nil {
					return nil, err
				}

				if pointer, err = ring.WriteCoeffsTo(pointer, N, level, rotationkey.evakeyRotColRight[i].evakey[j][x][1].Coeffs, data); err != nil {
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
func (rotationkey *RotationKeys) UnMarshalBinary(data []byte) (err error) {

	N := uint64(1 << data[0])
	level := uint64(data[1])
	decomposition := uint64(data[2])
	bitDecomp := uint64(data[3])
	mappingRow := uint64(data[4])
	mappingColL := make([]uint64, binary.BigEndian.Uint32(data[5:9]))
	mappingColR := make([]uint64, binary.BigEndian.Uint32(data[9:13]))

	rotationkey.bitDecomp = uint64(bitDecomp)

	rotationkey.evakeyRotColLeft = make(map[uint64]*SwitchingKey)
	//rotationkey.evakeyRotColRight = make(map[uint64][][][2]*ring.Poly)

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

		rotationkey.evakeyRotRows = new(SwitchingKey)
		rotationkey.evakeyRotRows.bitDecomp = bitDecomp
		rotationkey.evakeyRotRows.evakey = make([][][2]*ring.Poly, decomposition)

		for j := uint64(0); j < decomposition; j++ {

			bitLog = uint64(data[pointer])
			pointer++

			rotationkey.evakeyRotRows.evakey[j] = make([][2]*ring.Poly, bitLog)

			for x := uint64(0); x < bitLog; x++ {

				rotationkey.evakeyRotRows.evakey[j][x][0] = new(ring.Poly)
				rotationkey.evakeyRotRows.evakey[j][x][0].Coeffs = make([][]uint64, level)
				pointer, _ = ring.DecodeCoeffsNew(pointer, N, level, rotationkey.evakeyRotRows.evakey[j][x][0].Coeffs, data)

				rotationkey.evakeyRotRows.evakey[j][x][1] = new(ring.Poly)
				rotationkey.evakeyRotRows.evakey[j][x][1].Coeffs = make([][]uint64, level)
				pointer, _ = ring.DecodeCoeffsNew(pointer, N, level, rotationkey.evakeyRotRows.evakey[j][x][1].Coeffs, data)
			}
		}
	}

	if len(mappingColL) > 0 {

		rotationkey.evakeyRotColLeft = make(map[uint64]*SwitchingKey)

		for _, i := range mappingColL {

			rotationkey.evakeyRotColLeft[i] = new(SwitchingKey)
			rotationkey.evakeyRotColLeft[i].bitDecomp = bitDecomp
			rotationkey.evakeyRotColLeft[i].evakey = make([][][2]*ring.Poly, decomposition)

			for j := uint64(0); j < decomposition; j++ {

				bitLog = uint64(data[pointer])
				pointer++

				rotationkey.evakeyRotColLeft[i].evakey[j] = make([][2]*ring.Poly, bitLog)

				for x := uint64(0); x < bitLog; x++ {

					rotationkey.evakeyRotColLeft[i].evakey[j][x][0] = new(ring.Poly)
					rotationkey.evakeyRotColLeft[i].evakey[j][x][0].Coeffs = make([][]uint64, level)
					pointer, _ = ring.DecodeCoeffsNew(pointer, N, level, rotationkey.evakeyRotColLeft[i].evakey[j][x][0].Coeffs, data)

					rotationkey.evakeyRotColLeft[i].evakey[j][x][1] = new(ring.Poly)
					rotationkey.evakeyRotColLeft[i].evakey[j][x][1].Coeffs = make([][]uint64, level)
					pointer, _ = ring.DecodeCoeffsNew(pointer, N, level, rotationkey.evakeyRotColLeft[i].evakey[j][x][1].Coeffs, data)
				}
			}
		}
	}

	if len(mappingColR) > 0 {

		rotationkey.evakeyRotColRight = make(map[uint64]*SwitchingKey)

		for _, i := range mappingColR {

			rotationkey.evakeyRotColRight[i] = new(SwitchingKey)
			rotationkey.evakeyRotColRight[i].bitDecomp = bitDecomp
			rotationkey.evakeyRotColRight[i].evakey = make([][][2]*ring.Poly, decomposition)

			for j := uint64(0); j < decomposition; j++ {

				bitLog = uint64(data[pointer])
				pointer++

				rotationkey.evakeyRotColRight[i].evakey[j] = make([][2]*ring.Poly, bitLog)

				for x := uint64(0); x < bitLog; x++ {

					rotationkey.evakeyRotColRight[i].evakey[j][x][0] = new(ring.Poly)
					rotationkey.evakeyRotColRight[i].evakey[j][x][0].Coeffs = make([][]uint64, level)
					pointer, _ = ring.DecodeCoeffsNew(pointer, N, level, rotationkey.evakeyRotColRight[i].evakey[j][x][0].Coeffs, data)

					rotationkey.evakeyRotColRight[i].evakey[j][x][1] = new(ring.Poly)
					rotationkey.evakeyRotColRight[i].evakey[j][x][1].Coeffs = make([][]uint64, level)
					pointer, _ = ring.DecodeCoeffsNew(pointer, N, level, rotationkey.evakeyRotColRight[i].evakey[j][x][1].Coeffs, data)
				}
			}
		}
	}

	return
}
