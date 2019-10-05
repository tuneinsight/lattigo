package bfv

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

// Get returns the switching key backing slice
func (swk *SwitchingKey) Get() [][][2]*ring.Poly {
	return swk.evakey
}

// NewKeyGenerator creates a new KeyGenerator, from which the secret and public keys, as well as the evaluation,
// rotation and switching keys can be generated.
func (bfvcontext *BfvContext) NewKeyGenerator() (keygen *KeyGenerator) {
	keygen = new(KeyGenerator)
	keygen.bfvcontext = bfvcontext
	keygen.context = bfvcontext.contextQ
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

// MarshalBinary encodes a secret-key on a byte slice. The total size in byte is 1 + N/4.
func (sk *SecretKey) MarshalBinary() ([]byte, error) {

	var err error

	N := uint64(len(sk.sk.Coeffs[0]))
	numberModuli := uint64(len(sk.sk.Coeffs))

	if numberModuli > 0xFF {
		return nil, errors.New("error : max degree uint16 overflow")
	}

	data := make([]byte, 2+((N*numberModuli)<<3))

	data[0] = uint8((bits.Len64(uint64(N)) - 1))
	data[1] = uint8(numberModuli)

	pointer := uint64(2)

	if pointer, err = ring.WriteCoeffsTo(pointer, N, numberModuli, sk.sk.Coeffs, data); err != nil {
		return nil, err
	}

	return data, nil
}

// UnMarshalBinary decode a previously marshaled secret-key on the target secret-key.
// The target secret-key must be of the appropriate format, it can be created with the methode NewSecretKeyEmpty().
func (sk *SecretKey) UnmarshalBinary(data []byte) error {

	N := uint64(1 << data[0])
	numberModuli := uint64(data[1])

	pointer := uint64(2)

	if ((uint64(len(data)) - pointer) >> 3) != (N * numberModuli) {
		return errors.New("error : invalid secret-key encoding")
	}

	ring.DecodeCoeffs(pointer, N, numberModuli, sk.sk.Coeffs, data)

	return nil
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

// MarshalBinary encodes a public-key on a byte slice. The total size is 2 + 16 * N * numberModuliQ.
func (pk *PublicKey) MarshalBinary() ([]byte, error) {

	var err error

	N := uint64(len(pk.pk[0].Coeffs[0]))
	numberModuli := uint64(len(pk.pk[0].Coeffs))

	if numberModuli > 0xFF {
		return nil, errors.New("error : max degree uint16 overflow")
	}

	data := make([]byte, 2+((N*numberModuli)<<4))

	data[0] = uint8((bits.Len64(uint64(N)) - 1))
	data[1] = uint8(numberModuli)

	pointer := uint64(2)

	if pointer, err = ring.WriteCoeffsTo(pointer, N, numberModuli, pk.pk[0].Coeffs, data); err != nil {
		return nil, err
	}

	if pointer, err = ring.WriteCoeffsTo(pointer, N, numberModuli, pk.pk[1].Coeffs, data); err != nil {
		return nil, err
	}

	return data, nil
}

// UnMarshalBinary decodes a previously marshaled public-key on the target public-key.
// The target public-key must have the appropriate format and size, it can be created with
// the methode NewPublicKeyEmpty().
func (pk *PublicKey) UnmarshalBinary(data []byte) error {

	N := uint64(1 << data[0])
	numberModuli := uint64(data[1])

	pointer := uint64(2)

	if ((uint64(len(data)) - pointer) >> 4) != (N * numberModuli) {
		return errors.New("error : invalid publickey encoding")
	}

	pointer, _ = ring.DecodeCoeffs(pointer, N, numberModuli, pk.pk[0].Coeffs, data)
	pointer, _ = ring.DecodeCoeffs(pointer, N, numberModuli, pk.pk[1].Coeffs, data)

	return nil
}

// NewKeyPair generates a new secret-key with distribution [1/3, 1/3, 1/3] and a corresponding public-key.
func (keygen *KeyGenerator) NewKeyPair() (sk *SecretKey, pk *PublicKey) {
	sk = keygen.NewSecretKey()
	return sk, keygen.NewPublicKey(sk)
}

// NewRelinKey generates a new evaluation key from the provided secret-key. It will be used to relinearize a ciphertext (encrypted under a public-key generated from the provided secret-key)
// of degree > 1 to a ciphertext of degree 1. Max degree is the maximum degree of the ciphertext allowed to relinearize. Bitdecomp is the power of two binary decomposition of the key.
// A higher bigdecomp will induce smaller keys, faster key-switching, but at the cost of more noise.
func (keygen *KeyGenerator) NewRelinKey(sk *SecretKey, maxDegree, bitDecomp uint64) (newEvakey *EvaluationKey) {

	newEvakey = new(EvaluationKey)
	newEvakey.evakey = make([]*SwitchingKey, maxDegree)
	keygen.polypool.Copy(sk.Get())
	for i := uint64(0); i < maxDegree; i++ {
		keygen.context.MulCoeffsMontgomery(keygen.polypool, sk.Get(), keygen.polypool)
		newEvakey.evakey[i] = newswitchintkey(keygen.bfvcontext, keygen.polypool, sk.Get(), bitDecomp)
	}
	keygen.polypool.Zero()

	return
}

func (bfvcontext *BfvContext) NewRelinKey(maxDegree, bitDecomp uint64) (evakey *EvaluationKey) {
	evakey = new(EvaluationKey)

	if bitDecomp > bfvcontext.maxBit || bitDecomp == 0 {
		bitDecomp = bfvcontext.maxBit
	}

	context := bfvcontext.contextQ

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
func (keygen *KeyGenerator) NewSwitchingKey(sk_input, sk_output *SecretKey, bitDecomp uint64) (newevakey *SwitchingKey) {

	keygen.context.Sub(sk_input.Get(), sk_output.Get(), keygen.polypool)
	newevakey = newswitchintkey(keygen.bfvcontext, keygen.polypool, sk_output.Get(), bitDecomp)
	keygen.polypool.Zero()

	return
}

func (keygen *KeyGenerator) NewSwitchingKeyEmpty(bitDecomp uint64) (evakey *SwitchingKey) {
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
func (keygen *KeyGenerator) NewRotationKeys(sk *SecretKey, bitDecomp uint64, rotLeft []uint64, rotRight []uint64, row bool) (rotKey *RotationKeys) {

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

	return rotKey
}

func (keygen *KeyGenerator) NewRotationKeysEmpty() (rotKey *RotationKeys) {

	rotKey = new(RotationKeys)
	rotKey.bfvcontext = keygen.bfvcontext

	return rotKey

}

// Newrotationkeys generates a new struct of rotationkeys storing the keys of all the left and right powers of two rotations. The provided secret-key must be the secret-key used to generate the public-key under
// which the ciphertexts to rotate are encrypted under. rows is a boolean value indicatig if the keys for the row rotation have to be generated. Bitdecomp is the power of two binary decomposition of the key.
// A higher bigdecomp will induce smaller keys, faster key-switching, but at the cost of more noise.
func (keygen *KeyGenerator) NewRotationKeysPow2(sk *SecretKey, bitDecomp uint64, row bool) (rotKey *RotationKeys) {

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
func genrotkey(keygen *KeyGenerator, sk *ring.Poly, gen, bitDecomp uint64) (switchkey *SwitchingKey) {

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

// MarshalBinary encodes an evaluation key on a byte slice. The total size depends on each modulus size and the bit decomp.
// It will approximately be 5 + maxDegree * numberModuli * ( 1 + 2 * 8 * N * numberModuli * logQi/bitDecomp).
func (evaluationkey *EvaluationKey) MarshalBinary() ([]byte, error) {

	var err error

	N := uint64(len(evaluationkey.evakey[0].evakey[0][0][0].Coeffs[0]))
	numberModuli := uint64(len(evaluationkey.evakey[0].evakey[0][0][0].Coeffs))
	decomposition := numberModuli
	bitDecomp := evaluationkey.evakey[0].bitDecomp

	maxDegree := uint64(len(evaluationkey.evakey))

	if numberModuli > 0xFF {
		return nil, errors.New("cannot marshal evaluationkey -> max number moduli uint16 overflow")
	}

	if decomposition > 0xFF {
		return nil, errors.New("cannot marshal evaluationkey -> max decomposition uint16 overflow")
	}

	if bitDecomp > 0xFF {
		return nil, errors.New("cannot marshal evaluationkey -> max bitDecomp uint16 overflow")
	}

	if maxDegree > 0xFF {
		return nil, errors.New("cannot marshal evaluationkey -> max degree uint16 overflow")
	}

	var dataLen uint64
	dataLen = 5
	for i := uint64(0); i < maxDegree; i++ {
		for j := uint64(0); j < decomposition; j++ {
			dataLen += 1                                                                         //Information about the size of the bitdecomposition
			dataLen += 2 * 8 * N * numberModuli * uint64(len(evaluationkey.evakey[i].evakey[j])) // nb coefficients * 8
		}
	}

	data := make([]byte, dataLen)

	data[0] = uint8(bits.Len64(uint64(N)) - 1)
	data[1] = uint8(numberModuli)
	data[2] = uint8(decomposition)
	data[3] = uint8(bitDecomp)
	data[4] = uint8(maxDegree)

	pointer := uint64(5)

	var bitLog uint8
	for i := uint64(0); i < maxDegree; i++ {
		for j := uint64(0); j < decomposition; j++ {
			bitLog = uint8(len(evaluationkey.evakey[i].evakey[j]))
			data[pointer] = bitLog
			pointer += 1
			for x := uint8(0); x < bitLog; x++ {
				if pointer, err = ring.WriteCoeffsTo(pointer, N, numberModuli, evaluationkey.evakey[i].evakey[j][x][0].Coeffs, data); err != nil {
					return nil, err
				}

				if pointer, err = ring.WriteCoeffsTo(pointer, N, numberModuli, evaluationkey.evakey[i].evakey[j][x][1].Coeffs, data); err != nil {
					return nil, err
				}
			}
		}
	}

	return data, nil
}

// UnMarshalBinary decodes a previously marshaled evaluation-key on the target evaluation-key. The target evaluation-key
// must have the appropriate format and size, it can be created with the methode NewRelinKeyEmpty(uint64, uint64).
func (evaluationkey *EvaluationKey) UnmarshalBinary(data []byte) error {

	N := uint64(1 << data[0])
	numberModuli := uint64(data[1])
	decomposition := uint64(data[2])
	bitDecomp := uint64(data[3])
	maxDegree := uint64(data[4])

	pointer := uint64(5)
	var bitLog uint64
	for i := uint64(0); i < maxDegree; i++ {

		evaluationkey.evakey[i].bitDecomp = bitDecomp

		for j := uint64(0); j < decomposition; j++ {

			bitLog = uint64(data[pointer])
			pointer += 1

			for x := uint64(0); x < bitLog; x++ {

				if uint64(len(evaluationkey.evakey[i].evakey[j][x][0].Coeffs)) != numberModuli {
					return errors.New("cannot unmarshal evaluation-key -> receiver (numberModuli does not match data)")
				}

				if uint64(len(evaluationkey.evakey[i].evakey[j][x][1].Coeffs)) != numberModuli {
					return errors.New("cannot unmarshal evaluation-key -> receiver (numberModuli does not match data)")
				}

				pointer, _ = ring.DecodeCoeffs(pointer, N, numberModuli, evaluationkey.evakey[i].evakey[j][x][0].Coeffs, data)
				pointer, _ = ring.DecodeCoeffs(pointer, N, numberModuli, evaluationkey.evakey[i].evakey[j][x][1].Coeffs, data)
			}
		}
	}

	return nil
}

// MarshalBinary encodes an switching-key on a byte slice. The total size in byte will be approximately 5 + numberModuli * ( 1 + 2 * 8 * N * numberModuli * logQi/bitDecomp).
func (switchkey *SwitchingKey) MarshalBinary() ([]byte, error) {

	var err error

	N := uint64(len(switchkey.evakey[0][0][0].Coeffs[0]))
	level := uint64(len(switchkey.evakey[0][0][0].Coeffs))
	decomposition := level
	bitDecomp := switchkey.bitDecomp

	if level > 0xFF {
		return nil, errors.New("cannot marshal switching-key -> max number modulie uint8 overflow")
	}

	if decomposition > 0xFF {
		return nil, errors.New("cannot marshal switching-key -> max decomposition uint8 overflow")
	}

	if bitDecomp > 0xFF {
		return nil, errors.New("cannot marshal switching-key -> max bitDecomp uint8 overflow")
	}

	var dataLen uint64
	dataLen = 4

	for j := uint64(0); j < decomposition; j++ {
		dataLen += 1                                                                    //Information about the size of the bitdecomposition
		dataLen += 2 * 8 * N * level * decomposition * uint64(len(switchkey.evakey[j])) // nb coefficients * 8
	}

	data := make([]byte, dataLen)

	data[0] = uint8(bits.Len64(uint64(N)) - 1)
	data[1] = uint8(level)
	data[2] = uint8(decomposition)
	data[3] = uint8(bitDecomp)

	pointer := uint64(4)

	var bitLog uint8

	for j := uint64(0); j < decomposition; j++ {
		bitLog = uint8(len(switchkey.evakey[j]))
		data[pointer] = bitLog
		pointer += 1
		for x := uint8(0); x < bitLog; x++ {
			if pointer, err = ring.WriteCoeffsTo(pointer, N, level, switchkey.evakey[j][x][0].Coeffs, data); err != nil {
				return nil, err
			}

			if pointer, err = ring.WriteCoeffsTo(pointer, N, level, switchkey.evakey[j][x][1].Coeffs, data); err != nil {
				return nil, err
			}
		}
	}

	return data, nil
}

// UnMarshalBinary decode a previously marshaled switching-key on the target switching-key.
// The target switching-key must have the appropriate format and size, it can be created with the methode NewSwitchingKeyEmpty(uint64).
func (switchkey *SwitchingKey) UnmarshalBinary(data []byte) error {

	N := uint64(1 << data[0])
	level := uint64(data[1])
	decomposition := uint64(data[2])
	bitDecomp := uint64(data[3])

	pointer := uint64(4)
	var bitLog uint64

	switchkey.bitDecomp = bitDecomp

	for j := uint64(0); j < decomposition; j++ {

		bitLog = uint64(data[pointer])
		pointer += 1

		for x := uint64(0); x < bitLog; x++ {

			if uint64(len(switchkey.evakey[j][x][0].Coeffs)) != decomposition {
				return errors.New("cannot unmarshal switching-key -> receiver (numberModuli does not match data)")
			}

			if uint64(len(switchkey.evakey[j][x][1].Coeffs)) != decomposition {
				return errors.New("cannot unmarshal switching-key -> receiver (numberModuli does not match data)")
			}

			pointer, _ = ring.DecodeCoeffs(pointer, N, level, switchkey.evakey[j][x][0].Coeffs, data)
			pointer, _ = ring.DecodeCoeffs(pointer, N, level, switchkey.evakey[j][x][1].Coeffs, data)
		}
	}

	return nil
}

// MarshalBinary encodes a rotationkeys structure on a byte slice. The total size in byte is approximately
// 5 + 4*(nb left rot + num right rot) + (nb left rot + num right rot + 1 (if rotate row)) * numberModuli * ( 1 + 2 * 8 * N * numberModuli * logQi/bitDecomp).
func (rotationkey *RotationKeys) MarshalBinary() ([]byte, error) {

	var err error

	N := uint64(rotationkey.bfvcontext.n)
	numberModuli := uint64(len(rotationkey.bfvcontext.contextQ.Modulus))
	decomposition := numberModuli
	bitDecomp := rotationkey.bitDecomp
	mappingRow := 0
	mappingColL := []uint64{}
	mappingColR := []uint64{}

	if numberModuli > 0xFF {
		return nil, errors.New("error : max number modulies uint16 overflow")
	}

	if decomposition > 0xFF {
		return nil, errors.New("error : max decomposition uint16 overflow")
	}

	if bitDecomp > 0xFF {
		return nil, errors.New("error : max bitDecomp uint16 overflow")
	}

	var dataLen uint64
	dataLen = 13

	for i := uint64(1); i < rotationkey.bfvcontext.n>>1; i++ {
		if rotationkey.evakey_rot_col_L[i] != nil {

			mappingColL = append(mappingColL, i)

			for j := uint64(0); j < decomposition; j++ {
				dataLen += 1                                                                                                 //Information about the size of the bitdecomposition
				dataLen += 2 * 8 * N * numberModuli * decomposition * uint64(len(rotationkey.evakey_rot_col_L[i].evakey[j])) // nb coefficients * 8
			}
		}

		if rotationkey.evakey_rot_col_L[i] != nil {

			mappingColR = append(mappingColR, i)

			for j := uint64(0); j < decomposition; j++ {
				dataLen += 1                                                                                                 //Information about the size of the bitdecomposition
				dataLen += 2 * 8 * N * numberModuli * decomposition * uint64(len(rotationkey.evakey_rot_col_L[i].evakey[j])) // nb coefficients * 8
			}
		}
	}

	if rotationkey.evakey_rot_row != nil {
		mappingRow = 1
		for j := uint64(0); j < decomposition; j++ {
			dataLen += 1                                                                                            //Information about the size of the bitdecomposition
			dataLen += 2 * 8 * N * numberModuli * decomposition * uint64(len(rotationkey.evakey_rot_row.evakey[j])) // nb coefficients * 8
		}
	}

	dataLen += uint64(len(mappingColL)+len(mappingColR)) << 2 // size needed to encode what rotation are present

	data := make([]byte, dataLen)

	data[0] = uint8(bits.Len64(uint64(N)) - 1)
	data[1] = uint8(numberModuli)
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
				if pointer, err = ring.WriteCoeffsTo(pointer, N, numberModuli, rotationkey.evakey_rot_row.evakey[j][x][0].Coeffs, data); err != nil {
					return nil, err
				}

				if pointer, err = ring.WriteCoeffsTo(pointer, N, numberModuli, rotationkey.evakey_rot_row.evakey[j][x][1].Coeffs, data); err != nil {
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
				if pointer, err = ring.WriteCoeffsTo(pointer, N, numberModuli, rotationkey.evakey_rot_col_L[i].evakey[j][x][0].Coeffs, data); err != nil {
					return nil, err
				}

				if pointer, err = ring.WriteCoeffsTo(pointer, N, numberModuli, rotationkey.evakey_rot_col_L[i].evakey[j][x][1].Coeffs, data); err != nil {
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
				if pointer, err = ring.WriteCoeffsTo(pointer, N, numberModuli, rotationkey.evakey_rot_col_R[i].evakey[j][x][0].Coeffs, data); err != nil {
					return nil, err
				}

				if pointer, err = ring.WriteCoeffsTo(pointer, N, numberModuli, rotationkey.evakey_rot_col_R[i].evakey[j][x][1].Coeffs, data); err != nil {
					return nil, err
				}
			}
		}
	}

	return data, nil
}

// UnMarshalBinary decodes a previously marshaled rotation-keys on the target rotation-keys. In contrary to all
// the other structures, the unmarshaling for rotationkeys only need an empty receiver, as it is not possible to
// create receiver of the correct format and size without knowing all the content of the marshaled rotationkeys. The memory
// will be allocated on the fly.
func (rotationkey *RotationKeys) UnmarshalBinary(data []byte) error {

	N := uint64(1 << data[0])
	numberModuli := uint64(data[1])
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
				rotationkey.evakey_rot_row.evakey[j][x][0].Coeffs = make([][]uint64, numberModuli)
				pointer, _ = ring.DecodeCoeffsNew(pointer, N, numberModuli, rotationkey.evakey_rot_row.evakey[j][x][0].Coeffs, data)

				rotationkey.evakey_rot_row.evakey[j][x][1] = new(ring.Poly)
				rotationkey.evakey_rot_row.evakey[j][x][1].Coeffs = make([][]uint64, numberModuli)
				pointer, _ = ring.DecodeCoeffsNew(pointer, N, numberModuli, rotationkey.evakey_rot_row.evakey[j][x][1].Coeffs, data)
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
					rotationkey.evakey_rot_col_L[i].evakey[j][x][0].Coeffs = make([][]uint64, numberModuli)
					pointer, _ = ring.DecodeCoeffsNew(pointer, N, numberModuli, rotationkey.evakey_rot_col_L[i].evakey[j][x][0].Coeffs, data)

					rotationkey.evakey_rot_col_L[i].evakey[j][x][1] = new(ring.Poly)
					rotationkey.evakey_rot_col_L[i].evakey[j][x][1].Coeffs = make([][]uint64, numberModuli)
					pointer, _ = ring.DecodeCoeffsNew(pointer, N, numberModuli, rotationkey.evakey_rot_col_L[i].evakey[j][x][1].Coeffs, data)
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
					rotationkey.evakey_rot_col_R[i].evakey[j][x][0].Coeffs = make([][]uint64, numberModuli)
					pointer, _ = ring.DecodeCoeffsNew(pointer, N, numberModuli, rotationkey.evakey_rot_col_R[i].evakey[j][x][0].Coeffs, data)

					rotationkey.evakey_rot_col_R[i].evakey[j][x][1] = new(ring.Poly)
					rotationkey.evakey_rot_col_R[i].evakey[j][x][1].Coeffs = make([][]uint64, numberModuli)
					pointer, _ = ring.DecodeCoeffsNew(pointer, N, numberModuli, rotationkey.evakey_rot_col_R[i].evakey[j][x][1].Coeffs, data)
				}
			}
		}
	}

	return nil
}
