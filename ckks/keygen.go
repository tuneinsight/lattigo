package ckks

import (
	"encoding/binary"
	"github.com/ldsec/lattigo/ring"
)

// Keygenerator is a structure that stores the elements required to create new keys,
// as well as a small memory pool for intermediate values.
type KeyGenerator struct {
	ckkscontext *Context
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
	permuteNTTRightIndex     map[uint64][]uint64
	permuteNTTLeftIndex      map[uint64][]uint64
	permuteNTTConjugateIndex []uint64

	evakeyRotColLeft  map[uint64]*SwitchingKey
	evakeyRotColRight map[uint64]*SwitchingKey
	evakeyConjugate   *SwitchingKey
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
func (ckkscontext *Context) NewKeyGenerator() (keygen *KeyGenerator) {
	keygen = new(KeyGenerator)
	keygen.ckkscontext = ckkscontext
	keygen.context = ckkscontext.contextKeys
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
	sk.sk = keygen.ckkscontext.contextKeys.SampleTernaryMontgomeryNTTNew(p)
	return sk
}

func (keygen *KeyGenerator) NewSecretKeySparse(hw uint64) (sk *SecretKey) {
	sk = new(SecretKey)
	sk.sk = keygen.ckkscontext.contextKeys.SampleTernarySparseMontgomeryNTTNew(hw)
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

// MarshalBinary encodes a secret-key on a byte slice. The total size in byte is 1 + N/4.
func (sk *SecretKey) MarshalBinary() (data []byte, err error) {

	data = make([]byte, sk.GetDataLen(true))

	if _, err = sk.sk.WriteTo(data); err != nil {
		return nil, err
	}

	return data, nil
}

func (sk *SecretKey) GetDataLen(WithMetadata bool) (dataLen uint64) {
	return sk.sk.GetDataLen(WithMetadata)
}

// UnMarshalBinary decode a previously marshaled secret-key on the target secret-key.
// The target secret-key must be of the appropriate format, it can be created with the methode NewSecretKeyEmpty().
func (sk *SecretKey) UnmarshalBinary(data []byte) (err error) {

	sk.sk = new(ring.Poly)

	if _, err = sk.sk.DecodePolyNew(data); err != nil {
		return err
	}

	return nil
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

// MarshalBinary encodes a public-key on a byte slice. The total size is 2 + 16 * N * numberModuliQ.
func (pk *PublicKey) MarshalBinary() (data []byte, err error) {

	dataLen := pk.GetDataLen(true)

	data = make([]byte, dataLen)

	var pointer, inc uint64

	if inc, err = pk.pk[0].WriteTo(data[pointer:]); err != nil {
		return nil, err
	}

	if _, err = pk.pk[1].WriteTo(data[pointer+inc:]); err != nil {
		return nil, err
	}

	return data, err

}

func (pk *PublicKey) GetDataLen(WithMetadata bool) (dataLen uint64) {

	for _, el := range pk.pk {
		dataLen += el.GetDataLen(WithMetadata)
	}

	return
}

// UnMarshalBinary decodes a previously marshaled public-key on the target public-key.
func (pk *PublicKey) UnmarshalBinary(data []byte) (err error) {

	var pointer, inc uint64

	pk.pk[0] = new(ring.Poly)
	pk.pk[1] = new(ring.Poly)

	if inc, err = pk.pk[0].DecodePolyNew(data[pointer:]); err != nil {
		return err
	}

	if _, err = pk.pk[1].DecodePolyNew(data[pointer+inc:]); err != nil {
		return err
	}

	return nil
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

func (ckkscontext *Context) NewRelinKeyEmpty() (evakey *EvaluationKey) {
	evakey = new(EvaluationKey)
	evakey.evakey = new(SwitchingKey)

	// delta_sk = sk_input - sk_output = GaloisEnd(sk_output, rotation) - sk_output
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

func (evk *EvaluationKey) Set(rlk [][2]*ring.Poly) {

	evk.evakey = new(SwitchingKey)
	evk.evakey.evakey = make([][2]*ring.Poly, len(rlk))
	for j := range rlk {
		evk.evakey.evakey[j][0] = rlk[j][0].CopyNew()
		evk.evakey.evakey[j][1] = rlk[j][1].CopyNew()
	}
}

// MarshalBinary encodes an evaluation key on a byte slice.
func (evaluationkey *EvaluationKey) MarshalBinary() (data []byte, err error) {

	var pointer uint64

	dataLen := evaluationkey.evakey.GetDataLen(true)

	data = make([]byte, dataLen)

	if _, err = evaluationkey.evakey.encode(pointer, data); err != nil {
		return nil, err
	}

	return data, nil
}

// UnMarshalBinary decodes a previously marshaled evaluation-key on the target evaluation-key.
func (evaluationkey *EvaluationKey) UnmarshalBinary(data []byte) (err error) {
	evaluationkey.evakey = new(SwitchingKey)
	if _, err = evaluationkey.evakey.decode(data); err != nil {
		return err
	}
	return nil
}

func (evaluationkey *EvaluationKey) GetDataLen(WithMetadata bool) (dataLen uint64) {
	return evaluationkey.evakey.GetDataLen(WithMetadata)
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
	evakey.evakey = make([][2]*ring.Poly, keygen.ckkscontext.beta)

	for i := uint64(0); i < keygen.ckkscontext.beta; i++ {
		evakey.evakey[i][0] = keygen.context.NewPoly()
		evakey.evakey[i][1] = keygen.context.NewPoly()
	}

	return
}

func (keygen *KeyGenerator) newSwitchingKey(sk_in, sk_out *ring.Poly) (switchingkey *SwitchingKey) {

	switchingkey = new(SwitchingKey)

	context := keygen.ckkscontext.contextKeys

	// Computes P * sk_in
	for _, pj := range keygen.ckkscontext.specialprimes {
		context.MulScalar(sk_in, pj, sk_in)
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
			if index >= keygen.ckkscontext.levels-1 {
				break
			}
		}

		// (sk_in * P) * (q_star * q_tild) - a * sk_out + e mod QP
		context.MulCoeffsMontgomeryAndSub(switchingkey.evakey[i][1], sk_out, switchingkey.evakey[i][0])
	}

	return
}

// MarshalBinary encodes an switching-key on a byte slice.
func (switchkey *SwitchingKey) MarshalBinary() (data []byte, err error) {

	data = make([]byte, switchkey.GetDataLen(true))

	if _, err = switchkey.encode(0, data); err != nil {
		return nil, err
	}

	return data, nil
}

// UnMarshalBinary decode a previously marshaled switching-key on the target switching-key.
// The target switching-key must have the appropriate format and size, it can be created with the methode NewSwitchingKeyEmpty(uint64).
func (switchkey *SwitchingKey) UnmarshalBinary(data []byte) (err error) {

	if _, err = switchkey.decode(data); err != nil {
		return err
	}

	return nil
}

func (switchkey *SwitchingKey) GetDataLen(WithMetadata bool) (dataLen uint64) {

	if WithMetadata {
		dataLen++
	}

	for j := uint64(0); j < uint64(len(switchkey.evakey)); j++ {
		dataLen += switchkey.evakey[j][0].GetDataLen(WithMetadata)
		dataLen += switchkey.evakey[j][1].GetDataLen(WithMetadata)
	}

	return
}

func (switchkey *SwitchingKey) encode(pointer uint64, data []byte) (uint64, error) {

	var err error

	var inc uint64

	data[pointer] = uint8(len(switchkey.evakey))

	pointer++

	for j := uint64(0); j < uint64(len(switchkey.evakey)); j++ {

		if inc, err = switchkey.evakey[j][0].WriteTo(data[pointer:]); err != nil {
			return pointer, err
		}

		pointer += inc

		if inc, err = switchkey.evakey[j][1].WriteTo(data[pointer:]); err != nil {
			return pointer, err
		}

		pointer += inc
	}

	return pointer, nil
}

// UnMarshalBinary decode a previously marshaled switching-key on the target switching-key.
// The target switching-key must have the appropriate format and size, it can be created with the methode NewSwitchingKeyEmpty(uint64).
func (switchkey *SwitchingKey) decode(data []byte) (pointer uint64, err error) {

	decomposition := uint64(data[0])

	pointer = uint64(1)

	switchkey.evakey = make([][2]*ring.Poly, decomposition)

	var inc uint64

	for j := uint64(0); j < decomposition; j++ {

		switchkey.evakey[j][0] = new(ring.Poly)
		if inc, err = switchkey.evakey[j][0].DecodePolyNew(data[pointer:]); err != nil {
			panic(err)
		}
		pointer += inc

		switchkey.evakey[j][1] = new(ring.Poly)
		if inc, err = switchkey.evakey[j][1].DecodePolyNew(data[pointer:]); err != nil {
			panic(err)
		}
		pointer += inc

	}

	return pointer, nil
}

// NewRotationKeys generates a new instance of rotationkeys, with the provided rotation to the left, right and conjugation if asked.
// Here bitdecomp plays a role in the added noise if the scale of the input is smaller than the maximum size between the modulies.
func (ckkscontext *Context) NewRotationKeys() (rotKey *RotationKeys) {
	rotKey = new(RotationKeys)
	return
}

// Newrotationkeys generates a new struct of rotationkeys storing the keys for the specified rotations. The provided secret-key must be the secret-key used to generate the public-key under
// which the ciphertexts to rotate are encrypted under. Bitdecomp is the power of two binary decomposition of the key. A higher bigdecomp will induce smaller keys, faster key-switching,
// but at the cost of more noise. rotLeft and rotRight must be a slice of uint64 rotations, row is a boolean value indicating if the key for the row rotation must be generated.
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

// NewRotationkeysPow2 generates a new rotation key with all the power of two rotation to the left and right, as well as the conjugation
// key if asked. Here bitdecomp plays a role in the added noise if the scale of the input is smaller than the maximum size between the modulies.
func (keygen *KeyGenerator) NewRotationKeysPow2(sk_output *SecretKey) (rotKey *RotationKeys) {

	rotKey = new(RotationKeys)

	rotKey.evakeyRotColLeft = make(map[uint64]*SwitchingKey)
	rotKey.evakeyRotColRight = make(map[uint64]*SwitchingKey)

	rotKey.permuteNTTLeftIndex = make(map[uint64][]uint64)
	rotKey.permuteNTTRightIndex = make(map[uint64][]uint64)

	for n := uint64(1); n < keygen.ckkscontext.n>>1; n <<= 1 {

		rotKey.permuteNTTLeftIndex[n] = ring.PermuteNTTIndex(keygen.ckkscontext.galElRotColLeft[n], 1<<keygen.ckkscontext.logN)
		rotKey.permuteNTTRightIndex[n] = ring.PermuteNTTIndex(keygen.ckkscontext.galElRotColRight[n], 1<<keygen.ckkscontext.logN)

		rotKey.evakeyRotColLeft[n] = keygen.genrotKey(sk_output.Get(), keygen.ckkscontext.galElRotColLeft[n])
		rotKey.evakeyRotColRight[n] = keygen.genrotKey(sk_output.Get(), keygen.ckkscontext.galElRotColRight[n])
	}

	rotKey.permuteNTTConjugateIndex = ring.PermuteNTTIndex(keygen.ckkscontext.galElRotRow, 1<<keygen.ckkscontext.logN)
	rotKey.evakeyConjugate = keygen.genrotKey(sk_output.Get(), keygen.ckkscontext.galElRotRow)
	return
}

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

func (keygen *KeyGenerator) genrotKey(sk_output *ring.Poly, gen uint64) (switchingkey *SwitchingKey) {

	ring.PermuteNTT(sk_output, gen, keygen.polypool)
	keygen.context.Sub(keygen.polypool, sk_output, keygen.polypool)
	switchingkey = keygen.newSwitchingKey(keygen.polypool, sk_output)
	keygen.polypool.Zero()

	return
}

// MarshalBinary encodes a rotationkeys structure on a byte slice.
func (rotationkey *RotationKeys) MarshalBinary() ([]byte, error) {

	mappingColL := []uint64{}
	mappingColR := []uint64{}

	var dataLen uint64
	dataLen = 0

	for i := range rotationkey.evakeyRotColLeft {

		mappingColL = append(mappingColL, i)

		dataLen += 4 + rotationkey.evakeyRotColLeft[i].GetDataLen(true)
	}

	for i := range rotationkey.evakeyRotColRight {

		mappingColR = append(mappingColR, i)

		dataLen += 4 + rotationkey.evakeyRotColRight[i].GetDataLen(true)
	}

	if rotationkey.evakeyConjugate != nil {

		dataLen += 4 + rotationkey.evakeyConjugate.GetDataLen(true)
	}

	data := make([]byte, dataLen)

	pointer := uint64(0)

	for _, i := range mappingColL {

		binary.BigEndian.PutUint32(data[pointer:pointer+4], uint32(i))
		data[pointer] = uint8(RotationLeft)
		pointer += 4

		pointer, _ = rotationkey.evakeyRotColLeft[i].encode(pointer, data)
	}

	for _, i := range mappingColR {

		binary.BigEndian.PutUint32(data[pointer:pointer+4], uint32(i))
		data[pointer] = uint8(RotationRight)
		pointer += 4

		pointer, _ = rotationkey.evakeyRotColRight[i].encode(pointer, data)
	}

	if rotationkey.evakeyConjugate != nil {

		data[pointer] = uint8(Conjugate)
		pointer += 4

		_, _ = rotationkey.evakeyConjugate.encode(pointer, data)
	}

	return data, nil
}

func (rotationkey *RotationKeys) UnmarshalBinary(data []byte) (err error) {

	var rotationType int
	var rotationNumber uint64

	pointer := uint64(0)
	var inc uint64

	dataLen := len(data)

	for dataLen > 0 {

		rotationType = int(data[pointer])
		rotationNumber = (uint64(data[pointer+1]) << 16) | (uint64(data[pointer+2]) << 8) | (uint64(data[pointer+3]))

		pointer += 4

		if rotationType == RotationLeft {

			if rotationkey.evakeyRotColLeft == nil {
				rotationkey.evakeyRotColLeft = make(map[uint64]*SwitchingKey)
			}

			rotationkey.evakeyRotColLeft[rotationNumber] = new(SwitchingKey)
			if inc, err = rotationkey.evakeyRotColLeft[rotationNumber].decode(data[pointer:]); err != nil {
				return err
			}

		} else if rotationType == RotationRight {

			if rotationkey.evakeyRotColRight == nil {
				rotationkey.evakeyRotColRight = make(map[uint64]*SwitchingKey)
			}

			rotationkey.evakeyRotColRight[rotationNumber] = new(SwitchingKey)
			if inc, err = rotationkey.evakeyRotColRight[rotationNumber].decode(data[pointer:]); err != nil {
				return err
			}

		} else if rotationType == Conjugate {

			rotationkey.evakeyConjugate = new(SwitchingKey)
			if inc, err = rotationkey.evakeyConjugate.decode(data[pointer:]); err != nil {
				return err
			}

		} else {

			return err
		}

		pointer += inc

		dataLen -= int(4 + inc)
	}

	return nil
}
