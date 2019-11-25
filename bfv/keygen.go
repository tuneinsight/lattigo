package bfv

import (
	"encoding/binary"
	"github.com/ldsec/lattigo/ring"
)

// KeyGenerator is a structure that stores the elements required to create new keys,
// as well as a small memory pool for intermediate values.
type KeyGenerator struct {
	context  *Context
	polypool *ring.Poly
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
	evakeyRotColLeft  map[uint64]*SwitchingKey
	evakeyRotColRight map[uint64]*SwitchingKey
	evakeyRotRow      *SwitchingKey
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
func (context *Context) NewKeyGenerator() (keygen *KeyGenerator) {
	keygen = new(KeyGenerator)
	keygen.context = context
	keygen.polypool = keygen.context.ContextKeys().NewPoly()
	return
}

// Newsecretkey creates a new SecretKey with the distribution [1/3, 1/3, 1/3]
func (keygen *KeyGenerator) NewSecretKey() (sk *SecretKey) {
	return keygen.NewSecretkeyWithDistrib(1.0 / 3)
}

// Newsecretkey creates a new SecretKey with the distribution [(1-p)/2, p, (1-p)/2]
func (keygen *KeyGenerator) NewSecretkeyWithDistrib(p float64) (sk *SecretKey) {
	sk = new(SecretKey)
	sk.sk = keygen.context.contextKeys.SampleTernaryMontgomeryNTTNew(p)
	return sk
}

// NewSecretKeyEmpty creates a new SecretKey with all coeffcients set to zero, ready to received a marshaled SecretKey.
func (keygen *KeyGenerator) NewSecretKeyEmpty() *SecretKey {
	sk := new(SecretKey)
	sk.sk = keygen.context.ContextKeys().NewPoly()
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

// Newpublickey generates a new publickkey from the provided secret-key
func (keygen *KeyGenerator) NewPublicKey(sk *SecretKey) (pk *PublicKey) {

	pk = new(PublicKey)

	ringContext := keygen.context.ContextKeys()

	//pk[0] = [-(a*s + e)]
	//pk[1] = [a]
	pk.pk[0] = keygen.context.gaussianSampler.SampleNTTNew()
	pk.pk[1] = ringContext.NewUniformPoly()

	ringContext.MulCoeffsMontgomeryAndAdd(sk.sk, pk.pk[1], pk.pk[0])
	ringContext.Neg(pk.pk[0], pk.pk[0])

	return pk
}

func (context *Context) NewPublicKey() (pk *PublicKey) {
	pk = new(PublicKey)
	pk.pk[0] = context.contextKeys.NewPoly()
	pk.pk[1] = context.contextKeys.NewPoly()
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

	ringContext := keygen.context.contextKeys

	for _, pj := range keygen.context.specialprimes {
		ringContext.MulScalar(keygen.polypool, pj, keygen.polypool)
	}

	for i := uint64(0); i < maxDegree; i++ {
		ringContext.MulCoeffsMontgomery(keygen.polypool, sk.Get(), keygen.polypool)
		newEvakey.evakey[i] = newswitchintkey(keygen.context, keygen.polypool, sk.Get())
	}

	keygen.polypool.Zero()

	return newEvakey
}

func (context *Context) NewRelinKeyEmpty(maxDegree uint64) (evakey *EvaluationKey) {

	evakey = new(EvaluationKey)

	ringContext := context.contextKeys

	// delta_sk = sk_input - sk_output = GaloisEnd(sk_output, rotation) - sk_output

	evakey.evakey = make([]*SwitchingKey, maxDegree)

	for w := uint64(0); w < maxDegree; w++ {

		evakey.evakey[w] = new(SwitchingKey)
		evakey.evakey[w].evakey = make([][2]*ring.Poly, context.beta)

		for i := uint64(0); i < context.beta; i++ {
			evakey.evakey[w].evakey[i][0] = ringContext.NewPoly()
			evakey.evakey[w].evakey[i][1] = ringContext.NewPoly()

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

// MarshalBinary encodes an evaluation key on a byte slice.
func (evaluationkey *EvaluationKey) MarshalBinary() (data []byte, err error) {

	var pointer uint64

	dataLen := evaluationkey.GetDataLen(true)

	data = make([]byte, dataLen)

	data[0] = uint8(len(evaluationkey.evakey))

	pointer++

	for _, evakey := range evaluationkey.evakey {

		if pointer, err = evakey.encode(pointer, data); err != nil {
			return nil, err
		}
	}

	return data, nil
}

// UnMarshalBinary decodes a previously marshaled evaluation-key on the target evaluation-key.
func (evaluationkey *EvaluationKey) UnmarshalBinary(data []byte) (err error) {

	deg := uint64(data[0])

	evaluationkey.evakey = make([]*SwitchingKey, deg)

	pointer := uint64(1)
	var inc uint64
	for i := uint64(0); i < deg; i++ {
		evaluationkey.evakey[i] = new(SwitchingKey)
		if inc, err = evaluationkey.evakey[i].decode(data[pointer:]); err != nil {
			return err
		}
		pointer += inc
	}

	return nil
}

func (evaluationkey *EvaluationKey) GetDataLen(WithMetadata bool) (dataLen uint64) {

	if WithMetadata {
		dataLen++
	}

	for _, evakey := range evaluationkey.evakey {
		dataLen += evakey.GetDataLen(WithMetadata)
	}

	return
}

// Newswitchintkey generates a new key-switching key, that will allow to re-encrypt under the output-key a ciphertext encrypted under the input-key. Bitdecomp
// is the power of two binary decomposition of the key. A higher bigdecomp will induce smaller keys, faster key-switching, but at the cost of more noise.
func (keygen *KeyGenerator) NewSwitchingKey(sk_input, sk_output *SecretKey) (newevakey *SwitchingKey) {

	ringContext := keygen.context.contextKeys

	ringContext.Sub(sk_input.Get(), sk_output.Get(), keygen.polypool)

	for _, pj := range keygen.context.specialprimes {
		ringContext.MulScalar(keygen.polypool, pj, keygen.polypool)
	}

	newevakey = newswitchintkey(keygen.context, keygen.polypool, sk_output.Get())
	keygen.polypool.Zero()

	return
}

func (keygen *KeyGenerator) NewSwitchingKeyEmpty() (evakey *SwitchingKey) {
	evakey = new(SwitchingKey)

	ringContext := keygen.context.contextKeys

	// delta_sk = sk_input - sk_output = GaloisEnd(sk_output, rotation) - sk_output
	evakey.evakey = make([][2]*ring.Poly, keygen.context.beta)

	for i := uint64(0); i < keygen.context.beta; i++ {
		evakey.evakey[i][0] = ringContext.NewPoly()
		evakey.evakey[i][1] = ringContext.NewPoly()
	}

	return
}

// newswitchintkey is a generic methode to generate key-switching keys used in the evaluation, key-switching and rotation-keys generation.
func newswitchintkey(context *Context, sk_in, sk_out *ring.Poly) (switchkey *SwitchingKey) {

	switchkey = new(SwitchingKey)

	ringContext := context.contextKeys

	var index uint64

	// delta_sk = sk_input - sk_output = GaloisEnd(sk_output, rotation) - sk_output

	switchkey.evakey = make([][2]*ring.Poly, context.beta)

	for i := uint64(0); i < context.beta; i++ {

		// e
		switchkey.evakey[i][0] = context.gaussianSampler.SampleNTTNew()
		ringContext.MForm(switchkey.evakey[i][0], switchkey.evakey[i][0])
		// a
		switchkey.evakey[i][1] = ringContext.NewUniformPoly()

		// e + sk_in * (qiBarre*qiStar) * 2^w
		// (qiBarre*qiStar)%qi = 1, else 0

		for j := uint64(0); j < context.alpha; j++ {

			index = i*context.alpha + j

			qi := ringContext.Modulus[index]
			p0tmp := sk_in.Coeffs[index]
			p1tmp := switchkey.evakey[i][0].Coeffs[index]

			for w := uint64(0); w < ringContext.N; w++ {
				p1tmp[w] = ring.CRed(p1tmp[w]+p0tmp[w], qi)
			}

			// Handles the case where nb pj does not divides nb qi
			if index >= uint64(len(ringContext.Modulus)-1) {
				break
			}

		}

		// sk_in * (qiBarre*qiStar) * 2^w - a*sk + e
		ringContext.MulCoeffsMontgomeryAndSub(switchkey.evakey[i][1], sk_out, switchkey.evakey[i][0])
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

		if inc, err = switchkey.evakey[j][0].WriteTo(data[pointer : pointer+switchkey.evakey[j][0].GetDataLen(true)]); err != nil {
			return pointer, err
		}

		pointer += inc

		if inc, err = switchkey.evakey[j][1].WriteTo(data[pointer : pointer+switchkey.evakey[j][1].GetDataLen(true)]); err != nil {
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

func (context *Context) NewRotationKeys() (rotKey *RotationKeys) {
	rotKey = new(RotationKeys)
	return
}

// Newrotationkeys generates a new struct of rotationkeys storing the keys for the specified rotations. The provided secret-key must be the secret-key used to generate the public-key under
// which the ciphertexts to rotate are encrypted under. Bitdecomp is the power of two binary decomposition of the key. A higher bigdecomp will induce smaller keys, faster key-switching,
// but at the cost of more noise. rotLeft and rotRight must be a slice of uint64 rotations, row is a boolean value indicating if the key for the row rotation must be generated.
func (keygen *KeyGenerator) GenRot(rotType Rotation, sk *SecretKey, k uint64, rotKey *RotationKeys) {

	k &= ((keygen.context.n >> 1) - 1)

	switch rotType {
	case RotationLeft:
		if rotKey.evakeyRotColLeft == nil {
			rotKey.evakeyRotColLeft = make(map[uint64]*SwitchingKey)
		}
		if rotKey.evakeyRotColLeft[k] == nil && k != 0 {
			rotKey.evakeyRotColLeft[k] = genrotkey(keygen, sk.Get(), keygen.context.galElRotColLeft[k])
		}
	case RotationRight:
		if rotKey.evakeyRotColRight == nil {
			rotKey.evakeyRotColRight = make(map[uint64]*SwitchingKey)
		}
		if rotKey.evakeyRotColRight[k] == nil && k != 0 {
			rotKey.evakeyRotColRight[k] = genrotkey(keygen, sk.Get(), keygen.context.galElRotColRight[k])
		}
	case RotationRow:
		rotKey.evakeyRotRow = genrotkey(keygen, sk.Get(), keygen.context.galElRotRow)
	}
}

// Newrotationkeys generates a new struct of rotationkeys storing the keys of all the left and right powers of two rotations. The provided secret-key must be the secret-key used to generate the public-key under
// which the ciphertexts to rotate are encrypted under. rows is a boolean value indicatig if the keys for the row rotation have to be generated. Bitdecomp is the power of two binary decomposition of the key.
// A higher bigdecomp will induce smaller keys, faster key-switching, but at the cost of more noise.
func (keygen *KeyGenerator) NewRotationKeysPow2(sk *SecretKey) (rotKey *RotationKeys) {

	rotKey = new(RotationKeys)

	rotKey.evakeyRotColLeft = make(map[uint64]*SwitchingKey)
	rotKey.evakeyRotColRight = make(map[uint64]*SwitchingKey)

	for n := uint64(1); n < keygen.context.n>>1; n <<= 1 {

		rotKey.evakeyRotColLeft[n] = genrotkey(keygen, sk.Get(), keygen.context.galElRotColLeft[n])
		rotKey.evakeyRotColRight[n] = genrotkey(keygen, sk.Get(), keygen.context.galElRotColRight[n])
	}

	rotKey.evakeyRotRow = genrotkey(keygen, sk.Get(), keygen.context.galElRotRow)

	return
}

func (rotKey *RotationKeys) SetRotKey(rotType Rotation, k uint64, evakey [][2]*ring.Poly) {
	switch rotType {
	case RotationLeft:
		if rotKey.evakeyRotColLeft == nil {
			rotKey.evakeyRotColLeft = make(map[uint64]*SwitchingKey)
		}
		if rotKey.evakeyRotColLeft[k] == nil && k != 0 {
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
		if rotKey.evakeyRotColRight[k] == nil && k != 0 {
			rotKey.evakeyRotColRight[k] = new(SwitchingKey)
			rotKey.evakeyRotColRight[k].evakey = make([][2]*ring.Poly, len(evakey))
			for j := range evakey {
				rotKey.evakeyRotColRight[k].evakey[j][0] = evakey[j][0].CopyNew()
				rotKey.evakeyRotColRight[k].evakey[j][1] = evakey[j][1].CopyNew()
			}
		}
	case RotationRow:
		if rotKey.evakeyRotRow == nil {
			rotKey.evakeyRotRow = new(SwitchingKey)
			rotKey.evakeyRotRow.evakey = make([][2]*ring.Poly, len(evakey))
			for j := range evakey {
				rotKey.evakeyRotRow.evakey[j][0] = evakey[j][0].CopyNew()
				rotKey.evakeyRotRow.evakey[j][1] = evakey[j][1].CopyNew()
			}
		}
	}
}

// genrotkey is a methode used in the rotation-keys generation.
func genrotkey(keygen *KeyGenerator, sk *ring.Poly, gen uint64) (switchkey *SwitchingKey) {

	ringContext := keygen.context.contextKeys

	ring.PermuteNTT(sk, gen, keygen.polypool)
	ringContext.Sub(keygen.polypool, sk, keygen.polypool)

	for _, pj := range keygen.context.specialprimes {
		ringContext.MulScalar(keygen.polypool, pj, keygen.polypool)
	}

	switchkey = newswitchintkey(keygen.context, keygen.polypool, sk)
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

	if rotationkey.evakeyRotRow != nil {

		dataLen += 4 + rotationkey.evakeyRotRow.GetDataLen(true)
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

	if rotationkey.evakeyRotRow != nil {

		data[pointer] = uint8(RotationRow)
		pointer += 4

		_, _ = rotationkey.evakeyRotRow.encode(pointer, data)
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

		} else if rotationType == RotationRow {

			rotationkey.evakeyRotRow = new(SwitchingKey)
			if inc, err = rotationkey.evakeyRotRow.decode(data[pointer:]); err != nil {
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
