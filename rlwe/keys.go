package rlwe

import (
	"encoding/binary"

	"github.com/ldsec/lattigo/v2/ring"
)

// SecretKey is a type for generic RLWE secret keys.
type SecretKey struct {
	Value *ring.Poly
}

// PublicKey is a type for generic RLWE public keys.
type PublicKey struct {
	Value [2]*ring.Poly
}

// SwitchingKey is a type for generic RLWE public switching keys.
type SwitchingKey struct {
	Value [][2]*ring.Poly
}

// RelinearizationKey is a type for generic RLWE public relinearization keys. It stores a slice with a
// switching key per relinearizable degree. The switching key at index i is used to relinearize a degree
// i+2 ciphertexts back to a degree i + 1 one.
type RelinearizationKey struct {
	Keys []*SwitchingKey
}

// RotationKeySet is a type for storing generic RLWE public rotation keys. It stores a map indexed by the
// galois element defining the automorphism.
type RotationKeySet struct {
	Keys map[uint64]*SwitchingKey
}

// NewSecretKey generates a new SecretKey with zero values.
func NewSecretKey(ringDegree, moduliCount int) *SecretKey {

	sk := new(SecretKey)
	sk.Value = ring.NewPoly(ringDegree, moduliCount)
	return sk
}

// NewPublicKey returns a new PublicKey with zero values.
func NewPublicKey(ringDegree, moduliCount int) (pk *PublicKey) {
	return &PublicKey{Value: [2]*ring.Poly{ring.NewPoly(ringDegree, moduliCount), ring.NewPoly(ringDegree, moduliCount)}}
}

// NewRotationKeySet returns a new RotationKeySet with pre-allocated switching keys for each distinct galoisElement value.
func NewRotationKeySet(galoisElement []uint64, ringDegree, moduliCount, decompSize int) (rotKey *RotationKeySet) {
	rotKey = new(RotationKeySet)
	rotKey.Keys = make(map[uint64]*SwitchingKey, len(galoisElement))
	for _, galEl := range galoisElement {
		rotKey.Keys[galEl] = NewSwitchingKey(ringDegree, moduliCount, decompSize)
	}
	return
}

// GetRotationKey return the rotation key for the given galois element or nil if such key is not in the set. The
// second argument is true  iff the first one is non-nil.
func (rtks *RotationKeySet) GetRotationKey(galoisEl uint64) (*SwitchingKey, bool) {
	rotKey, inSet := rtks.Keys[galoisEl]
	return rotKey, inSet
}

// NewSwitchingKey returns a new public switching key with pre-allocated zero-value
func NewSwitchingKey(ringDegree, moduliCount, decompSize int) *SwitchingKey {

	swk := new(SwitchingKey)

	swk.Value = make([][2]*ring.Poly, int(decompSize))

	for i := 0; i < decompSize; i++ {
		swk.Value[i][0] = ring.NewPoly(ringDegree, moduliCount)
		swk.Value[i][1] = ring.NewPoly(ringDegree, moduliCount)
	}

	return swk
}

// NewRelinKey creates a new EvaluationKey with zero values.
func NewRelinKey(maxRelinDegree, ringDegree, moduliCount, decompSize int) (evakey *RelinearizationKey) {

	evakey = new(RelinearizationKey)

	evakey.Keys = make([]*SwitchingKey, maxRelinDegree)

	for d := 0; d < maxRelinDegree; d++ {
		evakey.Keys[d] = NewSwitchingKey(ringDegree, moduliCount, decompSize)
	}

	return
}

// GetDataLen returns the length in bytes of the target SecretKey.
func (sk *SecretKey) GetDataLen(WithMetadata bool) (dataLen int) {
	return sk.Value.GetDataLen(WithMetadata)
}

// MarshalBinary encodes a secret key in a byte slice.
func (sk *SecretKey) MarshalBinary() (data []byte, err error) {

	data = make([]byte, sk.GetDataLen(true))

	if _, err = sk.Value.WriteTo(data); err != nil {
		return nil, err
	}

	return data, nil
}

// UnmarshalBinary decodes a previously marshaled SecretKey in the target SecretKey.
func (sk *SecretKey) UnmarshalBinary(data []byte) (err error) {

	sk.Value = new(ring.Poly)

	if _, err = sk.Value.DecodePolyNew(data); err != nil {
		return err
	}

	return nil
}

// GetDataLen returns the length in bytes of the target PublicKey.
func (pk *PublicKey) GetDataLen(WithMetadata bool) (dataLen int) {

	for _, el := range pk.Value {
		dataLen += el.GetDataLen(WithMetadata)
	}

	return
}

// MarshalBinary encodes a PublicKey in a byte slice.
func (pk *PublicKey) MarshalBinary() (data []byte, err error) {

	dataLen := pk.GetDataLen(true)

	data = make([]byte, dataLen)

	var pointer, inc int

	if inc, err = pk.Value[0].WriteTo(data[pointer:]); err != nil {
		return nil, err
	}

	if _, err = pk.Value[1].WriteTo(data[pointer+inc:]); err != nil {
		return nil, err
	}

	return data, err

}

// UnmarshalBinary decodes a previously marshaled PublicKey in the target PublicKey.
func (pk *PublicKey) UnmarshalBinary(data []byte) (err error) {

	var pointer, inc int

	pk.Value[0] = new(ring.Poly)
	pk.Value[1] = new(ring.Poly)

	if inc, err = pk.Value[0].DecodePolyNew(data[pointer:]); err != nil {
		return err
	}

	if _, err = pk.Value[1].DecodePolyNew(data[pointer+inc:]); err != nil {
		return err
	}

	return nil
}

// GetDataLen returns the length in bytes of the target EvaluationKey.
func (evaluationkey *RelinearizationKey) GetDataLen(WithMetadata bool) (dataLen int) {

	if WithMetadata {
		dataLen++
	}

	for _, evakey := range evaluationkey.Keys {
		dataLen += (*SwitchingKey)(evakey).GetDataLen(WithMetadata)
	}

	return
}

// MarshalBinary encodes an EvaluationKey key in a byte slice.
func (evaluationkey *RelinearizationKey) MarshalBinary() (data []byte, err error) {

	var pointer int

	dataLen := evaluationkey.GetDataLen(true)

	data = make([]byte, dataLen)

	data[0] = uint8(len(evaluationkey.Keys))

	pointer++

	for _, evakey := range evaluationkey.Keys {

		if pointer, err = (*SwitchingKey)(evakey).encode(pointer, data); err != nil {
			return nil, err
		}
	}

	return data, nil
}

// UnmarshalBinary decodes a previously marshaled EvaluationKey in the target EvaluationKey.
func (evaluationkey *RelinearizationKey) UnmarshalBinary(data []byte) (err error) {

	deg := int(data[0])

	evaluationkey.Keys = make([]*SwitchingKey, deg)

	pointer := int(1)
	var inc int
	for i := 0; i < deg; i++ {
		evaluationkey.Keys[i] = new(SwitchingKey)
		if inc, err = evaluationkey.Keys[i].decode(data[pointer:]); err != nil {
			return err
		}
		pointer += inc
	}

	return nil
}

// GetDataLen returns the length in bytes of the target SwitchingKey.
func (switchkey *SwitchingKey) GetDataLen(WithMetadata bool) (dataLen int) {

	if WithMetadata {
		dataLen++
	}

	for j := uint64(0); j < uint64(len(switchkey.Value)); j++ {
		dataLen += switchkey.Value[j][0].GetDataLen(WithMetadata)
		dataLen += switchkey.Value[j][1].GetDataLen(WithMetadata)
	}

	return
}

// MarshalBinary encodes an SwitchingKey in a byte slice.
func (switchkey *SwitchingKey) MarshalBinary() (data []byte, err error) {

	data = make([]byte, switchkey.GetDataLen(true))

	if _, err = switchkey.encode(0, data); err != nil {
		return nil, err
	}

	return data, nil
}

// UnmarshalBinary decode a previously marshaled SwitchingKey in the target SwitchingKey.
func (switchkey *SwitchingKey) UnmarshalBinary(data []byte) (err error) {

	if _, err = switchkey.decode(data); err != nil {
		return err
	}

	return nil
}

func (switchkey *SwitchingKey) encode(pointer int, data []byte) (int, error) {

	var err error
	var inc int

	data[pointer] = uint8(len(switchkey.Value))

	pointer++

	for j := 0; j < len(switchkey.Value); j++ {

		if inc, err = switchkey.Value[j][0].WriteTo(data[pointer : pointer+switchkey.Value[j][0].GetDataLen(true)]); err != nil {
			return pointer, err
		}

		pointer += inc

		if inc, err = switchkey.Value[j][1].WriteTo(data[pointer : pointer+switchkey.Value[j][1].GetDataLen(true)]); err != nil {
			return pointer, err
		}

		pointer += inc
	}

	return pointer, nil
}

func (switchkey *SwitchingKey) decode(data []byte) (pointer int, err error) {

	decomposition := int(data[0])

	pointer = 1

	switchkey.Value = make([][2]*ring.Poly, decomposition)

	var inc int

	for j := 0; j < decomposition; j++ {

		switchkey.Value[j][0] = new(ring.Poly)
		if inc, err = switchkey.Value[j][0].DecodePolyNew(data[pointer:]); err != nil {
			return pointer, err
		}
		pointer += inc

		switchkey.Value[j][1] = new(ring.Poly)
		if inc, err = switchkey.Value[j][1].DecodePolyNew(data[pointer:]); err != nil {
			return pointer, err
		}
		pointer += inc

	}

	return pointer, nil
}

// GetDataLen returns the length in bytes of the target RotationKeys.
func (rtks *RotationKeySet) GetDataLen(WithMetaData bool) (dataLen int) {
	for _, k := range rtks.Keys {
		if WithMetaData {
			dataLen += 4
		}
		dataLen += k.GetDataLen(WithMetaData)
	}
	return
}

// MarshalBinary encodes a RotationKeys struct in a byte slice.
func (rtks *RotationKeySet) MarshalBinary() (data []byte, err error) {

	data = make([]byte, rtks.GetDataLen(true))

	pointer := int(0)

	for galEL, key := range rtks.Keys {

		binary.BigEndian.PutUint32(data[pointer:pointer+4], uint32(galEL))
		pointer += 4

		if pointer, err = key.encode(pointer, data); err != nil {
			return nil, err
		}
	}

	return data, nil
}

// UnmarshalBinary decodes a previously marshaled RotationKeys in the target RotationKeys.
func (rtks *RotationKeySet) UnmarshalBinary(data []byte) (err error) {

	rtks.Keys = make(map[uint64]*SwitchingKey)

	for len(data) > 0 {

		galEl := uint64(binary.BigEndian.Uint32(data))
		data = data[4:]

		swk := new(SwitchingKey)
		var inc int
		if inc, err = swk.decode(data); err != nil {
			return err
		}
		data = data[inc:]
		rtks.Keys[galEl] = swk

	}

	return nil
}
