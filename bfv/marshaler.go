package bfv

import (
	"encoding/binary"

	"github.com/ldsec/lattigo/v2/ring"
)

// MarshalBinary encodes a Ciphertext in a byte slice.
func (ciphertext *Ciphertext) MarshalBinary() (data []byte, err error) {

	data = make([]byte, ciphertext.GetDataLen(true))

	data[0] = uint8(len(ciphertext.value))

	var pointer, inc uint64

	pointer = 1

	for _, el := range ciphertext.value {

		if inc, err = el.WriteTo(data[pointer:]); err != nil {
			return nil, err
		}

		pointer += inc
	}

	return data, nil
}

// UnmarshalBinary decodes a previously marshaled Ciphertext in the target Ciphertext.
func (ciphertext *Ciphertext) UnmarshalBinary(data []byte) (err error) {

	ciphertext.Element = new(Element)

	ciphertext.value = make([]*ring.Poly, uint8(data[0]))

	var pointer, inc uint64
	pointer = 1

	for i := range ciphertext.value {

		ciphertext.value[i] = new(ring.Poly)

		if inc, err = ciphertext.value[i].DecodePolyNew(data[pointer:]); err != nil {
			return err
		}

		pointer += inc
	}

	return nil
}

// GetDataLen returns the length in bytes of the target Ciphertext.
func (ciphertext *Ciphertext) GetDataLen(WithMetaData bool) (dataLen uint64) {
	if WithMetaData {
		dataLen++
	}

	for _, el := range ciphertext.value {
		dataLen += el.GetDataLen(WithMetaData)
	}

	return dataLen
}

// GetDataLen returns the length in bytes of the target SecretKey.
func (sk *SecretKey) GetDataLen(WithMetadata bool) (dataLen uint64) {
	return sk.sk.GetDataLen(WithMetadata)
}

// MarshalBinary encodes a secret key in a byte slice.
func (sk *SecretKey) MarshalBinary() (data []byte, err error) {

	data = make([]byte, sk.GetDataLen(true))

	if _, err = sk.sk.WriteTo(data); err != nil {
		return nil, err
	}

	return data, nil
}

// UnmarshalBinary decodes a previously marshaled SecretKey in the target SecretKey.
func (sk *SecretKey) UnmarshalBinary(data []byte) (err error) {

	sk.sk = new(ring.Poly)

	if _, err = sk.sk.DecodePolyNew(data); err != nil {
		return err
	}

	return nil
}

// GetDataLen returns the length in bytes of the target PublicKey.
func (pk *PublicKey) GetDataLen(WithMetadata bool) (dataLen uint64) {

	for _, el := range pk.pk {
		dataLen += el.GetDataLen(WithMetadata)
	}

	return
}

// MarshalBinary encodes a PublicKey in a byte slice.
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

// UnmarshalBinary decodes a previously marshaled PublicKey in the target PublicKey.
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

// GetDataLen returns the length in bytes of the target EvaluationKey.
func (evaluationkey *RelinearizationKey) GetDataLen(WithMetadata bool) (dataLen uint64) {

	if WithMetadata {
		dataLen++
	}

	for _, evakey := range evaluationkey.keys {
		dataLen += evakey.GetDataLen(WithMetadata)
	}

	return
}

// MarshalBinary encodes an EvaluationKey key in a byte slice.
func (evaluationkey *RelinearizationKey) MarshalBinary() (data []byte, err error) {

	var pointer uint64

	dataLen := evaluationkey.GetDataLen(true)

	data = make([]byte, dataLen)

	data[0] = uint8(len(evaluationkey.keys))

	pointer++

	for _, evakey := range evaluationkey.keys {

		if pointer, err = evakey.encode(pointer, data); err != nil {
			return nil, err
		}
	}

	return data, nil
}

// UnmarshalBinary decodes a previously marshaled EvaluationKey in the target EvaluationKey.
func (evaluationkey *RelinearizationKey) UnmarshalBinary(data []byte) (err error) {

	deg := uint64(data[0])

	evaluationkey.keys = make([]*SwitchingKey, deg)

	pointer := uint64(1)
	var inc uint64
	for i := uint64(0); i < deg; i++ {
		evaluationkey.keys[i] = new(SwitchingKey)
		if inc, err = evaluationkey.keys[i].decode(data[pointer:]); err != nil {
			return err
		}
		pointer += inc
	}

	return nil
}

// GetDataLen returns the length in bytes of the target SwitchingKey.
func (switchkey *SwitchingKey) GetDataLen(WithMetadata bool) (dataLen uint64) {

	if WithMetadata {
		dataLen++
	}

	for j := uint64(0); j < uint64(len(switchkey.key)); j++ {
		dataLen += switchkey.key[j][0].GetDataLen(WithMetadata)
		dataLen += switchkey.key[j][1].GetDataLen(WithMetadata)
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

func (switchkey *SwitchingKey) encode(pointer uint64, data []byte) (uint64, error) {

	var err error

	var inc uint64

	data[pointer] = uint8(len(switchkey.key))

	pointer++

	for j := uint64(0); j < uint64(len(switchkey.key)); j++ {

		if inc, err = switchkey.key[j][0].WriteTo(data[pointer : pointer+switchkey.key[j][0].GetDataLen(true)]); err != nil {
			return pointer, err
		}

		pointer += inc

		if inc, err = switchkey.key[j][1].WriteTo(data[pointer : pointer+switchkey.key[j][1].GetDataLen(true)]); err != nil {
			return pointer, err
		}

		pointer += inc
	}

	return pointer, nil
}

func (switchkey *SwitchingKey) decode(data []byte) (pointer uint64, err error) {

	decomposition := uint64(data[0])

	pointer = uint64(1)

	switchkey.key = make([][2]*ring.Poly, decomposition)

	var inc uint64

	for j := uint64(0); j < decomposition; j++ {

		switchkey.key[j][0] = new(ring.Poly)
		if inc, err = switchkey.key[j][0].DecodePolyNew(data[pointer:]); err != nil {
			return pointer, err
		}
		pointer += inc

		switchkey.key[j][1] = new(ring.Poly)
		if inc, err = switchkey.key[j][1].DecodePolyNew(data[pointer:]); err != nil {
			return pointer, err
		}
		pointer += inc

	}

	return pointer, nil
}

// GetDataLen returns the length in bytes of the target RotationKeys.
func (rotationkey *RotationKeySet) GetDataLen(WithMetaData bool) (dataLen uint64) {
	for _, k := range rotationkey.keys {
		if WithMetaData {
			dataLen += 4
		}
		dataLen += k.GetDataLen(WithMetaData)
	}
	return
}

// MarshalBinary encodes a RotationKeys struct in a byte slice.
func (rotationkey *RotationKeySet) MarshalBinary() (data []byte, err error) {

	data = make([]byte, rotationkey.GetDataLen(true))

	pointer := uint64(0)

	for galEL, key := range rotationkey.keys {

		binary.BigEndian.PutUint32(data[pointer:pointer+4], uint32(galEL))
		pointer += 4

		if pointer, err = key.encode(pointer, data); err != nil {
			return nil, err
		}
	}

	return data, nil
}

// UnmarshalBinary decodes a previously marshaled RotationKeys in the target RotationKeys.
func (rotationkey *RotationKeySet) UnmarshalBinary(data []byte) (err error) {

	if rotationkey.keys == nil {
		rotationkey.keys = make(map[uint64]*SwitchingKey)
	} else {
		rotationkey.Delete()
	}

	for len(data) > 0 {

		galEl := uint64(binary.BigEndian.Uint32(data))
		data = data[4:]

		swk := new(SwitchingKey)
		var inc uint64
		if inc, err = swk.decode(data); err != nil {
			return err
		}
		data = data[inc:]
		rotationkey.keys[galEl] = swk

	}

	return nil
}
