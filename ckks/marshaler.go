package ckks

import (
	"encoding/binary"
	"errors"
	"math"

	"github.com/ldsec/lattigo/v2/ring"
)

// GetDataLen returns the length in bytes of the target Ciphertext.
func (ciphertext *Ciphertext) GetDataLen(WithMetaData bool) (dataLen uint64) {
	// MetaData is :
	// 1 byte : Degree
	// 9 byte : Scale
	// 1 byte : isNTT
	if WithMetaData {
		dataLen += 11
	}

	for _, el := range ciphertext.value {
		dataLen += el.GetDataLen(WithMetaData)
	}

	return dataLen
}

// MarshalBinary encodes a Ciphertext on a byte slice. The total size
// in byte is 4 + 8* N * numberModuliQ * (degree + 1).
func (ciphertext *Ciphertext) MarshalBinary() (data []byte, err error) {

	data = make([]byte, ciphertext.GetDataLen(true))

	data[0] = uint8(ciphertext.Degree() + 1)

	binary.LittleEndian.PutUint64(data[1:9], math.Float64bits(ciphertext.Scale()))

	if ciphertext.isNTT {
		data[10] = 1
	}

	var pointer, inc uint64

	pointer = 11

	for _, el := range ciphertext.value {

		if inc, err = el.WriteTo(data[pointer:]); err != nil {
			return nil, err
		}

		pointer += inc
	}

	return data, nil
}

// UnmarshalBinary decodes a previously marshaled Ciphertext on the target Ciphertext.
// The target Ciphertext must be of the appropriate format and size, it can be created with the
// method NewCiphertext(uint64).
func (ciphertext *Ciphertext) UnmarshalBinary(data []byte) (err error) {
	if len(data) < 11 { // cf. ciphertext.GetDataLen()
		return errors.New("too small bytearray")
	}

	ciphertext.Element = new(Element)

	ciphertext.value = make([]*ring.Poly, uint8(data[0]))

	ciphertext.scale = math.Float64frombits(binary.LittleEndian.Uint64(data[1:9]))

	if uint8(data[10]) == 1 {
		ciphertext.isNTT = true
	}

	var pointer, inc uint64
	pointer = 11

	for i := range ciphertext.value {

		ciphertext.value[i] = new(ring.Poly)

		if inc, err = ciphertext.value[i].DecodePolyNew(data[pointer:]); err != nil {
			return err
		}

		pointer += inc
	}

	if pointer != uint64(len(data)) {
		return errors.New("remaining unparsed data")
	}

	return nil
}

// GetDataLen returns the length in bytes of the target SecretKey.
func (sk *SecretKey) GetDataLen(WithMetaData bool) (dataLen uint64) {
	return sk.sk.GetDataLen(WithMetaData)
}

// MarshalBinary encodes a SecretKey in a byte slice.
func (sk *SecretKey) MarshalBinary() (data []byte, err error) {

	data = make([]byte, sk.GetDataLen(true))

	if _, err = sk.sk.WriteTo(data); err != nil {
		return nil, err
	}

	return data, nil
}

// UnmarshalBinary decodes a previously marshaled SecretKey on the target secret-key.
func (sk *SecretKey) UnmarshalBinary(data []byte) (err error) {

	sk.sk = new(ring.Poly)

	if _, err = sk.sk.DecodePolyNew(data); err != nil {
		return err
	}

	return nil
}

// GetDataLen returns the length in bytes of the target PublicKey.
func (pk *PublicKey) GetDataLen(WithMetaData bool) (dataLen uint64) {

	for _, el := range pk.pk {
		dataLen += el.GetDataLen(WithMetaData)
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
func (evaluationkey *EvaluationKey) GetDataLen(WithMetaData bool) (dataLen uint64) {
	return evaluationkey.evakey.GetDataLen(WithMetaData)
}

// MarshalBinary encodes an evaluation key in a byte slice.
func (evaluationkey *EvaluationKey) MarshalBinary() (data []byte, err error) {

	var pointer uint64

	dataLen := evaluationkey.evakey.GetDataLen(true)

	data = make([]byte, dataLen)

	if _, err = evaluationkey.evakey.encode(pointer, data); err != nil {
		return nil, err
	}

	return data, nil
}

// UnmarshalBinary decodes a previously marshaled evaluation-key in the target evaluation-key.
func (evaluationkey *EvaluationKey) UnmarshalBinary(data []byte) (err error) {
	evaluationkey.evakey = new(SwitchingKey)
	if _, err = evaluationkey.evakey.decode(data); err != nil {
		return err
	}
	return nil
}

// GetDataLen returns the length in bytes of the target SwitchingKey.
func (switchkey *SwitchingKey) GetDataLen(WithMetaData bool) (dataLen uint64) {

	if WithMetaData {
		dataLen++
	}

	for j := uint64(0); j < uint64(len(switchkey.key)); j++ {
		dataLen += switchkey.key[j][0].GetDataLen(WithMetaData)
		dataLen += switchkey.key[j][1].GetDataLen(WithMetaData)
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

		if inc, err = switchkey.key[j][0].WriteTo(data[pointer:]); err != nil {
			return pointer, err
		}

		pointer += inc

		if inc, err = switchkey.key[j][1].WriteTo(data[pointer:]); err != nil {
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
func (rotationkey *RotationKeys) GetDataLen(WithMetaData bool) (dataLen uint64) {
	for _, key := range rotationkey.keys {
		if WithMetaData {
			dataLen += 4
		}
		dataLen += key.GetDataLen(WithMetaData)
	}
	return
}

// MarshalBinary encodes a RotationKeys struct in a byte slice.
func (rotationkey *RotationKeys) MarshalBinary() (data []byte, err error) {

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
func (rotationkey *RotationKeys) UnmarshalBinary(data []byte) (err error) {

	if rotationkey.keys == nil {
		rotationkey.keys = make(map[uint64]*SwitchingKey)
		rotationkey.permuteNTTIndex = make(map[uint64][]uint64)
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
		rotationkey.keys[galEl] = swk
		rotationkey.permuteNTTIndex[galEl] = ring.PermuteNTTIndex(galEl, uint64(swk.key[0][0].GetDegree()))

		data = data[inc:]
	}

	return nil
}
