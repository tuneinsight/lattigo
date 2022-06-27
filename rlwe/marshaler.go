package rlwe

import (
	"encoding/binary"
	"errors"

	"github.com/tuneinsight/lattigo/v3/ring"
)

// GetDataLen returns the length in bytes of the target Ciphertext.
func (ciphertext *Ciphertext) GetDataLen(WithMetaData bool) (dataLen int) {
	// MetaData is :
	// 1 byte : Degree
	if WithMetaData {
		dataLen++
	}

	for _, el := range ciphertext.Value {
		dataLen += el.GetDataLen(WithMetaData)
	}

	return dataLen
}

// MarshalBinary encodes a Ciphertext on a byte slice. The total size
// in byte is 4 + 8* N * numberModuliQ * (degree + 1).
func (ciphertext *Ciphertext) MarshalBinary() (data []byte, err error) {

	data = make([]byte, ciphertext.GetDataLen(true))

	data[0] = uint8(ciphertext.Degree() + 1)

	var pointer, inc int

	pointer = 1

	for _, el := range ciphertext.Value {

		if inc, err = el.WriteTo(data[pointer:]); err != nil {
			return nil, err
		}

		pointer += inc
	}

	return data, nil
}

// UnmarshalBinary decodes a previously marshaled Ciphertext on the target Ciphertext.
func (ciphertext *Ciphertext) UnmarshalBinary(data []byte) (err error) {
	if len(data) < 10 { // cf. ciphertext.GetDataLen()
		return errors.New("too small bytearray")
	}

	ciphertext.Value = make([]*ring.Poly, uint8(data[0]))

	var pointer, inc int
	pointer = 1

	for i := range ciphertext.Value {

		ciphertext.Value[i] = new(ring.Poly)

		if inc, err = ciphertext.Value[i].DecodePolyNew(data[pointer:]); err != nil {
			return err
		}

		pointer += inc
	}

	if pointer != len(data) {
		return errors.New("remaining unparsed data")
	}

	return nil
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
	return
}

// UnmarshalBinary decodes a previously marshaled SecretKey in the target SecretKey.
func (sk *SecretKey) UnmarshalBinary(data []byte) (err error) {
	_, err = sk.Value.DecodePolyNew(data)
	return
}

// GetDataLen returns the length in bytes of the target PublicKey.
func (pk *PublicKey) GetDataLen(WithMetadata bool) (dataLen int) {
	return pk.Value[0].GetDataLen(WithMetadata) + pk.Value[1].GetDataLen(WithMetadata)
}

// MarshalBinary encodes a PublicKey in a byte slice.
func (pk *PublicKey) MarshalBinary() (data []byte, err error) {
	data = make([]byte, pk.GetDataLen(true))
	var inc, pt int
	if inc, err = pk.Value[0].WriteTo(data[pt:]); err != nil {
		return nil, err
	}
	pt += inc

	if _, err = pk.Value[1].WriteTo(data[pt:]); err != nil {
		return nil, err
	}

	return
}

// UnmarshalBinary decodes a previously marshaled PublicKey in the target PublicKey.
func (pk *PublicKey) UnmarshalBinary(data []byte) (err error) {

	var pt, inc int
	if inc, err = pk.Value[0].DecodePolyNew(data[pt:]); err != nil {
		return
	}
	pt += inc

	if _, err = pk.Value[1].DecodePolyNew(data[pt:]); err != nil {
		return
	}

	return
}

// GetDataLen returns the length in bytes of the target EvaluationKey.
func (rlk *RelinearizationKey) GetDataLen(WithMetadata bool) (dataLen int) {

	if WithMetadata {
		dataLen++
	}

	for _, evakey := range rlk.Keys {
		dataLen += (*SwitchingKey)(evakey).GetDataLen(WithMetadata)
	}

	return
}

// MarshalBinary encodes an EvaluationKey key in a byte slice.
func (rlk *RelinearizationKey) MarshalBinary() (data []byte, err error) {

	var pointer int

	dataLen := rlk.GetDataLen(true)

	data = make([]byte, dataLen)

	data[0] = uint8(len(rlk.Keys))

	pointer++

	for _, evakey := range rlk.Keys {

		if pointer, err = (*SwitchingKey)(evakey).encode(pointer, data); err != nil {
			return nil, err
		}
	}

	return data, nil
}

// UnmarshalBinary decodes a previously marshaled EvaluationKey in the target EvaluationKey.
func (rlk *RelinearizationKey) UnmarshalBinary(data []byte) (err error) {

	deg := int(data[0])

	rlk.Keys = make([]*SwitchingKey, deg)

	pointer := 1
	var inc int
	for i := 0; i < deg; i++ {
		rlk.Keys[i] = new(SwitchingKey)
		if inc, err = rlk.Keys[i].decode(data[pointer:]); err != nil {
			return err
		}
		pointer += inc
	}

	return nil
}

// GetDataLen returns the length in bytes of the target SwitchingKey.
func (swk *SwitchingKey) GetDataLen(WithMetadata bool) (dataLen int) {

	if WithMetadata {
		dataLen++
	}

	for j := uint64(0); j < uint64(len(swk.Value)); j++ {
		dataLen += swk.Value[j][0].GetDataLen(WithMetadata)
		dataLen += swk.Value[j][1].GetDataLen(WithMetadata)
	}

	return
}

// MarshalBinary encodes an SwitchingKey in a byte slice.
func (swk *SwitchingKey) MarshalBinary() (data []byte, err error) {

	data = make([]byte, swk.GetDataLen(true))

	if _, err = swk.encode(0, data); err != nil {
		return nil, err
	}

	return data, nil
}

// UnmarshalBinary decode a previously marshaled SwitchingKey in the target SwitchingKey.
func (swk *SwitchingKey) UnmarshalBinary(data []byte) (err error) {

	if _, err = swk.decode(data); err != nil {
		return err
	}

	return nil
}

func (swk *SwitchingKey) encode(pointer int, data []byte) (int, error) {

	var err error
	var inc int

	data[pointer] = uint8(len(swk.Value))

	pointer++

	for j := 0; j < len(swk.Value); j++ {

		if inc, err = swk.Value[j][0].WriteTo(data[pointer : pointer+swk.Value[j][0].GetDataLen(true)]); err != nil {
			return pointer, err
		}

		pointer += inc

		if inc, err = swk.Value[j][1].WriteTo(data[pointer : pointer+swk.Value[j][1].GetDataLen(true)]); err != nil {
			return pointer, err
		}

		pointer += inc
	}

	return pointer, nil
}

func (swk *SwitchingKey) decode(data []byte) (pointer int, err error) {

	decomposition := int(data[0])

	pointer = 1

	swk.Value = make([][2]PolyQP, decomposition)

	var inc int

	for j := 0; j < decomposition; j++ {

		swk.Value[j][0].Q = new(ring.Poly)
		if inc, err = swk.Value[j][0].DecodePolyNew(data[pointer:]); err != nil {
			return
		}
		pointer += inc

		swk.Value[j][1].P = new(ring.Poly)
		if inc, err = swk.Value[j][1].DecodePolyNew(data[pointer:]); err != nil {
			return
		}
		pointer += inc
	}

	return
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
