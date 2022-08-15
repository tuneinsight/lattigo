package rlwe

import (
	"encoding/binary"
	"errors"

	"github.com/tuneinsight/lattigo/v3/ring"
)

// GetDataLen returns the length in bytes of the target Ciphertext.
func (el *Ciphertext) GetDataLen(WithMetaData bool) (dataLen int) {
	// MetaData is :
	// 1 byte : Degree
	if WithMetaData {
		dataLen++
	}

	for _, el := range el.Value {
		dataLen += el.GetDataLen64(WithMetaData)
	}

	return dataLen
}

// MarshalBinary encodes a Ciphertext on a byte slice. The total size
// in byte is 4 + 8* N * numberModuliQ * (degree + 1).
func (el *Ciphertext) MarshalBinary() (data []byte, err error) {

	data = make([]byte, el.GetDataLen(true))

	data[0] = uint8(el.Degree() + 1)

	var pointer, inc int

	pointer = 1

	for _, el := range el.Value {

		if inc, err = el.WriteTo64(data[pointer:]); err != nil {
			return nil, err
		}

		pointer += inc
	}

	return data, nil
}

// UnmarshalBinary decodes a previously marshaled Ciphertext on the target Ciphertext.
func (el *Ciphertext) UnmarshalBinary(data []byte) (err error) {
	if len(data) < 10 { // cf. Ciphertext.GetDataLen()
		return errors.New("too small bytearray")
	}

	el.Value = make([]*ring.Poly, uint8(data[0]))

	var pointer, inc int
	pointer = 1

	for i := range el.Value {

		el.Value[i] = new(ring.Poly)

		if inc, err = el.Value[i].DecodePoly64(data[pointer:]); err != nil {
			return err
		}

		pointer += inc
	}

	if pointer != len(data) {
		return errors.New("remaining unparsed data")
	}

	return nil
}

// GetDataLen64 returns the length in bytes of the target SecretKey.
// Assumes that each coefficient uses 8 bytes.
func (sk *SecretKey) GetDataLen64(WithMetadata bool) (dataLen int) {
	return sk.Value.GetDataLen64(WithMetadata)
}

// MarshalBinary encodes a secret key in a byte slice.
func (sk *SecretKey) MarshalBinary() (data []byte, err error) {
	data = make([]byte, sk.GetDataLen64(true))
	if _, err = sk.Value.WriteTo64(data); err != nil {
		return nil, err
	}
	return
}

// UnmarshalBinary decodes a previously marshaled SecretKey in the target SecretKey.
func (sk *SecretKey) UnmarshalBinary(data []byte) (err error) {
	_, err = sk.Value.DecodePoly64(data)
	return
}

// GetDataLen64 returns the length in bytes of the target PublicKey.
func (pk *PublicKey) GetDataLen64(WithMetadata bool) (dataLen int) {
	return pk.CiphertextQP.GetDataLen64(WithMetadata)
}

// MarshalBinary encodes a PublicKey in a byte slice.
func (pk *PublicKey) MarshalBinary() (data []byte, err error) {
	data = make([]byte, pk.GetDataLen64(true))
	_, err = pk.CiphertextQP.WriteTo64(data)
	return
}

// UnmarshalBinary decodes a previously marshaled PublicKey in the target PublicKey.
func (pk *PublicKey) UnmarshalBinary(data []byte) (err error) {
	_, err = pk.CiphertextQP.Decode64(data)
	return
}

// MarshalBinary encodes the target SwitchingKey on a slice of bytes.
func (swk *SwitchingKey) MarshalBinary() (data []byte, err error) {
	return swk.GadgetCiphertext.MarshalBinary()
}

// UnmarshalBinary decodes a slice of bytes on the target SwitchingKey.
func (swk *SwitchingKey) UnmarshalBinary(data []byte) (err error) {
	return swk.GadgetCiphertext.UnmarshalBinary(data)
}

// GetDataLen returns the length in bytes of the target EvaluationKey.
func (rlk *RelinearizationKey) GetDataLen(WithMetadata bool) (dataLen int) {

	if WithMetadata {
		dataLen++
	}

	return dataLen + len(rlk.Keys)*rlk.Keys[0].GetDataLen(WithMetadata)
}

// MarshalBinary encodes an EvaluationKey key in a byte slice.
func (rlk *RelinearizationKey) MarshalBinary() (data []byte, err error) {

	var pointer int

	dataLen := rlk.GetDataLen(true)

	data = make([]byte, dataLen)

	data[0] = uint8(len(rlk.Keys))

	pointer++

	for _, evakey := range rlk.Keys {

		if pointer, err = evakey.Encode(pointer, data); err != nil {
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
		if inc, err = rlk.Keys[i].Decode(data[pointer:]); err != nil {
			return err
		}
		pointer += inc
	}

	return nil
}

// GetDataLen returns the length in bytes of the target RotationKeys.
func (rtks *RotationKeySet) GetDataLen(WithMetaData bool) (dataLen int) {
	for _, k := range rtks.Keys {
		if WithMetaData {
			dataLen += 8
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

		binary.BigEndian.PutUint64(data[pointer:pointer+8], galEL)
		pointer += 8

		if pointer, err = key.Encode(pointer, data); err != nil {
			return nil, err
		}
	}

	return data, nil
}

// UnmarshalBinary decodes a previously marshaled RotationKeys in the target RotationKeys.
func (rtks *RotationKeySet) UnmarshalBinary(data []byte) (err error) {

	rtks.Keys = make(map[uint64]*SwitchingKey)

	for len(data) > 0 {

		galEl := binary.BigEndian.Uint64(data)
		data = data[8:]

		swk := new(SwitchingKey)
		var inc int
		if inc, err = swk.Decode(data); err != nil {
			return err
		}
		data = data[inc:]
		rtks.Keys[galEl] = swk

	}

	return nil
}
