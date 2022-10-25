package rlwe

import (
	"encoding/binary"
)

// GetDataLen64 returns the length in bytes of the target SecretKey.
// Assumes that each coefficient uses 8 bytes.
func (sk *SecretKey) GetDataLen64(WithMetadata bool) (dataLen int) {
	return sk.Value.GetDataLen64(WithMetadata)
}

// MarshalBinary encodes a secret key in a byte slice.
func (sk *SecretKey) MarshalBinary() (data []byte, err error) {
	data = make([]byte, sk.GetDataLen64(true))
	if _, err = sk.Value.Encode64(data); err != nil {
		return nil, err
	}
	return
}

// UnmarshalBinary decodes a previously marshaled SecretKey in the target SecretKey.
func (sk *SecretKey) UnmarshalBinary(data []byte) (err error) {
	_, err = sk.Value.Decode64(data)
	return
}

// GetDataLen64 returns the length in bytes of the target PublicKey.
func (pk *PublicKey) GetDataLen64(WithMetadata bool) (dataLen int) {
	return pk.Value[0].GetDataLen64(WithMetadata) + pk.Value[1].GetDataLen64(WithMetadata)
}

// MarshalBinary encodes a PublicKey in a byte slice.
func (pk *PublicKey) MarshalBinary() (data []byte, err error) {
	data = make([]byte, pk.GetDataLen64(true))
	var inc, pt int
	if inc, err = pk.Value[0].Encode64(data[pt:]); err != nil {
		return nil, err
	}
	pt += inc

	if _, err = pk.Value[1].Encode64(data[pt:]); err != nil {
		return nil, err
	}

	return
}

// UnmarshalBinary decodes a previously marshaled PublicKey in the target PublicKey.
func (pk *PublicKey) UnmarshalBinary(data []byte) (err error) {

	var pt, inc int
	if inc, err = pk.Value[0].Decode64(data[pt:]); err != nil {
		return
	}
	pt += inc

	if _, err = pk.Value[1].Decode64(data[pt:]); err != nil {
		return
	}

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
