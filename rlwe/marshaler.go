package rlwe

import (
	"encoding/binary"
)

// MarshalBinarySize returns the length in bytes of the target SecretKey.
func (sk *SecretKey) MarshalBinarySize() (dataLen int) {
	return sk.Poly.MarshalBinarySize64() + sk.MetaData.MarshalBinarySize()
}

// MarshalBinary encodes a secret key in a byte slice.
func (sk *SecretKey) MarshalBinary() (data []byte, err error) {
	data = make([]byte, sk.MarshalBinarySize())

	var ptr int
	if ptr, err = sk.MetaData.Encode64(data[ptr:]); err != nil {
		return nil, err
	}

	if _, err = sk.Poly.Encode64(data[ptr:]); err != nil {
		return nil, err
	}

	return
}

// UnmarshalBinary decodes a previously marshaled SecretKey in the target SecretKey.
func (sk *SecretKey) UnmarshalBinary(data []byte) (err error) {

	var ptr int
	if ptr, err = sk.MetaData.Decode64(data[ptr:]); err != nil {
		return
	}

	if _, err = sk.Poly.Decode64(data[ptr:]); err != nil {
		return
	}

	return
}

// MarshalBinarySize returns the length in bytes of the target PublicKey.
func (pk *PublicKey) MarshalBinarySize() (dataLen int) {
	return pk.Value[0].MarshalBinarySize64() + pk.Value[1].MarshalBinarySize64() + pk.MetaData.MarshalBinarySize()
}

// MarshalBinary encodes a PublicKey in a byte slice.
func (pk *PublicKey) MarshalBinary() (data []byte, err error) {
	data = make([]byte, pk.MarshalBinarySize())
	var inc, ptr int

	if inc, err = pk.MetaData.Encode64(data[ptr:]); err != nil {
		return nil, err
	}

	ptr += inc

	if inc, err = pk.Value[0].Encode64(data[ptr:]); err != nil {
		return nil, err
	}

	ptr += inc

	if _, err = pk.Value[1].Encode64(data[ptr:]); err != nil {
		return nil, err
	}

	return
}

// UnmarshalBinary decodes a previously marshaled PublicKey in the target PublicKey.
func (pk *PublicKey) UnmarshalBinary(data []byte) (err error) {

	var ptr, inc int

	if inc, err = pk.MetaData.Decode64(data[ptr:]); err != nil {
		return
	}

	ptr += inc

	if inc, err = pk.Value[0].Decode64(data[ptr:]); err != nil {
		return
	}

	ptr += inc

	if _, err = pk.Value[1].Decode64(data[ptr:]); err != nil {
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

// MarshalBinarySize returns the length in bytes of the target EvaluationKey.
func (rlk *RelinearizationKey) MarshalBinarySize() (dataLen int) {
	return 1 + len(rlk.Keys)*rlk.Keys[0].MarshalBinarySize()
}

// MarshalBinary encodes an EvaluationKey key in a byte slice.
func (rlk *RelinearizationKey) MarshalBinary() (data []byte, err error) {

	data = make([]byte, rlk.MarshalBinarySize())

	var ptr int

	data[0] = uint8(len(rlk.Keys))

	ptr++

	var inc int
	for _, evakey := range rlk.Keys {

		if inc, err = evakey.Encode(data[ptr:]); err != nil {
			return nil, err
		}
		ptr += inc
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

// MarshalBinarySize returns the length in bytes of the target RotationKeys.
func (rtks *RotationKeySet) MarshalBinarySize() (dataLen int) {
	for _, k := range rtks.Keys {
		dataLen += 8 + k.MarshalBinarySize()
	}
	return
}

// MarshalBinary encodes a RotationKeys struct in a byte slice.
func (rtks *RotationKeySet) MarshalBinary() (data []byte, err error) {

	data = make([]byte, rtks.MarshalBinarySize())

	ptr := int(0)

	var inc int
	for galEL, key := range rtks.Keys {

		binary.BigEndian.PutUint64(data[ptr:], galEL)
		ptr += 8

		if inc, err = key.Encode(data[ptr:]); err != nil {
			return nil, err
		}

		ptr += inc
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
