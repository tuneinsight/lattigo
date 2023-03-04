package rlwe

import (
	"encoding/binary"
)

// MarshalBinarySize returns the length in bytes of the target SecretKey.
func (sk *SecretKey) MarshalBinarySize() (dataLen int) {
	return sk.Value.MarshalBinarySize64()
}

// MarshalBinary encodes a secret key in a byte slice.
func (sk *SecretKey) MarshalBinary() (data []byte, err error) {
	data = make([]byte, sk.MarshalBinarySize())

	if _, err = sk.Value.Encode64(data); err != nil {
		return nil, err
	}

	return
}

// UnmarshalBinary decodes a previously marshaled SecretKey in the target SecretKey.
func (sk *SecretKey) UnmarshalBinary(data []byte) (err error) {

	if _, err = sk.Value.Decode64(data); err != nil {
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

// MarshalBinary encodes the target EvaluationKey on a slice of bytes.
func (evk *EvaluationKey) MarshalBinary() (data []byte, err error) {
	return evk.GadgetCiphertext.MarshalBinary()
}

// UnmarshalBinary decodes a slice of bytes on the target EvaluationKey.
func (evk *EvaluationKey) UnmarshalBinary(data []byte) (err error) {
	return evk.GadgetCiphertext.UnmarshalBinary(data)
}

func (rlk *RelinearizationKey) MarshalBinary() (data []byte, err error) {
	return rlk.GadgetCiphertext.MarshalBinary()
}

func (rlk *RelinearizationKey) UnmarshalBinary(data []byte) (err error) {
	if rlk.EvaluationKey == nil {
		rlk.EvaluationKey = &EvaluationKey{}
	}

	return rlk.GadgetCiphertext.UnmarshalBinary(data)
}

func (gk *GaloisKey) MarshalBinary() (data []byte, err error) {
	data = make([]byte, gk.EvaluationKey.MarshalBinarySize()+16)
	binary.LittleEndian.PutUint64(data[0:], gk.GaloisElement)
	binary.LittleEndian.PutUint64(data[8:], gk.NthRoot)
	_, err = gk.EvaluationKey.GadgetCiphertext.Encode(data[16:])
	return
}

func (gk *GaloisKey) UnmarshalBinary(data []byte) (err error) {
	gk.GaloisElement = binary.LittleEndian.Uint64(data[0:])
	gk.NthRoot = binary.LittleEndian.Uint64(data[8:])

	if gk.EvaluationKey == nil {
		gk.EvaluationKey = &EvaluationKey{}
	}

	_, err = gk.EvaluationKey.GadgetCiphertext.Decode(data[16:])
	return
}
