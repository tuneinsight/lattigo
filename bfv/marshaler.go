package bfv

import (
	"github.com/ldsec/lattigo/v2/rlwe"
)

// MarshalBinary encodes a Ciphertext in a byte slice.
func (ciphertext *Ciphertext) MarshalBinary() (data []byte, err error) {
	return ciphertext.Ciphertext.MarshalBinary()
}

// UnmarshalBinary decodes a previously marshaled Ciphertext in the target Ciphertext.
func (ciphertext *Ciphertext) UnmarshalBinary(data []byte) (err error) {
	ciphertext.Ciphertext = new(rlwe.Ciphertext)
	return ciphertext.Ciphertext.UnmarshalBinary(data)
}

// GetDataLen returns the length in bytes of the target Ciphertext.
func (ciphertext *Ciphertext) GetDataLen(WithMetaData bool) (dataLen int) {
	return ciphertext.Ciphertext.GetDataLen(WithMetaData)
}
