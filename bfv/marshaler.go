package bfv

import (
	"github.com/ldsec/lattigo/v2/ring"
)

// MarshalBinary encodes a Ciphertext in a byte slice.
func (ciphertext *Ciphertext) MarshalBinary() (data []byte, err error) {

	data = make([]byte, ciphertext.GetDataLen(true))

	data[0] = uint8(len(ciphertext.value))

	var pointer, inc int

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

	var pointer, inc int
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
func (ciphertext *Ciphertext) GetDataLen(WithMetaData bool) (dataLen int) {
	if WithMetaData {
		dataLen++
	}

	for _, el := range ciphertext.value {
		dataLen += el.GetDataLen(WithMetaData)
	}

	return dataLen
}
