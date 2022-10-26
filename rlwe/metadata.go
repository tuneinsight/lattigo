package rlwe

import (
	"encoding/binary"
	"fmt"
)

// MetaData is a struct storing metadata.
type MetaData struct {
	Scale
	IsNTT        bool
	IsMontgomery bool
	ErrorScale   uint64
}

// Equal returns true if two MetaData structs are identical.
func (m *MetaData) Equal(other MetaData) (res bool) {
	res = m.Scale.Cmp(other.Scale) == 0
	res = res && m.IsNTT == other.IsNTT
	res = res && m.IsMontgomery == other.IsMontgomery
	return res && m.ErrorScale == other.ErrorScale
}

// MarshalBinarySize returns the length in bytes of the target MetaData.
func (m *MetaData) MarshalBinarySize() int {
	return 10 + m.Scale.MarshalBinarySize()
}

// MarshalBinary encodes a MetaData on a byte slice.
func (m *MetaData) MarshalBinary() (data []byte, err error) {
	data = make([]byte, m.MarshalBinarySize())
	_, err = m.Encode64(data)
	return
}

// UnmarshalBinary decodes a previously marshaled MetaData on the target MetaData.
func (m *MetaData) UnmarshalBinary(data []byte) (err error) {
	_, err = m.Decode64(data)
	return
}

// Encode64 encodes the target ciphertext on a byte array, using 8 bytes per coefficient.
// It returns the number of written bytes, and the corresponding error, if it occurred.
func (m *MetaData) Encode64(data []byte) (ptr int, err error) {

	if len(data) < m.MarshalBinarySize() {
		return 0, fmt.Errorf("Encode64: len(data) is too small")
	}

	if err = m.Scale.Encode(data[ptr:]); err != nil {
		return
	}

	ptr += m.Scale.MarshalBinarySize()

	if m.IsNTT {
		data[ptr] = 1
	}

	ptr++

	if m.IsMontgomery {
		data[ptr] = 1
	}

	ptr++

	binary.LittleEndian.PutUint64(data[ptr:], m.ErrorScale)
	ptr += 8

	return
}

// Decode64 decodes a slice of bytes in the target MetaData and returns the number of bytes decoded.
// The method will first try to write on the buffer. If this step fails, either because the buffer isn't
// allocated or because it is of the wrong size, the method will allocate the correct buffer.
// Assumes that each coefficient is encoded on 8 bytes.
func (m *MetaData) Decode64(data []byte) (ptr int, err error) {

	if len(data) < m.MarshalBinarySize() {
		return 0, fmt.Errorf("Decode64: len(data) is too small")
	}

	if err = m.Scale.Decode(data[ptr:]); err != nil {
		return
	}

	ptr += m.Scale.MarshalBinarySize()

	m.IsNTT = data[ptr] == 1
	ptr++

	m.IsMontgomery = data[ptr] == 1
	ptr++

	m.ErrorScale = binary.LittleEndian.Uint64(data[ptr:])
	ptr += 8

	return
}
