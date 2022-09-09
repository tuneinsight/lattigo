package rlwe

import (
	"github.com/tuneinsight/lattigo/v3/rlwe/ringqp"
)

// PlaintextQP is a generic type for RLWE plaintext in R_qp.
type PlaintextQP struct {
	Value ringqp.Poly
}

func (el *PlaintextQP) El() CiphertextQP {
	return CiphertextQP{Value: []ringqp.Poly{el.Value}}
}

func (el *PlaintextQP) Level() (int, int) {
	return el.Value.Q.Level(), el.Value.P.Level()
}

func (el *PlaintextQP) Degree() int {
	return 0
}

// GetDataLen64 returns the number of bytes that encoding the target
// would require.
func (el *PlaintextQP) GetDataLen64(WithMetaData bool) (dataLen int) {
	return el.Value.GetDataLen64(WithMetaData)
}

// WriteTo64 encodes the coefficients of the target PlaintextQP on a slice of pre-allocated bytes
// and returns the pointer (index) of the cursor in the byte slice.
// Each coefficient is encoded on 8 bytes.
func (el *PlaintextQP) WriteTo64(data []byte) (ptr int, err error) {
	return el.Value.WriteTo64(data)
}

// Decode64 decodes the input bytes on the target PlaintextQP.
// and returns the pointer (index) of the cursor in the byte slice.
// Assumes that each coefficient is encode on 8 bytes.
func (el *PlaintextQP) Decode64(data []byte) (ptr int, err error) {
	return el.Value.DecodePoly64(data)
}

// CiphertextQP is a generic type for RLWE ciphertext in R_qp.
type CiphertextQP struct {
	Value []ringqp.Poly
}

func (el *CiphertextQP) Level() (int, int) {
	return el.Value[0].Q.Level(), el.Value[0].P.Level()
}

func (el *CiphertextQP) Degree() int {
	return len(el.Value) - 1
}

func (el *CiphertextQP) El() CiphertextQP {
	return *el
}

func (el *CiphertextQP) CopyNew() *CiphertextQP {
	copy := &CiphertextQP{}
	copy.Value = make([]ringqp.Poly, el.Degree()+1)
	for i := range copy.Value {
		copy.Value[i] = el.Value[i].CopyNew()
	}
	return copy
}

// GetDataLen64 returns the number of bytes that encoding the target
// would require.
func (el *CiphertextQP) GetDataLen64(WithMetaData bool) (dataLen int) {
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

// WriteTo64 encodes the coefficients of the target CiphertextQP on a slice of pre-allocated bytes
// and returns the pointer (index) of the cursor in the byte slice.
// Each coefficient is encoded on 8 bytes.
func (el *CiphertextQP) WriteTo64(data []byte) (ptr int, err error) {

	data[0] = uint8(el.Degree())
	ptr++

	var inc int
	for i := range el.Value {
		if inc, err = el.Value[i].WriteTo64(data[ptr:]); err != nil {
			return
		}

		ptr += inc
	}

	return
}

// Decode64 decodes the input bytes on the target CiphertextQP.
// and returns the pointer (index) of the cursor in the byte slice.
// Assumes that each coefficient is encode on 8 bytes.
func (el *CiphertextQP) Decode64(data []byte) (ptr int, err error) {

	if pad := int(uint8(data[0])) + 1 - len(el.Value); pad > 0 {
		el.Value = append(el.Value, make([]ringqp.Poly, pad)...)
	} else {
		el.Value = el.Value[:uint8(data[0])+1]
	}

	ptr++

	var inc int
	for i := range el.Value {

		if inc, err = el.Value[i].DecodePoly64(data[ptr:]); err != nil {
			return
		}

		ptr += inc
	}

	return
}
