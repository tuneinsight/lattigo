package rlwe

import (
	"fmt"

	"github.com/tuneinsight/lattigo/v4/rlwe/ringqp"
)

// CiphertextQP is a generic type for RLWE ciphertexts in R_qp.
// It contains no MetaData.
type CiphertextQP struct {
	MetaData
	Value [2]ringqp.Poly
}

// NewCiphertextQP allocates a new CiphertextQP.
func NewCiphertextQP(params Parameters, levelQ, levelP int) CiphertextQP {
	ringQ := params.RingQ().AtLevel(levelQ)
	ringP := params.RingQ().AtLevel(levelP)

	return CiphertextQP{
		Value: [2]ringqp.Poly{
			{
				Q: ringQ.NewPoly(),
				P: ringP.NewPoly(),
			},
			{
				Q: ringQ.NewPoly(),
				P: ringP.NewPoly(),
			},
		},
		MetaData: MetaData{
			IsNTT: params.DefaultNTTFlag(),
		},
	}
}

// MarshalBinarySize returns the length in bytes of the target CiphertextQP.
func (ct *CiphertextQP) MarshalBinarySize() int {
	return ct.MetaData.MarshalBinarySize() + 2*ct.Value[0].MarshalBinarySize64()
}

// Encode64 encodes the target CiphertextQP on a byte array, using 8 bytes per coefficient.
// It returns the number of written bytes, and the corresponding error, if it occurred.
func (ct *CiphertextQP) Encode64(data []byte) (ptr int, err error) {

	if len(data) < ct.MarshalBinarySize() {
		return 0, fmt.Errorf("Encode64: len(data) is too small")
	}

	if ptr, err = ct.MetaData.Encode64(data); err != nil {
		return
	}

	var inc int

	if inc, err = ct.Value[0].Encode64(data[ptr:]); err != nil {
		return
	}

	ptr += inc

	if inc, err = ct.Value[1].Encode64(data[ptr:]); err != nil {
		return
	}

	ptr += inc

	return
}

// Decode64 decodes a slice of bytes in the target CiphertextQP and returns the number of bytes decoded.
// The method will first try to write on the buffer. If this step fails, either because the buffer isn't
// allocated or because it has the wrong size, the method will allocate the correct buffer.
// Assumes that each coefficient is encoded on 8 bytes.
func (ct *CiphertextQP) Decode64(data []byte) (ptr int, err error) {

	if ptr, err = ct.MetaData.Decode64(data); err != nil {
		return
	}

	var inc int

	if inc, err = ct.Value[0].Decode64(data[ptr:]); err != nil {
		return
	}

	ptr += inc

	if inc, err = ct.Value[1].Decode64(data[ptr:]); err != nil {
		return
	}

	ptr += inc

	return
}

// CopyNew returns a copy of the target CiphertextQP.
func (ct *CiphertextQP) CopyNew() *CiphertextQP {
	return &CiphertextQP{Value: [2]ringqp.Poly{ct.Value[0].CopyNew(), ct.Value[1].CopyNew()}, MetaData: ct.MetaData}
}
