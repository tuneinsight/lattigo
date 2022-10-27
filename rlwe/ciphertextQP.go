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

func (ct *CiphertextQP) MarshalBinarySize() int {
	return ct.MetaData.MarshalBinarySize() + 2*ct.Value[0].MarshalBinarySize64()
}

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
