// Package gadget implements the R-LWE gadget ciphertexts. A gadget ciphertext is a matrix of ciphertexts encrypting plaintexts
// decomposed in the RNS and power of two basis.
package gadget

import (
	"github.com/tuneinsight/lattigo/v3/rlwe/ringqp"
)

// Ciphertext is a struct for storing an encrypted
// plaintext times the gadget power matrix.
type Ciphertext struct {
	Value [][][2]ringqp.Poly
}

// NewCiphertextNTT returns a new Ciphertext key with pre-allocated zero-value
func NewCiphertextNTT(levelQ, levelP, decompRNS, decompBIT int, ringQP ringqp.Ring) (ct *Ciphertext) {
	ct = new(Ciphertext)
	ct.Value = make([][][2]ringqp.Poly, decompRNS)
	for i := 0; i < decompRNS; i++ {
		ct.Value[i] = make([][2]ringqp.Poly, decompBIT)
		for j := 0; j < decompBIT; j++ {
			ct.Value[i][j][0] = ringQP.NewPolyLvl(levelQ, levelP)
			ct.Value[i][j][1] = ringQP.NewPolyLvl(levelQ, levelP)

			ct.Value[i][j][0].Q.IsNTT = true
			ct.Value[i][j][1].Q.IsNTT = true

			if levelP != -1 {
				ct.Value[i][j][0].P.IsNTT = true
				ct.Value[i][j][1].P.IsNTT = true
			}
		}
	}

	return ct
}

// LevelQ returns the level of the modulus Q of the target Ciphertext.
func (ct *Ciphertext) LevelQ() int {
	return ct.Value[0][0][0].Q.Level()
}

// LevelP returns the level of the modulus P of the target Ciphertext.
func (ct *Ciphertext) LevelP() int {
	if ct.Value[0][0][0].P != nil {
		return ct.Value[0][0][0].P.Level()
	}

	return -1
}

// Equals checks two Ciphertexts for equality.
func (ct *Ciphertext) Equals(other *Ciphertext) bool {
	if ct == other {
		return true
	}
	if (ct == nil) != (other == nil) {
		return false
	}
	if len(ct.Value) != len(other.Value) {
		return false
	}

	if len(ct.Value[0]) != len(other.Value[0]) {
		return false
	}

	for i := range ct.Value {
		for j, pol := range ct.Value[i] {
			if !pol[0].Equals(other.Value[i][j][0]) && !pol[1].Equals(other.Value[i][j][1]) {
				return false
			}
		}
	}
	return true
}

// CopyNew creates a deep copy of the receiver Ciphertext and returns it.
func (ct *Ciphertext) CopyNew() (ctCopy *Ciphertext) {
	if ct == nil || len(ct.Value) == 0 {
		return nil
	}
	Value := make([][][2]ringqp.Poly, len(ct.Value))
	for i := range ct.Value {
		for j, el := range ct.Value[i] {
			Value[i][j] = [2]ringqp.Poly{el[0].CopyNew(), el[1].CopyNew()}
		}
	}
	return &Ciphertext{Value: Value}
}

// GetDataLen returns the length in bytes of the target Ciphertext.
func (ct *Ciphertext) GetDataLen(WithMetadata bool) (dataLen int) {

	if WithMetadata {
		dataLen += 2
	}

	for i := range ct.Value {
		for _, el := range ct.Value[i] {
			dataLen += el[0].GetDataLen64(WithMetadata)
			dataLen += el[1].GetDataLen64(WithMetadata)
		}
	}

	return
}

// MarshalBinary encodes the target Ciphertext on a slice of bytes.
func (ct *Ciphertext) MarshalBinary() (data []byte, err error) {
	data = make([]byte, ct.GetDataLen(true))
	if _, err = ct.Encode(0, data); err != nil {
		return
	}

	return
}

// UnmarshalBinary decode a slice of bytes on the target Ciphertext.
func (ct *Ciphertext) UnmarshalBinary(data []byte) (err error) {
	if _, err = ct.Decode(data); err != nil {
		return
	}

	return
}

// Encode encodes the target ciphertext on a pre-allocated slice of bytes.
func (ct *Ciphertext) Encode(pointer int, data []byte) (int, error) {

	var err error
	var inc int

	data[pointer] = uint8(len(ct.Value))
	pointer++
	data[pointer] = uint8(len(ct.Value[0]))
	pointer++

	for i := range ct.Value {
		for _, el := range ct.Value[i] {

			if inc, err = el[0].WriteTo64(data[pointer:]); err != nil {
				return pointer, err
			}
			pointer += inc

			if inc, err = el[1].WriteTo64(data[pointer:]); err != nil {
				return pointer, err
			}
			pointer += inc
		}
	}

	return pointer, nil
}

// Decode decodes a slice of bytes on the target ciphertext.
func (ct *Ciphertext) Decode(data []byte) (pointer int, err error) {

	decompRNS := int(data[0])
	decompBIT := int(data[1])

	pointer = 2

	ct.Value = make([][][2]ringqp.Poly, decompRNS)

	var inc int

	for i := range ct.Value {

		ct.Value[i] = make([][2]ringqp.Poly, decompBIT)

		for j := range ct.Value[i] {

			if inc, err = ct.Value[i][j][0].DecodePoly64(data[pointer:]); err != nil {
				return
			}
			pointer += inc

			if inc, err = ct.Value[i][j][1].DecodePoly64(data[pointer:]); err != nil {
				return
			}
			pointer += inc
		}
	}

	return
}
