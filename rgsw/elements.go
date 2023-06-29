package rgsw

import (
	"github.com/tuneinsight/lattigo/v4/rlwe"
)

// Ciphertext is a generic type for RGSW ciphertext.
type Ciphertext struct {
	Value [2]rlwe.GadgetCiphertext
}

// LevelQ returns the level of the modulus Q of the target.
func (ct *Ciphertext) LevelQ() int {
	return ct.Value[0].LevelQ()
}

// LevelP returns the level of the modulus P of the target.
func (ct *Ciphertext) LevelP() int {
	return ct.Value[0].LevelP()
}

// NewCiphertext allocates a new RGSW ciphertext in the NTT domain.
func NewCiphertext(params rlwe.Parameters, levelQ, levelP, BaseTwoDecomposition int) (ct *Ciphertext) {
	return &Ciphertext{
		Value: [2]rlwe.GadgetCiphertext{
			*rlwe.NewGadgetCiphertext(params, 1, levelQ, levelP, BaseTwoDecomposition),
			*rlwe.NewGadgetCiphertext(params, 1, levelQ, levelP, BaseTwoDecomposition),
		},
	}
}

// Plaintext stores an RGSW plaintext value.
type Plaintext rlwe.GadgetPlaintext

// NewPlaintext creates a new RGSW plaintext from value, which can be either uint64, int64 or *ring.Poly.
// Plaintext is returned in the NTT and Mongtomery domain.
func NewPlaintext(params rlwe.Parameters, value interface{}, levelQ, levelP, BaseTwoDecomposition int) (pt *Plaintext) {
	return &Plaintext{Value: rlwe.NewGadgetPlaintext(params, value, levelQ, levelP, BaseTwoDecomposition).Value}
}
