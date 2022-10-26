package bgv

import (
	"github.com/tuneinsight/lattigo/v4/rlwe"
)

// NewCiphertext creates a new Ciphertext parameterized by degree and level.
func NewCiphertext(params Parameters, degree, level int) (ct *rlwe.Ciphertext) {
	return rlwe.NewCiphertextNTT(params.Parameters, degree, level)
}

// NewPlaintext creates a new Plaintext parameterized by degree and level.
func NewPlaintext(params Parameters, level int) (pt *rlwe.Plaintext) {
	return rlwe.NewPlaintextNTT(params.Parameters, level)
}
