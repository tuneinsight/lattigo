package bfv

import (
	"github.com/ldsec/lattigo/ring"
)

// Ciphertext is a *ring.Poly array representing a polynomial of degree > 0 where coefficients are in R_Q.
type Ciphertext struct {
	*bfvElement
}

// NewCiphertext creates a new ciphertext parameterized by degree, level and scale.
func NewCiphertext(params *Parameters, degree uint64) (ciphertext *Ciphertext) {
	ciphertext = &Ciphertext{&bfvElement{}}

	ciphertext.value = make([]*ring.Poly, degree+1)
	for i := uint64(0); i < degree+1; i++ {
		ciphertext.value[i] = ring.NewPoly(1<<params.LogN, uint64(len(params.Q1)))
	}

	ciphertext.isNTT = true

	return ciphertext
}

// NewCiphertextRandom generates a new uniformely distributed ciphertext of degree, level and scale.
func NewCiphertextRandom(params *Parameters, degree uint64) (ciphertext *Ciphertext) {
	ciphertext = &Ciphertext{&bfvElement{}}

	ciphertext.value = make([]*ring.Poly, degree+1)
	for i := uint64(0); i < degree+1; i++ {
		ciphertext.value[i] = ring.NewPolyUniform(1<<params.LogN, uint64(len(params.Q1)))
	}

	ciphertext.isNTT = true

	return ciphertext
}
