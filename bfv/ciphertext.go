package bfv

import (
	"github.com/ldsec/lattigo/v2/rlwe"
	"github.com/ldsec/lattigo/v2/utils"
)

// Ciphertext is a *ring.Poly array representing a polynomial of degree > 0 with coefficients in R_Q.
type Ciphertext struct {
	*rlwe.Element
}

// NewCiphertext creates a new ciphertext parameterized by degree, level and scale.
func NewCiphertext(params *Parameters, degree uint64) (ciphertext *Ciphertext) {
	return &Ciphertext{rlwe.NewElement(params.RLWEParameters(), degree)}
}

// NewCiphertextRandom generates a new uniformly distributed ciphertext of degree, level and scale.
func NewCiphertextRandom(prng utils.PRNG, params *Parameters, degree uint64) (ciphertext *Ciphertext) {
	ciphertext = &Ciphertext{rlwe.NewElement(params.RLWEParameters(), degree)}
	rlwe.PopulateElementRandom(prng, params.RLWEParameters(), (*rlwe.Element)(ciphertext.Element))
	return
}
