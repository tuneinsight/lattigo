package bfv

import "github.com/ldsec/lattigo/utils"

// Ciphertext is a *ring.Poly array representing a polynomial of degree > 0 with coefficients in R_Q.
type Ciphertext struct {
	*bfvElement
}

// NewCiphertext creates a new ciphertext parameterized by degree, level and scale.
func NewCiphertext(params *Parameters, degree uint64) (ciphertext *Ciphertext) {

	if !params.isValid {
		panic("cannot NewCiphertext: params not valid (check if they were generated properly)")
	}

	return &Ciphertext{newBfvElement(params, degree)}
}

// NewCiphertextRandom generates a new uniformly distributed ciphertext of degree, level and scale.
func NewCiphertextRandom(prng utils.PRNG, params *Parameters, degree uint64) (ciphertext *Ciphertext) {

	if !params.isValid {
		panic("cannot NewCiphertextRandom: params not valid (check if they were generated properly)")
	}

	return &Ciphertext{newBfvElementRandom(prng, params, degree)}
}
