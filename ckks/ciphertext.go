package ckks

import (
	"github.com/ldsec/lattigo/ring"
	"github.com/ldsec/lattigo/utils"
)

// Ciphertext is *ring.Poly array representing a polynomial of degree > 0 with coefficients in R_Q.
type Ciphertext struct {
	*CkksElement
}

// NewCiphertext creates a new Ciphertext parameterized by degree, level and scale.
func NewCiphertext(params *Parameters, degree uint64, level uint64, scale float64) (ciphertext *Ciphertext) {

	ciphertext = &Ciphertext{&CkksElement{}}

	ciphertext.value = make([]*ring.Poly, degree+1)
	for i := uint64(0); i < degree+1; i++ {
		ciphertext.value[i] = ring.NewPoly(params.n, level+1)
	}

	ciphertext.scale = scale
	ciphertext.isNTT = true

	return ciphertext
}

// NewCiphertextRandom generates a new uniformly distributed Ciphertext of degree, level and scale.
func NewCiphertextRandom(prng utils.PRNG, params *Parameters, degree, level uint64, scale float64) (ciphertext *Ciphertext) {

	context, _ := ring.NewRing(params.n, params.qi[:level])

	sampler := ring.NewUniformSampler(prng, context)
	ciphertext = NewCiphertext(params, degree, level, scale)
	for i := uint64(0); i < degree+1; i++ {
		sampler.Read(ciphertext.value[i])
	}

	return ciphertext
}
