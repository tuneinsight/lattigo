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

	if !params.isValid {
		panic("cannot NewCiphertext: parameters are invalid (check if the generation was done properly)")
	}

	ciphertext = &Ciphertext{&CkksElement{}}

	ciphertext.value = make([]*ring.Poly, degree+1)
	for i := uint64(0); i < degree+1; i++ {
		ciphertext.value[i] = ring.NewPoly(1<<params.LogN, level+1)
	}

	ciphertext.scale = scale
	ciphertext.isNTT = true

	return ciphertext
}

// NewCiphertextRandom generates a new uniformly distributed Ciphertext of degree, level and scale.
func NewCiphertextRandom(prng utils.PRNG, params *Parameters, degree, level uint64, scale float64) (ciphertext *Ciphertext) {

	if !params.isValid {
		panic("cannot NewCiphertextRandom: parameters are invalid (check if the generation was done properly)")
	}
	context, err := ring.NewContextWithParams(params.N, params.Qi[:level])
	if err != nil {
		panic(err)
	}
	sampler := ring.NewUniformSampler(prng, context)
	ciphertext = NewCiphertext(params, degree, level, scale)
	for i := uint64(0); i < degree+1; i++ {
		sampler.Read(ciphertext.value[i])
	}

	return ciphertext
}
