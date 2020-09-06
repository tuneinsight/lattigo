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
		ciphertext.value[i] = ring.NewPoly(params.N(), level+1)
	}

	ciphertext.scale = scale
	ciphertext.isNTT = true

	return ciphertext
}

// NewCiphertext creates a new Ciphertext parameterized by degree, level and scale.
func NewCiphertextQP(params *Parameters, degree uint64, level uint64, scale float64) (ciphertextQ, ciphertextP *Ciphertext) {

	ciphertextQ = &Ciphertext{&CkksElement{}}

	ciphertextQ.value = make([]*ring.Poly, degree+1)
	for i := uint64(0); i < degree+1; i++ {
		ciphertextQ.value[i] = ring.NewPoly(params.N(), level+1)
	}

	ciphertextQ.scale = scale
	ciphertextQ.isNTT = true

	ciphertextP = &Ciphertext{&CkksElement{}}

	ciphertextP.value = make([]*ring.Poly, degree+1)
	for i := uint64(0); i < degree+1; i++ {
		ciphertextP.value[i] = ring.NewPoly(params.N(), params.PiCount())
	}

	ciphertextP.scale = scale
	ciphertextP.isNTT = true

	return ciphertextQ, ciphertextP
}

// NewCiphertextRandom generates a new uniformly distributed Ciphertext of degree, level and scale.
func NewCiphertextRandom(prng utils.PRNG, params *Parameters, degree, level uint64, scale float64) (ciphertext *Ciphertext) {

	context, err := ring.NewContextWithParams(params.N(), params.qi[:level+1])
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
