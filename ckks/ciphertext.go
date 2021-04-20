package ckks

import (
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/utils"
)

// Ciphertext is *ring.Poly array representing a polynomial of degree > 0 with coefficients in R_Q.
type Ciphertext struct {
	*Element
}

// NewCiphertext creates a new Ciphertext parameterized by degree, level and scale.
func NewCiphertext(params Parameters, degree uint64, level uint64, scale float64) (ciphertext *Ciphertext) {

	ciphertext = &Ciphertext{NewElement(params, degree, level, scale)}
	ciphertext.IsNTT = true // TODO are Ciphertexts are NTT by default ?

	return ciphertext
}

// NewCiphertextRandom generates a new uniformly distributed Ciphertext of degree, level and scale.
func NewCiphertextRandom(prng utils.PRNG, params Parameters, degree, level uint64, scale float64) (ciphertext *Ciphertext) {

	ringQ, err := ring.NewRing(params.N(), params.Q()[:level+1])
	if err != nil {
		panic(err)
	}

	sampler := ring.NewUniformSampler(prng, ringQ)
	ciphertext = NewCiphertext(params, degree, level, scale)
	for i := uint64(0); i < degree+1; i++ {
		sampler.Read(ciphertext.Value[i])
	}

	return ciphertext
}

func (ct *Ciphertext) Copy(ctp *Ciphertext) {
	ct.Element.Copy(ctp.Element)
}

func (ct *Ciphertext) CopyNew() *Ciphertext {
	return &Ciphertext{ct.Element.CopyNew()}
}
