package ckks

import (
	"github.com/ldsec/lattigo/ring"
)

// Ciphertext is a BigPoly of degree > 0.
type Ciphertext struct {
	*ckksElement
}

// NewCiphertext creates a new ciphertext parameterized by degree, level and scale.
func (ckkscontext *CkksContext) NewCiphertext(degree uint64, level uint64, scale float64) *Ciphertext {
	ciphertext := &Ciphertext{&ckksElement{}}

	ciphertext.value = make([]*ring.Poly, degree+1)
	for i := uint64(0); i < degree+1; i++ {
		ciphertext.value[i] = ckkscontext.contextLevel[level].NewPoly()
	}

	ciphertext.scale = scale
	ciphertext.currentModulus = ring.Copy(ckkscontext.contextLevel[level].ModulusBigint)
	ciphertext.isNTT = true

	return ciphertext
}

// NewRandoMCiphertext generates a new uniformely distributed ciphertext of degree, level and scale.
func (ckkscontext *CkksContext) NewRandomCiphertext(degree, level uint64, scale float64) (ciphertext *Ciphertext) {
	ciphertext = &Ciphertext{&ckksElement{}}

	ciphertext.value = make([]*ring.Poly, degree+1)
	for i := uint64(0); i < degree+1; i++ {
		ciphertext.value[i] = ckkscontext.contextLevel[level].NewUniformPoly()
	}

	ciphertext.scale = scale
	ciphertext.currentModulus = ring.Copy(ckkscontext.contextLevel[level].ModulusBigint)
	ciphertext.isNTT = true

	return ciphertext
}
