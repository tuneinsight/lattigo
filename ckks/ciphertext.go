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
		ciphertext.value[i] = ckkscontext.contextQ.NewPolyLvl(level)
	}

	ciphertext.scale = scale
	ciphertext.currentModulus = ring.Copy(ckkscontext.bigintChain[level])
	ciphertext.isNTT = true

	return ciphertext
}

// NewRandoMCiphertext generates a new uniformely distributed ciphertext of degree, level and scale.
func (ckkscontext *CkksContext) NewRandomCiphertext(degree, level uint64, scale float64) (ciphertext *Ciphertext) {
	ciphertext = &Ciphertext{&ckksElement{}}

	ciphertext.value = make([]*ring.Poly, degree+1)
	for i := uint64(0); i < degree+1; i++ {
		ciphertext.value[i] = ckkscontext.contextQ.NewUniformPolyLvl(level)
	}

	ciphertext.scale = scale
	ciphertext.currentModulus = ring.Copy(ckkscontext.bigintChain[level])
	ciphertext.isNTT = true

	return ciphertext
}
