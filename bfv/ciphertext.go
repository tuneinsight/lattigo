package bfv

import (
	"github.com/ldsec/lattigo/ring"
)

// Ciphertext is a *ring.Poly array representing a polynomial of degree > 0 where coefficients are in R_Q.
type Ciphertext struct {
	*bfvElement
}

// NewCiphertext creates a new empty ciphertext of degree degree.
func (bfvcontext *BfvContext) NewCiphertext(degree uint64) *Ciphertext {
	ciphertext := &Ciphertext{&bfvElement{}}
	ciphertext.value = make([]*ring.Poly, degree+1)
	for i := uint64(0); i < degree+1; i++ {
		ciphertext.value[i] = bfvcontext.contextQ.NewPoly()
	}
	ciphertext.isNTT = false

	return ciphertext
}

// NewRandomCiphertext creates a new ciphertext with uniform coefficients.
func (bfvcontext *BfvContext) NewRandomCiphertext(degree uint64) *Ciphertext {
	ciphertext := &Ciphertext{&bfvElement{}}
	ciphertext.value = make([]*ring.Poly, degree+1)
	for i := uint64(0); i < degree+1; i++ {
		ciphertext.value[i] = bfvcontext.contextQ.NewUniformPoly()
	}
	ciphertext.isNTT = false

	return ciphertext
}
