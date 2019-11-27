package bfv

import (
	"github.com/ldsec/lattigo/ring"
)

// Ciphertext is a *ring.Poly array representing a polynomial of degree > 0 where coefficients are in R_Q.
type Ciphertext struct {
	*bfvElement
}

// NewCiphertext creates a new Ciphertext structure.
func NewCiphertextStruct() (ciphertext *Ciphertext) {
	return &Ciphertext{&bfvElement{}}
}

func NewRingContext(params *Parameters) *ring.Context {
	ringCtx := ring.NewContext()
	ringCtx.SetParameters(params.N, params.Qi)

	err := ringCtx.GenNTTParams()
	if err != nil {
		panic(err)
	}

	return ringCtx
}

// NewCiphertext creates a new empty ciphertext of degree degree.
func NewCiphertext(degree uint64, ringCtx *ring.Context) *Ciphertext {
	ciphertext := &Ciphertext{&bfvElement{}}
	ciphertext.value = make([]*ring.Poly, degree+1)
	for i := uint64(0); i < degree+1; i++ {
		ciphertext.value[i] = ringCtx.NewPoly()
	}
	ciphertext.isNTT = false

	return ciphertext
}

// NewRandomCiphertext creates a new ciphertext with uniform coefficients.
func (context *Context) NewRandomCiphertext(degree uint64) *Ciphertext {
	ciphertext := &Ciphertext{&bfvElement{}}
	ciphertext.value = make([]*ring.Poly, degree+1)
	for i := uint64(0); i < degree+1; i++ {
		ciphertext.value[i] = context.contextQ.NewUniformPoly()
	}
	ciphertext.isNTT = false

	return ciphertext
}
