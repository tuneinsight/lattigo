package bfv

import (
	"github.com/ldsec/lattigo/ring"
)

// Plaintext is a bfvElement with only one Poly.
type Plaintext struct {
	*bfvElement
	value *ring.Poly
}

// NewPlaintext creates and allocates a new plaintext.
func NewPlaintext(params *Parameters) *Plaintext {

	plaintext := &Plaintext{newBfvElement(params, 0), nil}
	plaintext.value = plaintext.bfvElement.value[0]
	plaintext.isNTT = false
	return plaintext
}
