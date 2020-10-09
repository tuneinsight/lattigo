package bfv

import (
	"github.com/ldsec/lattigo/v2/ring"
)

// Plaintext is a Element with only one Poly.
type Plaintext struct {
	*Element
	value *ring.Poly
}

// NewPlaintext creates and allocates a new plaintext.
func NewPlaintext(params *Parameters) *Plaintext {

	plaintext := &Plaintext{NewElement(params, 0), nil}
	plaintext.value = plaintext.Element.value[0]
	plaintext.isNTT = false
	return plaintext
}
