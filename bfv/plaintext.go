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
func NewPlaintextZQ(params *Parameters) *Plaintext {

	plaintext := &Plaintext{NewElement(params, 0, true), nil}
	plaintext.value = plaintext.Element.value[0]
	plaintext.isNTT = false
	return plaintext
}

// NewPlaintext creates and allocates a new plaintext.
func NewPlaintextZT(params *Parameters) *Plaintext {

	plaintext := &Plaintext{NewElement(params, 0, false), nil}
	plaintext.value = plaintext.Element.value[0]
	plaintext.isNTT = false
	return plaintext
}
