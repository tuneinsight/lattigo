package ckks

import (
	"github.com/ldsec/lattigo/v2/ring"
)

// Plaintext is is a Element with only one Poly.
type Plaintext struct {
	*Element
	value *ring.Poly
}

// NewPlaintext creates a new Plaintext of level level and scale scale.
func NewPlaintext(params *Parameters, level int, scale float64) *Plaintext {

	plaintext := &Plaintext{Element: NewElement(params, 0, level, scale)}
	plaintext.value = plaintext.Element.Value[0]
	plaintext.Element.Element.IsNTT = true

	return plaintext
}
