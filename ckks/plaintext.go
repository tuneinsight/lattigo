package ckks

import (
	"github.com/ldsec/lattigo/ring"
)

// Plaintext is is a CkksElement with only one Poly.
type Plaintext struct {
	*CkksElement
	value *ring.Poly
}

// NewPlaintext creates a new Plaintext of level level and scale scale.
func NewPlaintext(params *Parameters, level uint64, scale float64) *Plaintext {

	plaintext := &Plaintext{&CkksElement{}, nil}

	plaintext.CkksElement.value = []*ring.Poly{ring.NewPoly(params.N(), level+1)}

	plaintext.value = plaintext.CkksElement.value[0]

	plaintext.scale = scale
	plaintext.isNTT = true

	return plaintext
}
