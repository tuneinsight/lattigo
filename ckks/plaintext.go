package ckks

import (
	"github.com/ldsec/lattigo/ring"
)

// Plaintext is BigPoly of degree 0.
type Plaintext struct {
	*ckksElement
	value *ring.Poly
}

// NewPlaintext creates a new plaintext of level level and scale scale.
func NewPlaintext(params *Parameters, level uint64, scale float64) *Plaintext {
	plaintext := &Plaintext{&ckksElement{}, nil}

	plaintext.ckksElement.value = []*ring.Poly{ring.NewPoly(1<<params.LogN, level+1)}

	plaintext.value = plaintext.ckksElement.value[0]

	plaintext.scale = scale
	plaintext.isNTT = true

	return plaintext
}
