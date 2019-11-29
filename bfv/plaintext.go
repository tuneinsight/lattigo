package bfv

import (
	"github.com/ldsec/lattigo/ring"
)

// Plaintext is a BigPoly of degree 0.
type Plaintext struct {
	*bfvElement
	value *ring.Poly
}

// NewPlaintext creates a new plaintext from the target context.
func NewPlaintext(params *Parameters) *Plaintext {
	plaintext := &Plaintext{&bfvElement{}, nil}

	plaintext.bfvElement.value = []*ring.Poly{ring.NewPoly(uint64(1<<params.LogN), uint64(len(params.Q1)))}
	plaintext.value = plaintext.bfvElement.value[0]
	plaintext.isNTT = false

	return plaintext
}
