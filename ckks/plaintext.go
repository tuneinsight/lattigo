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
func (ckkscontext *CkksContext) NewPlaintext(level uint64, scale float64) *Plaintext {
	plaintext := &Plaintext{&ckksElement{}, nil}

	plaintext.ckksElement.value = []*ring.Poly{ckkscontext.contextQ.NewPolyLvl(level)}
	plaintext.value = plaintext.ckksElement.value[0]

	plaintext.scale = scale
	plaintext.currentModulus = ring.Copy(ckkscontext.bigintChain[level])
	plaintext.isNTT = true

	return plaintext
}