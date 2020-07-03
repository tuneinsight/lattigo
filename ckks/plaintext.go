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

	if !params.isValid {
		panic("cannot NewPlaintext: parameters are invalid (check if the generation was done properly)")
	}

	plaintext := &Plaintext{&CkksElement{}, nil}

	plaintext.CkksElement.value = []*ring.Poly{ring.NewPoly(1<<params.LogN, level+1)}

	plaintext.value = plaintext.CkksElement.value[0]

	plaintext.scale = scale
	plaintext.isNTT = true

	return plaintext
}

// NewPlaintext creates a new Plaintext of level level and scale scale.
func NewPlaintextQP(params *Parameters, level uint64, scale float64) (*Plaintext, *Plaintext) {

	if !params.isValid {
		panic("cannot NewPlaintext: parameters are invalid (check if the generation was done properly)")
	}

	plaintextQ := &Plaintext{&CkksElement{}, nil}
	plaintextQ.CkksElement.value = []*ring.Poly{ring.NewPoly(1<<params.LogN, level+1)}
	plaintextQ.value = plaintextQ.CkksElement.value[0]
	plaintextQ.scale = scale
	plaintextQ.isNTT = true

	plaintextP := &Plaintext{&CkksElement{}, nil}
	plaintextP.CkksElement.value = []*ring.Poly{ring.NewPoly(1<<params.LogN, uint64(len(params.Pi)))}
	plaintextP.value = plaintextP.CkksElement.value[0]
	plaintextP.scale = scale
	plaintextP.isNTT = true

	return plaintextQ, plaintextP
}
