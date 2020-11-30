package bfv

import (
	"github.com/ldsec/lattigo/v2/ring"
)

// Plaintext is a Element with only one Poly.
type Plaintext struct {
	*Element
	value *ring.Poly
}

// NewPlaintextZQ creates and allocates a new plaintext in ZQ (multiple moduli of Q).
// The plaintext will be in ZQ and scaled by Q/t.
func NewPlaintextZQ(params *Parameters) *Plaintext {
	plaintext := &Plaintext{NewElement(params, 0, true), nil}
	plaintext.value = plaintext.Element.value[0]
	plaintext.isNTT = false
	plaintext.scaled = true
	return plaintext
}

// NewPlaintextZT creates and allocates a new plaintext in ZT (single modulus T).
// The plaintext will be in Zt.
func NewPlaintextZT(params *Parameters) *Plaintext {

	plaintext := &Plaintext{NewElement(params, 0, false), nil}
	plaintext.value = plaintext.Element.value[0]
	plaintext.isNTT = false
	plaintext.scaled = false
	return plaintext
}

// NewPlaintextMul creates and allocates a new plaintext optimized for ciphertext x plaintext multiplication.
// The plaintext will be in the NTT and Montgomery domain of ZQ and not scaled by Q/t.
func NewPlaintextMul(params *Parameters) *Plaintext {
	plaintext := &Plaintext{NewElement(params, 0, true), nil}
	plaintext.value = plaintext.Element.value[0]
	plaintext.isNTT = true
	plaintext.scaled = false
	return plaintext
}
