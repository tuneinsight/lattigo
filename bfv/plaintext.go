package bfv

import (
	"github.com/ldsec/lattigo/v2/ring"
)

// Plaintext is a Element with only one Poly.
type Plaintext struct {
	*Element
	value *ring.Poly
}

// NewPlaintext creates and allocates a new plaintext in ZQ (multiple moduli of Q).
// This is a generic all-purpose type of plaintext. It will work with all operations,
// but will likely not result in an optimized circuit.
// See NewPlaintextZQ, NewPlaintextZT and NewPlaintextMul for additional informations.
// The plaintext will be in ZQ and scaled by Q/t.
// Faster encryption
// Slower encoding and larger plaintext size
// Fast ct + pt
// Very slow ct * pt
func NewPlaintext(params *Parameters) *Plaintext {
	plaintext := &Plaintext{newElePTZQ(params), nil}
	plaintext.value = plaintext.Element.value[0]
	return plaintext
}

// NewPlaintextZQ creates and allocates a new plaintext in ZQ (multiple moduli of Q).
// The plaintext will be in ZQ and scaled by Q/t.
// Faster encryption
// Slower encoding and larger plaintext size
// Fast ct + pt
// Very slow ct * pt
func NewPlaintextZQ(params *Parameters) *Plaintext {
	plaintext := &Plaintext{newElePTZQ(params), nil}
	plaintext.value = plaintext.Element.value[0]
	return plaintext
}

// NewPlaintextZT creates and allocates a new plaintext in ZT (single modulus T).
// The plaintext will be in Zt.
// Faster encoding + smaller plaintext size
// Slower encryption
// Slow ct + pt
// Slow ct * pt
func NewPlaintextZT(params *Parameters) *Plaintext {

	plaintext := &Plaintext{newElePTZT(params), nil}
	plaintext.value = plaintext.Element.value[0]
	return plaintext
}

// NewPlaintextMul creates and allocates a new plaintext optimized for ciphertext x plaintext multiplication.
// The plaintext will be in the NTT and Montgomery domain of ZQ and not scaled by Q/t.
// Cannot do ct + pt
// Fast ct * pt
func NewPlaintextMul(params *Parameters) *Plaintext {
	plaintext := &Plaintext{newElePTMul(params), nil}
	plaintext.value = plaintext.Element.value[0]
	return plaintext
}
