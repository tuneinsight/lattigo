// Package bfv implements a RNS-accelerated Fan-Vercauteren version of Brakerski's scale invariant homomorphic encryption scheme. It provides modular arithmetic over the integers.
package bfv

import (
	"github.com/tuneinsight/lattigo/v4/rlwe"
)

// NewPlaintext creates and allocates a new plaintext in RingQ (multiple moduli of Q).
// The plaintext will be in RingQ and scaled by Q/t.
// Slower encoding and larger plaintext size
func NewPlaintext(params Parameters, level int) (pt *rlwe.Plaintext) {
	return rlwe.NewPlaintext(params.Parameters, level)
}

// PlaintextRingT represents a plaintext element in R_t.
// This is the most compact representation of a plaintext, but performing operations have the extra-cost of performing
// the scaling up by Q/t. See bfv/encoder.go for more information on plaintext types.
type PlaintextRingT struct {
	*rlwe.Plaintext
}

// PlaintextMul represents a plaintext element in R_q, in NTT and Montgomery form, but without scale up by Q/t.
// A PlaintextMul is a special-purpose plaintext for efficient Ciphertext x Plaintext multiplication. However,
// other operations on plaintexts are not supported. See bfv/encoder.go for more information on plaintext types.
type PlaintextMul struct {
	*rlwe.Plaintext
}

// NewPlaintextRingT creates and allocates a new plaintext in RingT (single modulus T).
// The plaintext will be in RingT.
func NewPlaintextRingT(params Parameters) *PlaintextRingT {
	return &PlaintextRingT{rlwe.NewPlaintext(params.Parameters, 0)}
}

// NewPlaintextMul creates and allocates a new plaintext optimized for Ciphertext x Plaintext multiplication.
// The Plaintext is allocated with level+1 moduli.
// The Plaintext will be in the NTT and Montgomery domain of RingQ and not scaled by Q/t.
func NewPlaintextMul(params Parameters, level int) *PlaintextMul {
	return &PlaintextMul{rlwe.NewPlaintext(params.Parameters, level)}
}

func NewCiphertext(params Parameters, degree, level int) (ct *rlwe.Ciphertext) {
	return rlwe.NewCiphertext(params.Parameters, degree, level)
}

func NewEncryptor(params Parameters, key interface{}) rlwe.Encryptor {
	return rlwe.NewEncryptor(params.Parameters, key)
}

func NewDecryptor(params Parameters, key *rlwe.SecretKey) rlwe.Decryptor {
	return rlwe.NewDecryptor(params.Parameters, key)
}

func NewKeyGenerator(params Parameters) rlwe.KeyGenerator {
	return rlwe.NewKeyGenerator(params.Parameters)
}

func NewPRNGEncryptor(params Parameters, key *rlwe.SecretKey) rlwe.PRNGEncryptor {
	return rlwe.NewPRNGEncryptor(params.Parameters, key)
}
