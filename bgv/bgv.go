// Package bgv implements a unified RNS-accelerated version of the Fan-Vercauteren version of the Brakerski's scale invariant homomorphic encryption scheme (BFV) and Brakerski-Gentry-Vaikuntanathan (BGV) homomorphic encryption scheme. It provides modular arithmetic over the integers.
package bgv

import (
	"github.com/tuneinsight/lattigo/v4/rlwe"
)

// NewPlaintext allocates a new rlwe.Plaintext.
//
// inputs:
// - params: bfv.Parameters
// - level: the level of the plaintext
//
// output: a newly allocated rlwe.Plaintext at the specified level.
func NewPlaintext(params Parameters, level int) (pt *rlwe.Plaintext) {
	return rlwe.NewPlaintext(params, level)
}

// NewCiphertext allocates a new rlwe.Ciphertext.
//
// inputs:
// - params: bfv.Parameters
// - degree: the degree of the ciphertext
// - level: the level of the Ciphertext
//
// output: a newly allocated rlwe.Ciphertext of the specified degree and level.
func NewCiphertext(params Parameters, degree, level int) (ct *rlwe.Ciphertext) {
	return rlwe.NewCiphertext(params, degree, level)
}

// NewEncryptor instantiates a new rlwe.Encryptor.
//
// inputs:
// - params: bfv.Parameters
// - key: *rlwe.SecretKey or *rlwe.PublicKey
//
// output: an rlwe.Encryptor instantiated with the provided key.
func NewEncryptor(params Parameters, key interface{}) rlwe.Encryptor {
	return rlwe.NewEncryptor(params, key)
}

// NewPRNGEncryptor instantiates a new rlwe.PRNGEncryptor.
//
// inputs:
// - params: bfv.Parameters
// - key: *rlwe.SecretKey
//
// output: an rlwe.PRNGEncryptor instantiated with the provided key.
func NewPRNGEncryptor(params Parameters, key *rlwe.SecretKey) rlwe.PRNGEncryptor {
	return rlwe.NewPRNGEncryptor(params, key)
}

// NewDecryptor instantiates a new rlwe.Decryptor.
//
// inputs:
// - params: bfv.Parameters
// - key: *rlwe.SecretKey
//
// output: an rlwe.Decryptor instantiated with the provided key.
func NewDecryptor(params Parameters, key *rlwe.SecretKey) rlwe.Decryptor {
	return rlwe.NewDecryptor(params, key)
}

// NewKeyGenerator instantiates a new rlwe.KeyGenerator.
//
// inputs:
// - params: bfv.Parameters
//
// output: an rlwe.KeyGenerator.
func NewKeyGenerator(params Parameters) *rlwe.KeyGenerator {
	return rlwe.NewKeyGenerator(params)
}
