// Package ckks implements a RNS-accelerated version of the Homomorphic Encryption for Arithmetic for Approximate Numbers
// (HEAAN, a.k.a. CKKS) scheme. It provides approximate arithmetic over the complex numbers.package ckks
package ckks

import (
	"github.com/tuneinsight/lattigo/v4/rlwe"
)

// NewPlaintext allocates a new rlwe.Plaintext.
//
// inputs:
// - params: an rlwe.ParametersInterface interface
// - level: the level of the plaintext
//
// output: a newly allocated rlwe.Plaintext at the specified level.
//
// Note: the user can update the field `MetaData` to set a specific scaling factor,
// plaintext dimensions (if applicable) or encoding domain, before encoding values
// on the created plaintext.
func NewPlaintext(params rlwe.ParametersInterface, level int) (pt *rlwe.Plaintext) {
	return rlwe.NewPlaintext(params, level)
}

// NewCiphertext allocates a new rlwe.Ciphertext.
//
// inputs:
// - params: an rlwe.ParametersInterface interface
// - degree: the degree of the ciphertext
// - level: the level of the Ciphertext
//
// output: a newly allocated rlwe.Ciphertext of the specified degree and level.
func NewCiphertext(params rlwe.ParametersInterface, degree, level int) (ct *rlwe.Ciphertext) {
	return rlwe.NewCiphertext(params, degree, level)
}

// NewEncryptor instantiates a new rlwe.Encryptor.
//
// inputs:
// - params: an rlwe.ParametersInterface interface
// - key: *rlwe.SecretKey or *rlwe.PublicKey
//
// output: an rlwe.Encryptor instantiated with the provided key.
func NewEncryptor[T *rlwe.SecretKey | *rlwe.PublicKey](params rlwe.ParametersInterface, key T) (rlwe.EncryptorInterface, error) {
	return rlwe.NewEncryptor(params, key)
}

// NewPRNGEncryptor instantiates a new rlwe.PRNGEncryptor.
//
// inputs:
// - params: an rlwe.ParametersInterface interface
// - key: *rlwe.SecretKey
//
// output: an rlwe.PRNGEncryptor instantiated with the provided key.
func NewPRNGEncryptor(params rlwe.ParametersInterface, key *rlwe.SecretKey) (rlwe.PRNGEncryptorInterface, error) {
	return rlwe.NewPRNGEncryptor(params, key)
}

// NewDecryptor instantiates a new rlwe.Decryptor.
//
// inputs:
// - params: an rlwe.ParametersInterface interface
// - key: *rlwe.SecretKey
//
// output: an rlwe.Decryptor instantiated with the provided key.
func NewDecryptor(params rlwe.ParametersInterface, key *rlwe.SecretKey) (*rlwe.Decryptor, error) {
	return rlwe.NewDecryptor(params, key)
}

// NewKeyGenerator instantiates a new rlwe.KeyGenerator.
//
// inputs:
// - params: an rlwe.ParametersInterface interface
//
// output: an rlwe.KeyGenerator.
func NewKeyGenerator(params rlwe.ParametersInterface) *rlwe.KeyGenerator {
	return rlwe.NewKeyGenerator(params)
}
