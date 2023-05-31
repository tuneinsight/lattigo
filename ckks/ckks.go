// Package ckks implements a RNS-accelerated version of the Homomorphic Encryption for Arithmetic for Approximate Numbers
// (HEAAN, a.k.a. CKKS) scheme. It provides approximate arithmetic over the complex numbers.package ckks
package ckks

import (
	"github.com/tuneinsight/lattigo/v4/rlwe"
)

func NewPlaintext(params rlwe.ParametersInterface, level int) (pt *rlwe.Plaintext) {
	return rlwe.NewPlaintext(params, level)
}

func NewCiphertext(params rlwe.ParametersInterface, degree, level int) (ct *rlwe.Ciphertext) {
	return rlwe.NewCiphertext(params, degree, level)
}

func NewEncryptor(params rlwe.ParametersInterface, key interface{}) rlwe.Encryptor {
	return rlwe.NewEncryptor(params, key)
}

func NewDecryptor(params rlwe.ParametersInterface, key *rlwe.SecretKey) rlwe.Decryptor {
	return rlwe.NewDecryptor(params, key)
}

func NewKeyGenerator(params rlwe.ParametersInterface) *rlwe.KeyGenerator {
	return rlwe.NewKeyGenerator(params)
}

func NewPRNGEncryptor(params rlwe.ParametersInterface, key *rlwe.SecretKey) rlwe.PRNGEncryptor {
	return rlwe.NewPRNGEncryptor(params, key)
}
