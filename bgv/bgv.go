// Package bgv implements a RNS-accelerated BGV homomorphic encryption scheme. It provides modular arithmetic over the integers.
package bgv

import (
	"github.com/tuneinsight/lattigo/v4/rlwe"
)

func NewPlaintext(params Parameters, level int) (pt *rlwe.Plaintext) {
	return rlwe.NewPlaintext(params.Parameters, level)
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
