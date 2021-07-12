package ckks

import (
	"github.com/ldsec/lattigo/v2/rlwe"
)

// Decryptor is an interface wrapping a rlwe.Decryptor.
type Decryptor interface {
	DecryptNew(ciphertext *Ciphertext) (plaintext *Plaintext)
	Decrypt(ciphertext *Ciphertext, plaintext *Plaintext)
}

type decryptor struct {
	rlwe.Decryptor
}

// NewDecryptor instantiates a new rlwe.Decryptor wrapped for CKKS ciphertexts and plaintexts.
func NewDecryptor(params Parameters, sk *rlwe.SecretKey) Decryptor {
	return &decryptor{Decryptor: rlwe.NewDecryptor(params.Parameters, sk)}
}

// DecryptNew calls rlwe.Decryptor.DecryptNTTNew.
func (dec *decryptor) DecryptNew(ciphertext *Ciphertext) (plaintext *Plaintext) {
	return &Plaintext{Plaintext: dec.Decryptor.DecryptNTTNew(&rlwe.Ciphertext{Value: ciphertext.Value}), Scale: ciphertext.Scale}
}

// Decrypt calls rlwe.Decryptor.Decrypt.
func (dec *decryptor) Decrypt(ciphertext *Ciphertext, plaintext *Plaintext) {
	dec.Decryptor.Decrypt(&rlwe.Ciphertext{Value: ciphertext.Value}, &rlwe.Plaintext{Value: plaintext.Value})
	plaintext.Scale = ciphertext.Scale
}
