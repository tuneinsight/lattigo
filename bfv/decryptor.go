package bfv

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

// NewDecryptor instantiates a new rlwe.Decryptor wrapped for BFV ciphertexts and plaintexts.
func NewDecryptor(params Parameters, sk *rlwe.SecretKey) Decryptor {
	return &decryptor{Decryptor: rlwe.NewDecryptor(params.Parameters, sk)}
}

// DecryptNew calls rlwe.Decryptor.DecryptNew.
func (dec *decryptor) DecryptNew(ciphertext *Ciphertext) (plaintext *Plaintext) {
	return &Plaintext{Plaintext: dec.Decryptor.DecryptNew(&rlwe.Element{Value: ciphertext.Value})}
}

// Decrypt calls rlwe.Decryptor.Decrypt.
func (dec *decryptor) Decrypt(ciphertext *Ciphertext, plaintext *Plaintext) {
	dec.Decryptor.Decrypt(&rlwe.Element{Value: ciphertext.Value}, &rlwe.Plaintext{Value: plaintext.Value})
}
