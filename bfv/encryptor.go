package bfv

import (
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/rlwe"
)

// Encryptor an inferface wrapping an rlwe.Encryptor.
type Encryptor interface {
	Encrypt(plaintext *Plaintext, ciphertext *Ciphertext)
	EncryptNew(plaintext *Plaintext) *Ciphertext
	EncryptFromCRP(plaintext *Plaintext, ciphertext *Ciphertext, crp *ring.Poly)
	EncryptFromCRPNew(plaintext *Plaintext, crp *ring.Poly) *Ciphertext
}

type encryptor struct {
	rlwe.Encryptor
	params Parameters
}

// NewEncryptor instanciate a new wrapped rlwe.Encryptor for BFV ciphertexts and plaintexts.
func NewEncryptor(params Parameters, key interface{}) Encryptor {
	return &encryptor{rlwe.NewEncryptor(params.Parameters, key), params}
}

func NewFastEncryptor(params Parameters, key *rlwe.PublicKey) Encryptor {
	return &encryptor{rlwe.NewFastEncryptor(params.Parameters, key), params}
}

// EncryptNew calls rlwe.Encryptor.EncryptNew.
func (enc *encryptor) EncryptNew(plaintext *Plaintext) *Ciphertext {
	ct := NewCiphertext(enc.params, 1)
	enc.Encryptor.Encrypt(plaintext.Plaintext, ct.Ciphertext)
	return ct
}

// Encrypt calls rlwe.Encryptor.Encrypt.
func (enc *encryptor) Encrypt(plaintext *Plaintext, ciphertext *Ciphertext) {
	enc.Encryptor.Encrypt(&rlwe.Plaintext{Value: plaintext.Value}, &rlwe.Ciphertext{Value: ciphertext.Value})
}

// EncryptFromCRP calls rlwe.Encryptor.EncryptFromCRP.
func (enc *encryptor) EncryptFromCRP(plaintext *Plaintext, ciphertext *Ciphertext, crp *ring.Poly) {
	enc.Encryptor.EncryptFromCRP(&rlwe.Plaintext{Value: plaintext.Value}, &rlwe.Ciphertext{Value: ciphertext.Value}, crp)
}

// EncryptFromCRPNew calls rlwe.Encryptor.EncryptFromCRPNew.
func (enc *encryptor) EncryptFromCRPNew(plaintext *Plaintext, crp *ring.Poly) *Ciphertext {
	ct := NewCiphertext(enc.params, 1)
	enc.Encryptor.EncryptFromCRP(&rlwe.Plaintext{Value: plaintext.Value}, ct.Ciphertext, crp)
	return ct
}
