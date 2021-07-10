package ckks

import (
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/rlwe"
)

// Encryptor an interface wrapping an rlwe.Encryptor.
type Encryptor interface {
	EncryptNew(plaintext *Plaintext) *Ciphertext
	Encrypt(plaintext *Plaintext, ciphertext *Ciphertext)
	EncryptFastNew(plaintext *Plaintext) *Ciphertext
	EncryptFast(plaintext *Plaintext, ciphertext *Ciphertext)
	EncryptFromCRP(plaintext *Plaintext, ciphertext *Ciphertext, crp *ring.Poly)
	EncryptFromCRPNew(plaintext *Plaintext, crp *ring.Poly) *Ciphertext
}

type encryptor struct {
	rlwe.Encryptor
}

// NewEncryptor instanciate a new wrapped rlwe.Encryptor for CKKS ciphertexts and plaintexts.
func NewEncryptor(params Parameters, key interface{}) Encryptor {
	return &encryptor{rlwe.NewEncryptor(params.Parameters, key)}
}

// EncryptNew calls rlwe.Encryptor.EncryptNew.
func (enc *encryptor) EncryptNew(plaintext *Plaintext) *Ciphertext {
	return &Ciphertext{Element: enc.Encryptor.EncryptNTTNew(&rlwe.Plaintext{Value: plaintext.Value}), Scale: plaintext.Scale}
}

// Encrypt calls rlwe.Encryptor.Encrypt.
func (enc *encryptor) Encrypt(plaintext *Plaintext, ciphertext *Ciphertext) {
	enc.Encryptor.Encrypt(&rlwe.Plaintext{Value: plaintext.Value}, &rlwe.Element{Value: ciphertext.Value})
	ciphertext.Scale = plaintext.Scale
}

// EncryptFastNew calls rlwe.Encryptor.EncryptFastNew.
func (enc *encryptor) EncryptFastNew(plaintext *Plaintext) *Ciphertext {
	return &Ciphertext{Element: enc.Encryptor.EncryptFastNTTNew(&rlwe.Plaintext{Value: plaintext.Value}), Scale: plaintext.Scale}
}

// EncryptFast calls rlwe.Encryptor.EncryptFast.
func (enc *encryptor) EncryptFast(plaintext *Plaintext, ciphertext *Ciphertext) {
	enc.Encryptor.EncryptFast(&rlwe.Plaintext{Value: plaintext.Value}, &rlwe.Element{Value: ciphertext.Value})
	ciphertext.Scale = plaintext.Scale
}

// EncryptFromCRP calls rlwe.Encryptor.EncryptFromCRP.
func (enc *encryptor) EncryptFromCRP(plaintext *Plaintext, ciphertext *Ciphertext, crp *ring.Poly) {
	enc.Encryptor.EncryptFromCRP(&rlwe.Plaintext{Value: plaintext.Value}, &rlwe.Element{Value: ciphertext.Value}, crp)
	ciphertext.Scale = plaintext.Scale
}

// EncryptFromCRPNew calls rlwe.Encryptor.EncryptFromCRPNew.
func (enc *encryptor) EncryptFromCRPNew(plaintext *Plaintext, crp *ring.Poly) *Ciphertext {
	return &Ciphertext{Element: enc.Encryptor.EncryptFromCRPNTTNew(&rlwe.Plaintext{Value: plaintext.Value}, crp), Scale: plaintext.Scale}
}
