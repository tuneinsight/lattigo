package bfv

import (
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/rlwe"
)

// Encryptor an inferface wrapping an rlwe.Encryptor.
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

// NewEncryptor instanciate a new wrapped rlwe.Encryptor for BFV ciphertexts and plaintexts.
func NewEncryptor(params Parameters, key interface{}) Encryptor {
	return &encryptor{rlwe.NewEncryptor(params.Parameters, key)}
}

// EncryptNew calls rlwe.Encryptor.EncryptNew.
func (enc *encryptor) EncryptNew(plaintext *Plaintext) *Ciphertext {
	return &Ciphertext{Element: enc.Encryptor.EncryptNew(&rlwe.Plaintext{Value: plaintext.Value})}
}

// Encrypt calls rlwe.Encryptor.Encrypt.
func (enc *encryptor) Encrypt(plaintext *Plaintext, ciphertext *Ciphertext) {
	enc.Encryptor.Encrypt(&rlwe.Plaintext{Value: plaintext.Value}, &rlwe.Element{Value: ciphertext.Value})
}

// EncryptFastNew calls rlwe.Encryptor.EncryptFastNew.
func (enc *encryptor) EncryptFastNew(plaintext *Plaintext) *Ciphertext {
	return &Ciphertext{Element: enc.Encryptor.EncryptFastNew(&rlwe.Plaintext{Value: plaintext.Value})}
}

// EncryptFast calls rlwe.Encryptor.EncryptFast.
func (enc *encryptor) EncryptFast(plaintext *Plaintext, ciphertext *Ciphertext) {
	enc.Encryptor.EncryptFast(&rlwe.Plaintext{Value: plaintext.Value}, &rlwe.Element{Value: ciphertext.Value})
}

// EncryptFromCRP calls rlwe.Encryptor.EncryptFromCRP.
func (enc *encryptor) EncryptFromCRP(plaintext *Plaintext, ciphertext *Ciphertext, crp *ring.Poly) {
	enc.Encryptor.EncryptFromCRP(&rlwe.Plaintext{Value: plaintext.Value}, &rlwe.Element{Value: ciphertext.Value}, crp)
}

// EncryptFromCRPNew calls rlwe.Encryptor.EncryptFromCRPNew.
func (enc *encryptor) EncryptFromCRPNew(plaintext *Plaintext, crp *ring.Poly) *Ciphertext {
	return &Ciphertext{Element: enc.Encryptor.EncryptFromCRPNew(&rlwe.Plaintext{Value: plaintext.Value}, crp)}
}
