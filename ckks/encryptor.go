package ckks

import (
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/rlwe"
)

// Encryptor an interface wrapping an rlwe.Encryptor.
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

// NewEncryptor instanciate a new wrapped rlwe.Encryptor for CKKS ciphertexts and plaintexts.
func NewEncryptor(params Parameters, key interface{}) Encryptor {
	return &encryptor{rlwe.NewEncryptor(params.Parameters, key), params}
}

func NewFastEncryptor(params Parameters, key *rlwe.PublicKey) Encryptor {
	return &encryptor{rlwe.NewFastEncryptor(params.Parameters, key), params}
}

// EncryptNew calls rlwe.Encryptor.EncryptNew.
func (enc *encryptor) EncryptNew(plaintext *Plaintext) *Ciphertext {
	ct := NewCiphertext(enc.params, 1, plaintext.Level(), plaintext.Scale)
	enc.Encryptor.Encrypt(&rlwe.Plaintext{Value: plaintext.Value}, ct.Element)
	return ct
}

// Encrypt calls rlwe.Encryptor.Encrypt.
func (enc *encryptor) Encrypt(plaintext *Plaintext, ciphertext *Ciphertext) {
	enc.Encryptor.Encrypt(&rlwe.Plaintext{Value: plaintext.Value}, &rlwe.Element{Value: ciphertext.Value})
	ciphertext.Scale = plaintext.Scale
}

// EncryptFromCRP calls rlwe.Encryptor.EncryptFromCRP.
func (enc *encryptor) EncryptFromCRP(plaintext *Plaintext, ciphertext *Ciphertext, crp *ring.Poly) {
	enc.Encryptor.EncryptFromCRP(&rlwe.Plaintext{Value: plaintext.Value}, &rlwe.Element{Value: ciphertext.Value}, crp)
	ciphertext.Scale = plaintext.Scale
}

// EncryptFromCRPNew calls rlwe.Encryptor.EncryptFromCRPNew.
func (enc *encryptor) EncryptFromCRPNew(plaintext *Plaintext, crp *ring.Poly) *Ciphertext {
	ct := NewCiphertext(enc.params, 1, plaintext.Level(), plaintext.Scale)
	enc.Encryptor.EncryptFromCRP(&rlwe.Plaintext{Value: plaintext.Value}, ct.Element, crp)
	return &Ciphertext{Element: ct.Element, Scale: plaintext.Scale}
}
