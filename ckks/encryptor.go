package ckks

import (
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/rlwe"
)

// Encryptor an encryption interface for the CKKS scheme.
type Encryptor interface {
	Encrypt(plaintext *Plaintext, ciphertext *Ciphertext)
	EncryptNew(plaintext *Plaintext) *Ciphertext
	EncryptFromCRP(plaintext *Plaintext, crp *ring.Poly, ciphertext *Ciphertext)
	EncryptFromCRPNew(plaintext *Plaintext, crp *ring.Poly) *Ciphertext
}

type encryptor struct {
	rlwe.Encryptor
	params Parameters
}

// NewEncryptor instantiates a new Encryptor for the CKKS scheme. The key argument can
// be either a *rlwe.PublicKey or a *rlwe.SecretKey.
func NewEncryptor(params Parameters, key interface{}) Encryptor {
	return &encryptor{rlwe.NewEncryptor(params.Parameters, key), params}
}

// NewFastEncryptor instantiates a new Encryptor for the CKKS scheme.
// This encryptor's Encrypt method first encrypts zero in Q and then adds the plaintext.
// This method is faster than the normal encryptor but result in a noisier ciphertext.
func NewFastEncryptor(params Parameters, key *rlwe.PublicKey) Encryptor {
	return &encryptor{rlwe.NewFastEncryptor(params.Parameters, key), params}
}

// Encrypt encrypts the input plaintext and write the result on ciphertext. The encryption
// algorithm depends on how the receiver encryptor was initialized (see NewEncryptor
// and NewFastEncryptor).
// The level of the output ciphertext is min(plaintext.Level(), ciphertext.Level()).
func (enc *encryptor) Encrypt(plaintext *Plaintext, ciphertext *Ciphertext) {
	enc.Encryptor.Encrypt(plaintext.Plaintext, ciphertext.Ciphertext)
	ciphertext.Scale = plaintext.Scale
}

// EncryptNew encrypts the input plaintext returns the result as a newly allocated ciphertext.
// The encryption algorithm depends on how the receiver encryptor was initialized (see
// NewEncryptor and NewFastEncryptor).
// The level of the output ciphertext is min(plaintext.Level(), ciphertext.Level()).
func (enc *encryptor) EncryptNew(plaintext *Plaintext) (ciphertext *Ciphertext) {
	ciphertext = NewCiphertext(enc.params, 1, plaintext.Level(), plaintext.Scale)
	enc.Encryptor.Encrypt(plaintext.Plaintext, ciphertext.Ciphertext)
	return
}

// EncryptFromCRP encrypts the input plaintext and writes the result in ciphertext.
// The passed crp is always treated as being in the NTT domain and the level of the output ciphertext is
// min(plaintext.Level(), ciphertext.Level()).
func (enc *encryptor) EncryptFromCRP(plaintext *Plaintext, crp *ring.Poly, ciphertext *Ciphertext) {
	enc.Encryptor.EncryptFromCRP(plaintext.Plaintext, crp, ciphertext.Ciphertext)
	ciphertext.Scale = plaintext.Scale
}

// EncryptFromCRPNew encrypts the input plaintext and returns the result as a newly allocated ciphertext.
// The passed crp is always treated as being in the NTT domain and the level of the output ciphertext is
// min(plaintext.Level(), ciphertext.Level()).
func (enc *encryptor) EncryptFromCRPNew(plaintext *Plaintext, crp *ring.Poly) (ciphertext *Ciphertext) {
	ciphertext = NewCiphertext(enc.params, 1, plaintext.Level(), plaintext.Scale)
	enc.Encryptor.EncryptFromCRP(plaintext.Plaintext, crp, ciphertext.Ciphertext)
	return
}
