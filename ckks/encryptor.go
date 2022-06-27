package ckks

import (
	"github.com/tuneinsight/lattigo/v3/ring"
	"github.com/tuneinsight/lattigo/v3/rlwe"
)

// Encryptor an encryption interface for the CKKS scheme.
type Encryptor interface {
	Encrypt(plaintext *Plaintext, ciphertext *Ciphertext)
	EncryptNew(plaintext *Plaintext) *Ciphertext
	EncryptFromCRP(plaintext *Plaintext, crp *ring.Poly, ciphertext *Ciphertext)
	EncryptFromCRPNew(plaintext *Plaintext, crp *ring.Poly) *Ciphertext
	ShallowCopy() Encryptor
	WithKey(key interface{}) Encryptor
}

type encryptor struct {
	rlwe.Encryptor
	params Parameters
}

// NewEncryptor instantiates a new Encryptor for the CKKS scheme. The key argument can
// be *rlwe.PublicKey, *rlwe.SecretKey or nil.
func NewEncryptor(params Parameters, key interface{}) Encryptor {
	return &encryptor{rlwe.NewEncryptor(params.Parameters, key), params}
}

// Encrypt encrypts the input plaintext and write the result on ciphertext.
// The level of the output ciphertext is min(plaintext.Level(), ciphertext.Level()).
func (enc *encryptor) Encrypt(plaintext *Plaintext, ciphertext *Ciphertext) {
	enc.Encryptor.Encrypt(plaintext.Plaintext, ciphertext.Ciphertext)
	ciphertext.Scale = plaintext.Scale
}

// EncryptNew encrypts the input plaintext returns the result as a newly allocated ciphertext.
// The level of the output ciphertext is min(plaintext.Level(), ciphertext.Level()).
func (enc *encryptor) EncryptNew(plaintext *Plaintext) (ciphertext *Ciphertext) {
	ciphertext = NewCiphertext(enc.params, 1, plaintext.Level(), plaintext.Scale)
	enc.Encryptor.Encrypt(plaintext.Plaintext, ciphertext.Ciphertext)
	return
}

// EncryptFromCRP encrypts the input plaintext and writes the result in ciphertext.
// This method of encryption only works if the encryptor has been instantiated with
// a secret key.
// The passed crp is always treated as being in the NTT domain and the level of the output ciphertext is
// min(plaintext.Level(), ciphertext.Level()).
func (enc *encryptor) EncryptFromCRP(plaintext *Plaintext, crp *ring.Poly, ciphertext *Ciphertext) {
	enc.Encryptor.EncryptFromCRP(plaintext.Plaintext, crp, ciphertext.Ciphertext)
	ciphertext.Scale = plaintext.Scale
}

// EncryptFromCRPNew encrypts the input plaintext and returns the result as a newly allocated ciphertext.
// This method of encryption only works if the encryptor has been instantiated with
// a secret key.
// The passed crp is always treated as being in the NTT domain and the level of the output ciphertext is
// min(plaintext.Level(), ciphertext.Level()).
func (enc *encryptor) EncryptFromCRPNew(plaintext *Plaintext, crp *ring.Poly) (ciphertext *Ciphertext) {
	ciphertext = NewCiphertext(enc.params, 1, plaintext.Level(), plaintext.Scale)
	enc.Encryptor.EncryptFromCRP(plaintext.Plaintext, crp, ciphertext.Ciphertext)
	return
}

// ShallowCopy creates a shallow copy of this encryptor in which all the read-only data-structures are
// shared with the receiver and the temporary buffers are reallocated. The receiver and the returned
// Encryptors can be used concurrently.
func (enc *encryptor) ShallowCopy() Encryptor {
	return &encryptor{enc.Encryptor.ShallowCopy(), enc.params}
}

// WithKey creates a shallow copy of this encryptor with a new key in which all the read-only data-structures are
// shared with the receiver and the temporary buffers are reallocated. The receiver and the returned
// Encryptors can be used concurrently.
// Key can be *rlwe.PublicKey or *rlwe.SecretKey.
func (enc *encryptor) WithKey(key interface{}) Encryptor {
	return &encryptor{enc.Encryptor.WithKey(key), enc.params}
}
