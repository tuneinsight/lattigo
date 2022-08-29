package bgv

import (
	"github.com/tuneinsight/lattigo/v3/rlwe"
)

// Encryptor is an interface wrapping a rlwe.Encryptor, and modified
// to accomodate for the BGV encryption by scaling the encryptions
// of zero by the plaintext modulus before adding the plaintext.
type Encryptor interface {
	Encrypt(pt *Plaintext, ct *Ciphertext)
	EncryptNew(pt *Plaintext) *Ciphertext
	ShallowCopy() Encryptor
	WithKey(key interface{}) Encryptor
}

type encryptor struct {
	rlwe.Encryptor
	params Parameters
}

// NewEncryptor instantiates a new Encryptor for the BGV scheme. The key argument can
// be *rlwe.PublicKey, *rlwe.SecretKey or nil.
func NewEncryptor(params Parameters, key interface{}) Encryptor {
	return &encryptor{rlwe.NewEncryptor(params.Parameters, key), params}
}

// Encrypt encrypts the input plaintext and write the result on ciphertext.
func (enc *encryptor) Encrypt(pt *Plaintext, ct *Ciphertext) {
	enc.Encryptor.EncryptZero(ct.Ciphertext)
	ringQ := enc.params.RingQ()
	level := ct.Level()

	// The returned encryption of zero (as well as the public key) are
	// regular R-LWE encryptions of zero: (-as + e, a).
	// BGV operates with LSB plaintext, so we must scale the encryption of
	// zero by T: (-as + e, a) * T = (-bs + e*T, b) before adding the plaintext.
	ringQ.MulScalarLvl(level, ct.Value[0], enc.params.T(), ct.Value[0])
	ringQ.MulScalarLvl(level, ct.Value[1], enc.params.T(), ct.Value[1])
	ringQ.AddLvl(level, ct.Value[0], pt.Value, ct.Value[0])
	ct.scale = pt.scale
}

// EncryptNew encrypts the input plaintext returns the result as a newly allocated ct.
func (enc *encryptor) EncryptNew(pt *Plaintext) (ct *Ciphertext) {
	ct = NewCiphertext(enc.params, 1, pt.Level(), pt.scale)
	enc.Encrypt(pt, ct)
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
