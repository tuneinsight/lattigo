package ckks

import (
	"github.com/tuneinsight/lattigo/v4/rlwe"
)

// Decryptor is an interface wrapping a rlwe.Decryptor.
type Decryptor interface {
	Decrypt(ct *rlwe.Ciphertext, pt *rlwe.Plaintext)
	DecryptNew(ct *rlwe.Ciphertext) (pt *rlwe.Plaintext)
	ShallowCopy() Decryptor
	WithKey(sk *rlwe.SecretKey) Decryptor
}
type decryptor struct {
	rlwe.Decryptor
	params Parameters
}

// NewDecryptor instantiates a Decryptor for the CKKS scheme.
func NewDecryptor(params Parameters, sk *rlwe.SecretKey) Decryptor {
	return &decryptor{rlwe.NewDecryptor(params.Parameters, sk), params}
}

// DecryptNew decrypts the ciphertext and returns the result in a newly allocated Plaintext.
func (dec *decryptor) DecryptNew(ct *rlwe.Ciphertext) (pt *rlwe.Plaintext) {
	pt = rlwe.NewPlaintextNTT(dec.params.Parameters, ct.Level())
	dec.Decrypt(ct, pt)
	return
}

// Decrypt decrypts the ciphertext and writes the result in pt.
func (dec *decryptor) Decrypt(ct *rlwe.Ciphertext, pt *rlwe.Plaintext) {
	dec.Decryptor.Decrypt(ct, pt)
}

// ShallowCopy creates a shallow copy of Decryptor in which all the read-only data-structures are
// shared with the receiver and the temporary buffers are reallocated. The receiver and the returned
// Decryptor can be used concurrently.
func (dec *decryptor) ShallowCopy() Decryptor {
	return &decryptor{dec.Decryptor.ShallowCopy(), dec.params}
}

// WithKey creates a shallow copy of Decryptor with a new decryption key, in which all the
// read-only data-structures are shared with the receiver and the temporary buffers
// are reallocated. The receiver and the returned Decryptor can be used concurrently.
func (dec *decryptor) WithKey(sk *rlwe.SecretKey) Decryptor {
	return &decryptor{dec.Decryptor.WithKey(sk), dec.params}
}
