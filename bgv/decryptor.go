package bgv

import (
	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/rlwe"
)

// Decryptor is an interface wrapping a rlwe.Decryptor.
type Decryptor interface {
	DecryptNew(ct *rlwe.Ciphertext) (pt *rlwe.Plaintext)
	Decrypt(ct *rlwe.Ciphertext, pt *rlwe.Plaintext)
	ShallowCopy() Decryptor
	WithKey(sk *rlwe.SecretKey) Decryptor
}

type decryptor struct {
	rlwe.Decryptor
	params Parameters
	buffQ  *ring.Poly
}

// NewDecryptor instantiates a Decryptor for the BGV scheme.
func NewDecryptor(params Parameters, sk *rlwe.SecretKey) Decryptor {
	buffQ := params.RingQ().NewPoly()
	buffQ.IsNTT = true
	return &decryptor{rlwe.NewDecryptor(params.Parameters, sk), params, buffQ}
}

// Decrypt decrypts the ciphertext and writes the result in pt.
func (dec *decryptor) Decrypt(ct *rlwe.Ciphertext, pt *rlwe.Plaintext) {
	dec.Decryptor.Decrypt(ct, pt)
}

// DecryptNew decrypts the ciphertext and returns the result in a newly allocated Plaintext.
func (dec *decryptor) DecryptNew(ct *rlwe.Ciphertext) (pt *rlwe.Plaintext) {
	pt = NewPlaintext(dec.params, ct.Level())
	dec.Decryptor.Decrypt(ct, pt)
	return
}

// ShallowCopy creates a shallow copy of Decryptor in which all the read-only data-structures are
// shared with the receiver and the temporary buffers are reallocated. The receiver and the returned
// Decryptor can be used concurrently.
func (dec *decryptor) ShallowCopy() Decryptor {
	buffQ := dec.params.RingQ().NewPoly()
	buffQ.IsNTT = true
	return &decryptor{dec.Decryptor.ShallowCopy(), dec.params, buffQ}
}

// WithKey creates a shallow copy of Decryptor with a new decryption key, in which all the
// read-only data-structures are shared with the receiver and the temporary buffers
// are reallocated. The receiver and the returned Decryptor can be used concurrently.
func (dec *decryptor) WithKey(sk *rlwe.SecretKey) Decryptor {
	return &decryptor{dec.Decryptor.WithKey(sk), dec.params, dec.buffQ}
}
