package rlwe

import (
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/utils"
)

// Decryptor is an interface generic RLWE encryption.
type Decryptor interface {
	// Decrypt decrypts the ciphertext and write the result in ptOut.
	// The level of the output plaintext is min(ciphertext.Level(), plaintext.Level())
	// Output domain will match plaintext.Value.IsNTT value.
	Decrypt(ciphertext *Ciphertext, plaintext *Plaintext)

	ShallowCopy() Decryptor
	WithKey(sk *SecretKey) Decryptor
}

// decryptor is a structure used to decrypt ciphertext. It stores the secret-key.
type decryptor struct {
	ringQ *ring.Ring
	pool  *ring.Poly
	sk    *SecretKey
}

// ShallowCopy creates a shallow copy of Decryptor in which all the read-only data-structures are
// shared with the receiver and the temporary buffers are reallocated. The receiver and the returned
// Decryptor can be used concurrently.
func (d *decryptor) ShallowCopy() Decryptor {
	return &decryptor{
		ringQ: d.ringQ,
		pool:  d.ringQ.NewPoly(),
		sk:    d.sk,
	}
}

// WithKey creates a shallow copy of Decryptor with a new decryption key, in which all the
// read-only data-structures are shared with the receiver and the temporary buffers
// are reallocated. The receiver and the returned Decryptor can be used concurrently.
func (d *decryptor) WithKey(sk *SecretKey) Decryptor {
	return &decryptor{
		ringQ: d.ringQ,
		pool:  d.ringQ.NewPoly(),
		sk:    sk,
	}
}

// NewDecryptor instantiates a new generic RLWE Decryptor.
func NewDecryptor(params Parameters, sk *SecretKey) Decryptor {

	if sk.Value.Q.Degree() != params.N() {
		panic("secret_key is invalid for the provided parameters")
	}

	return &decryptor{
		ringQ: params.RingQ(),
		pool:  params.RingQ().NewPoly(),
		sk:    sk,
	}
}

// Decrypt decrypts the ciphertext and write the result in ptOut.
// The level of the output plaintext is min(ciphertext.Level(), plaintext.Level())
// Output domain will match plaintext.Value.IsNTT value.
func (d *decryptor) Decrypt(ciphertext *Ciphertext, plaintext *Plaintext) {

	ringQ := d.ringQ

	level := utils.MinInt(ciphertext.Level(), plaintext.Level())

	plaintext.Value.Coeffs = plaintext.Value.Coeffs[:level+1]

	if ciphertext.Value[0].IsNTT {
		ring.CopyValuesLvl(level, ciphertext.Value[ciphertext.Degree()], plaintext.Value)
	} else {
		ringQ.NTTLazyLvl(level, ciphertext.Value[ciphertext.Degree()], plaintext.Value)
	}

	for i := ciphertext.Degree(); i > 0; i-- {

		ringQ.MulCoeffsMontgomeryLvl(level, plaintext.Value, d.sk.Value.Q, plaintext.Value)

		if !ciphertext.Value[0].IsNTT {
			ringQ.NTTLazyLvl(level, ciphertext.Value[i-1], d.pool)
			ringQ.AddLvl(level, plaintext.Value, d.pool, plaintext.Value)
		} else {
			ringQ.AddLvl(level, plaintext.Value, ciphertext.Value[i-1], plaintext.Value)
		}

		if i&7 == 7 {
			ringQ.ReduceLvl(level, plaintext.Value, plaintext.Value)
		}
	}

	if (ciphertext.Degree())&7 != 7 {
		ringQ.ReduceLvl(level, plaintext.Value, plaintext.Value)
	}

	if !plaintext.Value.IsNTT {
		ringQ.InvNTTLvl(level, plaintext.Value, plaintext.Value)
	}
}
