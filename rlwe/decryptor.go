package rlwe

import (
	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/utils"
)

// Decryptor is an interface generic RLWE encryption.
type Decryptor interface {
	Decrypt(ciphertext *Ciphertext, plaintext *Plaintext)
	ShallowCopy() Decryptor
	WithKey(sk *SecretKey) Decryptor
}

// decryptor is a structure used to decrypt ciphertext. It stores the secret-key.
type decryptor struct {
	ringQ *ring.Ring
	buff  *ring.Poly
	sk    *SecretKey
}

// NewDecryptor instantiates a new generic RLWE Decryptor.
func NewDecryptor(params Parameters, sk *SecretKey) Decryptor {

	if sk.Q.N() != params.N() {
		panic("secret_key is invalid for the provided parameters")
	}

	return &decryptor{
		ringQ: params.RingQ(),
		buff:  params.RingQ().NewPoly(),
		sk:    sk,
	}
}

// Decrypt decrypts the ciphertext and writes the result in ptOut.
// The level of the output plaintext is min(ciphertext.Level(), plaintext.Level())
// Output pt MetaData will match the input ct MetaData.
func (d *decryptor) Decrypt(ct *Ciphertext, pt *Plaintext) {

	ringQ := d.ringQ

	level := utils.MinInt(ct.Level(), pt.Level())

	pt.Value.Resize(level)

	pt.MetaData = ct.MetaData

	if ct.IsNTT {
		ring.CopyLvl(level, ct.Value[ct.Degree()], pt.Value)
	} else {
		ringQ.NTTLazyLvl(level, ct.Value[ct.Degree()], pt.Value)
	}

	for i := ct.Degree(); i > 0; i-- {

		ringQ.MulCoeffsMontgomeryLvl(level, pt.Value, d.sk.Q, pt.Value)

		if !ct.IsNTT {
			ringQ.NTTLazyLvl(level, ct.Value[i-1], d.buff)
			ringQ.AddLvl(level, pt.Value, d.buff, pt.Value)
		} else {
			ringQ.AddLvl(level, pt.Value, ct.Value[i-1], pt.Value)
		}

		if i&7 == 7 {
			ringQ.ReduceLvl(level, pt.Value, pt.Value)
		}
	}

	if (ct.Degree())&7 != 7 {
		ringQ.ReduceLvl(level, pt.Value, pt.Value)
	}
}

// ShallowCopy creates a shallow copy of Decryptor in which all the read-only data-structures are
// shared with the receiver and the temporary buffers are reallocated. The receiver and the returned
// Decryptor can be used concurrently.
func (d *decryptor) ShallowCopy() Decryptor {
	return &decryptor{
		ringQ: d.ringQ,
		buff:  d.ringQ.NewPoly(),
		sk:    d.sk,
	}
}

// WithKey creates a shallow copy of Decryptor with a new decryption key, in which all the
// read-only data-structures are shared with the receiver and the temporary buffers
// are reallocated. The receiver and the returned Decryptor can be used concurrently.
func (d *decryptor) WithKey(sk *SecretKey) Decryptor {
	return &decryptor{
		ringQ: d.ringQ,
		buff:  d.ringQ.NewPoly(),
		sk:    sk,
	}
}
