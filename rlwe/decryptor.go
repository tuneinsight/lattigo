package rlwe

import (
	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/utils"
)

// Decryptor is an RLWE decryption interface.
type Decryptor interface {
	Decrypt(ct *Ciphertext, pt *Plaintext)
	DecryptNew(ct *Ciphertext) (pt *Plaintext)
	ShallowCopy() Decryptor
	WithKey(sk *SecretKey) Decryptor
}

// decryptor is a structure used to decrypt Ciphertext. It stores the secret-key.
type decryptor struct {
	params Parameters
	ringQ  *ring.Ring
	buff   *ring.Poly
	sk     *SecretKey
}

// NewDecryptor instantiates a new generic RLWE Decryptor.
func NewDecryptor(params Parameters, sk *SecretKey) Decryptor {

	if sk.Value.Q.N() != params.N() {
		panic("cannot NewDecryptor: secret_key is invalid for the provided parameters")
	}

	return &decryptor{
		params: params,
		ringQ:  params.RingQ(),
		buff:   params.RingQ().NewPoly(),
		sk:     sk,
	}
}

// Decrypt decrypts the Ciphertext and returns the result in a new Plaintext.
// Output pt MetaData will match the input ct MetaData.
func (d *decryptor) DecryptNew(ct *Ciphertext) (pt *Plaintext) {
	pt = NewPlaintext(d.params, ct.Level())
	d.Decrypt(ct, pt)
	return
}

// Decrypt decrypts the Ciphertext and writes the result in pt.
// The level of the output Plaintext is min(ct.Level(), pt.Level())
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

		ringQ.MulCoeffsMontgomeryLvl(level, pt.Value, d.sk.Value.Q, pt.Value)

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

	if !ct.IsNTT {
		ringQ.InvNTTLvl(level, pt.Value, pt.Value)
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
