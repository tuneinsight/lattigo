package rlwe

import (
	"fmt"

	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/utils"
	"github.com/tuneinsight/lattigo/v6/utils/structs"
)

// Decryptor is a structure used to decrypt [Ciphertext]. It stores the secret-key.
type Decryptor struct {
	params    Parameters
	ringQ     *ring.Ring
	buffQPool structs.BufferPool[*ring.Poly]
	sk        *SecretKey
}

// NewDecryptor instantiates a new generic RLWE [Decryptor].
func NewDecryptor(params ParameterProvider, sk *SecretKey) *Decryptor {

	p := params.GetRLWEParameters()

	if sk.Value.Q.N() != p.N() {
		panic(fmt.Errorf("cannot NewDecryptor: secret_key ring degree does not match parameters ring degree"))
	}

	return &Decryptor{
		params:    *p,
		ringQ:     p.RingQ(),
		buffQPool: p.ringQ.NewBuffFromUintPool(),
		sk:        sk,
	}
}

// GetRLWEParameters returns the underlying [Parameters].
func (d Decryptor) GetRLWEParameters() *Parameters {
	return &d.params
}

// DecryptNew decrypts the [Ciphertext] and returns the result in a new [Plaintext].
// Output pt [MetaData] will match the input ct [MetaData].
func (d Decryptor) DecryptNew(ct *Ciphertext) (pt *Plaintext) {
	pt = NewPlaintext(d.params, ct.Level())
	d.Decrypt(ct, pt)
	return
}

// Decrypt decrypts the [Ciphertext] and writes the result in pt.
// The level of the output [Plaintext] is min(ct.Level(), pt.Level())
// Output pt [MetaData] will match the input ct [MetaData].
func (d Decryptor) Decrypt(ct *Ciphertext, pt *Plaintext) {

	level := utils.Min(ct.Level(), pt.Level())

	ringQ := d.ringQ.AtLevel(level)

	pt.Resize(0, level)

	*pt.MetaData = *ct.MetaData

	if ct.IsNTT {
		pt.Value.CopyLvl(level, ct.Value[ct.Degree()])
	} else {
		ringQ.NTTLazy(ct.Value[ct.Degree()], pt.Value)
	}

	for i := ct.Degree(); i > 0; i-- {

		ringQ.MulCoeffsMontgomery(pt.Value, d.sk.Value.Q, pt.Value)

		if !ct.IsNTT {
			buff := d.buffQPool.Get()
			ringQ.NTTLazy(ct.Value[i-1], *buff)
			ringQ.Add(pt.Value, *buff, pt.Value)
			d.buffQPool.Put(buff)
		} else {
			ringQ.Add(pt.Value, ct.Value[i-1], pt.Value)
		}

		if i&7 == 7 {
			ringQ.Reduce(pt.Value, pt.Value)
		}
	}

	if (ct.Degree())&7 != 7 {
		ringQ.Reduce(pt.Value, pt.Value)
	}

	if !ct.IsNTT {
		ringQ.INTT(pt.Value, pt.Value)
	}
}

// ShallowCopy creates a shallow copy of [Decryptor] in which all the read-only data-structures are
// shared with the receiver and the temporary buffers are reallocated. The receiver and the returned
// [Decryptor] can be used concurrently.
func (d Decryptor) ShallowCopy() *Decryptor {
	return &Decryptor{
		params:    d.params,
		ringQ:     d.ringQ,
		buffQPool: d.buffQPool,
		sk:        d.sk,
	}
}

// WithKey creates a shallow copy of [Decryptor] with a new decryption key, in which all the
// data-structures are shared with the receiver.
// The receiver and the returned [Decryptor] can be used concurrently.
func (d Decryptor) WithKey(sk *SecretKey) *Decryptor {
	return &Decryptor{
		params:    d.params,
		ringQ:     d.ringQ,
		buffQPool: d.buffQPool,
		sk:        sk,
	}
}
