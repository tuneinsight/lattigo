package rlwe

import (
	"math"
	"math/big"

	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/utils"
)

// Decryptor is an interface generic RLWE encryption.
type Decryptor interface {
	Decrypt(ct *Ciphertext, pt *Plaintext)
	DecryptNew(ct *Ciphertext) (pt *Plaintext)
	ShallowCopy() Decryptor
	WithKey(sk *SecretKey) Decryptor
}

// decryptor is a structure used to decrypt ciphertext. It stores the secret-key.
type decryptor struct {
	params Parameters
	ringQ  *ring.Ring
	buff   *ring.Poly
	sk     *SecretKey
}

// NewDecryptor instantiates a new generic RLWE Decryptor.
func NewDecryptor(params Parameters, sk *SecretKey) Decryptor {

	if sk.Q.N() != params.N() {
		panic("secret_key is invalid for the provided parameters")
	}

	return &decryptor{
		params: params,
		ringQ:  params.RingQ(),
		buff:   params.RingQ().NewPoly(),
		sk:     sk,
	}
}

// Decrypt decrypts the ciphertext and returns the result in a new plaintext.
// Output pt MetaData will match the input ct MetaData.
func (d *decryptor) DecryptNew(ct *Ciphertext) (pt *Plaintext) {
	pt = NewPlaintext(d.params, ct.Level())
	d.Decrypt(ct, pt)
	return
}

// Decrypt decrypts the ciphertext and writes the result in pt.
// The level of the output plaintext is min(ct.Level(), pt.Level())
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

// Norm returns the log2 of the standard deviation, minimum and maximum absolute norm of
// the decrypted ciphertext, before the decoding (i.e. including the error).
func Norm(ct *Ciphertext, dec Decryptor) (std, min, max float64) {

	params := dec.(*decryptor).params

	coeffsBigint := make([]*big.Int, params.N())
	for i := range coeffsBigint {
		coeffsBigint[i] = new(big.Int)
	}

	pt := NewPlaintext(params, ct.Level())

	dec.Decrypt(ct, pt)

	if pt.IsNTT {
		params.RingQ().InvNTTLvl(ct.Level(), pt.Value, pt.Value)
	}

	params.RingQ().PolyToBigintCenteredLvl(ct.Level(), pt.Value, 1, coeffsBigint)

	return normStats(coeffsBigint)
}

func normStats(vec []*big.Int) (float64, float64, float64) {

	vecfloat := make([]*big.Float, len(vec))
	minErr := new(big.Float).SetFloat64(0)
	maxErr := new(big.Float).SetFloat64(0)
	tmp := new(big.Float)
	minErr.SetInt(vec[0])
	minErr.Abs(minErr)
	for i := range vec {
		vecfloat[i] = new(big.Float)
		vecfloat[i].SetInt(vec[i])

		tmp.Abs(vecfloat[i])

		if minErr.Cmp(tmp) == 1 {
			minErr.Set(tmp)
		}

		if maxErr.Cmp(tmp) == -1 {
			maxErr.Set(tmp)
		}
	}

	n := new(big.Float).SetFloat64(float64(len(vec)))

	mean := new(big.Float).SetFloat64(0)

	for _, c := range vecfloat {
		mean.Add(mean, c)
	}

	mean.Quo(mean, n)

	err := new(big.Float).SetFloat64(0)
	for _, c := range vecfloat {
		tmp.Sub(c, mean)
		tmp.Mul(tmp, tmp)
		err.Add(err, tmp)
	}

	err.Quo(err, n)
	err.Sqrt(err)

	x, _ := err.Float64()
	y, _ := minErr.Float64()
	z, _ := maxErr.Float64()

	return math.Log2(x), math.Log2(y), math.Log2(z)
}
