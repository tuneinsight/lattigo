package bgv

import (
	"math"
	"math/big"

	"github.com/tuneinsight/lattigo/v3/ring"
	"github.com/tuneinsight/lattigo/v3/rlwe"
)

// Decryptor is an interface wrapping a rlwe.Decryptor.
type Decryptor interface {
	DecryptNew(ciphertext *Ciphertext) (plaintext *Plaintext)
	Decrypt(ciphertext *Ciphertext, plaintext *Plaintext)
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

// Decrypt decrypts the ciphertext and write the result in ptOut.
func (dec *decryptor) Decrypt(ct *Ciphertext, ptOut *Plaintext) {
	dec.Decryptor.Decrypt(ct.Ciphertext, ptOut.Plaintext)
	ptOut.Scale = ct.Scale
}

// DecryptNew decrypts the ciphertext and returns the result in a newly allocated Plaintext.
func (dec *decryptor) DecryptNew(ct *Ciphertext) (ptOut *Plaintext) {
	pt := NewPlaintext(dec.params, ct.Level(), ct.Scale)
	dec.Decryptor.Decrypt(ct.Ciphertext, pt.Plaintext)
	return pt
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

func errorStats(vec []*big.Int) (float64, float64, float64) {

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
