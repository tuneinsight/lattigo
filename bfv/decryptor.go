package bfv

import (
	"github.com/ldsec/lattigo/v2/ring"
)

// Decryptor is an interface for decryptors
type Decryptor interface {
	// DecryptNew decrypts the input ciphertext and returns the result on a new
	// plaintext.
	DecryptNew(ciphertext *Ciphertext) *Plaintext

	// Decrypt decrypts the input ciphertext and returns the result on the
	// provided receiver plaintext. Acts accordingly depending on the plaintext type.
	Decrypt(ciphertext *Ciphertext, plaintext *Plaintext)
}

// decryptor is a structure used to decrypt ciphertexts. It stores the secret-key.
type decryptor struct {
	params   *Parameters
	ringQ    *ring.Ring
	sk       *SecretKey
	scaler   ring.Scaler
	polypool [2]*ring.Poly
}

// NewDecryptor creates a new Decryptor from the parameters with the secret-key
// given as input.
func NewDecryptor(params *Parameters, sk *SecretKey) Decryptor {

	var ringQ *ring.Ring
	var err error
	if ringQ, err = ring.NewRing(params.N(), params.qi); err != nil {
		panic(err)
	}

	return &decryptor{
		params:   params.Copy(),
		ringQ:    ringQ,
		sk:       sk,
		scaler:   ring.NewRNSScaler(params.t, ringQ),
		polypool: [2]*ring.Poly{ringQ.NewPoly(), ringQ.NewPoly()},
	}
}

func (decryptor *decryptor) DecryptNew(ciphertext *Ciphertext) *Plaintext {
	p := NewPlaintext(decryptor.params)
	decryptor.Decrypt(ciphertext, p)
	return p
}

func (decryptor *decryptor) Decrypt(ciphertext *Ciphertext, p *Plaintext) {

	ringQ := decryptor.ringQ
	tmp := decryptor.polypool[0]
	accumulator := decryptor.polypool[1]

	ringQ.NTTLazy(ciphertext.value[ciphertext.Degree()], accumulator)

	for i := uint64(ciphertext.Degree()); i > 0; i-- {
		ringQ.MulCoeffsMontgomery(accumulator, decryptor.sk.sk, accumulator)
		ringQ.NTTLazy(ciphertext.value[i-1], tmp)
		ringQ.Add(accumulator, tmp, accumulator)

		if i&3 == 3 {
			ringQ.Reduce(accumulator, accumulator)
		}
	}

	if (ciphertext.Degree())&3 != 3 {
		ringQ.Reduce(accumulator, accumulator)
	}

	ringQ.InvNTT(accumulator, p.value)
}
