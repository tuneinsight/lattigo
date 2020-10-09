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
	// provided receiver plaintext.
	Decrypt(ciphertext *Ciphertext, plaintext *Plaintext)
}

// decryptor is a structure used to decrypt ciphertexts. It stores the secret-key.
type decryptor struct {
	params   *Parameters
	ringQ    *ring.Ring
	sk       *SecretKey
	polypool *ring.Poly
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
		polypool: ringQ.NewPoly(),
	}
}

func (decryptor *decryptor) DecryptNew(ciphertext *Ciphertext) *Plaintext {
	plaintext := NewPlaintext(decryptor.params)

	decryptor.Decrypt(ciphertext, plaintext)

	return plaintext
}

func (decryptor *decryptor) Decrypt(ciphertext *Ciphertext, plaintext *Plaintext) {
	ringQ := decryptor.ringQ

	ringQ.NTT(ciphertext.value[ciphertext.Degree()], plaintext.value)

	for i := uint64(ciphertext.Degree()); i > 0; i-- {
		ringQ.MulCoeffsMontgomery(plaintext.value, decryptor.sk.sk, plaintext.value)
		ringQ.NTT(ciphertext.value[i-1], decryptor.polypool)
		ringQ.Add(plaintext.value, decryptor.polypool, plaintext.value)

		if i&7 == 7 {
			ringQ.Reduce(plaintext.value, plaintext.value)
		}
	}

	if (ciphertext.Degree())&7 != 7 {
		ringQ.Reduce(plaintext.value, plaintext.value)
	}

	ringQ.InvNTT(plaintext.value, plaintext.value)
}
