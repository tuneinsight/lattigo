package bfv

import (
	"github.com/ldsec/lattigo/ring"
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
	params     *Parameters
	bfvContext *bfvContext
	sk         *SecretKey
	polypool   *ring.Poly
}

// NewDecryptor creates a new Decryptor from the parameters with the secret-key
// given as input.
func NewDecryptor(params *Parameters, sk *SecretKey) Decryptor {

	ctx := newBFVContext(params)

	return &decryptor{
		params:     params.Copy(),
		bfvContext: ctx,
		sk:         sk,
		polypool:   ctx.ringQ.NewPoly(),
	}
}

func (decryptor *decryptor) DecryptNew(ciphertext *Ciphertext) *Plaintext {
	plaintext := NewPlaintext(decryptor.params)

	decryptor.Decrypt(ciphertext, plaintext)

	return plaintext
}

func (decryptor *decryptor) Decrypt(ciphertext *Ciphertext, plaintext *Plaintext) {
	ringContext := decryptor.bfvContext.ringQ

	ringContext.NTT(ciphertext.value[ciphertext.Degree()], plaintext.value)

	for i := uint64(ciphertext.Degree()); i > 0; i-- {
		ringContext.MulCoeffsMontgomery(plaintext.value, decryptor.sk.sk, plaintext.value)
		ringContext.NTT(ciphertext.value[i-1], decryptor.polypool)
		ringContext.Add(plaintext.value, decryptor.polypool, plaintext.value)

		if i&7 == 7 {
			ringContext.Reduce(plaintext.value, plaintext.value)
		}
	}

	if (ciphertext.Degree())&7 != 7 {
		ringContext.Reduce(plaintext.value, plaintext.value)
	}

	ringContext.InvNTT(plaintext.value, plaintext.value)
}
