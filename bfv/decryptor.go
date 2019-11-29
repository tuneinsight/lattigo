package bfv

import (
	"github.com/ldsec/lattigo/ring"
)

// Decryptor is a structure used to decrypt ciphertext. It stores the secret-key.
type Decryptor struct {
	params     *Parameters
	bfvContext *Context
	sk         *SecretKey
	polypool   *ring.Poly
}

// NewDecryptor creates a new Decryptor from the target context with the secret-key given as input.
func NewDecryptor(params *Parameters, sk *SecretKey) (decryptor *Decryptor) {

	if sk.sk.GetDegree() != int(params.N) {
		panic("error : secret_key degree must match context degree")
	}

	decryptor = new(Decryptor)

	decryptor.params = params.Copy()

	decryptor.bfvContext = NewContext(params)

	decryptor.sk = sk

	decryptor.polypool = decryptor.bfvContext.contextQ1.NewPoly()

	return decryptor
}

// DecryptNew decrypts the input ciphertext and returns the result on a new plaintext.
func (decryptor *Decryptor) DecryptNew(ciphertext *Ciphertext) (plaintext *Plaintext) {

	plaintext = NewPlaintextFromParams(decryptor.params)

	decryptor.Decrypt(ciphertext, plaintext)

	return
}

// Decrypt decrypts the input ciphertext and returns the result on the provided receiver plaintext.
func (decryptor *Decryptor) Decrypt(ciphertext *Ciphertext, plaintext *Plaintext) {

	ringContext := decryptor.bfvContext.contextQ1

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
