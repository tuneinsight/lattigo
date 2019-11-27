package bfv

import (
	"github.com/ldsec/lattigo/ring"
)

type decryptorContext struct {
	// Polynomial degree
	n uint64

	// Polynomial contexts
	contextQ *ring.Context
}

func newDecryptorContext(params *Parameters) *decryptorContext {
	n := params.N

	contextQ := ring.NewContext()
	contextQ.SetParameters(n, params.Qi)
	if err := contextQ.GenNTTParams(); err != nil {
		panic(err)
	}

	return &decryptorContext{
		n:        n,
		contextQ: contextQ,
	}
}

// Decryptor is a structure used to decrypt ciphertext. It stores the secret-key.
type Decryptor struct {
	context  *decryptorContext
	sk       *SecretKey
	polypool *ring.Poly
}

// NewDecryptor creates a new Decryptor from the target context with the secret-key given as input.
func NewDecryptor(sk *SecretKey, params *Parameters) (decryptor *Decryptor) {
	context := newDecryptorContext(params)

	if sk.sk.GetDegree() != int(context.n) {
		panic("error : secret_key degree must match context degree")
	}

	decryptor = new(Decryptor)

	decryptor.context = context

	decryptor.sk = sk

	decryptor.polypool = context.contextQ.NewPoly()

	return decryptor
}

// DecryptNew decrypts the input ciphertext and returns the result on a new plaintext.
func (decryptor *Decryptor) DecryptNew(ciphertext *Ciphertext) (plaintext *Plaintext) {

	plaintext = NewPlaintext(decryptor.context.contextQ)

	decryptor.Decrypt(ciphertext, plaintext)

	return
}

// Decrypt decrypts the input ciphertext and returns the result on the provided receiver plaintext.
func (decryptor *Decryptor) Decrypt(ciphertext *Ciphertext, plaintext *Plaintext) {

	ringContext := decryptor.context.contextQ

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
