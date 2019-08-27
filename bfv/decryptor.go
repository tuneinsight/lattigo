package bfv

import (
	"errors"
	"github.com/lca1/lattigo/ring"
)

// Decryptor is a structure used to decrypt ciphertext. It stores the secret-key.
type Decryptor struct {
	bfvcontext *BfvContext
	sk         *SecretKey
	polypool   *ring.Poly
}

// NewDecryptor creates a new Decryptor from the target bfvcontext with the secret-key given as input.
func (bfvcontext *BfvContext) NewDecryptor(sk *SecretKey) (decryptor *Decryptor, err error) {

	if sk.sk.GetDegree() != int(bfvcontext.n) {
		return nil, errors.New("error : secret_key degree must match context degree")
	}

	decryptor = new(Decryptor)

	decryptor.bfvcontext = bfvcontext

	decryptor.sk = sk

	decryptor.polypool = bfvcontext.contextQ.NewPoly()

	return decryptor, nil
}

// DecryptNew decrypts the input ciphertext and returns the result on a new plaintext.
func (decryptor *Decryptor) DecryptNew(ciphertext *Ciphertext) (plaintext *Plaintext, err error) {

	plaintext = decryptor.bfvcontext.NewPlaintext()

	if err = decryptor.Decrypt(ciphertext, plaintext); err != nil {
		return nil, err
	}

	return plaintext, nil
}

// Decrypt decrypts the input ciphertext and returns the value on the provided receiver plaintext.
func (decryptor *Decryptor) Decrypt(ciphertext *Ciphertext, plaintext *Plaintext) (err error) {

	decryptor.bfvcontext.contextQ.NTT(ciphertext.value[ciphertext.Degree()], plaintext.value[0])

	for i := uint64(ciphertext.Degree()); i > 0; i-- {
		decryptor.bfvcontext.contextQ.MulCoeffsMontgomery(plaintext.value[0], decryptor.sk.sk, plaintext.value[0])
		decryptor.bfvcontext.contextQ.NTT(ciphertext.value[i-1], decryptor.polypool)
		decryptor.bfvcontext.contextQ.Add(plaintext.value[0], decryptor.polypool, plaintext.value[0])

		if i&7 == 7 {
			decryptor.bfvcontext.contextQ.Reduce(plaintext.value[0], plaintext.value[0])
		}
	}

	if (ciphertext.Degree())&7 != 7 {
		decryptor.bfvcontext.contextQ.Reduce(plaintext.value[0], plaintext.value[0])
	}

	decryptor.bfvcontext.contextQ.InvNTT(plaintext.value[0], plaintext.value[0])

	return nil
}
