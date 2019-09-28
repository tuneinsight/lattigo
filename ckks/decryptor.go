package ckks

import (
	"errors"
)

// Decryptor is a structure used to decrypt ciphertext. It stores the secret-key.
type Decryptor struct {
	ckkscontext *CkksContext
	sk          *SecretKey
}

// NewDecryptor instanciates a new decryptor that will be able to decrypt ciphertext
// encrypted under the provided secret-key.
func (ckkscontext *CkksContext) NewDecryptor(sk *SecretKey) (*Decryptor, error) {

	if sk.sk.GetDegree() != int(ckkscontext.n) {
		return nil, errors.New("error : secret_key degree must match context degree")
	}

	decryptor := new(Decryptor)

	decryptor.ckkscontext = ckkscontext

	decryptor.sk = sk

	return decryptor, nil
}

// DecryptNew decrypts the ciphertext and returns a newly created plaintext.
// A Horner methode is used for evaluating the decryption.
func (decryptor *Decryptor) DecryptNew(ciphertext *Ciphertext) (plaintext *Plaintext, err error) {

	plaintext = decryptor.ckkscontext.NewPlaintext(ciphertext.Level(), ciphertext.Scale())

	if err = decryptor.Decrypt(ciphertext, plaintext); err != nil {
		return nil, err
	}

	return plaintext, nil
}

// Decrypt decrypts the ciphertext and returns the result on the provided receiver plaintext.
// A Horner methode is used for evaluating the decryption.
func (decryptor *Decryptor) Decrypt(ciphertext *Ciphertext, plaintext *Plaintext) (err error) {

	level := ciphertext.Level()

	plaintext.SetScale(ciphertext.Scale())
	plaintext.currentModulus.SetBigInt(ciphertext.currentModulus)

	plaintext.value.Copy(ciphertext.value[ciphertext.Degree()])

	for i := uint64(ciphertext.Degree()); i > 0; i-- {

		decryptor.ckkscontext.contextLevel[level].MulCoeffsMontgomery(plaintext.value, decryptor.sk.sk, plaintext.value)
		decryptor.ckkscontext.contextLevel[level].Add(plaintext.value, ciphertext.value[i-1], plaintext.value)

		if i&7 == 7 {
			decryptor.ckkscontext.contextLevel[level].Reduce(plaintext.value, plaintext.value)
		}
	}

	if (ciphertext.Degree())&7 != 7 {
		decryptor.ckkscontext.contextLevel[level].Reduce(plaintext.value, plaintext.value)
	}

	return nil
}
