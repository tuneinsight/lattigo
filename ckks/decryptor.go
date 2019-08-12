package ckks

import (
	"errors"
)

type Decryptor struct {
	ckkscontext *CkksContext
	sk          *SecretKey
}

// NewDecryptor instanciates a new decryptor that will be able to decrypt ciphertext
// encrypted under the provided key.
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

	if err = ciphertext.value[ciphertext.Degree()].Copy(plaintext.value[0]); err != nil {
		return err
	}

	for i := uint64(ciphertext.Degree()); i > 0; i-- {

		decryptor.ckkscontext.contextLevel[level].MulCoeffsMontgomery(plaintext.value[0], decryptor.sk.sk, plaintext.value[0])
		decryptor.ckkscontext.contextLevel[level].Add(plaintext.value[0], ciphertext.value[i-1], plaintext.value[0])

		if i&7 == 7 {
			decryptor.ckkscontext.contextLevel[level].Reduce(plaintext.value[0], plaintext.value[0])
		}
	}

	if (ciphertext.Degree())&7 != 7 {
		decryptor.ckkscontext.contextLevel[level].Reduce(plaintext.value[0], plaintext.value[0])
	}

	return nil
}
