package ckks

import (
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/rlwe"
	"github.com/ldsec/lattigo/v2/utils"
)

// Decryptor is an interface for decrypting Ciphertexts. A Decryptor stores the secret-key.
type Decryptor interface {
	// DecryptNew decrypts the ciphertext and returns a newly created
	// plaintext. A Horner method is used for evaluating the decryption.
	// The level of the output plaintext is ciphertext.Level().
	DecryptNew(ciphertext *Ciphertext) (plaintext *Plaintext)

	// Decrypt decrypts the ciphertext and returns the result on the provided
	// receiver plaintext. A Horner method is used for evaluating the
	// decryption.
	// The level of the output plaintext is min(ciphertext.Level(), plaintext.Level())
	Decrypt(ciphertext *Ciphertext, plaintext *Plaintext)
}

// decryptor is a structure used to decrypt ciphertext. It stores the secret-key.
type decryptor struct {
	params Parameters
	ringQ  *ring.Ring
	sk     *rlwe.SecretKey
}

// NewDecryptor instantiates a new Decryptor that will be able to decrypt ciphertexts
// encrypted under the provided secret-key.
func NewDecryptor(params Parameters, sk *rlwe.SecretKey) Decryptor {

	if sk.Value.Degree() != params.N() {
		panic("secret_key is invalid for the provided parameters")
	}

	return &decryptor{
		params: params,
		ringQ:  params.RingQ(),
		sk:     sk,
	}
}

func (decryptor *decryptor) DecryptNew(ciphertext *Ciphertext) (plaintext *Plaintext) {

	plaintext = NewPlaintext(decryptor.params, ciphertext.Level(), ciphertext.Scale)

	decryptor.Decrypt(ciphertext, plaintext)

	return plaintext
}

func (decryptor *decryptor) Decrypt(ciphertext *Ciphertext, plaintext *Plaintext) {

	level := utils.MinInt(ciphertext.Level(), plaintext.Level())

	plaintext.Scale = ciphertext.Scale

	ring.CopyValuesLvl(level, ciphertext.Value[ciphertext.Degree()], plaintext.Value)

	plaintext.Value.Coeffs = plaintext.Value.Coeffs[:ciphertext.Level()+1]

	for i := ciphertext.Degree(); i > 0; i-- {

		decryptor.ringQ.MulCoeffsMontgomeryLvl(level, plaintext.Value, decryptor.sk.Value, plaintext.Value)
		decryptor.ringQ.AddLvl(level, plaintext.Value, ciphertext.Value[i-1], plaintext.Value)

		if i&7 == 7 {
			decryptor.ringQ.ReduceLvl(level, plaintext.Value, plaintext.Value)
		}
	}

	if (ciphertext.Degree())&7 != 7 {
		decryptor.ringQ.ReduceLvl(level, plaintext.Value, plaintext.Value)
	}

	plaintext.Value.Coeffs = plaintext.Value.Coeffs[:level+1]
}
