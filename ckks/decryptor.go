package ckks

// Decryptor is a structure used to decrypt ciphertext. It stores the secret-key.
type Decryptor struct {
	params      *Parameters
	ckksContext *Context
	sk          *SecretKey
}

// NewDecryptor instanciates a new decryptor that will be able to decrypt ciphertext
// encrypted under the provided secret-key.
func NewDecryptor(params *Parameters, sk *SecretKey) *Decryptor {

	if sk.sk.GetDegree() != int(1<<params.LogN) {
		panic("secret_key degree must match context degree")
	}

	decryptor := new(Decryptor)

	decryptor.params = params.Copy()

	decryptor.ckksContext = newContext(params)

	decryptor.sk = sk

	return decryptor
}

// DecryptNew decrypts the ciphertext and returns a newly created plaintext.
// A Horner methode is used for evaluating the decryption.
func (decryptor *Decryptor) DecryptNew(ciphertext *Ciphertext) (plaintext *Plaintext) {

	plaintext = NewPlaintext(decryptor.params, ciphertext.Level(), ciphertext.Scale())

	decryptor.Decrypt(ciphertext, plaintext)

	return
}

// Decrypt decrypts the ciphertext and returns the result on the provided receiver plaintext.
// A Horner methode is used for evaluating the decryption.
func (decryptor *Decryptor) Decrypt(ciphertext *Ciphertext, plaintext *Plaintext) {

	context := decryptor.ckksContext.contextQ

	level := ciphertext.Level()

	plaintext.SetScale(ciphertext.Scale())

	context.CopyLvl(level, ciphertext.value[ciphertext.Degree()], plaintext.value)

	plaintext.value.Coeffs = plaintext.value.Coeffs[:ciphertext.Level()+1]

	for i := uint64(ciphertext.Degree()); i > 0; i-- {

		context.MulCoeffsMontgomeryLvl(level, plaintext.value, decryptor.sk.sk, plaintext.value)
		context.AddLvl(level, plaintext.value, ciphertext.value[i-1], plaintext.value)

		if i&7 == 7 {
			context.ReduceLvl(level, plaintext.value, plaintext.value)
		}
	}

	if (ciphertext.Degree())&7 != 7 {
		context.ReduceLvl(level, plaintext.value, plaintext.value)
	}
}
