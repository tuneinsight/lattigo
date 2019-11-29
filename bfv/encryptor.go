package bfv

import (
	"github.com/ldsec/lattigo/ring"
)

// Encryptor is a structure holding the parameters needed to encrypt plaintexts.
type Encryptor struct {
	params     *Parameters
	bfvContext *Context
	pk         *PublicKey
	sk         *SecretKey
	polypool   [3]*ring.Poly

	baseconverter *ring.FastBasisExtender
}

// NewEncryptorFromPk creates a new Encryptor with the provided public-key.
// This encryptor can be used to encrypt plaintexts, using the stored key.
func NewEncryptorFromPk(params *Parameters, pk *PublicKey) *Encryptor {
	return newEncryptor(params, pk, nil)
}

// NewEncryptorFromSk creates a new Encryptor with the provided secret-key.
// This encryptor can be used to encrypt plaintexts, using the stored key.
func NewEncryptorFromSk(params *Parameters, sk *SecretKey) *Encryptor {
	return newEncryptor(params, nil, sk)
}

func newEncryptor(params *Parameters, pk *PublicKey, sk *SecretKey) (encryptor *Encryptor) {

	if pk != nil && (uint64(pk.pk[0].GetDegree()) != uint64(1<<params.LogN) || uint64(pk.pk[1].GetDegree()) != uint64(1<<params.LogN)) {
		panic("error : pk ring degree doesn't match bfvcontext ring degree")
	}

	if sk != nil && uint64(sk.sk.GetDegree()) != uint64(1<<params.LogN) {
		panic("error : sk ring degree doesn't match bfvcontext ring degree")
	}

	encryptor = new(Encryptor)
	encryptor.params = params.Copy()
	encryptor.bfvContext = NewContext(params)
	encryptor.pk = pk
	encryptor.sk = sk

	encryptor.polypool[0] = encryptor.bfvContext.contextQ1P.NewPoly()
	encryptor.polypool[1] = encryptor.bfvContext.contextQ1P.NewPoly()
	encryptor.polypool[2] = encryptor.bfvContext.contextQ1P.NewPoly()

	encryptor.baseconverter = ring.NewFastBasisExtender(encryptor.bfvContext.contextQ1.Modulus, encryptor.bfvContext.contextP.Modulus)

	return
}

// EncryptNew encrypts the input plaintext using the stored key and returns
// the result on a newly created ciphertext.
//
// encrypt with pk : ciphertext = [pk[0]*u + m + e_0, pk[1]*u + e_1]
// encrypt with sk : ciphertext = [-a*sk + m + e, a]
func (encryptor *Encryptor) EncryptNew(plaintext *Plaintext) (ciphertext *Ciphertext) {

	ciphertext = NewCiphertextFromParams(encryptor.params, 1)
	encryptor.Encrypt(plaintext, ciphertext)
	return
}

// Encrypt encrypts the input plaintext using the stored key, and returns the result
// on the receiver ciphertext.
//
// encrypt with pk : ciphertext = [pk[0]*u + m + e_0, pk[1]*u + e_1]
// encrypt with sk : ciphertext = [-a*sk + m + e, a]
func (encryptor *Encryptor) Encrypt(plaintext *Plaintext, ciphertext *Ciphertext) {

	if encryptor.sk != nil {

		encryptfromsk(encryptor, plaintext, ciphertext)

	} else if encryptor.pk != nil {

		encryptfrompk(encryptor, plaintext, ciphertext)

	} else {

		panic("cannot encrypt -> public-key and/or secret-key has not been set")
	}
}

func encryptfrompk(encryptor *Encryptor, plaintext *Plaintext, ciphertext *Ciphertext) {

	ringContext := encryptor.bfvContext.contextQ1P

	// u
	ringContext.SampleTernaryMontgomeryNTT(encryptor.polypool[2], 0.5)

	// ct[0] = pk[0]*u
	// ct[1] = pk[1]*u
	ringContext.MulCoeffsMontgomery(encryptor.polypool[2], encryptor.pk.pk[0], encryptor.polypool[0])
	ringContext.MulCoeffsMontgomery(encryptor.polypool[2], encryptor.pk.pk[1], encryptor.polypool[1])

	ringContext.InvNTT(encryptor.polypool[0], encryptor.polypool[0])
	ringContext.InvNTT(encryptor.polypool[1], encryptor.polypool[1])

	// ct[0] = pk[0]*u + e0
	encryptor.bfvContext.gaussianSampler.Sample(encryptor.polypool[2])
	ringContext.Add(encryptor.polypool[0], encryptor.polypool[2], encryptor.polypool[0])

	// ct[1] = pk[1]*u + e1
	encryptor.bfvContext.gaussianSampler.Sample(encryptor.polypool[2])
	ringContext.Add(encryptor.polypool[1], encryptor.polypool[2], encryptor.polypool[1])

	// We rescal the encryption of zero by the special prime, dividing the error by this prime
	encryptor.baseconverter.ModDown(ringContext, encryptor.bfvContext.rescaleParamsKeys, uint64(len(plaintext.Value()[0].Coeffs))-1, encryptor.polypool[0], ciphertext.value[0], encryptor.polypool[2])
	encryptor.baseconverter.ModDown(ringContext, encryptor.bfvContext.rescaleParamsKeys, uint64(len(plaintext.Value()[0].Coeffs))-1, encryptor.polypool[1], ciphertext.value[1], encryptor.polypool[2])

	ringContext = encryptor.bfvContext.contextQ1

	// ct[0] = pk[0]*u + e0 + m
	// ct[1] = pk[1]*u + e1
	ringContext.Add(ciphertext.value[0], plaintext.value, ciphertext.value[0])

}

func encryptfromsk(encryptor *Encryptor, plaintext *Plaintext, ciphertext *Ciphertext) {

	ringContext := encryptor.bfvContext.contextQ1P

	// ct = [(-a*s + e)/P , a/P]
	ringContext.UniformPoly(encryptor.polypool[1])
	ringContext.MulCoeffsMontgomeryAndSub(encryptor.polypool[1], encryptor.sk.sk, encryptor.polypool[0])

	// We rescal the encryption of zero by the special prime, dividing the error by this prime
	ringContext.InvNTT(encryptor.polypool[0], encryptor.polypool[0])
	ringContext.InvNTT(encryptor.polypool[1], encryptor.polypool[1])

	encryptor.bfvContext.gaussianSampler.SampleAndAdd(encryptor.polypool[0])

	encryptor.baseconverter.ModDown(ringContext, encryptor.bfvContext.rescaleParamsKeys, uint64(len(plaintext.Value()[0].Coeffs))-1, encryptor.polypool[0], ciphertext.value[0], encryptor.polypool[2])
	encryptor.baseconverter.ModDown(ringContext, encryptor.bfvContext.rescaleParamsKeys, uint64(len(plaintext.Value()[0].Coeffs))-1, encryptor.polypool[1], ciphertext.value[1], encryptor.polypool[2])

	ringContext = encryptor.bfvContext.contextQ1

	// ct = [-a*s + m + e , a]
	ringContext.Add(ciphertext.value[0], plaintext.value, ciphertext.value[0])

}
