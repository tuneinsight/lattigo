package ckks

import (
	"github.com/ldsec/lattigo/ring"
)

// Encryptor is a struct used to encrypt Plaintexts. It stores the public-key and/or secret-key.
type Encryptor struct {
	params      *Parameters
	ckksContext *Context
	pk          *PublicKey
	sk          *SecretKey
	polypool    [3]*ring.Poly

	rescalepool []uint64

	baseconverter *ring.FastBasisExtender
}

// NewEncryptorFromPk creates a new Encryptor with the provided public-key.
// This Encryptor can be used to encrypt Plaintexts, using the stored key.
func NewEncryptorFromPk(params *Parameters, pk *PublicKey) *Encryptor {
	return newEncryptor(params, pk, nil)
}

// NewEncryptorFromSk creates a new Encryptor with the provided secret-key.
// This Encryptor can be used to encrypt Plaintexts, using the stored key.
func NewEncryptorFromSk(params *Parameters, sk *SecretKey) *Encryptor {
	return newEncryptor(params, nil, sk)
}

// NewEncryptor creates a new Encryptor with the input public-key and/or secret-key.
// This Encryptor can be used to encrypt Plaintexts, using the stored keys.
func newEncryptor(params *Parameters, pk *PublicKey, sk *SecretKey) (encryptor *Encryptor) {

	if !params.isValid {
		panic("cannot newEncryptor: parameters are invalid (check if the generation was done properly)")
	}

	if pk != nil && (uint64(pk.pk[0].GetDegree()) != uint64(1<<params.LogN) || uint64(pk.pk[1].GetDegree()) != uint64(1<<params.LogN)) {
		panic("cannot newEncrpytor: pk ring degree does not match params ring degree")
	}

	if sk != nil && uint64(sk.sk.GetDegree()) != uint64(1<<params.LogN) {
		panic("cannot newEncryptor: sk ring degree does not match params ring degree")
	}

	encryptor = new(Encryptor)
	encryptor.params = params.Copy()
	encryptor.ckksContext = newContext(params)
	encryptor.pk = pk
	encryptor.sk = sk

	encryptor.polypool[0] = encryptor.ckksContext.contextQP.NewPoly()
	encryptor.polypool[1] = encryptor.ckksContext.contextQP.NewPoly()
	encryptor.polypool[2] = encryptor.ckksContext.contextQP.NewPoly()

	encryptor.rescalepool = make([]uint64, encryptor.ckksContext.n)

	encryptor.baseconverter = ring.NewFastBasisExtender(encryptor.ckksContext.contextQ, encryptor.ckksContext.contextP)

	return encryptor
}

// EncryptNew encrypts the input Plaintext using the stored key and returns
// the result on a newly created Ciphertext.
//
// encrypt with pk: ciphertext = [pk[0]*u + m + e_0, pk[1]*u + e_1]
// encrypt with sk: ciphertext = [-a*sk + m + e, a]
func (encryptor *Encryptor) EncryptNew(plaintext *Plaintext) (ciphertext *Ciphertext) {

	ciphertext = NewCiphertext(encryptor.params, 1, plaintext.Level(), plaintext.Scale())
	encryptor.Encrypt(plaintext, ciphertext)
	return
}

// Encrypt encrypts the input Plaintext using the stored key, and returns the result
// on the receiver Ciphertext.
//
// encrypt with pk: ciphertext = [pk[0]*u + m + e_0, pk[1]*u + e_1]
// encrypt with sk: ciphertext = [-a*sk + m + e, a]
func (encryptor *Encryptor) Encrypt(plaintext *Plaintext, ciphertext *Ciphertext) {

	if plaintext.Level() != encryptor.ckksContext.levels-1 {
		panic("cannot Encrypt: Plaintext not at maximum level")
	}

	if encryptor.sk != nil {

		encryptfromsk(encryptor, plaintext, ciphertext)

	} else if encryptor.pk != nil {

		encryptfrompk(encryptor, plaintext, ciphertext)

	} else {

		panic("cannot Encrypt: public-key and/or secret-key has not been set")
	}
}

func encryptfrompk(encryptor *Encryptor, plaintext *Plaintext, ciphertext *Ciphertext) {

	// We sample a R-LWE instance (encryption of zero) over the keys context (ciphertext context + special prime)

	contextQP := encryptor.ckksContext.contextQP
	contextQ := encryptor.ckksContext.contextQ

	encryptor.ckksContext.contextQP.SampleTernaryMontgomeryNTT(encryptor.polypool[2], 0.5)

	// ct0 = u*pk0
	contextQP.MulCoeffsMontgomery(encryptor.polypool[2], encryptor.pk.pk[0], encryptor.polypool[0])
	// ct1 = u*pk1
	contextQP.MulCoeffsMontgomery(encryptor.polypool[2], encryptor.pk.pk[1], encryptor.polypool[1])

	// 2*(#Q + #P) NTT
	contextQP.InvNTT(encryptor.polypool[0], encryptor.polypool[0])
	contextQP.InvNTT(encryptor.polypool[1], encryptor.polypool[1])

	// ct0 = u*pk0 + e0
	encryptor.ckksContext.gaussianSampler.SampleAndAdd(encryptor.polypool[0])
	// ct1 = u*pk1 + e1
	encryptor.ckksContext.gaussianSampler.SampleAndAdd(encryptor.polypool[1])

	// ct0 = (u*pk0 + e0)/P
	encryptor.baseconverter.ModDownPQ(plaintext.Level(), encryptor.polypool[0], ciphertext.value[0])

	// ct1 = (u*pk1 + e1)/P
	encryptor.baseconverter.ModDownPQ(plaintext.Level(), encryptor.polypool[1], ciphertext.value[1])

	// 2*#Q NTT
	contextQ.NTT(ciphertext.value[0], ciphertext.value[0])
	contextQ.NTT(ciphertext.value[1], ciphertext.value[1])

	// ct0 = (u*pk0 + e0)/P + m
	contextQ.Add(ciphertext.value[0], plaintext.value, ciphertext.value[0])

	ciphertext.isNTT = true
}

func encryptfromsk(encryptor *Encryptor, plaintext *Plaintext, ciphertext *Ciphertext) {
	contextQP := encryptor.ckksContext.contextQP
	contextQ := encryptor.ckksContext.contextQ

	// ct1 = a
	contextQP.UniformPoly(encryptor.polypool[1])

	// ct0 = -s*a
	contextQP.MulCoeffsMontgomery(encryptor.polypool[1], encryptor.sk.sk, encryptor.polypool[0])
	contextQP.Neg(encryptor.polypool[0], encryptor.polypool[0])

	// #Q + #P NTT
	contextQP.InvNTT(encryptor.polypool[0], encryptor.polypool[0])

	// ct0 = -s*a + e
	encryptor.ckksContext.gaussianSampler.SampleAndAdd(encryptor.polypool[0])

	// We rescale by the special prime, dividing the error by this prime
	// ct0 = (-s*a + e)/P
	encryptor.baseconverter.ModDownPQ(plaintext.Level(), encryptor.polypool[0], ciphertext.value[0])

	// #Q + #P NTT
	// ct1 = a/P
	encryptor.baseconverter.ModDownNTTPQ(plaintext.Level(), encryptor.polypool[1], ciphertext.value[1])

	// #Q NTT
	contextQ.NTT(ciphertext.value[0], ciphertext.value[0])

	// ct0 = -s*a + m + e
	contextQ.Add(ciphertext.value[0], plaintext.value, ciphertext.value[0])

	ciphertext.isNTT = true
}
