package ckks

import (
	"github.com/ldsec/lattigo/ring"
)

// Encryptor is a struct used to encrypt plaintext and storing the public-key and/or secret-key.
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
// This encryptor can be used to encrypt plaintexts, using the stored key.
func NewEncryptorFromPk(params *Parameters, pk *PublicKey) *Encryptor {
	return newEncryptor(params, pk, nil)
}

// NewEncryptorFromSk creates a new Encryptor with the provided secret-key.
// This encryptor can be used to encrypt plaintexts, using the stored key.
func NewEncryptorFromSk(params *Parameters, sk *SecretKey) *Encryptor {
	return newEncryptor(params, nil, sk)
}

// NewEncryptor creates a new Encryptor with the input public-key and/or secret-key.
// This encryptor can be used to encrypt plaintexts, using the stored keys.
func newEncryptor(params *Parameters, pk *PublicKey, sk *SecretKey) (encryptor *Encryptor) {

	if pk != nil && (uint64(pk.pk[0].GetDegree()) != uint64(1<<params.LogN) || uint64(pk.pk[1].GetDegree()) != uint64(1<<params.LogN)) {
		panic("pk ring degree doesn't match ckkscontext ring degree")
	}

	if sk != nil && uint64(sk.sk.GetDegree()) != uint64(1<<params.LogN) {
		panic("sk ring degree doesn't match ckkscontext ring degree")
	}

	encryptor = new(Encryptor)
	encryptor.params = params.Copy()
	encryptor.ckksContext = NewContext(params)
	encryptor.pk = pk
	encryptor.sk = sk

	encryptor.polypool[0] = encryptor.ckksContext.contextKeys.NewPoly()
	encryptor.polypool[1] = encryptor.ckksContext.contextKeys.NewPoly()
	encryptor.polypool[2] = encryptor.ckksContext.contextKeys.NewPoly()

	encryptor.rescalepool = make([]uint64, encryptor.ckksContext.n)

	encryptor.baseconverter = ring.NewFastBasisExtender(encryptor.ckksContext.contextQ.Modulus, encryptor.ckksContext.contextP.Modulus)

	return encryptor
}

// EncryptNew encrypts the input plaintext using the stored key and returns
// the result on a newly created ciphertext.
//
// encrypt with pk : ciphertext = [pk[0]*u + m + e_0, pk[1]*u + e_1]
// encrypt with sk : ciphertext = [-a*sk + m + e, a]
func (encryptor *Encryptor) EncryptNew(plaintext *Plaintext) (ciphertext *Ciphertext) {

	ciphertext = NewCiphertextFromParams(encryptor.params, 1, plaintext.Level(), plaintext.Scale())
	encryptor.Encrypt(plaintext, ciphertext)
	return
}

// Encrypt encrypts the input plaintext using the stored key, and returns the result
// on the reciver ciphertext.
//
// encrypt with pk : ciphertext = [pk[0]*u + m + e_0, pk[1]*u + e_1]
// encrypt with sk : ciphertext = [-a*sk + m + e, a]
func (encryptor *Encryptor) Encrypt(plaintext *Plaintext, ciphertext *Ciphertext) {

	if plaintext.Level() != encryptor.ckksContext.levels-1 {
		panic("cannot encrypt -> plaintext not at maximum level")
	}

	if encryptor.sk != nil {

		encryptfromsk(encryptor, plaintext, ciphertext)

	} else if encryptor.pk != nil {

		encryptfrompk(encryptor, plaintext, ciphertext)

	} else {

		panic("cannot encrypt -> public-key and/or secret-key has not been set")
	}
}

func encryptfrompk(encryptor *Encryptor, plaintext *Plaintext, ciphertext *Ciphertext) {

	// We sample a R-WLE instance (encryption of zero) over the keys context (ciphertext context + special prime)

	contextKeys := encryptor.ckksContext.contextKeys
	contextQ := encryptor.ckksContext.contextQ

	encryptor.ckksContext.contextKeys.SampleTernaryMontgomeryNTT(encryptor.polypool[2], 0.5)

	// ct0 = u*pk0
	contextKeys.MulCoeffsMontgomery(encryptor.polypool[2], encryptor.pk.pk[0], encryptor.polypool[0])
	// ct1 = u*pk1
	contextKeys.MulCoeffsMontgomery(encryptor.polypool[2], encryptor.pk.pk[1], encryptor.polypool[1])

	// 2*(#Q + #P) NTT
	contextKeys.InvNTT(encryptor.polypool[0], encryptor.polypool[0])
	contextKeys.InvNTT(encryptor.polypool[1], encryptor.polypool[1])

	// ct0 = u*pk0 + e0
	encryptor.ckksContext.gaussianSampler.SampleAndAdd(encryptor.polypool[0])
	// ct1 = u*pk1 + e1
	encryptor.ckksContext.gaussianSampler.SampleAndAdd(encryptor.polypool[1])

	// ct0 = (u*pk0 + e0)/P
	encryptor.baseconverter.ModDown(contextKeys, encryptor.ckksContext.rescaleParamsKeys, plaintext.Level(), encryptor.polypool[0], ciphertext.value[0], encryptor.polypool[2])

	// ct1 = (u*pk1 + e1)/P
	encryptor.baseconverter.ModDown(contextKeys, encryptor.ckksContext.rescaleParamsKeys, plaintext.Level(), encryptor.polypool[1], ciphertext.value[1], encryptor.polypool[2])

	// 2*#Q NTT
	contextQ.NTT(ciphertext.value[0], ciphertext.value[0])
	contextQ.NTT(ciphertext.value[1], ciphertext.value[1])

	// ct0 = (u*pk0 + e0)/P + m
	contextQ.Add(ciphertext.value[0], plaintext.value, ciphertext.value[0])

	ciphertext.isNTT = true
}

func encryptfromsk(encryptor *Encryptor, plaintext *Plaintext, ciphertext *Ciphertext) {
	contextKeys := encryptor.ckksContext.contextKeys
	contextP := encryptor.ckksContext.contextP
	contextQ := encryptor.ckksContext.contextQ

	// ct1 = a
	contextKeys.UniformPoly(encryptor.polypool[1])

	// ct0 = -s*a
	contextKeys.MulCoeffsMontgomery(encryptor.polypool[1], encryptor.sk.sk, encryptor.polypool[0])
	contextKeys.Neg(encryptor.polypool[0], encryptor.polypool[0])

	// #Q + #P NTT
	contextKeys.InvNTT(encryptor.polypool[0], encryptor.polypool[0])

	// ct0 = -s*a + e
	encryptor.ckksContext.gaussianSampler.SampleAndAdd(encryptor.polypool[0])

	// We rescal by the special prime, dividing the error by this prime
	// ct0 = (-s*a + e)/P
	encryptor.baseconverter.ModDown(contextKeys, encryptor.ckksContext.rescaleParamsKeys, plaintext.Level(), encryptor.polypool[0], ciphertext.value[0], encryptor.polypool[2])

	// #Q + #P NTT
	// ct1 = a/P
	encryptor.baseconverter.ModDownNTT(contextQ, contextP, encryptor.ckksContext.rescaleParamsKeys, plaintext.Level(), encryptor.polypool[1], ciphertext.value[1], encryptor.polypool[2])

	// #Q NTT
	contextQ.NTT(ciphertext.value[0], ciphertext.value[0])

	// ct0 = -s*a + m + e
	contextQ.Add(ciphertext.value[0], plaintext.value, ciphertext.value[0])

	ciphertext.isNTT = true
}
