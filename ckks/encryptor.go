package ckks

import (
	"github.com/ldsec/lattigo/ring"
)

// Encreyptor is a struct used to encrypt plaintext and storing the public-key and/or secret-key.
type Encryptor struct {
	ckkscontext *Context
	pk          *PublicKey
	sk          *SecretKey
	polypool    [3]*ring.Poly

	rescalepool []uint64

	baseconverter *ring.FastBasisExtender
}

// NewEncryptorFromPk creates a new Encryptor with the provided public-key.
// This encryptor can be used to encrypt plaintexts, using the stored key.
func (ckkscontext *Context) NewEncryptorFromPk(pk *PublicKey) *Encryptor {
	return ckkscontext.newEncryptor(pk, nil)
}

// NewEncryptorFromSk creates a new Encryptor with the provided secret-key.
// This encryptor can be used to encrypt plaintexts, using the stored key.
func (ckkscontext *Context) NewEncryptorFromSk(sk *SecretKey) *Encryptor {
	return ckkscontext.newEncryptor(nil, sk)
}

// NewEncryptor creates a new Encryptor with the input public-key and/or secret-key.
// This encryptor can be used to encrypt plaintexts, using the stored keys.
func (ckkscontext *Context) newEncryptor(pk *PublicKey, sk *SecretKey) (encryptor *Encryptor) {

	if pk != nil && (uint64(pk.pk[0].GetDegree()) != ckkscontext.n || uint64(pk.pk[1].GetDegree()) != ckkscontext.n) {
		panic("pk ring degree doesn't match ckkscontext ring degree")
	}

	if sk != nil && uint64(sk.sk.GetDegree()) != ckkscontext.n {
		panic("sk ring degree doesn't match ckkscontext ring degree")
	}

	encryptor = new(Encryptor)
	encryptor.ckkscontext = ckkscontext
	encryptor.pk = pk
	encryptor.sk = sk

	encryptor.polypool[0] = ckkscontext.contextKeys.NewPoly()
	encryptor.polypool[1] = ckkscontext.contextKeys.NewPoly()
	encryptor.polypool[2] = ckkscontext.contextKeys.NewPoly()

	encryptor.rescalepool = make([]uint64, ckkscontext.n)

	encryptor.baseconverter = ring.NewFastBasisExtender(ckkscontext.contextQ.Modulus, ckkscontext.specialprimes)

	return encryptor
}

// EncryptFromPkNew encrypts the input plaintext using the stored public-key and returns
// the result on a newly created ciphertext. It will encrypt the plaintext with the stored key, which can be
// private or public, a private-key encryption puts initial noise.
//
// encrypt with pk : ciphertext = [pk[0]*u + m + e_0, pk[1]*u + e_1]
// encrypt with sk : ciphertext = [-a*sk + m + e, a]
func (encryptor *Encryptor) EncryptNew(plaintext *Plaintext) (ciphertext *Ciphertext) {

	ciphertext = encryptor.ckkscontext.NewCiphertext(1, plaintext.Level(), plaintext.Scale())
	encryptor.Encrypt(plaintext, ciphertext)
	return
}

// EncryptFromPk encrypts the input plaintext using the stored public-key, and returns the result
// on the reciver ciphertext. It will encrypt the plaintext with the stored key, which can be
// private or public, a private-key encryption puts initial noise.
//
// encrypt with pk : ciphertext = [pk[0]*u + m + e_0, pk[1]*u + e_1]
// encrypt with sk : ciphertext = [-a*sk + m + e, a]
func (encryptor *Encryptor) Encrypt(plaintext *Plaintext, ciphertext *Ciphertext) {

	if plaintext.Level() != encryptor.ckkscontext.levels-1 {
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

	contextKeys := encryptor.ckkscontext.contextKeys
	contextQ := encryptor.ckkscontext.contextQ

	encryptor.ckkscontext.contextKeys.SampleTernaryMontgomeryNTT(encryptor.polypool[2], 0.5)

	// ct0 = u*pk0
	contextKeys.MulCoeffsMontgomery(encryptor.polypool[2], encryptor.pk.pk[0], encryptor.polypool[0])
	// ct1 = u*pk1
	contextKeys.MulCoeffsMontgomery(encryptor.polypool[2], encryptor.pk.pk[1], encryptor.polypool[1])

	// 2*(#Q + #P) NTT
	contextKeys.InvNTT(encryptor.polypool[0], encryptor.polypool[0])
	contextKeys.InvNTT(encryptor.polypool[1], encryptor.polypool[1])

	// ct0 = u*pk0 + e0
	encryptor.ckkscontext.gaussianSampler.SampleAndAdd(encryptor.polypool[0])
	// ct1 = u*pk1 + e1
	encryptor.ckkscontext.gaussianSampler.SampleAndAdd(encryptor.polypool[1])

	// ct0 = (u*pk0 + e0)/P
	encryptor.baseconverter.ModDown(contextKeys, encryptor.ckkscontext.rescaleParamsKeys, plaintext.Level(), encryptor.polypool[0], ciphertext.value[0], encryptor.polypool[2])

	// ct1 = (u*pk1 + e1)/P
	encryptor.baseconverter.ModDown(contextKeys, encryptor.ckkscontext.rescaleParamsKeys, plaintext.Level(), encryptor.polypool[1], ciphertext.value[1], encryptor.polypool[2])

	// 2*#Q NTT
	contextQ.NTT(ciphertext.value[0], ciphertext.value[0])
	contextQ.NTT(ciphertext.value[1], ciphertext.value[1])

	// ct0 = (u*pk0 + e0)/P + m
	contextQ.Add(ciphertext.value[0], plaintext.value, ciphertext.value[0])

	ciphertext.isNTT = true
}

func encryptfromsk(encryptor *Encryptor, plaintext *Plaintext, ciphertext *Ciphertext) {
	contextKeys := encryptor.ckkscontext.contextKeys
	contextP := encryptor.ckkscontext.contextP
	contextQ := encryptor.ckkscontext.contextQ

	// ct1 = a
	contextKeys.UniformPoly(encryptor.polypool[1])

	// ct0 = -s*a
	contextKeys.MulCoeffsMontgomery(encryptor.polypool[1], encryptor.sk.sk, encryptor.polypool[0])
	contextKeys.Neg(encryptor.polypool[0], encryptor.polypool[0])

	// #Q + #P NTT
	contextKeys.InvNTT(encryptor.polypool[0], encryptor.polypool[0])

	// ct0 = -s*a + e
	encryptor.ckkscontext.gaussianSampler.SampleAndAdd(encryptor.polypool[0])

	// We rescal by the special prime, dividing the error by this prime
	// ct0 = (-s*a + e)/P
	encryptor.baseconverter.ModDown(contextKeys, encryptor.ckkscontext.rescaleParamsKeys, plaintext.Level(), encryptor.polypool[0], ciphertext.value[0], encryptor.polypool[2])

	// #Q + #P NTT
	// ct1 = a/P
	encryptor.baseconverter.ModDownNTT(contextQ, contextP, encryptor.ckkscontext.rescaleParamsKeys, plaintext.Level(), encryptor.polypool[1], ciphertext.value[1], encryptor.polypool[2])

	// #Q NTT
	contextQ.NTT(ciphertext.value[0], ciphertext.value[0])

	// ct0 = -s*a + m + e
	contextQ.Add(ciphertext.value[0], plaintext.value, ciphertext.value[0])

	ciphertext.isNTT = true
}
