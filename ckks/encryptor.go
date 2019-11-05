package ckks

import (
	"errors"
	"github.com/ldsec/lattigo/ring"
)

// Encreyptor is a struct used to encrypt plaintext and storing the public-key and/or secret-key.
type Encryptor struct {
	ckkscontext *CkksContext
	pk          *PublicKey
	sk          *SecretKey
	polypool    [3]*ring.Poly

	rescalepool []uint64

	baseconverter *ring.FastBasisExtender
}

// NewEncryptorFromPk creates a new Encryptor with the provided public-key.
// This encryptor can be used to encrypt plaintexts, using the stored key.
func (ckkscontext *CkksContext) NewEncryptorFromPk(pk *PublicKey) (*Encryptor, error) {
	return ckkscontext.newEncryptor(pk, nil)
}

// NewEncryptorFromSk creates a new Encryptor with the provided secret-key.
// This encryptor can be used to encrypt plaintexts, using the stored key.
func (ckkscontext *CkksContext) NewEncryptorFromSk(sk *SecretKey) (*Encryptor, error) {
	return ckkscontext.newEncryptor(nil, sk)
}

// NewEncryptor creates a new Encryptor with the input public-key and/or secret-key.
// This encryptor can be used to encrypt plaintexts, using the stored keys.
func (ckkscontext *CkksContext) newEncryptor(pk *PublicKey, sk *SecretKey) (encryptor *Encryptor, err error) {

	if pk != nil && (uint64(pk.pk[0].GetDegree()) != ckkscontext.n || uint64(pk.pk[1].GetDegree()) != ckkscontext.n) {
		return nil, errors.New("error : pk ring degree doesn't match ckkscontext ring degree")
	}

	if sk != nil && uint64(sk.sk.GetDegree()) != ckkscontext.n {
		return nil, errors.New("error : sk ring degree doesn't match ckkscontext ring degree")
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

	return encryptor, nil
}

// EncryptFromPkNew encrypts the input plaintext using the stored public-key and returns
// the result on a newly created ciphertext. It will encrypt the plaintext with the stored key, which can be
// private or public, a private-key encryption puts initial noise.
//
// encrypt with pk : ciphertext = [pk[0]*u + m + e_0, pk[1]*u + e_1]
// encrypt with sk : ciphertext = [-a*sk + m + e, a]
func (encryptor *Encryptor) EncryptNew(plaintext *Plaintext) (ciphertext *Ciphertext, err error) {

	ciphertext = encryptor.ckkscontext.NewCiphertext(1, plaintext.Level(), plaintext.Scale())

	return ciphertext, encryptor.Encrypt(plaintext, ciphertext)
}

// EncryptFromPk encrypts the input plaintext using the stored public-key, and returns the result
// on the reciver ciphertext. It will encrypt the plaintext with the stored key, which can be
// private or public, a private-key encryption puts initial noise.
//
// encrypt with pk : ciphertext = [pk[0]*u + m + e_0, pk[1]*u + e_1]
// encrypt with sk : ciphertext = [-a*sk + m + e, a]
func (encryptor *Encryptor) Encrypt(plaintext *Plaintext, ciphertext *Ciphertext) (err error) {

	if plaintext.Level() != encryptor.ckkscontext.levels-1 {
		return errors.New("cannot encrypt -> plaintext not at maximum level")
	}

	if encryptor.sk != nil {

		encryptfromsk(encryptor, plaintext, ciphertext)

	} else if encryptor.pk != nil {

		encryptfrompk(encryptor, plaintext, ciphertext)

	} else {

		return errors.New("cannot encrypt -> public-key and/or secret-key has not been set")
	}

	return nil
}

func encryptfrompk(encryptor *Encryptor, plaintext *Plaintext, ciphertext *Ciphertext) {

	// We sample a R-WLE instance (encryption of zero) over the keys context (ciphertext context + special prime)

	contextKeys := encryptor.ckkscontext.contextKeys
	contextQ := encryptor.ckkscontext.contextQ
	contextP := encryptor.ckkscontext.contextP

	encryptor.ckkscontext.ternarySampler.SampleNTT(0.5, encryptor.polypool[2])

	contextKeys.MulCoeffsMontgomery(encryptor.polypool[2], encryptor.pk.pk[0], encryptor.polypool[0])
	contextKeys.MulCoeffsMontgomery(encryptor.polypool[2], encryptor.pk.pk[1], encryptor.polypool[1])

	encryptor.ckkscontext.gaussianSampler.SampleNTT(encryptor.polypool[2])
	contextKeys.Add(encryptor.polypool[0], encryptor.polypool[2], encryptor.polypool[0])

	encryptor.ckkscontext.gaussianSampler.SampleNTT(encryptor.polypool[2])
	contextKeys.Add(encryptor.polypool[1], encryptor.polypool[2], encryptor.polypool[1])

	// We rescal the encryption of zero by the special prime, dividing the error by this prime
	encryptor.baseconverter.ModDownNTT(contextQ, contextP, encryptor.ckkscontext.rescaleParamsKeys, plaintext.Level(), encryptor.polypool[0], ciphertext.value[0], encryptor.polypool[2])
	encryptor.baseconverter.ModDownNTT(contextQ, contextP, encryptor.ckkscontext.rescaleParamsKeys, plaintext.Level(), encryptor.polypool[1], ciphertext.value[1], encryptor.polypool[2])

	// We switch to the ciphertext context and add the message to the encryption of zero

	ciphertext.SetCurrentModulus(contextQ.ModulusBigint)
	contextQ.Add(ciphertext.value[0], plaintext.value, ciphertext.value[0])

	ciphertext.isNTT = true
}

func encryptfromsk(encryptor *Encryptor, plaintext *Plaintext, ciphertext *Ciphertext) {
	contextKeys := encryptor.ckkscontext.contextKeys
	contextQ := encryptor.ckkscontext.contextQ
	contextP := encryptor.ckkscontext.contextP

	// ct = [e, a]
	contextKeys.UniformPoly(encryptor.polypool[1])
	encryptor.ckkscontext.gaussianSampler.SampleNTT(encryptor.polypool[0])

	// ct = [-s*a + e, a]
	contextKeys.MulCoeffsMontgomeryAndSub(encryptor.polypool[1], encryptor.sk.sk, encryptor.polypool[0])

	// We rescal by the special prime, dividing the error by this prime
	encryptor.baseconverter.ModDownNTT(contextQ, contextP, encryptor.ckkscontext.rescaleParamsKeys, plaintext.Level(), encryptor.polypool[0], ciphertext.value[0], encryptor.polypool[2])
	encryptor.baseconverter.ModDownNTT(contextQ, contextP, encryptor.ckkscontext.rescaleParamsKeys, plaintext.Level(), encryptor.polypool[1], ciphertext.value[1], encryptor.polypool[2])
	// We switch to the ciphertext context and add the message
	// ct = [-s*a + m + e, a]

	ciphertext.SetCurrentModulus(contextQ.ModulusBigint)
	contextQ.Add(ciphertext.value[0], plaintext.value, ciphertext.value[0])

	ciphertext.isNTT = true
}
