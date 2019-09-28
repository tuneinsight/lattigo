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
	polypool    *ring.Poly
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
	encryptor.polypool = ckkscontext.contextLevel[ckkscontext.levels-1].NewPoly()

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

	context := encryptor.ckkscontext.contextLevel[plaintext.Level()]

	encryptor.ckkscontext.ternarySampler.SampleMontgomeryNTT(0.5, encryptor.polypool)

	context.MulCoeffsMontgomery(encryptor.polypool, encryptor.pk.pk[0], ciphertext.value[0])
	context.MulCoeffsMontgomery(encryptor.polypool, encryptor.pk.pk[1], ciphertext.value[1])

	encryptor.ckkscontext.gaussianSampler.SampleNTT(encryptor.polypool)
	context.Add(ciphertext.value[0], encryptor.polypool, ciphertext.value[0])

	encryptor.ckkscontext.gaussianSampler.SampleNTT(encryptor.polypool)
	context.Add(ciphertext.value[1], encryptor.polypool, ciphertext.value[1])

	context.Add(ciphertext.value[0], plaintext.value, ciphertext.value[0])

	encryptor.polypool.Zero()

	ciphertext.isNTT = true
}

func encryptfromsk(encryptor *Encryptor, plaintext *Plaintext, ciphertext *Ciphertext) {
	context := encryptor.ckkscontext.contextLevel[plaintext.Level()]

	// ct = [e, a]
	ciphertext.value[1] = context.NewUniformPoly()
	encryptor.ckkscontext.gaussianSampler.SampleNTT(ciphertext.value[0])

	// ct = [-s*a + e, a]
	context.MulCoeffsMontgomeryAndSub(ciphertext.value[1], encryptor.sk.sk, ciphertext.value[0])

	// ct = [-s*a + m + e, a]
	context.Add(ciphertext.value[0], plaintext.value, ciphertext.value[0])

}
