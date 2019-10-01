package bfv

import (
	"errors"
	"github.com/ldsec/lattigo/ring"
)

// Encryptor is a structure holding the parameters needed to encrypt plaintexts.
type Encryptor struct {
	bfvcontext *BfvContext
	pk         *PublicKey
	sk         *SecretKey
	polypool   *ring.Poly
}

// NewEncryptorFromPk creates a new Encryptor with the provided public-key.
// This encryptor can be used to encrypt plaintexts, using the stored key.
func (bfvcontext *BfvContext) NewEncryptorFromPk(pk *PublicKey) (*Encryptor, error) {
	return bfvcontext.newEncryptor(pk, nil)
}

// NewEncryptorFromSk creates a new Encryptor with the provided secret-key.
// This encryptor can be used to encrypt plaintexts, using the stored key.
func (bfvcontext *BfvContext) NewEncryptorFromSk(sk *SecretKey) (*Encryptor, error) {
	return bfvcontext.newEncryptor(nil, sk)
}

func (bfvcontext *BfvContext) newEncryptor(pk *PublicKey, sk *SecretKey) (encryptor *Encryptor, err error) {

	if pk != nil && (uint64(pk.pk[0].GetDegree()) != bfvcontext.n || uint64(pk.pk[1].GetDegree()) != bfvcontext.n) {
		return nil, errors.New("error : pk ring degree doesn't match bfvcontext ring degree")
	}

	if sk != nil && uint64(sk.sk.GetDegree()) != bfvcontext.n {
		return nil, errors.New("error : sk ring degree doesn't match bfvcontext ring degree")
	}

	encryptor = new(Encryptor)
	encryptor.bfvcontext = bfvcontext
	encryptor.pk = pk
	encryptor.sk = sk
	encryptor.polypool = bfvcontext.contextQ.NewPoly()

	return encryptor, nil
}

// EncryptFromPkNew encrypts the input plaintext using the stored public-key and returns
// the result on a newly created ciphertext. It will encrypt the plaintext with the stored key, which can be
// private or public, a private-key encryption puts initial noise.
//
// encrypt with pk : ciphertext = [pk[0]*u + m + e_0, pk[1]*u + e_1]
// encrypt with sk : ciphertext = [-a*sk + m + e, a]
func (encryptor *Encryptor) EncryptNew(plaintext *Plaintext) (ciphertext *Ciphertext, err error) {

	ciphertext = encryptor.bfvcontext.NewCiphertext(1)

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

	context := encryptor.bfvcontext.contextQ

	// u
	encryptor.bfvcontext.ternarySampler.SampleMontgomeryNTT(0.5, encryptor.polypool)

	// ct[0] = pk[0]*u
	// ct[1] = pk[1]*u
	context.MulCoeffsMontgomery(encryptor.polypool, encryptor.pk.pk[0], ciphertext.value[0])
	context.MulCoeffsMontgomery(encryptor.polypool, encryptor.pk.pk[1], ciphertext.value[1])

	context.InvNTT(ciphertext.value[0], ciphertext.value[0])
	context.InvNTT(ciphertext.value[1], ciphertext.value[1])

	// ct[0] = pk[0]*u + e0
	encryptor.bfvcontext.gaussianSampler.Sample(encryptor.polypool)
	context.Add(ciphertext.value[0], encryptor.polypool, ciphertext.value[0])

	// ct[1] = pk[1]*u + e1
	encryptor.bfvcontext.gaussianSampler.Sample(encryptor.polypool)
	context.Add(ciphertext.value[1], encryptor.polypool, ciphertext.value[1])

	// ct[0] = pk[0]*u + e0 + m
	// ct[1] = pk[1]*u + e1
	context.Add(ciphertext.value[0], plaintext.value, ciphertext.value[0])

	encryptor.polypool.Zero()
}

func encryptfromsk(encryptor *Encryptor, plaintext *Plaintext, ciphertext *Ciphertext) {

	context := encryptor.bfvcontext.contextQ

	// ct = [-a*s , a]
	ciphertext.value[1] = context.NewUniformPoly()
	context.MulCoeffsMontgomery(ciphertext.value[1], encryptor.sk.sk, ciphertext.value[0])
	context.Neg(ciphertext.value[0], ciphertext.value[0])
	context.InvNTT(ciphertext.value[0], ciphertext.value[0])
	context.InvNTT(ciphertext.value[1], ciphertext.value[1])

	// ct = [-a*s + e, a]
	encryptor.bfvcontext.gaussianSampler.Sample(encryptor.polypool)
	context.Add(ciphertext.value[0], encryptor.polypool, ciphertext.value[0])

	// ct = [-a*s + m + e , a]
	context.Add(ciphertext.value[0], plaintext.value, ciphertext.value[0])

	encryptor.polypool.Zero()
}
