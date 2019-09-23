package ckks

import (
	"errors"
	"github.com/ldsec/lattigo/ring"
	//"fmt"
)

type Encryptor struct {
	ckkscontext *CkksContext
	pk          *PublicKey // secret key
	polypool    *ring.Poly
}

// NewEncryptor instanciates a new Encryptor with the provided public key.
// This encryptor can be used to encrypt plaintexts, using its the stored public key.
func (ckkscontext *CkksContext) NewEncryptor(pk *PublicKey) (*Encryptor, error) {
	if uint64(pk.pk[0].GetDegree()+pk.pk[1].GetDegree())>>1 != ckkscontext.n {
		return nil, errors.New("error : pk ring degree doesn't match ckkscontext ring degree")
	}
	//TODO : check that pk degree matches ckksContext degree
	encryptor := new(Encryptor)
	encryptor.ckkscontext = ckkscontext
	encryptor.pk = pk
	encryptor.polypool = ckkscontext.contextLevel[ckkscontext.levels-1].NewPoly()
	return encryptor, nil
}

// EncryptNew encrypts the provided plaintext under the stored public key and returns
// the ciphertext on a newly created element.
func (encryptor *Encryptor) EncryptNew(plaintext *Plaintext) (*Ciphertext, error) {

	if uint64(plaintext.value[0].GetDegree()) != encryptor.ckkscontext.n {
		return nil, errors.New("error : plaintext ring degree doesn't match encryptor ckkscontext ring degree")
	}

	ciphertext := encryptor.ckkscontext.NewCiphertext(1, plaintext.Level(), plaintext.Scale())

	encrypt(encryptor, plaintext, ciphertext)

	return ciphertext, nil
}

// Encrypt encrypts the provided plaitnext under the stored public key, and returns the result
// on the provided reciver ciphertext.
func (encryptor *Encryptor) Encrypt(plaintext *Plaintext, ciphertext *Ciphertext) error {

	if uint64(plaintext.value[0].GetDegree()) != encryptor.ckkscontext.n {
		return errors.New("error : plaintext ring degree doesn't match encryptor ckkscontext ring degree")
	}

	if ciphertext.Degree() != 1 {
		return errors.New("error : invalide receiver -> ciphertext degree > 1")
	}

	plaintext.CopyParams(ciphertext)

	encrypt(encryptor, plaintext, ciphertext)

	return nil
}

func encrypt(encryptor *Encryptor, plaintext *Plaintext, ciphertext *Ciphertext) {

	context := encryptor.ckkscontext.contextLevel[plaintext.Level()]

	encryptor.ckkscontext.ternarySampler.SampleMontgomeryNTT(encryptor.polypool)

	context.MulCoeffsMontgomery(encryptor.polypool, encryptor.pk.pk[0], ciphertext.value[0])
	context.MulCoeffsMontgomery(encryptor.polypool, encryptor.pk.pk[1], ciphertext.value[1])

	context.Add(ciphertext.value[0], plaintext.value[0], ciphertext.value[0])

	encryptor.ckkscontext.gaussianSampler.SampleNTT(encryptor.polypool)
	context.Add(ciphertext.value[0], encryptor.polypool, ciphertext.value[0])

	encryptor.ckkscontext.gaussianSampler.SampleNTT(encryptor.polypool)
	context.Add(ciphertext.value[1], encryptor.polypool, ciphertext.value[1])

	encryptor.polypool.Zero()

	ciphertext.isNTT = true

}
