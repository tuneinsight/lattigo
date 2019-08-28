package bfv

import (
	"errors"
	"github.com/lca1/lattigo/ring"
)

// Encryptor is a structure holding the parameters needed to encrypt plaintexts.
type Encryptor struct {
	bfvcontext *BfvContext
	pk         *PublicKey // secret key
	polypool   *ring.Poly
}

// NewEncryptor creates a new Encryptor from the target bfvcontext, using the input public-key.
func (bfvcontext *BfvContext) NewEncryptor(pk *PublicKey) (*Encryptor, error) {
	if uint64(pk.pk[0].GetDegree()+pk.pk[1].GetDegree())>>1 != bfvcontext.n {
		return nil, errors.New("error : pk ring degree doesn't match bfvcontext ring degree")
	}
	//TODO : check that pk degree matches bfvContext degree
	encryptor := new(Encryptor)
	encryptor.bfvcontext = bfvcontext
	encryptor.pk = pk
	encryptor.polypool = bfvcontext.contextQ.NewPoly()
	return encryptor, nil
}

// EncryptNew encrypts the input plaintext and returns the result on a new ciphertext.
func (encryptor *Encryptor) EncryptNew(plaintext *Plaintext) (ciphertext *Ciphertext, err error) {

	if uint64(plaintext.value[0].GetDegree()) != encryptor.bfvcontext.n {
		return nil, errors.New("error : plaintext ring degree doesn't match encryptor bfvcontext ring degree")
	}

	ciphertext = encryptor.bfvcontext.NewCiphertext(1)

	encrypt(encryptor, plaintext, ciphertext)

	return ciphertext, nil
}

// Encrypt encrypts the input plaintext and returns the result on the receiver ciphertext.
func (encryptor *Encryptor) Encrypt(plaintext *Plaintext, ciphertext *Ciphertext) (err error) {

	if uint64(plaintext.value[0].GetDegree()) != encryptor.bfvcontext.n {
		return errors.New("error : plaintext ring degree doesn't match encryptor bfvcontext ring degree")
	}

	if ciphertext.Degree() != 1 {
		return errors.New("error : invalide receiver -> ciphertext degree > 1")
	}

	encrypt(encryptor, plaintext, ciphertext)

	return nil
}

// encrypt performes the bfv enryption.
func encrypt(encryptor *Encryptor, plaintext *Plaintext, ciphertext *Ciphertext) {

	context := encryptor.bfvcontext.contextQ

	// u
	encryptor.bfvcontext.ternarySampler.SampleMontgomeryNTT(encryptor.polypool)

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
	context.Add(ciphertext.value[0], plaintext.value[0], ciphertext.value[0])

	encryptor.polypool.Zero()
}
