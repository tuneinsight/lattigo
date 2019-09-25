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

// NewEncryptor creates a new Encryptor with the provided public-key and/or secret-key.
// This encryptor can be used to encrypt plaintexts, using the stored keys.
func (bfvcontext *BfvContext) NewEncryptor(pk *PublicKey, sk *SecretKey) (*Encryptor, error) {
	if uint64(pk.pk[0].GetDegree()+pk.pk[1].GetDegree())>>1 != bfvcontext.n {
		return nil, errors.New("error : pk ring degree doesn't match bfvcontext ring degree")
	}

	if sk != nil && uint64(sk.sk.GetDegree()) != bfvcontext.n {
		return nil, errors.New("error : sk ring degree doesn't match bfvcontext ring degree")
	}

	encryptor := new(Encryptor)
	encryptor.bfvcontext = bfvcontext
	encryptor.pk = pk
	encryptor.sk = sk
	encryptor.polypool = bfvcontext.contextQ.NewPoly()
	return encryptor, nil
}

// EncryptFromPkNew encrypts the input plaintext using the stored public-key and returns
// the result on a newly created ciphertext.
//
// ciphertext = [pk[0]*u + m + e_0, pk[1]*u + e_1]
func (encryptor *Encryptor) EncryptFromPkNew(plaintext *Plaintext) (ciphertext *Ciphertext, err error) {

	if encryptor.pk == nil {
		return nil, errors.New("cannot encrypt -> public-key has not been set")
	}

	if uint64(plaintext.value.GetDegree()) != encryptor.bfvcontext.n {
		return nil, errors.New("cannot encrypt -> plaintext ring degree doesn't match encryptor bfvcontext ring degree")
	}

	ciphertext = encryptor.bfvcontext.NewCiphertext(1)

	encryptfrompk(encryptor, plaintext, ciphertext)

	return ciphertext, nil
}

// EncryptFromPk encrypts the input plaintext using the stored public-key, and returns the result
// on the reciver ciphertext.
//
// ciphertext = [pk[0]*u + m + e_0, pk[1]*u + e_1]
func (encryptor *Encryptor) EncryptFromPk(plaintext *Plaintext, ciphertext *Ciphertext) (err error) {

	if encryptor.pk == nil {
		return errors.New("cannot encrypt -> public-key has not been set")
	}

	if uint64(plaintext.value.GetDegree()) != encryptor.bfvcontext.n {
		return errors.New("cannot encrypt -> plaintext ring degree doesn't match encryptor bfvcontext ring degree")
	}

	if ciphertext.Degree() != 1 {
		return errors.New("cannot encrypt -> invalide receiver -> ciphertext degree > 1")
	}

	encryptfrompk(encryptor, plaintext, ciphertext)

	return nil
}

// EncryptFromSkNew encrypts the input plaintext using the stored secret-key and returns
// the result on a newly created ciphertext. Encrypting with the secret-key introduces less noise than encrypting
// with the public-key.
//
// ciphertext = [-a*sk + m + e, a]
func (encryptor *Encryptor) EncryptFromSkNew(plaintext *Plaintext) (ciphertext *Ciphertext, err error) {

	if encryptor.sk == nil {
		return nil, errors.New("cannot encrypt -> secret-key has not been set")
	}

	if uint64(plaintext.value.GetDegree()) != encryptor.bfvcontext.n {
		return nil, errors.New("cannot encrypt -> plaintext ring degree doesn't match encryptor bfvcontext ring degree")
	}

	ciphertext = encryptor.bfvcontext.NewCiphertext(1)

	encryptfromsk(encryptor, plaintext, ciphertext)

	return ciphertext, nil
}

// EncryptFromSk encrypts the input plaintext using the stored secret-key, and returns the result
// on the reciever ciphertext. Encrypting with the secret-key introduces less noise than encrypting
// with the public-key.
//
// ciphertext = [-a*sk + m + e, a]
func (encryptor *Encryptor) EncryptFromSk(plaintext *Plaintext, ciphertext *Ciphertext) (err error) {

	if encryptor.sk == nil {
		return errors.New("cannot encrypt -> secret-key has not been set")
	}

	if uint64(plaintext.value.GetDegree()) != encryptor.bfvcontext.n {
		return errors.New("cannot encrypt -> plaintext ring degree doesn't match encryptor bfvcontext ring degree")
	}

	if ciphertext.Degree() != 1 {
		return errors.New("cannot encrypt -> invalide receiver -> ciphertext degree > 1")
	}

	encryptfromsk(encryptor, plaintext, ciphertext)

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
