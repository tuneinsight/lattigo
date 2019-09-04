package ckks

import (
	"errors"
	"github.com/lca1/lattigo/ring"
	//"fmt"
)

type Encryptor struct {
	ckkscontext *CkksContext
	pk          *PublicKey
	sk          *SecretKey
	polypool    *ring.Poly
}

// NewEncryptor instanciates a new Encryptor with the provided public key.
// This encryptor can be used to encrypt plaintexts, using its the stored public key.
func (ckkscontext *CkksContext) NewEncryptor(pk *PublicKey, sk *SecretKey) (encryptor *Encryptor, err error) {

	if uint64(pk.pk[0].GetDegree()+pk.pk[1].GetDegree())>>1 != ckkscontext.n {
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

// EncryptFromPkNew encrypts the provided plaintext using the stored public-key and returns
// the ciphertext on a newly created element.
func (encryptor *Encryptor) EncryptFromPkNew(plaintext *Plaintext) (ciphertext *Ciphertext, err error) {

	if encryptor.pk == nil {
		return nil, errors.New("cannot encrypt -> public-key hasn't been set")
	}

	if uint64(plaintext.value[0].GetDegree()) != encryptor.ckkscontext.n {
		return nil, errors.New("cannot encrypt -> plaintext ring degree doesn't match encryptor ckkscontext ring degree")
	}

	ciphertext = encryptor.ckkscontext.NewCiphertext(1, plaintext.Level(), plaintext.Scale())

	encryptfrompk(encryptor, plaintext, ciphertext)

	return ciphertext, nil
}

// EncryptFromPk encrypts the provided plaintext using the stored public-key, and returns the result
// on the provided reciver ciphertext.
func (encryptor *Encryptor) EncryptFromPk(plaintext *Plaintext, ciphertext *Ciphertext) (err error) {

	if encryptor.pk == nil {
		return errors.New("cannot encrypt -> public-key hasn't been set")
	}

	if uint64(plaintext.value[0].GetDegree()) != encryptor.ckkscontext.n {
		return errors.New("cannot encrypt -> plaintext ring degree doesn't match encryptor ckkscontext ring degree")
	}

	if ciphertext.Degree() != 1 {
		return errors.New("cannot encrypt -> invalide receiver -> ciphertext degree > 1")
	}

	plaintext.CopyParams(ciphertext)

	encryptfrompk(encryptor, plaintext, ciphertext)

	return nil
}

// EncryptFromSkNew encrypts the provided plaintext using the stored secret-key and returns
// the ciphertext on a newly created element.
func (encryptor *Encryptor) EncryptFromSkNew(plaintext *Plaintext) (ciphertext *Ciphertext, err error) {

	if encryptor.sk == nil {
		return nil, errors.New("cannot encrypt -> secret-key hasn't been set")
	}

	if uint64(plaintext.value[0].GetDegree()) != encryptor.ckkscontext.n {
		return nil, errors.New("cannot encrypt -> plaintext ring degree doesn't match encryptor ckkscontext ring degree")
	}

	ciphertext = encryptor.ckkscontext.NewCiphertext(1, plaintext.Level(), plaintext.Scale())

	encryptfromsk(encryptor, plaintext, ciphertext)

	return ciphertext, nil
}

// EncryptFromSk encrypts the provided plaintext using the stored secret-key, and returns the result
// on the provided reciver ciphertext.
func (encryptor *Encryptor) EncryptFromSk(plaintext *Plaintext, ciphertext *Ciphertext) (err error) {

	if encryptor.sk == nil {
		return errors.New("cannot encrypt -> secret-key hasn't been set")
	}

	if uint64(plaintext.value[0].GetDegree()) != encryptor.ckkscontext.n {
		return errors.New("cannot encrypt -> plaintext ring degree doesn't match encryptor ckkscontext ring degree")
	}

	if ciphertext.Degree() != 1 {
		return errors.New("cannot encrypt -> invalide receiver -> ciphertext degree > 1")
	}

	plaintext.CopyParams(ciphertext)

	encryptfromsk(encryptor, plaintext, ciphertext)

	return nil
}

func encryptfrompk(encryptor *Encryptor, plaintext *Plaintext, ciphertext *Ciphertext) {

	context := encryptor.ckkscontext.contextLevel[plaintext.Level()]

	encryptor.ckkscontext.ternarySampler.SampleMontgomeryNTT(encryptor.polypool)

	context.MulCoeffsMontgomery(encryptor.polypool, encryptor.pk.pk[0], ciphertext.value[0])
	context.MulCoeffsMontgomery(encryptor.polypool, encryptor.pk.pk[1], ciphertext.value[1])

	encryptor.ckkscontext.gaussianSampler.SampleNTT(encryptor.polypool)
	context.Add(ciphertext.value[0], encryptor.polypool, ciphertext.value[0])

	encryptor.ckkscontext.gaussianSampler.SampleNTT(encryptor.polypool)
	context.Add(ciphertext.value[1], encryptor.polypool, ciphertext.value[1])

	context.Add(ciphertext.value[0], plaintext.value[0], ciphertext.value[0])

	encryptor.polypool.Zero()

	ciphertext.isNTT = true
}

func encryptfromsk(encryptor *Encryptor, plaintext *Plaintext, ciphertext *Ciphertext) {
	context := encryptor.ckkscontext.contextLevel[plaintext.Level()]

	// ct = [e, a]
	ciphertext.value[1] = context.NewUniformPoly()
	encryptor.ckkscontext.gaussianSampler.SampleNTT(ciphertext.value[1])

	// ct = [-s*a + e, a]
	context.MulCoeffsMontgomeryAndSub(ciphertext.value[1], encryptor.sk.sk, ciphertext.value[0])

	// ct = [-s*a + m + e, a]
	context.Add(ciphertext.value[0], plaintext.value[0], ciphertext.value[0])

}
