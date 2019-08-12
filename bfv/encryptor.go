package bfv

import (
	"errors"
	"github.com/lca1/lattigo/ring"
)

type Encryptor struct {
	bfvcontext *BfvContext
	pk         *PublicKey // secret key
	polypool   *ring.Poly
}

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

func (encryptor *Encryptor) EncryptNew(plaintext *Plaintext) (*Ciphertext, error) {

	if uint64(plaintext.value[0].GetDegree()) != encryptor.bfvcontext.n {
		return nil, errors.New("error : plaintext ring degree doesn't match encryptor bfvcontext ring degree")
	}

	ciphertext := encryptor.bfvcontext.NewCiphertext(1)

	encrypt(encryptor, plaintext, ciphertext)

	return ciphertext, nil
}

func (encryptor *Encryptor) Encrypt(plaintext *Plaintext, ciphertext *Ciphertext) error {

	if uint64(plaintext.value[0].GetDegree()) != encryptor.bfvcontext.n {
		return errors.New("error : plaintext ring degree doesn't match encryptor bfvcontext ring degree")
	}

	if ciphertext.Degree() != 1 {
		return errors.New("error : invalide receiver -> ciphertext degree > 1")
	}

	encrypt(encryptor, plaintext, ciphertext)

	return nil
}

func encrypt(encryptor *Encryptor, plaintext *Plaintext, ciphertext *Ciphertext) {

	context := encryptor.bfvcontext.contextQ

	encryptor.bfvcontext.ternarySampler.SampleMontgomeryNTT(encryptor.polypool)

	context.MulCoeffsMontgomery(encryptor.polypool, encryptor.pk.pk[0], ciphertext.value[0])
	context.MulCoeffsMontgomery(encryptor.polypool, encryptor.pk.pk[1], ciphertext.value[1])

	context.InvNTT(ciphertext.value[0], ciphertext.value[0])
	context.InvNTT(ciphertext.value[1], ciphertext.value[1])

	context.Add(ciphertext.value[0], plaintext.value[0], ciphertext.value[0])

	encryptor.bfvcontext.gaussianSampler.Sample(encryptor.polypool)
	context.Add(ciphertext.value[0], encryptor.polypool, ciphertext.value[0])

	encryptor.bfvcontext.gaussianSampler.Sample(encryptor.polypool)
	context.Add(ciphertext.value[1], encryptor.polypool, ciphertext.value[1])

	encryptor.polypool.Zero()
}
