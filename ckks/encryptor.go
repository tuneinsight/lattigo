package ckks

import (
	"github.com/ldsec/lattigo/ring"
	"github.com/ldsec/lattigo/utils"
)

// Encryptor in an interface for encryptors
//
// encrypt with pk : ciphertext = [pk[0]*u + m + e_0, pk[1]*u + e_1]
// encrypt with sk : ciphertext = [-a*sk + m + e, a]
type Encryptor interface {
	// EncryptNew encrypts the input plaintext using the stored key and returns
	// the result on a newly created ciphertext. The encryption is done by first
	// encrypting zero in QP, dividing by P and then adding the plaintext.
	EncryptNew(plaintext *Plaintext) *Ciphertext

	// Encrypt encrypts the input plaintext using the stored key, and returns
	// the result on the receiver ciphertext. The encryption is done by first
	// encrypting zero in QP, dividing by P and then adding the plaintext.
	Encrypt(plaintext *Plaintext, ciphertext *Ciphertext)

	// EncryptFastNew encrypts the input plaintext using the stored key and returns
	// the result on a newly created ciphertext. The encryption is done by first
	// encrypting zero in Q and then adding the plaintext.
	EncryptFastNew(plaintext *Plaintext) *Ciphertext

	// EncryptFsat encrypts the input plaintext using the stored-key, and returns
	// the result onthe receiver ciphertext. The encryption is done by first
	// encrypting zero in Q and then adding the plaintext.
	EncryptFast(plaintext *Plaintext, ciphertext *Ciphertext)

	// EncryptFromCRPNew encrypts the input plaintext using the stored key and returns
	// the result on a newly created ciphertext. The encryption is done by first encrypting
	// zero in QP, using the provided polynomial as the uniform polynomial, dividing by P and
	// then adding the plaintext.
	EncryptFromCRPNew(plaintext *Plaintext, crp *ring.Poly) *Ciphertext

	// EncryptFromCRP encrypts the input plaintext using the stored key and returns
	// the result tge receiver ciphertext. The encryption is done by first encrypting
	// zero in QP, using the provided polynomial as the uniform polynomial, dividing by P and
	// then adding the plaintext.
	EncryptFromCRP(plaintext *Plaintext, ciphertetx *Ciphertext, crp *ring.Poly)

	// EncryptFromCRPNew encrypts the input plaintext using the stored key and returns
	// the result on a newly created ciphertext. The encryption is done by first encrypting
	// zero in Q, using the provided polynomial as the uniform polynomial, and
	// then adding the plaintext.
	EncryptFromCRPFastNew(plaintext *Plaintext, crp *ring.Poly) *Ciphertext

	// EncryptFromCRP encrypts the input plaintext using the stored key and returns
	// the result tge receiver ciphertext. The encryption is done by first encrypting
	// zero in Q, using the provided polynomial as the uniform polynomial, and
	// then adding the plaintext.
	EncryptFromCRPFast(plaintext *Plaintext, ciphertetx *Ciphertext, crp *ring.Poly)
}

// encryptor is a struct used to encrypt Plaintexts. It stores the public-key and/or secret-key.
type encryptor struct {
	params      *Parameters
	ckksContext *Context
	polypool    [3]*ring.Poly

	baseconverter     *ring.FastBasisExtender
	gaussianSamplerQ  *ring.GaussianSampler
	uniformSamplerQ   *ring.UniformSampler
	ternarySamplerQ   *ring.TernarySampler
	gaussianSamplerQP *ring.GaussianSampler
	uniformSamplerQP  *ring.UniformSampler
	ternarySamplerQP  *ring.TernarySampler
}

type pkEncryptor struct {
	encryptor
	pk *PublicKey
}

type skEncryptor struct {
	encryptor
	sk *SecretKey
}

// NewEncryptorFromPk creates a new Encryptor with the provided public-key.
// This Encryptor can be used to encrypt Plaintexts, using the stored key.
func NewEncryptorFromPk(params *Parameters, pk *PublicKey) Encryptor {
	enc := newEncryptor(params)

	if uint64(pk.pk[0].GetDegree()) != uint64(1<<params.LogN) || uint64(pk.pk[1].GetDegree()) != uint64(1<<params.LogN) {
		panic("cannot newEncrpytor: pk ring degree does not match params ring degree")
	}

	return &pkEncryptor{enc, pk}
}

// NewEncryptorFromSk creates a new Encryptor with the provided secret-key.
// This Encryptor can be used to encrypt Plaintexts, using the stored key.
func NewEncryptorFromSk(params *Parameters, sk *SecretKey) Encryptor {
	enc := newEncryptor(params)

	if uint64(sk.sk.GetDegree()) != uint64(1<<params.LogN) {
		panic("cannot newEncryptor: sk ring degree does not match params ring degree")
	}

	return &skEncryptor{enc, sk}
}

func newEncryptor(params *Parameters) encryptor {
	if !params.isValid {
		panic("cannot newEncryptor: parameters are invalid (check if the generation was done properly)")
	}

	ctx := newContext(params)
	qp := ctx.contextQP

	var baseconverter *ring.FastBasisExtender
	if len(params.Pi) != 0 {
		baseconverter = ring.NewFastBasisExtender(ctx.contextQ, ctx.contextP)
	}

	prng, err := utils.NewPRNG()
	if err != nil {
		panic(err)
	}
	gaussianSamplerQ := ring.NewGaussianSampler(prng, ctx.contextQ)
	uniformSamplerQ := ring.NewUniformSampler(prng, ctx.contextQ)
	gaussianSamplerQP := ring.NewGaussianSampler(prng, ctx.contextQP)
	uniformSamplerQP := ring.NewUniformSampler(prng, ctx.contextQP)

	return encryptor{
		params:            params.Copy(),
		ckksContext:       ctx,
		polypool:          [3]*ring.Poly{qp.NewPoly(), qp.NewPoly(), qp.NewPoly()},
		baseconverter:     baseconverter,
		gaussianSamplerQ:  gaussianSamplerQ,
		uniformSamplerQ:   uniformSamplerQ,
		ternarySamplerQ:   ring.NewTernarySampler(prng, ctx.contextQ, 0.5, true),
		gaussianSamplerQP: gaussianSamplerQP,
		uniformSamplerQP:  uniformSamplerQP,
		ternarySamplerQP:  ring.NewTernarySampler(prng, ctx.contextQP, 0.5, true),
	}
}

// EncryptNew encrypts the input Plaintext using the stored key and returns
// the result on a newly created Ciphertext.
//
// encrypt with pk: ciphertext = [pk[0]*u + m + e_0, pk[1]*u + e_1]
// encrypt with sk: ciphertext = [-a*sk + m + e, a]
func (encryptor *pkEncryptor) EncryptNew(plaintext *Plaintext) *Ciphertext {

	if encryptor.baseconverter == nil {
		panic("Cannot EncryptNew : modulus P is empty -> use instead EncryptFastNew")
	}

	ciphertext := NewCiphertext(encryptor.params, 1, plaintext.Level(), plaintext.Scale())
	encryptor.encrypt(plaintext, ciphertext, false)

	return ciphertext
}

func (encryptor *pkEncryptor) Encrypt(plaintext *Plaintext, ciphertext *Ciphertext) {

	if encryptor.baseconverter == nil {
		panic("Cannot Encrypt : modulus P is empty -> use instead EncryptFast")
	}

	encryptor.encrypt(plaintext, ciphertext, false)
}

func (encryptor *pkEncryptor) EncryptFastNew(plaintext *Plaintext) *Ciphertext {
	ciphertext := NewCiphertext(encryptor.params, 1, plaintext.Level(), plaintext.Scale())
	encryptor.encrypt(plaintext, ciphertext, true)

	return ciphertext
}

func (encryptor *pkEncryptor) EncryptFast(plaintext *Plaintext, ciphertext *Ciphertext) {
	encryptor.encrypt(plaintext, ciphertext, true)
}

func (encryptor *pkEncryptor) EncryptFromCRP(plaintext *Plaintext, ciphertext *Ciphertext, crp *ring.Poly) {
	panic("Cannot encrypt with CRP using an encryptor created with the public-key")
}

func (encryptor *pkEncryptor) EncryptFromCRPNew(plaintext *Plaintext, crp *ring.Poly) *Ciphertext {
	panic("Cannot encrypt with CRP using an encryptor created with the public-key")
}

func (encryptor *pkEncryptor) EncryptFromCRPFast(plaintext *Plaintext, ciphertext *Ciphertext, crp *ring.Poly) {
	panic("Cannot encrypt with CRP using an encryptor created with the public-key")
}

func (encryptor *pkEncryptor) EncryptFromCRPFastNew(plaintext *Plaintext, crp *ring.Poly) *Ciphertext {
	panic("Cannot encrypt with CRP using an encryptor created with the public-key")
}

// Encrypt encrypts the input Plaintext using the stored key, and returns the result
// on the receiver Ciphertext.
//
// encrypt with pk: ciphertext = [pk[0]*u + m + e_0, pk[1]*u + e_1]
// encrypt with sk: ciphertext = [-a*sk + m + e, a]
func (encryptor *pkEncryptor) encrypt(plaintext *Plaintext, ciphertext *Ciphertext, fast bool) {

	// We sample a R-WLE instance (encryption of zero) over the keys context (ciphertext context + special prime)

	contextQ := encryptor.ckksContext.contextQ

	if fast {

		level := uint64(len(contextQ.Modulus) - 1)

		encryptor.ternarySamplerQ.ReadNTT(encryptor.polypool[2])

		// ct0 = u*pk0
		contextQ.MulCoeffsMontgomery(encryptor.polypool[2], encryptor.pk.pk[0], ciphertext.value[0])
		// ct1 = u*pk1
		contextQ.MulCoeffsMontgomery(encryptor.polypool[2], encryptor.pk.pk[1], ciphertext.value[1])

		// ct1 = u*pk1 + e1
		encryptor.gaussianSamplerQ.SampleNTTLvl(level, encryptor.polypool[0], encryptor.params.Sigma, uint64(6*encryptor.params.Sigma))
		contextQ.Add(ciphertext.value[1], encryptor.polypool[0], ciphertext.value[1])

		if !plaintext.isNTT {

			// ct0 = u*pk0 + e0
			encryptor.gaussianSamplerQ.SampleLvl(level, encryptor.polypool[0], encryptor.params.Sigma, uint64(6*encryptor.params.Sigma))
			// ct0 = (u*pk0 + e0)/P + m
			contextQ.Add(encryptor.polypool[0], plaintext.value, encryptor.polypool[0])
			contextQ.NTT(encryptor.polypool[0], encryptor.polypool[0])
			contextQ.Add(ciphertext.value[0], encryptor.polypool[0], ciphertext.value[0])

		} else {
			// ct0 = u*pk0 + e0
			encryptor.gaussianSamplerQ.SampleNTTLvl(level, encryptor.polypool[0], encryptor.params.Sigma, uint64(6*encryptor.params.Sigma))
			contextQ.Add(ciphertext.value[0], encryptor.polypool[0], ciphertext.value[0])
			contextQ.Add(ciphertext.value[0], plaintext.value, ciphertext.value[0])
		}

	} else {

		contextQP := encryptor.ckksContext.contextQP

		level := uint64(len(contextQP.Modulus) - 1)

		encryptor.ternarySamplerQP.ReadNTT(encryptor.polypool[2])

		// ct0 = u*pk0
		contextQP.MulCoeffsMontgomery(encryptor.polypool[2], encryptor.pk.pk[0], encryptor.polypool[0])
		// ct1 = u*pk1
		contextQP.MulCoeffsMontgomery(encryptor.polypool[2], encryptor.pk.pk[1], encryptor.polypool[1])

		// 2*(#Q + #P) NTT
		contextQP.InvNTT(encryptor.polypool[0], encryptor.polypool[0])
		contextQP.InvNTT(encryptor.polypool[1], encryptor.polypool[1])

		// ct0 = u*pk0 + e0
		encryptor.gaussianSamplerQP.SampleAndAddLvl(level, encryptor.polypool[0], encryptor.params.Sigma, uint64(6*encryptor.params.Sigma))
		// ct1 = u*pk1 + e1
		encryptor.gaussianSamplerQP.SampleAndAddLvl(level, encryptor.polypool[1], encryptor.params.Sigma, uint64(6*encryptor.params.Sigma))

		// ct0 = (u*pk0 + e0)/P
		encryptor.baseconverter.ModDownPQ(plaintext.Level(), encryptor.polypool[0], ciphertext.value[0])

		// ct1 = (u*pk1 + e1)/P
		encryptor.baseconverter.ModDownPQ(plaintext.Level(), encryptor.polypool[1], ciphertext.value[1])

		if !plaintext.isNTT {
			contextQ.Add(ciphertext.value[0], plaintext.value, ciphertext.value[0])
		}

		// 2*#Q NTT
		contextQ.NTT(ciphertext.value[0], ciphertext.value[0])
		contextQ.NTT(ciphertext.value[1], ciphertext.value[1])

		if plaintext.isNTT {
			// ct0 = (u*pk0 + e0)/P + m
			contextQ.Add(ciphertext.value[0], plaintext.value, ciphertext.value[0])
		}
	}

	ciphertext.isNTT = true
}

func (encryptor *skEncryptor) EncryptNew(plaintext *Plaintext) *Ciphertext {

	if encryptor.baseconverter == nil {
		panic("Cannot EncryptNew : modulus P is empty -> use instead EncryptFastNew")
	}

	ciphertext := NewCiphertext(encryptor.params, 1, plaintext.Level(), plaintext.Scale())
	encryptor.Encrypt(plaintext, ciphertext)
	return ciphertext
}

func (encryptor *skEncryptor) Encrypt(plaintext *Plaintext, ciphertext *Ciphertext) {

	if encryptor.baseconverter == nil {
		panic("Cannot Encrypt : modulus P is empty -> use instead EncryptFast")
	}

	encryptor.encryptSample(plaintext, ciphertext, false)
}

func (encryptor *skEncryptor) EncryptFastNew(plaintext *Plaintext) *Ciphertext {
	ciphertext := NewCiphertext(encryptor.params, 1, plaintext.Level(), plaintext.Scale())
	encryptor.EncryptFast(plaintext, ciphertext)
	return ciphertext
}

func (encryptor *skEncryptor) EncryptFast(plaintext *Plaintext, ciphertext *Ciphertext) {
	encryptor.encryptSample(plaintext, ciphertext, true)
}

func (encryptor *skEncryptor) EncryptFromCRPNew(plaintext *Plaintext, crp *ring.Poly) *Ciphertext {

	if encryptor.baseconverter == nil {
		panic("Cannot EncryptFromCRPNew : modulus P is empty -> use instead EncryptFromCRPFastNew")
	}

	ciphertext := NewCiphertext(encryptor.params, 1, plaintext.Level(), plaintext.Scale())
	encryptor.EncryptFromCRP(plaintext, ciphertext, crp)
	return ciphertext
}

func (encryptor *skEncryptor) EncryptFromCRP(plaintext *Plaintext, ciphertext *Ciphertext, crp *ring.Poly) {

	if encryptor.baseconverter == nil {
		panic("Cannot EncryptFromCRP : modulus P is empty -> use instead EncryptFromCRPFast")
	}

	encryptor.encryptFromCRP(plaintext, ciphertext, crp, false)
}

func (encryptor *skEncryptor) EncryptFromCRPFastNew(plaintext *Plaintext, crp *ring.Poly) *Ciphertext {
	ciphertext := NewCiphertext(encryptor.params, 1, plaintext.Level(), plaintext.Scale())
	encryptor.EncryptFromCRPFast(plaintext, ciphertext, crp)
	return ciphertext
}

func (encryptor *skEncryptor) EncryptFromCRPFast(plaintext *Plaintext, ciphertext *Ciphertext, crp *ring.Poly) {
	encryptor.encryptFromCRP(plaintext, ciphertext, crp, true)

}

func (encryptor *skEncryptor) encryptSample(plaintext *Plaintext, ciphertext *Ciphertext, fast bool) {
	if fast {
		encryptor.uniformSamplerQ.Sample(ciphertext.value[1])
		encryptor.encrypt(plaintext, ciphertext, ciphertext.value[1], fast)
	} else {
		encryptor.uniformSamplerQP.Sample(encryptor.polypool[1])
		encryptor.encrypt(plaintext, ciphertext, encryptor.polypool[1], fast)
	}

}

func (encryptor *skEncryptor) encryptFromCRP(plaintext *Plaintext, ciphertext *Ciphertext, crp *ring.Poly, fast bool) {
	if fast {
		encryptor.ckksContext.contextQ.Copy(crp, ciphertext.value[1])
		encryptor.encrypt(plaintext, ciphertext, ciphertext.value[1], fast)
	} else {
		encryptor.ckksContext.contextQP.Copy(crp, encryptor.polypool[1])
		encryptor.encrypt(plaintext, ciphertext, encryptor.polypool[1], fast)
	}

}

func (encryptor *skEncryptor) encrypt(plaintext *Plaintext, ciphertext *Ciphertext, crp *ring.Poly, fast bool) {

	contextQ := encryptor.ckksContext.contextQ

	if fast {

		level := uint64(len(contextQ.Modulus) - 1)

		contextQ.MulCoeffsMontgomery(ciphertext.value[1], encryptor.sk.sk, ciphertext.value[0])
		contextQ.Neg(ciphertext.value[0], ciphertext.value[0])

		if plaintext.isNTT {
			encryptor.gaussianSamplerQ.SampleNTTLvl(level, encryptor.polypool[0], encryptor.params.Sigma, uint64(6*encryptor.params.Sigma))
			contextQ.Add(ciphertext.value[0], encryptor.polypool[0], ciphertext.value[0])
			contextQ.Add(ciphertext.value[0], plaintext.value, ciphertext.value[0])
		} else {
			encryptor.gaussianSamplerQ.SampleLvl(level, encryptor.polypool[0], encryptor.params.Sigma, uint64(6*encryptor.params.Sigma))
			contextQ.Add(encryptor.polypool[0], plaintext.value, encryptor.polypool[0])
			contextQ.NTT(encryptor.polypool[0], encryptor.polypool[0])
			contextQ.Add(ciphertext.value[0], encryptor.polypool[0], ciphertext.value[0])
		}

	} else {

		contextQP := encryptor.ckksContext.contextQP

		level := uint64(len(contextQP.Modulus) - 1)

		// ct0 = -s*a
		contextQP.MulCoeffsMontgomery(crp, encryptor.sk.sk, encryptor.polypool[0])
		contextQP.Neg(encryptor.polypool[0], encryptor.polypool[0])

		// #Q + #P NTT
		contextQP.InvNTT(encryptor.polypool[0], encryptor.polypool[0])

		// ct0 = -s*a + e
		encryptor.gaussianSamplerQP.SampleAndAddLvl(level, encryptor.polypool[0], encryptor.params.Sigma, uint64(6*encryptor.params.Sigma))

		// We rescale by the special prime, dividing the error by this prime
		// ct0 = (-s*a + e)/P
		encryptor.baseconverter.ModDownPQ(plaintext.Level(), encryptor.polypool[0], ciphertext.value[0])

		// #Q + #P NTT
		// ct1 = a/P
		encryptor.baseconverter.ModDownNTTPQ(plaintext.Level(), crp, ciphertext.value[1])

		if !plaintext.isNTT {
			contextQ.Add(ciphertext.value[0], plaintext.value, ciphertext.value[0])
		}

		// #Q NTT
		contextQ.NTT(ciphertext.value[0], ciphertext.value[0])

		if plaintext.isNTT {
			// ct0 = (u*pk0 + e0)/P + m
			contextQ.Add(ciphertext.value[0], plaintext.value, ciphertext.value[0])
		}
	}

	ciphertext.isNTT = true
}
