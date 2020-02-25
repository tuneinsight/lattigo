package bfv

import (
	"github.com/ldsec/lattigo/ring"
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

// encryptor is a structure that holds the parameters needed to encrypt plaintexts.
type encryptor struct {
	params     *Parameters
	bfvContext *bfvContext
	polypool   [3]*ring.Poly

	baseconverter *ring.FastBasisExtender
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
// This encryptor can be used to encrypt plaintexts, using the stored key.
func NewEncryptorFromPk(params *Parameters, pk *PublicKey) Encryptor {
	enc := newEncryptor(params)

	if uint64(pk.pk[0].GetDegree()) != uint64(1<<params.LogN) || uint64(pk.pk[1].GetDegree()) != uint64(1<<params.LogN) {
		panic("error: pk ring degree doesn't match bfvcontext ring degree")
	}

	return &pkEncryptor{enc, pk}
}

// NewEncryptorFromSk creates a new Encryptor with the provided secret-key.
// This encryptor can be used to encrypt plaintexts, using the stored key.
func NewEncryptorFromSk(params *Parameters, sk *SecretKey) Encryptor {
	enc := newEncryptor(params)

	if uint64(sk.sk.GetDegree()) != uint64(1<<params.LogN) {
		panic("error: sk ring degree doesn't match bfvcontext ring degree")
	}

	return &skEncryptor{enc, sk}
}

func newEncryptor(params *Parameters) encryptor {
	if !params.isValid {
		panic("cannot NewEncryptor: params not valid (check if they were generated properly)")
	}

	ctx := newBFVContext(params)
	qp := ctx.contextQP

	var baseconverter *ring.FastBasisExtender
	if len(params.Pi) != 0 {
		baseconverter = ring.NewFastBasisExtender(ctx.contextQ, ctx.contextP)
	}

	return encryptor{
		params:        params.Copy(),
		bfvContext:    ctx,
		polypool:      [3]*ring.Poly{qp.NewPoly(), qp.NewPoly(), qp.NewPoly()},
		baseconverter: baseconverter,
	}
}

func (encryptor *pkEncryptor) EncryptNew(plaintext *Plaintext) *Ciphertext {

	if encryptor.baseconverter == nil {
		panic("Cannot EncryptNew : modulus P is empty -> use instead EncryptFastNew")
	}

	ciphertext := NewCiphertext(encryptor.params, 1)
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
	ciphertext := NewCiphertext(encryptor.params, 1)
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

func (encryptor *pkEncryptor) encrypt(plaintext *Plaintext, ciphertext *Ciphertext, fast bool) {

	var ringContext *ring.Context

	if fast {
		ringContext = encryptor.bfvContext.contextQ

		ringContext.SampleTernaryMontgomeryNTT(encryptor.polypool[2], 0.5)

		ringContext.MulCoeffsMontgomery(encryptor.polypool[2], encryptor.pk.pk[0], encryptor.polypool[0])
		ringContext.MulCoeffsMontgomery(encryptor.polypool[2], encryptor.pk.pk[1], encryptor.polypool[1])

		ringContext.InvNTT(encryptor.polypool[0], encryptor.polypool[0])
		ringContext.InvNTT(encryptor.polypool[1], encryptor.polypool[1])

		// ct[0] = pk[0]*u + e0
		encryptor.bfvContext.gaussianSampler.Sample(encryptor.polypool[2])
		ringContext.Add(encryptor.polypool[0], encryptor.polypool[2], encryptor.polypool[0])

		// ct[1] = pk[1]*u + e1
		encryptor.bfvContext.gaussianSampler.Sample(encryptor.polypool[2])
		ringContext.Add(encryptor.polypool[1], encryptor.polypool[2], encryptor.polypool[1])

	} else {

		ringContext = encryptor.bfvContext.contextQP
		// u
		ringContext.SampleTernaryMontgomeryNTT(encryptor.polypool[2], 0.5)

		// ct[0] = pk[0]*u
		// ct[1] = pk[1]*u
		ringContext.MulCoeffsMontgomery(encryptor.polypool[2], encryptor.pk.pk[0], encryptor.polypool[0])
		ringContext.MulCoeffsMontgomery(encryptor.polypool[2], encryptor.pk.pk[1], encryptor.polypool[1])

		ringContext.InvNTT(encryptor.polypool[0], encryptor.polypool[0])
		ringContext.InvNTT(encryptor.polypool[1], encryptor.polypool[1])

		// ct[0] = pk[0]*u + e0
		encryptor.bfvContext.gaussianSampler.Sample(encryptor.polypool[2])
		ringContext.Add(encryptor.polypool[0], encryptor.polypool[2], encryptor.polypool[0])

		// ct[1] = pk[1]*u + e1
		encryptor.bfvContext.gaussianSampler.Sample(encryptor.polypool[2])
		ringContext.Add(encryptor.polypool[1], encryptor.polypool[2], encryptor.polypool[1])

		// We rescale the encryption of zero by the special prime, dividing the error by this prime
		encryptor.baseconverter.ModDownPQ(uint64(len(plaintext.Value()[0].Coeffs))-1, encryptor.polypool[0], ciphertext.value[0])
		encryptor.baseconverter.ModDownPQ(uint64(len(plaintext.Value()[0].Coeffs))-1, encryptor.polypool[1], ciphertext.value[1])
	}
	ringContext = encryptor.bfvContext.contextQ

	// ct[0] = pk[0]*u + e0 + m
	// ct[1] = pk[1]*u + e1
	ringContext.Add(ciphertext.value[0], plaintext.value, ciphertext.value[0])
}

func (encryptor *skEncryptor) EncryptNew(plaintext *Plaintext) *Ciphertext {

	if encryptor.baseconverter == nil {
		panic("Cannot EncryptNew : modulus P is empty -> use instead EncryptFastNew")
	}

	ciphertext := NewCiphertext(encryptor.params, 1)
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
	ciphertext := NewCiphertext(encryptor.params, 1)
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

	ciphertext := NewCiphertext(encryptor.params, 1)
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
	ciphertext := NewCiphertext(encryptor.params, 1)
	encryptor.EncryptFromCRPFast(plaintext, ciphertext, crp)
	return ciphertext
}

func (encryptor *skEncryptor) EncryptFromCRPFast(plaintext *Plaintext, ciphertext *Ciphertext, crp *ring.Poly) {
	encryptor.encryptFromCRP(plaintext, ciphertext, crp, true)

}

func (encryptor *skEncryptor) encryptSample(plaintext *Plaintext, ciphertext *Ciphertext, fast bool) {
	if fast {
		encryptor.bfvContext.contextQ.UniformPoly(encryptor.polypool[1])
	} else {
		encryptor.bfvContext.contextQP.UniformPoly(encryptor.polypool[1])
	}

	encryptor.encrypt(plaintext, ciphertext, encryptor.polypool[1], fast)
}

func (encryptor *skEncryptor) encryptFromCRP(plaintext *Plaintext, ciphertext *Ciphertext, crp *ring.Poly, fast bool) {
	if fast {
		encryptor.bfvContext.contextQ.Copy(crp, encryptor.polypool[1])
	} else {
		encryptor.bfvContext.contextQP.Copy(crp, encryptor.polypool[1])
	}

	encryptor.encrypt(plaintext, ciphertext, encryptor.polypool[1], fast)
}

func (encryptor *skEncryptor) encrypt(plaintext *Plaintext, ciphertext *Ciphertext, crp *ring.Poly, fast bool) {

	var ringContext *ring.Context

	if fast {

		ringContext = encryptor.bfvContext.contextQ

		ringContext.MulCoeffsMontgomery(crp, encryptor.sk.sk, ciphertext.value[0])
		ringContext.Neg(ciphertext.value[0], ciphertext.value[0])

		ringContext.InvNTT(ciphertext.value[0], ciphertext.value[0])
		ringContext.InvNTT(crp, ciphertext.value[1])

		ringContext.SampleGaussianAndAdd(ciphertext.value[0], encryptor.params.Sigma, uint64(6*encryptor.params.Sigma))

	} else {
		ringContext = encryptor.bfvContext.contextQP

		// ct = [(-a*s + e)/P , a/P]
		ringContext.MulCoeffsMontgomery(crp, encryptor.sk.sk, encryptor.polypool[0])
		ringContext.Neg(encryptor.polypool[0], encryptor.polypool[0])

		// We rescale the encryption of zero by the special prime, dividing the error by this prime
		ringContext.InvNTT(encryptor.polypool[0], encryptor.polypool[0])
		ringContext.InvNTT(crp, crp)

		encryptor.bfvContext.gaussianSampler.SampleAndAdd(encryptor.polypool[0])

		encryptor.baseconverter.ModDownPQ(uint64(len(plaintext.Value()[0].Coeffs))-1, encryptor.polypool[0], ciphertext.value[0])
		encryptor.baseconverter.ModDownPQ(uint64(len(plaintext.Value()[0].Coeffs))-1, crp, ciphertext.value[1])

		ringContext = encryptor.bfvContext.contextQ
	}

	ringContext = encryptor.bfvContext.contextQ

	// ct = [-a*s + m + e , a]
	ringContext.Add(ciphertext.value[0], plaintext.value, ciphertext.value[0])
}
