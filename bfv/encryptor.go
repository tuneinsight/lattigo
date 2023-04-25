package bfv

import (
	"github.com/tuneinsight/lattigo/ring"
)

// Encryptor is an interface for encryptors
//
// encrypt with pk: ciphertext = [pk[0]*u + m + e_0, pk[1]*u + e_1]
// encrypt with sk: ciphertext = [-a*sk + m + e, a]
type Encryptor interface {
	// EncryptNew encrypts the input plaintext using the stored key and returns
	// the result on a newly created ciphertext.
	EncryptNew(plaintext *Plaintext) *Ciphertext

	// Encrypt encrypts the input plaintext using the stored key, and returns
	// the result on the receiver ciphertext.
	Encrypt(plaintext *Plaintext, ciphertext *Ciphertext)
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

	return encryptor{
		params:        params.Copy(),
		bfvContext:    ctx,
		polypool:      [3]*ring.Poly{qp.NewPoly(), qp.NewPoly(), qp.NewPoly()},
		baseconverter: ring.NewFastBasisExtender(ctx.contextQ, ctx.contextP),
	}
}

func (encryptor *pkEncryptor) EncryptNew(plaintext *Plaintext) *Ciphertext {
	ciphertext := NewCiphertext(encryptor.params, 1)
	encryptor.Encrypt(plaintext, ciphertext)

	return ciphertext
}

func (encryptor *pkEncryptor) Encrypt(plaintext *Plaintext, ciphertext *Ciphertext) {
	ringContext := encryptor.bfvContext.contextQP

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

	ringContext = encryptor.bfvContext.contextQ

	// ct[0] = pk[0]*u + e0 + m
	// ct[1] = pk[1]*u + e1
	ringContext.Add(ciphertext.value[0], plaintext.value, ciphertext.value[0])
}

func (encryptor *skEncryptor) EncryptNew(plaintext *Plaintext) *Ciphertext {
	ciphertext := NewCiphertext(encryptor.params, 1)
	encryptor.Encrypt(plaintext, ciphertext)

	return ciphertext
}

func (encryptor *skEncryptor) Encrypt(plaintext *Plaintext, ciphertext *Ciphertext) {
	ringContext := encryptor.bfvContext.contextQP

	// ct = [(-a*s + e)/P , a/P]
	ringContext.UniformPoly(encryptor.polypool[1])
	ringContext.MulCoeffsMontgomeryAndSub(encryptor.polypool[1], encryptor.sk.sk, encryptor.polypool[0])

	// We rescale the encryption of zero by the special prime, dividing the error by this prime
	ringContext.InvNTT(encryptor.polypool[0], encryptor.polypool[0])
	ringContext.InvNTT(encryptor.polypool[1], encryptor.polypool[1])

	encryptor.bfvContext.gaussianSampler.SampleAndAdd(encryptor.polypool[0])

	encryptor.baseconverter.ModDownPQ(uint64(len(plaintext.Value()[0].Coeffs))-1, encryptor.polypool[0], ciphertext.value[0])
	encryptor.baseconverter.ModDownPQ(uint64(len(plaintext.Value()[0].Coeffs))-1, encryptor.polypool[1], ciphertext.value[1])

	ringContext = encryptor.bfvContext.contextQ

	// ct = [-a*s + m + e , a]
	ringContext.Add(ciphertext.value[0], plaintext.value, ciphertext.value[0])
}
