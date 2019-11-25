package bfv

import (
	"math/big"

	"github.com/ldsec/lattigo/ring"
)

type encryptorContext struct {
	// Polynomial degree
	n uint64

	// Ternary and Gaussian samplers
	gaussianSampler *ring.KYSampler

	// Polynomial contexts
	contextQ *ring.Context

	contextKeys       *ring.Context
	specialPrimes     []uint64
	rescaleParamsKeys []uint64 // (P^-1) mod each qi
}

func newEncryptorContext(params *Parameters) *encryptorContext {
	n := params.N

	contextQ := ring.NewContext()
	contextQ.SetParameters(n, params.Qi)

	contextKeys := ring.NewContext()
	contextKeys.SetParameters(n, append(params.Qi, params.KeySwitchPrimes...))
	err := contextKeys.GenNTTParams()
	if err != nil {
		panic(err)
	}

	gaussianSampler := contextKeys.NewKYSampler(params.Sigma, int(6*params.Sigma))

	specialPrimes := make([]uint64, len(params.KeySwitchPrimes))
	for i := range params.KeySwitchPrimes {
		specialPrimes[i] = params.KeySwitchPrimes[i]
	}

	rescaleParamsKeys := make([]uint64, len(params.Qi))
	PBig := ring.NewUint(1)
	for _, pj := range specialPrimes {
		PBig.Mul(PBig, ring.NewUint(pj))
	}

	tmp := new(big.Int)
	bredParams := contextQ.GetBredParams()
	for i, Qi := range params.Qi {
		tmp.Mod(PBig, ring.NewUint(Qi))
		rescaleParamsKeys[i] = ring.MForm(ring.ModExp(ring.BRedAdd(tmp.Uint64(), Qi, bredParams[i]), Qi-2, Qi), Qi, bredParams[i])
	}

	return &encryptorContext{
		n:                 n,
		gaussianSampler:   gaussianSampler,
		contextQ:          contextQ,
		contextKeys:       contextKeys,
		specialPrimes:     specialPrimes,
		rescaleParamsKeys: rescaleParamsKeys,
	}
}

// Encryptor is a structure holding the parameters needed to encrypt plaintexts.
type Encryptor struct {
	context  *encryptorContext
	pk       *PublicKey
	sk       *SecretKey
	polypool [3]*ring.Poly

	baseconverter *ring.FastBasisExtender
}

// NewEncryptorFromPk creates a new Encryptor with the provided public-key.
// This encryptor can be used to encrypt plaintexts, using the stored key.
func NewEncryptorFromPk(pk *PublicKey, params *Parameters) *Encryptor {
	return newEncryptor(pk, nil, params)
}

// NewEncryptorFromSk creates a new Encryptor with the provided secret-key.
// This encryptor can be used to encrypt plaintexts, using the stored key.
func NewEncryptorFromSk(sk *SecretKey, params *Parameters) *Encryptor {
	return newEncryptor(nil, sk, params)
}

func newEncryptor(pk *PublicKey, sk *SecretKey, params *Parameters) (encryptor *Encryptor) {
	context := newEncryptorContext(params)

	if pk != nil && (uint64(pk.pk[0].GetDegree()) != context.n || uint64(pk.pk[1].GetDegree()) != context.n) {
		panic("error : pk ring degree doesn't match bfvcontext ring degree")
	}

	if sk != nil && uint64(sk.sk.GetDegree()) != context.n {
		panic("error : sk ring degree doesn't match bfvcontext ring degree")
	}

	encryptor = new(Encryptor)
	encryptor.context = context
	encryptor.pk = pk
	encryptor.sk = sk

	encryptor.polypool[0] = context.contextKeys.NewPoly()
	encryptor.polypool[1] = context.contextKeys.NewPoly()
	encryptor.polypool[2] = context.contextKeys.NewPoly()

	encryptor.baseconverter = ring.NewFastBasisExtender(context.contextQ.Modulus, context.specialPrimes)

	return
}

// EncryptNew encrypts the input plaintext using the stored key and returns
// the result on a newly created ciphertext.
//
// encrypt with pk : ciphertext = [pk[0]*u + m + e_0, pk[1]*u + e_1]
// encrypt with sk : ciphertext = [-a*sk + m + e, a]
func (encryptor *Encryptor) EncryptNew(plaintext *Plaintext) (ciphertext *Ciphertext) {

	ciphertext = NewCiphertext(1, encryptor.context.contextQ)
	encryptor.Encrypt(plaintext, ciphertext)
	return
}

// Encrypt encrypts the input plaintext using the stored key, and returns the result
// on the receiver ciphertext.
//
// encrypt with pk : ciphertext = [pk[0]*u + m + e_0, pk[1]*u + e_1]
// encrypt with sk : ciphertext = [-a*sk + m + e, a]
func (encryptor *Encryptor) Encrypt(plaintext *Plaintext, ciphertext *Ciphertext) {

	if encryptor.sk != nil {

		encryptfromsk(encryptor, plaintext, ciphertext)

	} else if encryptor.pk != nil {

		encryptfrompk(encryptor, plaintext, ciphertext)

	} else {

		panic("cannot encrypt -> public-key and/or secret-key has not been set")
	}
}

func encryptfrompk(encryptor *Encryptor, plaintext *Plaintext, ciphertext *Ciphertext) {

	ringContext := encryptor.context.contextKeys

	// u
	ringContext.SampleTernaryMontgomeryNTT(encryptor.polypool[2], 0.5)

	// ct[0] = pk[0]*u
	// ct[1] = pk[1]*u
	ringContext.MulCoeffsMontgomery(encryptor.polypool[2], encryptor.pk.pk[0], encryptor.polypool[0])
	ringContext.MulCoeffsMontgomery(encryptor.polypool[2], encryptor.pk.pk[1], encryptor.polypool[1])

	ringContext.InvNTT(encryptor.polypool[0], encryptor.polypool[0])
	ringContext.InvNTT(encryptor.polypool[1], encryptor.polypool[1])

	// ct[0] = pk[0]*u + e0
	encryptor.context.gaussianSampler.Sample(encryptor.polypool[2])
	ringContext.Add(encryptor.polypool[0], encryptor.polypool[2], encryptor.polypool[0])

	// ct[1] = pk[1]*u + e1
	encryptor.context.gaussianSampler.Sample(encryptor.polypool[2])
	ringContext.Add(encryptor.polypool[1], encryptor.polypool[2], encryptor.polypool[1])

	// We rescal the encryption of zero by the special prime, dividing the error by this prime
	encryptor.baseconverter.ModDown(ringContext, encryptor.context.rescaleParamsKeys, uint64(len(plaintext.Value()[0].Coeffs))-1, encryptor.polypool[0], ciphertext.value[0], encryptor.polypool[2])
	encryptor.baseconverter.ModDown(ringContext, encryptor.context.rescaleParamsKeys, uint64(len(plaintext.Value()[0].Coeffs))-1, encryptor.polypool[1], ciphertext.value[1], encryptor.polypool[2])

	ringContext = encryptor.context.contextQ

	// ct[0] = pk[0]*u + e0 + m
	// ct[1] = pk[1]*u + e1
	ringContext.Add(ciphertext.value[0], plaintext.value, ciphertext.value[0])

}

func encryptfromsk(encryptor *Encryptor, plaintext *Plaintext, ciphertext *Ciphertext) {

	ringContext := encryptor.context.contextKeys

	// ct = [(-a*s + e)/P , a/P]
	ringContext.UniformPoly(encryptor.polypool[1])
	ringContext.MulCoeffsMontgomeryAndSub(encryptor.polypool[1], encryptor.sk.sk, encryptor.polypool[0])

	// We rescal the encryption of zero by the special prime, dividing the error by this prime
	ringContext.InvNTT(encryptor.polypool[0], encryptor.polypool[0])
	ringContext.InvNTT(encryptor.polypool[1], encryptor.polypool[1])

	encryptor.context.gaussianSampler.SampleAndAdd(encryptor.polypool[0])

	encryptor.baseconverter.ModDown(ringContext, encryptor.context.rescaleParamsKeys, uint64(len(plaintext.Value()[0].Coeffs))-1, encryptor.polypool[0], ciphertext.value[0], encryptor.polypool[2])
	encryptor.baseconverter.ModDown(ringContext, encryptor.context.rescaleParamsKeys, uint64(len(plaintext.Value()[0].Coeffs))-1, encryptor.polypool[1], ciphertext.value[1], encryptor.polypool[2])

	ringContext = encryptor.context.contextQ

	// ct = [-a*s + m + e , a]
	ringContext.Add(ciphertext.value[0], plaintext.value, ciphertext.value[0])

}
