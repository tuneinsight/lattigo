package ckks

import (
	"github.com/ldsec/lattigo/ring"
)

type encryptorContext struct {
	// Context parameters
	n uint64

	// Number of available levels
	levels uint64

	// Contexts
	specialPrimes []uint64
	contextQ      *ring.Context
	contextP      *ring.Context
	contextKeys   *ring.Context

	// Pre-computed values for the rescaling
	rescaleParamsKeys []uint64 // (P^-1) mod each qi

	// Samplers
	gaussianSampler *ring.KYSampler
}

func newEncryptorContext(params *Parameters) *encryptorContext {
	n := uint64(1 << uint64(params.LogN))
	levels := uint64(len(params.Modulichain))

	scalechain := make([]float64, len(params.Modulichain))

	// Extracts all the different primes bit size and maps their number
	primesbitlen := make(map[uint64]uint64)
	for i, qi := range params.Modulichain {

		primesbitlen[uint64(qi)]++

		if uint64(params.Modulichain[i]) > 60 {
			panic("provided moduli must be smaller than 61")
		}
	}

	for _, pj := range params.P {
		primesbitlen[uint64(pj)]++

		if uint64(pj) > 60 {
			panic("provided P must be smaller than 61")
		}
	}

	// For each bitsize, finds that many primes
	primes := make(map[uint64][]uint64)
	for key, value := range primesbitlen {
		primes[key] = GenerateCKKSPrimes(key, uint64(params.LogN), value)
	}

	// Assigns the primes to the ckks moduli chain
	moduli := make([]uint64, len(params.Modulichain))
	for i, qi := range params.Modulichain {
		moduli[i] = primes[uint64(params.Modulichain[i])][0]
		primes[uint64(qi)] = primes[uint64(qi)][1:]

		scalechain[i] = float64(moduli[i])
	}

	// Assigns the primes to the special primes list for the the keyscontext
	specialPrimes := make([]uint64, len(params.P))
	for i, pj := range params.P {
		specialPrimes[i] = primes[uint64(pj)][0]
		primes[uint64(pj)] = primes[uint64(pj)][1:]
	}

	// Contexts
	contextQ := ring.NewContext()
	contextQ.SetParameters(1<<params.LogN, moduli)

	err := contextQ.GenNTTParams()
	if err != nil {
		panic(err)
	}

	contextP := ring.NewContext()
	contextP.SetParameters(1<<params.LogN, specialPrimes)

	err = contextP.GenNTTParams()
	if err != nil {
		panic(err)
	}

	contextKeys := ring.NewContext()
	contextKeys.SetParameters(1<<params.LogN, append(moduli, specialPrimes...))

	err = contextKeys.GenNTTParams()
	if err != nil {
		panic(err)
	}

	var Qi uint64

	bredParams := contextQ.GetBredParams()

	rescaleParamsKeys := make([]uint64, levels)

	PBig := ring.NewUint(1)
	for _, pj := range specialPrimes {
		PBig.Mul(PBig, ring.NewUint(pj))
	}

	tmp := ring.NewUint(0)

	for i := uint64(0); i < levels; i++ {

		Qi = moduli[i]

		tmp.Mod(PBig, ring.NewUint(Qi))

		rescaleParamsKeys[i] = ring.MForm(ring.ModExp(ring.BRedAdd(tmp.Uint64(), Qi, bredParams[i]), Qi-2, Qi), Qi, bredParams[i])
	}

	gaussianSampler := contextKeys.NewKYSampler(params.Sigma, int(6*params.Sigma))

	return &encryptorContext{
		n:                 n,
		levels:            levels,
		specialPrimes:     specialPrimes,
		contextQ:          contextQ,
		contextP:          contextP,
		contextKeys:       contextKeys,
		rescaleParamsKeys: rescaleParamsKeys,
		gaussianSampler:   gaussianSampler,
	}
}

// Encryptor is a struct used to encrypt plaintext and storing the public-key and/or secret-key.
type Encryptor struct {
	context  *encryptorContext
	pk       *PublicKey
	sk       *SecretKey
	polypool [3]*ring.Poly

	rescalepool []uint64

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

// NewEncryptor creates a new Encryptor with the input public-key and/or secret-key.
// This encryptor can be used to encrypt plaintexts, using the stored keys.
func newEncryptor(pk *PublicKey, sk *SecretKey, params *Parameters) (encryptor *Encryptor) {
	context := newEncryptorContext(params)

	if pk != nil && (uint64(pk.pk[0].GetDegree()) != context.n || uint64(pk.pk[1].GetDegree()) != context.n) {
		panic("pk ring degree doesn't match ckkscontext ring degree")
	}

	if sk != nil && uint64(sk.sk.GetDegree()) != context.n {
		panic("sk ring degree doesn't match ckkscontext ring degree")
	}

	encryptor = new(Encryptor)
	encryptor.context = context
	encryptor.pk = pk
	encryptor.sk = sk

	encryptor.polypool[0] = context.contextKeys.NewPoly()
	encryptor.polypool[1] = context.contextKeys.NewPoly()
	encryptor.polypool[2] = context.contextKeys.NewPoly()

	encryptor.rescalepool = make([]uint64, context.n)

	encryptor.baseconverter = ring.NewFastBasisExtender(context.contextQ.Modulus, context.specialPrimes)

	return encryptor
}

// EncryptNew encrypts the input plaintext using the stored key and returns
// the result on a newly created ciphertext.
//
// encrypt with pk : ciphertext = [pk[0]*u + m + e_0, pk[1]*u + e_1]
// encrypt with sk : ciphertext = [-a*sk + m + e, a]
func (encryptor *Encryptor) EncryptNew(plaintext *Plaintext) (ciphertext *Ciphertext) {

	ciphertext = NewCiphertext(1, plaintext.Level(), plaintext.Scale(), encryptor.context.contextQ)
	encryptor.Encrypt(plaintext, ciphertext)
	return
}

// Encrypt encrypts the input plaintext using the stored key, and returns the result
// on the reciver ciphertext.
//
// encrypt with pk : ciphertext = [pk[0]*u + m + e_0, pk[1]*u + e_1]
// encrypt with sk : ciphertext = [-a*sk + m + e, a]
func (encryptor *Encryptor) Encrypt(plaintext *Plaintext, ciphertext *Ciphertext) {

	if plaintext.Level() != encryptor.context.levels-1 {
		panic("cannot encrypt -> plaintext not at maximum level")
	}

	if encryptor.sk != nil {

		encryptfromsk(encryptor, plaintext, ciphertext)

	} else if encryptor.pk != nil {

		encryptfrompk(encryptor, plaintext, ciphertext)

	} else {

		panic("cannot encrypt -> public-key and/or secret-key has not been set")
	}
}

func encryptfrompk(encryptor *Encryptor, plaintext *Plaintext, ciphertext *Ciphertext) {

	// We sample a R-WLE instance (encryption of zero) over the keys context (ciphertext context + special prime)

	contextKeys := encryptor.context.contextKeys
	contextQ := encryptor.context.contextQ

	encryptor.context.contextKeys.SampleTernaryMontgomeryNTT(encryptor.polypool[2], 0.5)

	// ct0 = u*pk0
	contextKeys.MulCoeffsMontgomery(encryptor.polypool[2], encryptor.pk.pk[0], encryptor.polypool[0])
	// ct1 = u*pk1
	contextKeys.MulCoeffsMontgomery(encryptor.polypool[2], encryptor.pk.pk[1], encryptor.polypool[1])

	// 2*(#Q + #P) NTT
	contextKeys.InvNTT(encryptor.polypool[0], encryptor.polypool[0])
	contextKeys.InvNTT(encryptor.polypool[1], encryptor.polypool[1])

	// ct0 = u*pk0 + e0
	encryptor.context.gaussianSampler.SampleAndAdd(encryptor.polypool[0])
	// ct1 = u*pk1 + e1
	encryptor.context.gaussianSampler.SampleAndAdd(encryptor.polypool[1])

	// ct0 = (u*pk0 + e0)/P
	encryptor.baseconverter.ModDown(contextKeys, encryptor.context.rescaleParamsKeys, plaintext.Level(), encryptor.polypool[0], ciphertext.value[0], encryptor.polypool[2])

	// ct1 = (u*pk1 + e1)/P
	encryptor.baseconverter.ModDown(contextKeys, encryptor.context.rescaleParamsKeys, plaintext.Level(), encryptor.polypool[1], ciphertext.value[1], encryptor.polypool[2])

	// 2*#Q NTT
	contextQ.NTT(ciphertext.value[0], ciphertext.value[0])
	contextQ.NTT(ciphertext.value[1], ciphertext.value[1])

	// ct0 = (u*pk0 + e0)/P + m
	contextQ.Add(ciphertext.value[0], plaintext.value, ciphertext.value[0])

	ciphertext.isNTT = true
}

func encryptfromsk(encryptor *Encryptor, plaintext *Plaintext, ciphertext *Ciphertext) {
	contextKeys := encryptor.context.contextKeys
	contextP := encryptor.context.contextP
	contextQ := encryptor.context.contextQ

	// ct1 = a
	contextKeys.UniformPoly(encryptor.polypool[1])

	// ct0 = -s*a
	contextKeys.MulCoeffsMontgomery(encryptor.polypool[1], encryptor.sk.sk, encryptor.polypool[0])
	contextKeys.Neg(encryptor.polypool[0], encryptor.polypool[0])

	// #Q + #P NTT
	contextKeys.InvNTT(encryptor.polypool[0], encryptor.polypool[0])

	// ct0 = -s*a + e
	encryptor.context.gaussianSampler.SampleAndAdd(encryptor.polypool[0])

	// We rescal by the special prime, dividing the error by this prime
	// ct0 = (-s*a + e)/P
	encryptor.baseconverter.ModDown(contextKeys, encryptor.context.rescaleParamsKeys, plaintext.Level(), encryptor.polypool[0], ciphertext.value[0], encryptor.polypool[2])

	// #Q + #P NTT
	// ct1 = a/P
	encryptor.baseconverter.ModDownNTT(contextQ, contextP, encryptor.context.rescaleParamsKeys, plaintext.Level(), encryptor.polypool[1], ciphertext.value[1], encryptor.polypool[2])

	// #Q NTT
	contextQ.NTT(ciphertext.value[0], ciphertext.value[0])

	// ct0 = -s*a + m + e
	contextQ.Add(ciphertext.value[0], plaintext.value, ciphertext.value[0])

	ciphertext.isNTT = true
}
