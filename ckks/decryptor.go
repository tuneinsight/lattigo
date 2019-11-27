package ckks

import (
	"github.com/ldsec/lattigo/ring"
)

type decryptorContext struct {
	// Context parameters
	n uint64

	// Contexts
	contextQ *ring.Context
}

func newDecryptorContext(params *Parameters) *decryptorContext {
	n := uint64(1 << uint64(params.LogN))

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

	return &decryptorContext{
		n:        n,
		contextQ: contextQ,
	}
}

// Decryptor is a structure used to decrypt ciphertext. It stores the secret-key.
type Decryptor struct {
	context *decryptorContext
	sk      *SecretKey
}

// NewDecryptor instanciates a new decryptor that will be able to decrypt ciphertext
// encrypted under the provided secret-key.
func NewDecryptor(sk *SecretKey, params *Parameters) *Decryptor {
	context := newDecryptorContext(params)

	if sk.sk.GetDegree() != int(context.n) {
		panic("secret_key degree must match context degree")
	}

	decryptor := new(Decryptor)

	decryptor.context = context

	decryptor.sk = sk

	return decryptor
}

// DecryptNew decrypts the ciphertext and returns a newly created plaintext.
// A Horner methode is used for evaluating the decryption.
func (decryptor *Decryptor) DecryptNew(ciphertext *Ciphertext) (plaintext *Plaintext) {

	plaintext = NewPlaintext(ciphertext.Level(), ciphertext.Scale(), decryptor.context.contextQ)

	decryptor.Decrypt(ciphertext, plaintext)

	return
}

// Decrypt decrypts the ciphertext and returns the result on the provided receiver plaintext.
// A Horner methode is used for evaluating the decryption.
func (decryptor *Decryptor) Decrypt(ciphertext *Ciphertext, plaintext *Plaintext) {

	context := decryptor.context.contextQ

	level := ciphertext.Level()

	plaintext.SetScale(ciphertext.Scale())

	context.CopyLvl(level, ciphertext.value[ciphertext.Degree()], plaintext.value)

	plaintext.value.Coeffs = plaintext.value.Coeffs[:ciphertext.Level()+1]

	for i := uint64(ciphertext.Degree()); i > 0; i-- {

		context.MulCoeffsMontgomeryLvl(level, plaintext.value, decryptor.sk.sk, plaintext.value)
		context.AddLvl(level, plaintext.value, ciphertext.value[i-1], plaintext.value)

		if i&7 == 7 {
			context.ReduceLvl(level, plaintext.value, plaintext.value)
		}
	}

	if (ciphertext.Degree())&7 != 7 {
		context.ReduceLvl(level, plaintext.value, plaintext.value)
	}
}
