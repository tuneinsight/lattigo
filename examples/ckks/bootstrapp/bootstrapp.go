package main

import (
	"fmt"
	"github.com/ldsec/lattigo/ckks"
	"math"
	"math/rand"
	"sort"
)

func randomFloat(min, max float64) float64 {
	return min + rand.Float64()*(max-min)
}

func randomComplex(min, max float64) complex128 {
	return complex(randomFloat(min, max), randomFloat(min, max))
}

func main() {

	var bootcontext *ckks.BootContext
	var kgen ckks.KeyGenerator
	var encoder ckks.Encoder
	var sk *ckks.SecretKey
	var pk *ckks.PublicKey
	var encryptor ckks.Encryptor
	var decryptor ckks.Decryptor
	var ciphertext *ckks.Ciphertext
	var plaintext *ckks.Plaintext

	// Bootstrapping parameters
	// Five sets of parameters (index 0 to 4) ensuring 128 bit of security
	// are avaliable in github.com/ldsec/lattigo/ckks/bootparams
	// LogSlots is hardcoded to 10 in the parameters, but can be changed from 1 to 15.
	bootparams := ckks.BootstrappParams[4]
	bootparams.Gen()
	params := &bootparams.Parameters

	fmt.Println()
	fmt.Printf("CKKS parameters : logN = %d, logSlots = %d, logQP = %d, levels = %d, scale= %f, sigma = %f \n", params.LogN(), params.LogSlots(), params.LogQP(), params.Levels(), params.Scale(), params.Sigma())

	// Scheme context and keys
	kgen = ckks.NewKeyGenerator(params)

	sk, pk = kgen.GenKeyPairSparse(128)

	encoder = ckks.NewEncoder(params)
	decryptor = ckks.NewDecryptor(params, sk)
	encryptor = ckks.NewEncryptorFromPk(params, pk)

	// Bootstrapp context (generates public evaluation keys and relevant parameters)
	// TODO : seperate the key-generation from the rest, since those evaluation keys only need to be generated once.
	// Make it so they can be given as input to the bootstrapp.
	fmt.Println()
	fmt.Println("Generating bootstrappign keys...")
	bootcontext = ckks.NewBootContext(bootparams)
	bootcontext.GenBootKeys(sk)
	fmt.Println("Done")

	// Generates a random plaintext
	valuesWant := make([]complex128, params.Slots())
	for i := range valuesWant {
		valuesWant[i] = randomComplex(-1, 1)
	}

	plaintext = ckks.NewPlaintext(params, params.MaxLevel(), params.Scale())
	encoder.Encode(plaintext, valuesWant, params.Slots())

	// Encrypts
	ciphertext = encryptor.EncryptNew(plaintext)

	// Decrypts, prints and compares with the plaintext values
	printDebug(ciphertext, valuesWant, decryptor, encoder)

	// Bootstrapps the ciphertext (homomorphic re-encryption)
	// Takes a ciphertext at level 0 (if not at level 0, then it will reduce it to level 0)
	// and returns a ciphertext at level MaxLevel - k, where k is the depth of the bootstrapping circuit.
	// CAUTION : the scale of the ciphertext MUST be equal to params.Scale (or very close)
	// To equalize the scale, the function evaluator.SetScale(ciphertext, parameters.Scale) can be used at the expanse of one level.
	fmt.Println()
	fmt.Println("Bootstrapping...")
	ciphertext = bootcontext.Bootstrapp(ciphertext)
	fmt.Println("Done")

	// Decrypts, prints and compares with the plaintext values
	printDebug(ciphertext, valuesWant, decryptor, encoder)
}

func printDebug(ciphertext *ckks.Ciphertext, valuesWant []complex128, decryptor ckks.Decryptor, encoder ckks.Encoder) {

	slots := uint64(len(valuesWant))

	valuesTest := encoder.Decode(decryptor.DecryptNew(ciphertext), slots)

	fmt.Println()
	fmt.Printf("ValuesTest : %6.10f %6.10f %6.10f %6.10f...\n", valuesTest[0], valuesTest[1], valuesTest[2], valuesTest[3])
	fmt.Printf("ValuesWant : %6.10f %6.10f %6.10f %6.10f...\n", valuesWant[0], valuesWant[1], valuesWant[2], valuesWant[3])
	verifyVector(valuesWant, valuesTest)
}

func verifyVector(valuesWant, valuesTest []complex128) (err error) {

	var deltaReal, deltaImag float64
	var delta, minprec, maxprec, meanprec, medianprec complex128

	diff := make([]complex128, len(valuesWant))

	minprec = complex(0, 0)
	maxprec = complex(1, 1)

	meanprec = complex(0, 0)

	for i := range valuesWant {

		delta = valuesTest[i] - valuesWant[i]
		deltaReal = math.Abs(real(delta))
		deltaImag = math.Abs(imag(delta))

		diff[i] += complex(deltaReal, deltaImag)

		meanprec += diff[i]

		if deltaReal > real(minprec) {
			minprec = complex(deltaReal, imag(minprec))
		}

		if deltaImag > imag(minprec) {
			minprec = complex(real(minprec), deltaImag)
		}

		if deltaReal < real(maxprec) {
			maxprec = complex(deltaReal, imag(maxprec))
		}

		if deltaImag < imag(maxprec) {
			maxprec = complex(real(maxprec), deltaImag)
		}
	}

	meanprec /= complex(float64(len(valuesWant)), 0)
	medianprec = calcmedian(diff)

	fmt.Println()
	fmt.Printf("Minimum precision : (%.2f, %.2f) bits \n", math.Log2(1/real(minprec)), math.Log2(1/imag(minprec)))
	fmt.Printf("Maximum precision : (%.2f, %.2f) bits \n", math.Log2(1/real(maxprec)), math.Log2(1/imag(maxprec)))
	fmt.Printf("Mean    precision : (%.2f, %.2f) bits \n", math.Log2(1/real(meanprec)), math.Log2(1/imag(meanprec)))
	fmt.Printf("Median  precision : (%.2f, %.2f) bits \n", math.Log2(1/real(medianprec)), math.Log2(1/imag(medianprec)))
	fmt.Println()

	return nil
}

func calcmedian(values []complex128) (median complex128) {

	tmp := make([]float64, len(values))

	for i := range values {
		tmp[i] = real(values[i])
	}

	sort.Float64s(tmp)

	for i := range values {
		values[i] = complex(tmp[i], imag(values[i]))
	}

	for i := range values {
		tmp[i] = imag(values[i])
	}

	sort.Float64s(tmp)

	for i := range values {
		values[i] = complex(real(values[i]), tmp[i])
	}

	index := len(values) / 2

	if len(values)&1 == 1 {
		return values[index]
	}

	return (values[index] + values[index+1]) / 2

}
