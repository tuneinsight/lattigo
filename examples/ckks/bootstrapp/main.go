package main

import (
	"fmt"
	"math"
	"math/rand"

	"github.com/ldsec/lattigo/ckks"
)

func randomFloat(min, max float64) float64 {
	return min + rand.Float64()*(max-min)
}

func randomComplex(min, max float64) complex128 {
	return complex(randomFloat(min, max), randomFloat(min, max))
}

func main() {

	var err error

	var btp *ckks.Bootstrapper
	var kgen ckks.KeyGenerator
	var encoder ckks.Encoder
	var sk *ckks.SecretKey
	var pk *ckks.PublicKey
	var encryptor ckks.Encryptor
	var decryptor ckks.Decryptor
	var plaintext *ckks.Plaintext

	// Bootstrapping parameters
	// Five sets of parameters (index 0 to 4) ensuring 128 bit of security
	// are avaliable in github.com/ldsec/lattigo/ckks/bootparams
	// LogSlots is hardcoded to 15 in the parameters, but can be changed from 1 to 15.
	// When changing logSlots make sure that the number of levels allocated to CtS and StC is
	// smaller or equal to logSlots.
	params := ckks.DefaultBootstrappSchemeParams[4]
	btpParams := ckks.DefaultBootstrappParams[4]

	fmt.Println()
	fmt.Printf("CKKS parameters : logN = %d, logSlots = %d, h = %d, logQP = %d, levels = %d, scale= 2^%f, sigma = %f \n", params.LogN(), params.LogSlots(), btpParams.H, params.LogQP(), params.Levels(), math.Log2(params.Scale()), params.Sigma())

	// Scheme context and keys
	kgen = ckks.NewKeyGenerator(params)

	sk, pk = kgen.GenKeyPairSparse(btpParams.H)

	encoder = ckks.NewEncoder(params)
	decryptor = ckks.NewDecryptor(params, sk)
	encryptor = ckks.NewEncryptorFromPk(params, pk)

	// Bootstrapp context (generates public evaluation keys and relevant parameters)
	// TODO : seperate the key-generation from the rest, since those evaluation keys only need to be generated once.
	// Make it so they can be given as input to the bootstrapp.
	fmt.Println()
	fmt.Println("Generating bootstrappign keys...")
	btpKey := kgen.GenBootstrappingKey(params.LogSlots(), btpParams, sk)
	if btp, err = ckks.NewBootstrapper(params, btpParams, btpKey); err != nil {
		panic(err)
	}
	fmt.Println("Done")

	// Generates a random plaintext
	valuesWant := make([]complex128, params.Slots())
	for i := range valuesWant {
		valuesWant[i] = randomComplex(-1, 1)
	}

	plaintext = ckks.NewPlaintext(params, params.MaxLevel(), params.Scale())
	encoder.Encode(plaintext, valuesWant, params.Slots())

	// Encrypts
	ciphertext1 := encryptor.EncryptNew(plaintext)

	// Decrypts, prints and compares with the plaintext values
	fmt.Println()
	fmt.Println("Precision of values vs. ciphertext")
	valuesTest1 := printDebug(params, ciphertext1, valuesWant, decryptor, encoder)

	// Bootstrapps the ciphertext (homomorphic re-encryption)
	// Takes a ciphertext at level 0 (if not at level 0, then it will reduce it to level 0)
	// and returns a ciphertext at level MaxLevel - k, where k is the depth of the bootstrapping circuit.
	// CAUTION : the scale of the ciphertext MUST be equal to params.Scale (or very close)
	// To equalize the scale, the function evaluator.SetScale(ciphertext, parameters.Scale) can be used at the expanse of one level.
	fmt.Println()
	fmt.Println("Bootstrapping...")
	ciphertext2 := btp.Bootstrapp(ciphertext1)
	fmt.Println("Done")

	// Decrypts, prints and compares with the plaintext values
	fmt.Println()
	fmt.Println("Precision of ciphertext vs. Bootstrapp(ciphertext)")
	printDebug(params, ciphertext2, valuesTest1, decryptor, encoder)
}

func printDebug(params *ckks.Parameters, ciphertext *ckks.Ciphertext, valuesWant []complex128, decryptor ckks.Decryptor, encoder ckks.Encoder) (valuesTest []complex128) {

	slots := uint64(len(valuesWant))

	valuesTest = encoder.Decode(decryptor.DecryptNew(ciphertext), slots)

	fmt.Println()
	fmt.Printf("Level : %d (logQ = %d)\n", ciphertext.Level(), params.LogQLvl(ciphertext.Level()))
	fmt.Printf("Scale : 2^%f\n", math.Log2(ciphertext.Scale()))
	fmt.Printf("ValuesTest : %6.10f %6.10f %6.10f %6.10f...\n", valuesTest[0], valuesTest[1], valuesTest[2], valuesTest[3])
	fmt.Printf("ValuesWant : %6.10f %6.10f %6.10f %6.10f...\n", valuesWant[0], valuesWant[1], valuesWant[2], valuesWant[3])

	precStats := ckks.GetPrecisionStats(params, nil, nil, valuesWant, valuesTest)

	fmt.Println(precStats.String())
	fmt.Println()

	return
}
