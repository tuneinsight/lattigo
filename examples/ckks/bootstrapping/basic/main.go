// Package main implements an example showcasing the basics of the bootstrapping for the CKKS scheme.
// The CKKS bootstrapping is a circuit that homomorphically re-encrypts a ciphertext at level zero to a ciphertext at a higher level, enabling further computations.
// Note that, unlike the BGV or BFV bootstrapping, the CKKS bootstrapping does not reduce the error in the ciphertext, but only enables further computations.
// Use the flag -short to run the examples fast but with insecure parameters.
package main

import (
	"flag"
	"fmt"
	"math"

	"github.com/tuneinsight/lattigo/v4/circuits/float/bootstrapper"
	"github.com/tuneinsight/lattigo/v4/ckks"
	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/rlwe"
	"github.com/tuneinsight/lattigo/v4/utils"
	"github.com/tuneinsight/lattigo/v4/utils/sampling"
)

var flagShort = flag.Bool("short", false, "run the example with a smaller and insecure ring degree.")

func main() {

	flag.Parse()

	LogN := 16

	if *flagShort {
		LogN -= 3
	}

	// First we define the residual CKKS parameters.
	// For this example, we have a logQ = 55 + 10*40 and logP = 3*61
	// These are the parameters that the regular circuit will use outside of the
	// circuit bootstrapping.
	// The bootstrapping circuit use its own ckks.Parameters which are automatically
	// parameterized given the residual parameters and the bootsrappping parameters.
	params, err := ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
		LogN:            LogN,                                              // Log2 of the ringdegree
		LogQ:            []int{55, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40}, // Log2 of the ciphertext prime moduli
		LogP:            []int{61, 61, 61},                                 // Log2 of the key-switch auxiliary prime moduli
		LogDefaultScale: 40,                                                // Log2 of the scale
		Xs:              ring.Ternary{H: 192},                              // Hamming weight of the secret
	})

	if err != nil {
		panic(err)
	}

	// Note that with H=192 and LogN=16 the bootstrapping parameters are at least 128-bit if their LogQP <= 1550.

	// For this first example, we do not specify any optional field of the bootstrapping parameters.
	// Thus we expect the bootstrapping to give a precision of 27.25 bits with H=192 (and 23.8 with H=N/2)
	// if the plaintext values are uniformly distributed in [-1, 1] for both the real and imaginary part.
	// See `/ckks/bootstrapping/parameters.go` for information about the optional fields.
	btpParametersLit := bootstrapper.ParametersLiteral{
		// We specify LogN to ensure that both the residual parameters and the bootstrapping parameters
		// have the same LogN
		LogN: utils.Pointy(params.LogN()),

		// We manually specify the number of auxiliary primes used by the evaluation keys of the bootstrapping
		// circuit, so that the security target of LogQP is met.
		NumberOfPi: utils.Pointy(4),
	}

	// Now we generate the updated ckks.ParametersLiteral that contain our residual moduli and the moduli for
	// the bootstrapping circuit, as well as the bootstrapping.Parameters that contain all the necessary information
	// of the bootstrapping circuit.
	btpParams, err := bootstrapper.NewParametersFromLiteral(params, btpParametersLit)
	if err != nil {
		panic(err)
	}

	if *flagShort {
		// Corrects the message ratio to take into account the smaller number of slots and keep the same precision
		btpParams.Mod1ParametersLiteral.LogMessageRatio += 3
	}

	// Here we print some information about the residual parameters and the bootstrapping parameters
	// We can notably check that the LogQP of the bootstrapping parameters is smaller than 1550, which ensures
	// 128-bit of security as explained above.
	fmt.Printf("Residual parameters: logN=%d, logSlots=%d, H=%d, sigma=%f, logQP=%f, levels=%d, scale=2^%d\n", params.LogN(), params.LogMaxSlots(), params.XsHammingWeight(), params.Xe(), params.LogQP(), params.MaxLevel(), params.LogDefaultScale())
	fmt.Printf("Bootstrapping parameters: logN=%d, logSlots=%d, H(%d; %d), sigma=%f, logQP=%f, levels=%d, scale=2^%d\n", btpParams.LogN(), btpParams.LogMaxSlots(), btpParams.XsHammingWeight(), btpParams.EphemeralSecretWeight, btpParams.Xe(), btpParams.LogQP(), btpParams.QCount(), btpParams.LogDefaultScale())

	// Scheme context and keys
	kgen := ckks.NewKeyGenerator(params)

	sk, pk := kgen.GenKeyPairNew()

	encoder := ckks.NewEncoder(params)
	decryptor := ckks.NewDecryptor(params, sk)
	encryptor := ckks.NewEncryptor(params, pk)

	fmt.Println()
	fmt.Println("Generating bootstrapping keys...")
	evk, err := btpParams.GenBootstrappingKeys(sk)
	if err != nil {
		panic(err)
	}
	fmt.Println("Done")

	var btp *bootstrapper.Bootstrapper
	if btp, err = bootstrapper.NewBootstrapper(params, btpParams, evk); err != nil {
		panic(err)
	}

	// Generate a random plaintext with values uniformely distributed in [-1, 1] for the real and imaginary part.
	valuesWant := make([]complex128, params.MaxSlots())
	for i := range valuesWant {
		valuesWant[i] = sampling.RandComplex128(-1, 1)
	}

	plaintext := ckks.NewPlaintext(params, params.MaxLevel())
	if err := encoder.Encode(valuesWant, plaintext); err != nil {
		panic(err)
	}

	// Encrypt
	ciphertext1, err := encryptor.EncryptNew(plaintext)
	if err != nil {
		panic(err)
	}

	// Decrypt, print and compare with the plaintext values
	fmt.Println()
	fmt.Println("Precision of values vs. ciphertext")
	valuesTest1 := printDebug(params, ciphertext1, valuesWant, decryptor, encoder)

	// Bootstrap the ciphertext (homomorphic re-encryption)
	// It takes a ciphertext at level 0 (if not at level 0, then it will reduce it to level 0)
	// and returns a ciphertext with the max level of `ckksParamsResidualLit`.
	// CAUTION: the scale of the ciphertext MUST be equal (or very close) to params.DefaultScale()
	// To equalize the scale, the function evaluator.SetScale(ciphertext, parameters.DefaultScale()) can be used at the expense of one level.
	// If the ciphertext is is at level one or greater when given to the bootstrapper, this equalization is automatically done.
	fmt.Println("Bootstrapping...")
	ciphertext2, err := btp.Bootstrap(ciphertext1)
	if err != nil {
		panic(err)
	}
	fmt.Println("Done")

	// Decrypt, print and compare with the plaintext values
	fmt.Println()
	fmt.Println("Precision of ciphertext vs. Bootstrap(ciphertext)")
	printDebug(params, ciphertext2, valuesTest1, decryptor, encoder)
}

func printDebug(params ckks.Parameters, ciphertext *rlwe.Ciphertext, valuesWant []complex128, decryptor *rlwe.Decryptor, encoder *ckks.Encoder) (valuesTest []complex128) {

	valuesTest = make([]complex128, ciphertext.Slots())

	if err := encoder.Decode(decryptor.DecryptNew(ciphertext), valuesTest); err != nil {
		panic(err)
	}

	fmt.Println()
	fmt.Printf("Level: %d (logQ = %d)\n", ciphertext.Level(), params.LogQLvl(ciphertext.Level()))

	fmt.Printf("Scale: 2^%f\n", math.Log2(ciphertext.Scale.Float64()))
	fmt.Printf("ValuesTest: %6.10f %6.10f %6.10f %6.10f...\n", valuesTest[0], valuesTest[1], valuesTest[2], valuesTest[3])
	fmt.Printf("ValuesWant: %6.10f %6.10f %6.10f %6.10f...\n", valuesWant[0], valuesWant[1], valuesWant[2], valuesWant[3])

	precStats := ckks.GetPrecisionStats(params, encoder, nil, valuesWant, valuesTest, nil, false)

	fmt.Println(precStats.String())
	fmt.Println()

	return
}
