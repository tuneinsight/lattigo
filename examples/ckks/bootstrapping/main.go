// Package main implements an example showcasing the basics of the bootstrapping for the CKKS scheme.
// The CKKS bootstrapping is a circuit that homomorphically re-encrypts a ciphertext at level zero to a ciphertext at a higher level, enabling further computations.
// Note that, unlike the BGV or BFV bootstrapping, the CKKS bootstrapping does not reduce the error in the ciphertext, but only enables further computations.
// Use the flag -short to run the examples fast but with insecure parameters.
package main

import (
	"flag"
	"fmt"
	"math"

	"github.com/tuneinsight/lattigo/v4/ckks"
	"github.com/tuneinsight/lattigo/v4/ckks/bootstrapping"
	"github.com/tuneinsight/lattigo/v4/ring/distribution"
	"github.com/tuneinsight/lattigo/v4/rlwe"
	"github.com/tuneinsight/lattigo/v4/utils/sampling"
)

var flagShort = flag.Bool("short", false, "run the example with a smaller and insecure ring degree.")

func main() {

	flag.Parse()

	// First we define the residual CKKS parameters. This is only a template that will be given
	// to the constructor along with the specificities of the bootstrapping circuit we choose, to
	// enable it to create the appropriate ckks.ParametersLiteral that enable the evaluation of the
	// bootstrapping circuit on top of the residual moduli that we defined.
	ckksParamsResidualLit := ckks.ParametersLiteral{
		LogN:     16,                                                // Log2 of the ringdegree
		LogQ:     []int{55, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40}, // Log2 of the ciphertext prime moduli
		LogP:     []int{61, 61, 61, 61},                             // Log2 of the key-switch auxiliary prime moduli
		LogScale: 40,                                                // Log2 of the scale
		Xs:       &distribution.Ternary{H: 192},                     // Hamming weight of the secret
	}

	LogSlots := ckksParamsResidualLit.LogN - 2

	if *flagShort {
		ckksParamsResidualLit.LogN -= 3
		LogSlots -= 3
	}

	// Note that with H=192 and LogN=16, parameters are at least 128-bit if LogQP <= 1550.
	// Our default parameters have an expected logQP of 55 + 10*40 + 4*61 = 699, meaning
	// that the depth of the bootstrapping shouldn't be larger than 1550-699 = 851.

	// For this first example, we do not specify any optional field of the bootstrapping
	// Thus we expect the bootstrapping to give a precision of 27.25 bits with H=192 (and 23.8 with H=N/2)
	// if the plaintext values are uniformly distributed in [-1, 1] for both the real and imaginary part.
	// See `/ckks/bootstrapping/parameters.go` for information about the optional fields.
	btpParametersLit := bootstrapping.ParametersLiteral{
		// Since a ciphertext with message m and LogSlots = x is equivalent to a ciphertext with message m|m and LogSlots = x+1
		// it is possible to run the bootstrapping on any ciphertext with LogSlots <= bootstrapping.LogSlots, however doing so
		// will increase the runtime, so it is recommanded to have the LogSlots of the ciphertext and bootstrapping parameters
		// be the same.
		LogSlots: &LogSlots,
	}

	// The default bootstrapping parameters consume 822 bits which is smaller than the maximum
	// allowed of 851 in our example, so the target security is easily met.
	// We can print and verify the expected bit consumption of bootstrapping parameters with:
	bits, err := btpParametersLit.BitConsumption(LogSlots)
	if err != nil {
		panic(err)
	}
	fmt.Printf("Bootstrapping depth (bits): %d\n", bits)

	// Now we generate the updated ckks.ParametersLiteral that contain our residual moduli and the moduli for
	// the bootstrapping circuit, as well as the bootstrapping.Parameters that contain all the necessary information
	// of the bootstrapping circuit.
	ckksParamsLit, btpParams, err := bootstrapping.NewParametersFromLiteral(ckksParamsResidualLit, btpParametersLit)
	if err != nil {
		panic(err)
	}

	if *flagShort {
		// Corrects the message ratio to take into account the smaller number of slots and keep the same precision
		btpParams.EvalModParameters.LogMessageRatio += 3
	}

	// This generate ckks.Parameters, with the NTT tables and other pre-computations from the ckks.ParametersLiteral (which is only a template).
	params, err := ckks.NewParametersFromLiteral(ckksParamsLit)
	if err != nil {
		panic(err)
	}

	// Here we print some information about the generated ckks.Parameters
	// We can notably check that the LogQP of the generated ckks.Parameters is equal to 699 + 822 = 1521.
	// Not that this value can be overestimated by one bit.
	fmt.Printf("CKKS parameters: logN=%d, logSlots=%d, H(%d; %d), sigma=%f, logQP=%f, levels=%d, scale=2^%f\n", params.LogN(), LogSlots, params.XsHammingWeight(), btpParams.EphemeralSecretWeight, params.Xe(), params.LogQP(), params.QCount(), math.Log2(params.DefaultScale().Float64()))

	// Scheme context and keys
	kgen := ckks.NewKeyGenerator(params)

	sk, pk := kgen.GenKeyPairNew()

	encoder := ckks.NewEncoder(params)
	decryptor := ckks.NewDecryptor(params, sk)
	encryptor := ckks.NewEncryptor(params, pk)

	fmt.Println()
	fmt.Println("Generating bootstrapping keys...")
	evk := bootstrapping.GenEvaluationKeySetNew(btpParams, params, sk)
	fmt.Println("Done")

	var btp *bootstrapping.Bootstrapper
	if btp, err = bootstrapping.NewBootstrapper(params, btpParams, evk); err != nil {
		panic(err)
	}

	// Generate a random plaintext with values uniformely distributed in [-1, 1] for the real and imaginary part.
	valuesWant := make([]complex128, 1<<LogSlots)
	for i := range valuesWant {
		valuesWant[i] = sampling.RandComplex128(-1, 1)
	}

	plaintext := ckks.NewPlaintext(params, params.MaxLevel())
	plaintext.LogSlots = [2]int{0, LogSlots}
	if err := encoder.Encode(valuesWant, plaintext); err != nil {
		panic(err)
	}

	// Encrypt
	ciphertext1 := encryptor.EncryptNew(plaintext)

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
	fmt.Println(ciphertext1.LogSlots)
	fmt.Println()
	fmt.Println("Bootstrapping...")
	ciphertext2 := btp.Bootstrap(ciphertext1)
	fmt.Println("Done")

	// Decrypt, print and compare with the plaintext values
	fmt.Println()
	fmt.Println("Precision of ciphertext vs. Bootstrap(ciphertext)")
	printDebug(params, ciphertext2, valuesTest1, decryptor, encoder)
}

func printDebug(params ckks.Parameters, ciphertext *rlwe.Ciphertext, valuesWant []complex128, decryptor rlwe.Decryptor, encoder *ckks.Encoder) (valuesTest []complex128) {

	valuesTest = make([]complex128, 1<<ciphertext.LogSlots[1])

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
