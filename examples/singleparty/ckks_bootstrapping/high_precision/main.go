// Package main implements an example showcasing high-precision bootstrapping for high-precision fixed-
// point approximate arithmetic over the reals/complexes.
// High-precision bootstrapping is achieved by bootstrapping the residual error iteratively until it
// becomes smaller than the initial ciphertext error. Assume that the bootstrapping circuit has `d` bits
// of precision and that a ciphertext encrypts a message `m` with `kd` bits of scaling factor. Then the procedure
// works as follow:
//
// 1) Input: ctIn[2^{kd} * m]
// 2) Bootstrap(ctIn[2^{kd} * m] / 2^{(k-1)d})*2^{(k-1)d} -> ctOut[2^{kd} * m + 2^{(k-1)d} * e0]
// 3) Bootstrap((ctOut - ctIn)/2^{(k-2)d}) * 2^{(k-2)d} -> ctIn = [2^{(k-1)d} * e0 + 2^{(k-2)d} * e1]
// 4) ctOut = ctOut - ctIn = [2^{kd} * m + 2^{(k-2)d} * e1]
// 4) We can repeat this process k-2 additional times to get a bootstrapping of k * d bits of precision.
//
// The method is described in details by Bae et al. in META-BTS: Bootstrapping Precision Beyond the Limit (https://eprint.iacr.org/2022/1167).
//
// This example assumes that the user is already familiar with the bootstrapping and its different steps.
// See the basic example `lattigo/single_party/applications/reals_bootstrapping/basics` for an introduction into the
// bootstrapping.
// Use the flag -short to run the examples fast but with insecure parameters.
package main

import (
	"flag"
	"fmt"
	"math"

	"github.com/tuneinsight/lattigo/v6/circuits/ckks/bootstrapping"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
	"github.com/tuneinsight/lattigo/v6/utils"
	"github.com/tuneinsight/lattigo/v6/utils/bignum"
	"github.com/tuneinsight/lattigo/v6/utils/sampling"
)

var flagShort = flag.Bool("short", false, "run the example with a smaller and insecure ring degree.")

func main() {

	flag.Parse()

	// Default LogN, which with the following defined parameters
	// provides a security of 128-bit.
	LogN := 16

	if *flagShort {
		LogN -= 3
	}

	//==============================
	//=== 1) RESIDUAL PARAMETERS ===
	//==============================

	// First we must define the residual parameters.
	// The residual parameters are the parameters used outside of the bootstrapping circuit.
	// For this example, we have a LogN=16, logQ = (55+45) + 5*(45+45) and logP = 3*61, so LogQP = 638.
	// With LogN=16, LogQP=638 and H=192, these parameters achieve well over 128-bit of security.
	params, err := ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
		LogN:            LogN,              // Log2 of the ring degree
		LogQ:            []int{60, 45},     // Log2 of the ciphertext prime moduli
		LogP:            []int{61, 61, 61}, // Log2 of the key-switch auxiliary prime moduli
		LogDefaultScale: 90,                // Log2 of the scale
		Xs:              ring.Ternary{H: 192},
	})

	if err != nil {
		panic(err)
	}

	prec := params.EncodingPrecision()

	//==========================================
	//=== 2) BOOTSTRAPPING PARAMETERSLITERAL ===
	//==========================================

	// The bootstrapping circuit use its own Parameters which will be automatically
	// instantiated given the residual parameters and the bootstrapping parameters.

	// !WARNING! The bootstrapping parameters are not ensure to be 128-bit secure, it is the
	// responsibility of the user to check that the meet the security requirement and tweak them if necessary.

	// Note that the default bootstrapping parameters use LogN=16 and a ternary secret with H=192 non-zero coefficients
	// which provides parameters which are at least 128-bit if their LogQP <= 1550.

	// For this first example, we do not specify any circuit specific optional field in the bootstrapping parameters literal.
	// Thus we expect the bootstrapping to give an average precision of 27.9 bits with H=192 (and 24.4 with H=N/2)
	// if the plaintext values are uniformly distributed in [-1, 1] for both the real and imaginary part.
	// See `circuits/bootstrapping/parameters_literal.go` for detailed information about the optional fields.
	btpParametersLit := bootstrapping.ParametersLiteral{
		// We specify LogN to ensure that both the residual parameters and the bootstrapping parameters
		// have the same LogN. This is not required, but we want it for this example.
		LogN: utils.Pointy(LogN),

		// In this example we need manually specify the number of auxiliary primes (i.e. #Pi) used by the
		// evaluation keys of the bootstrapping circuit, so that the size of LogQP  meets the security target.
		LogP: []int{61, 61, 61, 61},

		// Sets the IterationsParameters.
		// The default bootstrapping parameters have 27.9 bits of average precision and
		// ~25 bits of minimum precision, and the maximum precision that can be theoretically
		// achieved is LogScale - LogN/2.
		// Therefore we start with 27.9 bits and each can in theory increase the precision an additional 25 bits.
		// However, to achieve the best possible precision, we must carefully adjust each iteration by hand so
		// that the sum of all the minimum precision is as close as possible
		// to LogScale - LogN/2. Here 27.9+25+25+5 ~= 82.5 (for the insecure parameters with LogN=13, with
		// the secure parameters using LogN=16 achieve 82.5 - (16-13)/2 = 81 bits of precision).
		IterationsParameters: &bootstrapping.IterationsParameters{
			BootstrappingPrecision: []float64{25, 25, 5},
			ReservedPrimeBitSize:   28,
		},

		// In this example we manually specify the bootstrapping parameters' secret distribution.
		// This is not necessary, but we ensure here that they are the same as the residual parameters.
		Xs: params.Xs(),
	}

	//===================================
	//=== 3) BOOTSTRAPPING PARAMETERS ===
	//===================================

	// Now that the residual parameters and the bootstrapping parameters literals are defined, we can instantiate
	// the bootstrapping parameters.
	// The instantiated bootstrapping parameters store their own ckks.Parameter, which are the parameters of the
	// ring used by the bootstrapping circuit.
	// The bootstrapping parameters are a wrapper of ckks.Parameters, with additional information.
	// They therefore has the same API as the ckks.Parameters and we can use this API to print some information.
	btpParams, err := bootstrapping.NewParametersFromLiteral(params, btpParametersLit)
	if err != nil {
		panic(err)
	}

	if *flagShort {
		// Corrects the message ratio Q0/|m(X)| to take into account the smaller number of slots and keep the same precision
		btpParams.Mod1ParametersLiteral.LogMessageRatio += 16 - params.LogN()
	}

	// We print some information about the residual parameters.
	fmt.Printf("Residual parameters: logN=%d, logSlots=%d, H=%d, sigma=%f, logQP=%f, levels=%d, scale=2^%d\n",
		btpParams.ResidualParameters.LogN(),
		btpParams.ResidualParameters.LogMaxSlots(),
		btpParams.ResidualParameters.XsHammingWeight(),
		btpParams.ResidualParameters.Xe(), params.LogQP(),
		btpParams.ResidualParameters.MaxLevel(),
		btpParams.ResidualParameters.LogDefaultScale())

	// And some information about the bootstrapping parameters.
	// We can notably check that the LogQP of the bootstrapping parameters is smaller than 1550, which ensures
	// 128-bit of security as explained above.
	fmt.Printf("Bootstrapping parameters: logN=%d, logSlots=%d, H(%d; %d), sigma=%f, logQP=%f, levels=%d, scale=2^%d\n",
		btpParams.BootstrappingParameters.LogN(),
		btpParams.BootstrappingParameters.LogMaxSlots(),
		btpParams.BootstrappingParameters.XsHammingWeight(),
		btpParams.EphemeralSecretWeight,
		btpParams.BootstrappingParameters.Xe(),
		btpParams.BootstrappingParameters.LogQP(),
		btpParams.BootstrappingParameters.QCount(),
		btpParams.BootstrappingParameters.LogDefaultScale())

	//===========================
	//=== 4) KEYGEN & ENCRYPT ===
	//===========================

	// Now that both the residual and bootstrapping parameters are instantiated, we can
	// instantiate the usual necessary object to encode, encrypt and decrypt.

	// Scheme context and keys
	kgen := rlwe.NewKeyGenerator(params)

	sk, pk := kgen.GenKeyPairNew()

	encoder := ckks.NewEncoder(params)
	decryptor := rlwe.NewDecryptor(params, sk)
	encryptor := rlwe.NewEncryptor(params, pk)

	fmt.Println()
	fmt.Println("Generating bootstrapping evaluation keys...")
	evk, _, err := btpParams.GenEvaluationKeys(sk)
	if err != nil {
		panic(err)
	}
	fmt.Println("Done")

	//========================
	//=== 5) BOOTSTRAPPING ===
	//========================

	// Instantiates the bootstrapper
	var eval *bootstrapping.Evaluator
	if eval, err = bootstrapping.NewEvaluator(btpParams, evk); err != nil {
		panic(err)
	}

	// Generate a random plaintext with values uniformely distributed in [-1, 1] for the real and imaginary part.
	valuesWant := make([]*bignum.Complex, params.MaxSlots())
	for i := range valuesWant {
		valuesWant[i] = &bignum.Complex{
			bignum.NewFloat(sampling.RandFloat64(-1, 1), prec),
			bignum.NewFloat(sampling.RandFloat64(-1, 1), prec),
		}
	}

	// We encrypt at level=LevelsConsumedPerRescaling-1
	plaintext := ckks.NewPlaintext(params, params.LevelsConsumedPerRescaling()-1)
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
	// and returns a ciphertext with the max level of `floatParamsResidualLit`.
	// CAUTION: the scale of the ciphertext MUST be equal (or very close) to params.DefaultScale()
	// To equalize the scale, the function evaluator.SetScale(ciphertext, parameters.DefaultScale()) can be used at the expense of one level.
	// If the ciphertext is is at level one or greater when given to the bootstrapper, this equalization is automatically done.
	fmt.Println("Bootstrapping...")
	ciphertext2, err := eval.Bootstrap(ciphertext1)

	if err != nil {
		panic(err)
	}
	fmt.Println("Done")

	//==================
	//=== 6) DECRYPT ===
	//==================

	// Decrypt, print and compare with the plaintext values
	fmt.Println()
	fmt.Println("Precision of ciphertext vs. Bootstrap(ciphertext)")
	printDebug(params, ciphertext2, valuesTest1, decryptor, encoder)
}

func printDebug(params ckks.Parameters, ciphertext *rlwe.Ciphertext, valuesWant []*bignum.Complex, decryptor *rlwe.Decryptor, encoder *ckks.Encoder) (valuesTest []*bignum.Complex) {

	valuesTest = make([]*bignum.Complex, ciphertext.Slots())

	if err := encoder.Decode(decryptor.DecryptNew(ciphertext), valuesTest); err != nil {
		panic(err)
	}

	fmt.Println()
	fmt.Printf("Level: %d (logQ = %d)\n", ciphertext.Level(), params.LogQLvl(ciphertext.Level()))

	fmt.Printf("Scale: 2^%f\n", math.Log2(ciphertext.Scale.Float64()))
	fmt.Printf("ValuesTest: %6.27f %6.27f...\n", valuesTest[0], valuesTest[1])
	fmt.Printf("ValuesWant: %6.27f %6.27f...\n", valuesWant[0], valuesWant[1])

	precStats := ckks.GetPrecisionStats(params, encoder, nil, valuesWant, valuesTest, 0, false)

	fmt.Println(precStats.String())
	fmt.Println()

	return
}
