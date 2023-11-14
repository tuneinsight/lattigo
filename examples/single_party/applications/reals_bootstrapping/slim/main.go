// Package main implements an example showcasing slim for bootstrapping for fixed-point approximate
// arithmetic over the reals/complexes numbers.
// This re-ordering of the bootstrapping steps was first proposed for the BFV/BGV schemes by Chen and Han
// in Homomorphic Lower Digits Removal and Improved FHE Bootstrapping (https://eprint.iacr.org/2018/067).
// It was also used by Kim and Guyot in Optimized Privacy-Preserving CNN Inference With Fully Homomorphic
// Encryption (https://ieeexplore.ieee.org/document/10089847) to efficiently perform the convolution in
// the coefficient domain.
//
// This example assumes that the user is already familiar with the bootstrapping and its different steps.
// See the basic example `lattigo/examples/he/hefloat/bootstrapping/basic` for an introduction into the
// bootstrapping.
//
// The usual order of the bootstrapping operations is:
//
// 0) User defined circuit in the slots domain
// 1) ScaleDown: Scale the ciphertext to q0/|m|
// 2) ModUp: Raise modulus from q0 to qL
// 3) CoeffsToSlots: Homomorphic encoding
// 4) EvalMod: Homomorphic modular reduction
// 5) SlotsToCoeffs (and go back to 0): Homomorphic Decoding
//
// This example instantiates a custom order of the circuit evaluating:
//
// 0) User defined circuit in the slots domain
// 1) SlotsToCoeffs: Homomorphic Decoding
// 2) User defined circuit in the coeffs domain
// 3) ScaleDown: Scale the ciphertext to q0/|m|
// 4) ModUp: Raise modulus from q0 to qL
// 5) CoeffsToSlots: Homomorphic encoding
// 6) EvalMod (and to back to 0): Homomorphic modular reduction
//
// Use the flag -short to run the examples fast but with insecure parameters.
package main

import (
	"flag"
	"fmt"
	"math"
	"math/big"

	"github.com/tuneinsight/lattigo/v5/core/rlwe"
	"github.com/tuneinsight/lattigo/v5/he/hefloat"
	"github.com/tuneinsight/lattigo/v5/he/hefloat/bootstrapping"
	"github.com/tuneinsight/lattigo/v5/ring"
	"github.com/tuneinsight/lattigo/v5/utils"
	"github.com/tuneinsight/lattigo/v5/utils/sampling"
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

	//============================
	//=== 1) SCHEME PARAMETERS ===
	//============================

	// In this example, for a pratical purpose, the residual parameters and bootstrapping
	// parameters are the same. But in practice the residual parameters would not contain the
	// moduli for the CoeffsToSlots and EvalMod steps.
	// With LogN=16, LogQP=1221 and H=192, these parameters achieve well over 128-bit of security.
	// For the purpose of the example, only one prime is allocated to the circuit in the slots domain
	// and no prime is allocated to the circuit in the coeffs domain.

	LogDefaultScale := 40

	q0 := []int{55}                                    // 3) ScaleDown & 4) ModUp
	qiSlotsToCoeffs := []int{39, 39, 39}               // 1) SlotsToCoeffs
	qiCircuitSlots := []int{LogDefaultScale}           // 0) Circuit in the slot domain
	qiEvalMod := []int{60, 60, 60, 60, 60, 60, 60, 60} // 6) EvalMod
	qiCoeffsToSlots := []int{56, 56, 56, 56}           // 5) CoeffsToSlots

	LogQ := append(q0, qiSlotsToCoeffs...)
	LogQ = append(LogQ, qiCircuitSlots...)
	LogQ = append(LogQ, qiEvalMod...)
	LogQ = append(LogQ, qiCoeffsToSlots...)

	params, err := hefloat.NewParametersFromLiteral(hefloat.ParametersLiteral{
		LogN:            LogN,                      // Log2 of the ring degree
		LogQ:            LogQ,                      // Log2 of the ciphertext modulus
		LogP:            []int{61, 61, 61, 61, 61}, // Log2 of the key-switch auxiliary prime moduli
		LogDefaultScale: LogDefaultScale,           // Log2 of the scale
		Xs:              ring.Ternary{H: 192},
	})

	if err != nil {
		panic(err)
	}

	//====================================
	//=== 2) BOOTSTRAPPING PARAMETERS ===
	//====================================

	// CoeffsToSlots parameters (homomorphic encoding)
	CoeffsToSlotsParameters := hefloat.DFTMatrixLiteral{
		Type:         hefloat.HomomorphicEncode,
		Format:       hefloat.RepackImagAsReal, // Returns the real and imaginary part into separate ciphertexts
		LogSlots:     params.LogMaxSlots(),
		LevelStart:   params.MaxLevel(),
		LogBSGSRatio: 1,
		Levels:       []int{1, 1, 1, 1}, //qiCoeffsToSlots
	}

	// Parameters of the homomorphic modular reduction x mod 1
	Mod1ParametersLiteral := hefloat.Mod1ParametersLiteral{
		LogScale:        60,                  // Matches qiEvalMod
		Mod1Type:        hefloat.CosDiscrete, // Multi-interval Chebyshev interpolation
		Mod1Degree:      30,                  // Depth 5
		DoubleAngle:     3,                   // Depth 3
		K:               16,                  // With EphemeralSecretWeight = 32 and 2^{15} slots, ensures < 2^{-138.7} failure probability
		LogMessageRatio: 5,                   // q/|m| = 2^5
		Mod1InvDegree:   0,                   // Depth 0
		LevelStart:      params.MaxLevel() - len(CoeffsToSlotsParameters.Levels),
	}

	// Since we scale the values by 1/2^{LogMessageRatio} during CoeffsToSlots,
	// we must scale them back by 2^{LogMessageRatio} after EvalMod.
	// This is done by scaling the EvalMod polynomial coefficients by 2^{LogMessageRatio}.
	Mod1ParametersLiteral.Scaling = math.Exp2(-float64(Mod1ParametersLiteral.LogMessageRatio))

	// SlotsToCoeffs parameters (homomorphic decoding)
	SlotsToCoeffsParameters := hefloat.DFTMatrixLiteral{
		Type:         hefloat.HomomorphicDecode,
		LogSlots:     params.LogMaxSlots(),
		Scaling:      new(big.Float).SetFloat64(math.Exp2(float64(Mod1ParametersLiteral.LogMessageRatio))),
		LogBSGSRatio: 1,
		Levels:       []int{1, 1, 1}, // qiSlotsToCoeffs
	}

	SlotsToCoeffsParameters.LevelStart = len(SlotsToCoeffsParameters.Levels)

	// Custom bootstrapping.Parameters.
	// All fields are public and can be manually instantiated.
	btpParams := bootstrapping.Parameters{
		ResidualParameters:      params,
		BootstrappingParameters: params,
		SlotsToCoeffsParameters: SlotsToCoeffsParameters,
		Mod1ParametersLiteral:   Mod1ParametersLiteral,
		CoeffsToSlotsParameters: CoeffsToSlotsParameters,
		EphemeralSecretWeight:   32, // > 128bit secure for LogN=16 and LogQP = 115.
		CircuitOrder:            bootstrapping.DecodeThenModUp,
	}

	if *flagShort {
		// Corrects the message ratio Q0/|m(X)| to take into account the smaller number of slots and keep the same precision
		btpParams.Mod1ParametersLiteral.LogMessageRatio += 16 - params.LogN()
	}

	// We pring some information about the bootstrapping parameters (which are identical to the residual parameters in this example).
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
	//=== 3) KEYGEN & ENCRYPT ===
	//===========================

	// Now that both the residual and bootstrapping parameters are instantiated, we can
	// instantiate the usual necessary object to encode, encrypt and decrypt.

	// Scheme context and keys
	kgen := rlwe.NewKeyGenerator(params)

	sk, pk := kgen.GenKeyPairNew()

	encoder := hefloat.NewEncoder(params)
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
	//=== 4) BOOTSTRAPPING ===
	//========================

	// Instantiates the bootstrapper
	var eval *bootstrapping.Evaluator
	if eval, err = bootstrapping.NewEvaluator(btpParams, evk); err != nil {
		panic(err)
	}

	// Generate a random plaintext with values uniformely distributed in [-1, 1] for the real and imaginary part.
	valuesWant := make([]complex128, params.MaxSlots())
	for i := range valuesWant {
		valuesWant[i] = sampling.RandComplex128(-1, 1)
	}

	// We encrypt at level 0
	plaintext := hefloat.NewPlaintext(params, SlotsToCoeffsParameters.LevelStart)
	if err := encoder.Encode(valuesWant, plaintext); err != nil {
		panic(err)
	}

	// Encrypt
	ciphertext, err := encryptor.EncryptNew(plaintext)
	if err != nil {
		panic(err)
	}

	// Decrypt, print and compare with the plaintext values
	fmt.Println()
	fmt.Println("Precision of values vs. ciphertext")
	valuesTest := printDebug(params, ciphertext, valuesWant, decryptor, encoder)

	fmt.Println("Bootstrapping...")

	// Step 0: Some circuit in the slots domain

	// Step 1 : SlotsToCoeffs (Homomorphic decoding)
	if ciphertext, err = eval.SlotsToCoeffs(ciphertext, nil); err != nil {
		panic(err)
	}

	// Step 2: Some circuit in the coefficient domain
	// Note: the result of SlotsToCoeffs is naturaly given in bit-reversed order
	// In this example, we multiply by the monomial X^{N/2} (which is the imaginary
	// unit in the slots domain)
	if err = eval.Evaluator.Mul(ciphertext, 1i, ciphertext); err != nil {
		panic(err)
	}

	// Then we need to apply the same mapping to the reference values:

	// Maps C^{N/2} to R[X]/(X^N+1) (bit-reversed)
	utils.BitReverseInPlaceSlice(valuesTest, len(valuesTest))
	valuesTestFloat := make([]float64, ciphertext.Slots()*2)
	for i, j := 0, params.N()/2; i < params.N()/2; i, j = i+1, j+1 {
		valuesTestFloat[i] = real(valuesTest[i])
		valuesTestFloat[j] = imag(valuesTest[i])
	}

	// Multiplication by X^{N/2}
	utils.RotateSliceInPlace(valuesTestFloat, -params.N()/2)
	for i := 0; i < params.N()/2; i++ {
		valuesTestFloat[i] *= -1
	}

	// Maps R[X]/(X^N+1) to C^{N/2} (bit-reversed)
	for i, j := 0, params.N()/2; i < params.N()/2; i, j = i+1, j+1 {
		valuesTest[i] = complex(valuesTestFloat[i], valuesTestFloat[j])
	}
	utils.BitReverseInPlaceSlice(valuesTest, len(valuesTest))

	// Step 3: scale to q/|m|
	if ciphertext, _, err = eval.ScaleDown(ciphertext); err != nil {
		panic(err)
	}

	// Step 4 : Extend the basis from q to Q
	if ciphertext, err = eval.ModUp(ciphertext); err != nil {
		panic(err)
	}

	// Step 5 : CoeffsToSlots (Homomorphic encoding)
	// Note: expects the result to be given in bit-reversed order
	// Also, we need the homomorphic encoding to split the real and
	// imaginary parts into two pure real ciphertexts, because the
	// homomorphic modular reduction is only defined on the reals.
	// The `imag` ciphertext can be ignored if the original input
	// is purely real.
	var real, imag *rlwe.Ciphertext
	if real, imag, err = eval.CoeffsToSlots(ciphertext); err != nil {
		panic(err)
	}

	// Step 6 : EvalMod (Homomorphic modular reduction)
	if real, err = eval.EvalMod(real); err != nil {
		panic(err)
	}

	if imag, err = eval.EvalMod(imag); err != nil {
		panic(err)
	}

	// Recombines the real and imaginary part
	if err = eval.Evaluator.Mul(imag, 1i, imag); err != nil {
		panic(err)
	}

	if err = eval.Evaluator.Add(real, imag, ciphertext); err != nil {
		panic(err)
	}

	fmt.Println("Done")

	//==================
	//=== 5) DECRYPT ===
	//==================

	// Decrypt, print and compare with the plaintext values
	fmt.Println()
	fmt.Println("Precision of ciphertext vs. Bootstrap(ciphertext)")
	printDebug(params, ciphertext, valuesTest, decryptor, encoder)
}

func printDebug(params hefloat.Parameters, ciphertext *rlwe.Ciphertext, valuesWant []complex128, decryptor *rlwe.Decryptor, encoder *hefloat.Encoder) (valuesTest []complex128) {

	slots := ciphertext.Slots()

	if !ciphertext.IsBatched {
		slots *= 2
	}

	valuesTest = make([]complex128, slots)

	if err := encoder.Decode(decryptor.DecryptNew(ciphertext), valuesTest); err != nil {
		panic(err)
	}

	fmt.Println()
	fmt.Printf("Level: %d (logQ = %d)\n", ciphertext.Level(), params.LogQLvl(ciphertext.Level()))

	fmt.Printf("Scale: 2^%f\n", math.Log2(ciphertext.Scale.Float64()))
	fmt.Printf("ValuesTest: %10.14f %10.14f %10.14f %10.14f...\n", valuesTest[0], valuesTest[1], valuesTest[2], valuesTest[3])
	fmt.Printf("ValuesWant: %10.14f %10.14f %10.14f %10.14f...\n", valuesWant[0], valuesWant[1], valuesWant[2], valuesWant[3])

	precStats := hefloat.GetPrecisionStats(params, encoder, nil, valuesWant, valuesTest, 0, false)

	fmt.Println(precStats.String())
	fmt.Println()

	return
}
