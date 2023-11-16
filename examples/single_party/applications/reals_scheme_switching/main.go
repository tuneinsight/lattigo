// Package main showcases how lookup tables can complement fixed-point approximate
// homomorphic encryption to compute non-linear functions such as sign.
// The example starts by homomorphically decoding the ciphertext from the SIMD
// encoding to the coefficient encoding: IDFT(m(X)) -> m(X).
// It then evaluates a Lookup-Table (LUT) on each coefficient of m(X): m(X)[i] -> LUT(m(X)[i])
// and repacks each LUT(m(X)[i]) in a single RLWE ciphertext: Repack(LUT(m(X)[i])) -> LUT(m(X)).
// Finally, it homomorphically switches LUT(m(X)) back to the SIMD domain: LUT(m(X)) -> IDFT(LUT(m(X))).
package main

import (
	"flag"
	"fmt"
	"math/big"
	"time"

	"github.com/tuneinsight/lattigo/v5/core/rlwe"
	"github.com/tuneinsight/lattigo/v5/he/hebin"
	"github.com/tuneinsight/lattigo/v5/he/hefloat"
	"github.com/tuneinsight/lattigo/v5/ring"
	"github.com/tuneinsight/lattigo/v5/utils"
)

// ========================================
// Functions to evaluate with BlindRotation
// ========================================
func sign(x float64) (y float64) {
	if x > 0 {
		return 1
	} else if x < 0 {
		return -1
	} else {
		return 0
	}
}

func main() {

	var err error

	// Base ring degree
	LogN := 12

	// Q modulus Q
	Q := []uint64{0x800004001, 0x40002001} // 65.0000116961637 bits

	// P modulus P
	P := []uint64{0x4000026001} // 38.00000081692261 bits

	flagShort := flag.Bool("short", false, "runs the example with insecure parameters for fast testing")
	flag.Parse()

	if *flagShort {
		LogN = 6
	}

	LogSlots := 4
	slots := 1 << LogSlots

	// Starting RLWE params, size of these params
	// determine the complexity of the BlindRotation:
	// each BlindRotation takes ~N RGSW ciphertext-ciphertext mul.
	// LogN = 12 & LogQP = ~103 -> >128-bit secure.
	var paramsN12 hefloat.Parameters
	if paramsN12, err = hefloat.NewParametersFromLiteral(hefloat.ParametersLiteral{
		LogN:            LogN,
		Q:               Q,
		P:               P,
		LogDefaultScale: 32,
	}); err != nil {
		panic(err)
	}

	// BlindRotation RLWE params, N of these params determine
	// the test poly degree and therefore precision.
	// LogN = 11 & LogQP = ~54 -> 128-bit secure.
	var paramsN11 hefloat.Parameters
	if paramsN11, err = hefloat.NewParametersFromLiteral(hefloat.ParametersLiteral{
		LogN: LogN - 1,
		Q:    Q[:1],
		P:    []uint64{0x42001},
	}); err != nil {
		panic(err)
	}

	// Set the parameters for the blind rotation keys
	evkParams := rlwe.EvaluationKeyParameters{BaseTwoDecomposition: utils.Pointy(12)}

	// function interval
	a, b := -8.0, 8.0

	// Rescale inputs during Homomorphic Decoding by the normalization of the
	// test poly inputs and change of scale to ensure that upperbound on the homomorphic
	// decryption of LWE during the BlindRotation evaluation X^{dec(lwe)} is smaller than N
	// to avoid negacyclic wrapping of X^{dec(lwe)}.
	diffScale := float64(paramsN11.Q()[0]) / (4.0 * paramsN12.DefaultScale().Float64())
	normalization := 2.0 / (b - a) // all inputs are normalized before the BlindRotation evaluation.

	// SlotsToCoeffsParameters homomorphic encoding parameters
	var SlotsToCoeffsParameters = hefloat.DFTMatrixLiteral{
		Type:       hefloat.HomomorphicDecode,
		LogSlots:   LogSlots,
		Scaling:    new(big.Float).SetFloat64(normalization * diffScale),
		LevelStart: 1,        // starting level
		Levels:     []int{1}, // Decomposition levels of the encoding matrix (this will use one one matrix in one level)
	}

	// CoeffsToSlotsParameters homomorphic decoding parameters
	var CoeffsToSlotsParameters = hefloat.DFTMatrixLiteral{
		Type:       hefloat.HomomorphicEncode,
		LogSlots:   LogSlots,
		LevelStart: 1,        // starting level
		Levels:     []int{1}, // Decomposition levels of the encoding matrix (this will use one one matrix in one level)
	}

	fmt.Printf("Generating Test Poly... ")
	now := time.Now()
	// Generate test polynomial, provide function, outputscale, ring and interval.
	testPoly := hebin.InitTestPolynomial(sign, paramsN12.DefaultScale(), paramsN12.RingQ(), a, b)
	fmt.Printf("Done (%s)\n", time.Since(now))

	// Index of the test poly and repacking after evaluating the BlindRotation.
	testPolyMap := make(map[int]*ring.Poly) // Which slot to evaluate on the BlindRotation
	repackIndex := make(map[int]int)        // Where to repack slots after the BlindRotation
	gapN11 := paramsN11.N() / (2 * slots)
	gapN12 := paramsN12.N() / (2 * slots)

	for i := 0; i < slots; i++ {
		testPolyMap[i*gapN11] = &testPoly
		repackIndex[i*gapN11] = i * gapN12
	}

	kgenN12 := rlwe.NewKeyGenerator(paramsN12)
	skN12 := kgenN12.GenSecretKeyNew()
	encoderN12 := hefloat.NewEncoder(paramsN12)
	encryptorN12 := rlwe.NewEncryptor(paramsN12, skN12)
	decryptorN12 := rlwe.NewDecryptor(paramsN12, skN12)

	kgenN11 := rlwe.NewKeyGenerator(paramsN11)
	skN11 := kgenN11.GenSecretKeyNew()

	// EvaluationKey RLWEN12 -> RLWEN11
	evkN12ToN11 := rlwe.NewKeyGenerator(paramsN12).GenEvaluationKeyNew(skN12, skN11)

	fmt.Printf("Gen SlotsToCoeffs Matrices... ")
	now = time.Now()
	SlotsToCoeffsMatrix, err := hefloat.NewDFTMatrixFromLiteral(paramsN12, SlotsToCoeffsParameters, encoderN12)
	if err != nil {
		panic(err)
	}
	CoeffsToSlotsMatrix, err := hefloat.NewDFTMatrixFromLiteral(paramsN12, CoeffsToSlotsParameters, encoderN12)
	if err != nil {
		panic(err)
	}
	fmt.Printf("Done (%s)\n", time.Since(now))

	// GaloisKeys
	galEls := paramsN12.GaloisElementsForTrace(0)
	galEls = append(galEls, SlotsToCoeffsParameters.GaloisElements(paramsN12)...)
	galEls = append(galEls, CoeffsToSlotsParameters.GaloisElements(paramsN12)...)
	galEls = append(galEls, paramsN12.GaloisElementForComplexConjugation())

	evk := rlwe.NewMemEvaluationKeySet(nil, kgenN12.GenGaloisKeysNew(galEls, skN12)...)

	// BlindRotation Evaluator
	evalBR := hebin.NewEvaluator(paramsN12, paramsN11)

	// Evaluator
	eval := hefloat.NewEvaluator(paramsN12, evk)
	evalHDFT := hefloat.NewDFTEvaluator(paramsN12, eval)

	fmt.Printf("Encrypting bits of skLWE in RGSW... ")
	now = time.Now()
	blindRotateKey := hebin.GenEvaluationKeyNew(paramsN12, skN12, paramsN11, skN11, evkParams) // Generate RGSW(sk_i) for all coefficients of sk
	fmt.Printf("Done (%s)\n", time.Since(now))

	// Generates the starting plaintext values.
	interval := (b - a) / float64(slots)
	values := make([]float64, slots)
	for i := 0; i < slots; i++ {
		values[i] = a + float64(i)*interval
	}

	pt := hefloat.NewPlaintext(paramsN12, paramsN12.MaxLevel())
	pt.LogDimensions.Cols = LogSlots
	if err := encoderN12.Encode(values, pt); err != nil {
		panic(err)
	}

	ctN12, err := encryptorN12.EncryptNew(pt)
	if err != nil {
		panic(err)
	}

	fmt.Printf("Homomorphic Decoding... ")
	now = time.Now()

	// Homomorphic Decoding: [(a+bi), (c+di)] -> [a, c, b, d]
	ctN12, err = evalHDFT.SlotsToCoeffsNew(ctN12, nil, SlotsToCoeffsMatrix)
	if err != nil {
		panic(err)
	}
	ctN12.IsBatched = false

	// Key-Switch from LogN = 12 to LogN = 11
	ctN11 := hefloat.NewCiphertext(paramsN11, 1, paramsN11.MaxLevel())
	// key-switch to LWE degree
	if err := eval.ApplyEvaluationKey(ctN12, evkN12ToN11, ctN11); err != nil {
		panic(err)
	}
	fmt.Printf("Done (%s)\n", time.Since(now))

	fmt.Printf("Evaluating BlindRotations... ")
	now = time.Now()
	// Extracts & EvalBR(LWEs, indexTestPoly) on the fly -> Repack(LWEs, indexRepack) -> RLWE
	ctN12, err = evalBR.EvaluateAndRepack(ctN11, testPolyMap, repackIndex, blindRotateKey, evk)
	if err != nil {
		panic(err)
	}
	fmt.Printf("Done (%s)\n", time.Since(now))
	ctN12.IsBatched = false
	ctN12.LogDimensions = paramsN12.LogMaxDimensions()
	ctN12.Scale = paramsN12.DefaultScale()

	fmt.Printf("Homomorphic Encoding... ")
	now = time.Now()
	// Homomorphic Encoding: [BR(a), BR(c), BR(b), BR(d)] -> [(BR(a)+BR(b)i), (BR(c)+BR(d)i)]
	ctN12, _, err = evalHDFT.CoeffsToSlotsNew(ctN12, CoeffsToSlotsMatrix)
	if err != nil {
		panic(err)
	}
	fmt.Printf("Done (%s)\n", time.Since(now))

	res := make([]float64, slots)
	ctN12.IsBatched = true
	ctN12.LogDimensions.Cols = LogSlots
	if err := encoderN12.Decode(decryptorN12.DecryptNew(ctN12), res); err != nil {
		panic(err)
	}
	for i, v := range res {
		fmt.Printf("%7.4f -> %7.4f\n", values[i], v)
	}
}
