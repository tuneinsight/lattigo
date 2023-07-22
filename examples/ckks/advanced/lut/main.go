package main

import (
	"flag"
	"fmt"
	"math/big"
	"time"

	"github.com/tuneinsight/lattigo/v4/ckks"
	"github.com/tuneinsight/lattigo/v4/rgsw/lut"
	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/rlwe"
)

// This example showcases how lookup tables can complement the CKKS scheme to compute non-linear functions
// such as sign. The example starts by homomorphically decoding the CKKS ciphertext from the canonical embeding
// to the coefficient embeding. It then evaluates the Look-Up-Table (LUT) on each coefficient and repacks the
// outputs of each LUT in a single RLWE ciphertext. Finally, it homomorphically encodes the RLWE ciphertext back
// to the canonical embeding of the CKKS scheme.

// ==============================
// Functions to evaluate with LUT
// ==============================
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
	// determine the complexity of the LUT:
	// each LUT takes N RGSW ciphertext-ciphetext mul.
	// LogN = 12 & LogQP = ~103 -> >128-bit secure.
	var paramsN12 ckks.Parameters
	if paramsN12, err = ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
		LogN:              LogN,
		Q:                 Q,
		P:                 P,
		LogPlaintextScale: 32,
	}); err != nil {
		panic(err)
	}

	// LUT RLWE params, N of these params determine
	// the LUT poly and therefore precision.
	// LogN = 11 & LogQP = ~54 -> 128-bit secure.
	var paramsN11 ckks.Parameters
	if paramsN11, err = ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
		LogN: LogN - 1,
		Q:    Q[:1],
		P:    []uint64{0x42001},
	}); err != nil {
		panic(err)
	}

	Base2Decomposition := 12

	// LUT interval
	a, b := -8.0, 8.0

	// Rescale inputs during Homomorphic Decoding by the normalization of the
	// LUT inputs and change of scale to ensure that upperbound on the homomorphic
	// decryption of LWE during the LUT evaluation X^{dec(lwe)} is smaller than N
	// to avoid negacyclic wrapping of X^{dec(lwe)}.
	diffScale := float64(paramsN11.Q()[0]) / (4.0 * paramsN12.PlaintextScale().Float64())
	normalization := 2.0 / (b - a) // all inputs are normalized before the LUT evaluation.

	// SlotsToCoeffsParameters homomorphic encoding parameters
	var SlotsToCoeffsParameters = ckks.HomomorphicDFTMatrixLiteral{
		Type:       ckks.Decode,
		LogSlots:   LogSlots,
		Scaling:    new(big.Float).SetFloat64(normalization * diffScale),
		LevelStart: 1,        // starting level
		Levels:     []int{1}, // Decomposition levels of the encoding matrix (this will use one one matrix in one level)
	}

	// CoeffsToSlotsParameters homomorphic decoding parameters
	var CoeffsToSlotsParameters = ckks.HomomorphicDFTMatrixLiteral{
		Type:       ckks.Encode,
		LogSlots:   LogSlots,
		LevelStart: 1,        // starting level
		Levels:     []int{1}, // Decomposition levels of the encoding matrix (this will use one one matrix in one level)
	}

	fmt.Printf("Generating LUT... ")
	now := time.Now()
	// Generate LUT, provide function, outputscale, ring and interval.
	LUTPoly := lut.InitLUT(sign, paramsN12.PlaintextScale(), paramsN12.RingQ(), a, b)
	fmt.Printf("Done (%s)\n", time.Since(now))

	// Index of the LUT poly and repacking after evaluating the LUT.
	lutPolyMap := make(map[int]*ring.Poly) // Which slot to evaluate on the LUT
	repackIndex := make(map[int]int)       // Where to repack slots after the LUT
	gapN11 := paramsN11.N() / (2 * slots)
	gapN12 := paramsN12.N() / (2 * slots)

	for i := 0; i < slots; i++ {
		lutPolyMap[i*gapN11] = &LUTPoly
		repackIndex[i*gapN11] = i * gapN12
	}

	kgenN12 := ckks.NewKeyGenerator(paramsN12)
	skN12 := kgenN12.GenSecretKeyNew()
	encoderN12 := ckks.NewEncoder(paramsN12)
	encryptorN12, err := ckks.NewEncryptor(paramsN12, skN12)
	if err != nil {
		panic(err)
	}
	decryptorN12, err := ckks.NewDecryptor(paramsN12, skN12)
	if err != nil {
		panic(err)
	}

	kgenN11 := ckks.NewKeyGenerator(paramsN11)
	skN11 := kgenN11.GenSecretKeyNew()

	// EvaluationKey RLWEN12 -> RLWEN11
	evkN12ToN11, err := ckks.NewKeyGenerator(paramsN12).GenEvaluationKeyNew(skN12, skN11)
	if err != nil {
		panic(err)
	}

	fmt.Printf("Gen SlotsToCoeffs Matrices... ")
	now = time.Now()
	SlotsToCoeffsMatrix, err := ckks.NewHomomorphicDFTMatrixFromLiteral(SlotsToCoeffsParameters, encoderN12)
	if err != nil {
		panic(err)
	}
	CoeffsToSlotsMatrix, err := ckks.NewHomomorphicDFTMatrixFromLiteral(CoeffsToSlotsParameters, encoderN12)
	if err != nil {
		panic(err)
	}
	fmt.Printf("Done (%s)\n", time.Since(now))

	// GaloisKeys
	galEls := paramsN12.GaloisElementsForTrace(0)
	galEls = append(galEls, SlotsToCoeffsParameters.GaloisElements(paramsN12)...)
	galEls = append(galEls, CoeffsToSlotsParameters.GaloisElements(paramsN12)...)
	galEls = append(galEls, paramsN12.GaloisElementForComplexConjugation())

	gks, err := kgenN12.GenGaloisKeysNew(galEls, skN12)
	if err != nil {
		panic(err)
	}

	evk := rlwe.NewMemEvaluationKeySet(nil, gks...)

	// LUT Evaluator
	evalLUT := lut.NewEvaluator(paramsN12.Parameters, paramsN11.Parameters, Base2Decomposition, evk)

	// CKKS Evaluator
	evalCKKS := ckks.NewEvaluator(paramsN12, evk)

	fmt.Printf("Encrypting bits of skLWE in RGSW... ")
	now = time.Now()
	LUTKEY, err := lut.GenEvaluationKeyNew(paramsN12.Parameters, skN12, paramsN11.Parameters, skN11, Base2Decomposition) // Generate RGSW(sk_i) for all coefficients of sk
	if err != nil {
		panic(err)
	}
	fmt.Printf("Done (%s)\n", time.Since(now))

	// Generates the starting plaintext values.
	interval := (b - a) / float64(slots)
	values := make([]float64, slots)
	for i := 0; i < slots; i++ {
		values[i] = a + float64(i)*interval
	}

	pt := ckks.NewPlaintext(paramsN12, paramsN12.MaxLevel())
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
	ctN12, err = evalCKKS.SlotsToCoeffsNew(ctN12, nil, SlotsToCoeffsMatrix)
	if err != nil {
		panic(err)
	}
	ctN12.IsBatched = false

	// Key-Switch from LogN = 12 to LogN = 11
	ctN11 := ckks.NewCiphertext(paramsN11, 1, paramsN11.MaxLevel())
	// key-switch to LWE degree
	if err := evalCKKS.ApplyEvaluationKey(ctN12, evkN12ToN11, ctN11); err != nil {
		panic(err)
	}
	fmt.Printf("Done (%s)\n", time.Since(now))

	fmt.Printf("Evaluating LUT... ")
	now = time.Now()
	// Extracts & EvalLUT(LWEs, indexLUT) on the fly -> Repack(LWEs, indexRepack) -> RLWE
	ctN12, err = evalLUT.EvaluateAndRepack(ctN11, lutPolyMap, repackIndex, LUTKEY)
	if err != nil {
		panic(err)
	}
	fmt.Printf("Done (%s)\n", time.Since(now))
	ctN12.IsBatched = false
	ctN12.LogDimensions = paramsN12.PlaintextLogDimensions()
	ctN12.Scale = paramsN12.PlaintextScale()

	fmt.Println(ctN12.MetaData)

	fmt.Printf("Homomorphic Encoding... ")
	now = time.Now()
	// Homomorphic Encoding: [LUT(a), LUT(c), LUT(b), LUT(d)] -> [(LUT(a)+LUT(b)i), (LUT(c)+LUT(d)i)]
	ctN12, _, err = evalCKKS.CoeffsToSlotsNew(ctN12, CoeffsToSlotsMatrix)
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
