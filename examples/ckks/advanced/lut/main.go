package main

import (
	"flag"
	"fmt"
	"github.com/tuneinsight/lattigo/v3/ckks"
	ckksAdvanced "github.com/tuneinsight/lattigo/v3/ckks/advanced"
	"github.com/tuneinsight/lattigo/v3/ring"
	"github.com/tuneinsight/lattigo/v3/rlwe"
	"github.com/tuneinsight/lattigo/v3/rlwe/lut"
	"time"
)

// This examples showcases how lookup tables can complement the CKKS scheme to compute non-linear functions
// such as sign. The examples starts by homomorphically decoding the CKKS ciphertext from the canonical embeding
// to the coefficient embeding. It then evaluates the Look-Up-Table (LUT) on each coefficient and repacks the
// outputs of each LUT in a single RLWE ciphertext. Finally it homomorphically encodes the RLWE ciphertext back
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

	// Starting RLWE params, size of these params
	// determine the complexity of the LUT:
	// each LUT takes N RGSW ciphertext-ciphetext mul.
	// LogN = 12 & LogQP = ~103 -> >128-bit secure.
	var paramsN12 ckks.Parameters
	if paramsN12, err = ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
		ParametersLiteral: rlwe.ParametersLiteral{
			LogN:     LogN,
			Q:        Q,
			P:        P,
			LogBase2: 0,
			H:        0,
			Sigma:    rlwe.DefaultSigma,
			RingType: ring.Standard,
		},
		LogSlots:     4,
		DefaultScale: 1 << 32,
	}); err != nil {
		panic(err)
	}

	// Params for Key-switching N12 to N11.
	// LogN = 12 & LogQP = ~54 -> >>>128-bit secure.
	var paramsN12ToN11 ckks.Parameters
	if paramsN12ToN11, err = ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
		ParametersLiteral: rlwe.ParametersLiteral{
			LogN:     LogN,
			Q:        Q[:1],
			P:        []uint64{0x42001},
			LogBase2: 16,
			H:        0,
			Sigma:    rlwe.DefaultSigma,
			RingType: ring.Standard,
		},
		LogSlots:     0,
		DefaultScale: 0,
	}); err != nil {
		panic(err)
	}

	// LUT RLWE params, N of these params determine
	// the LUT poly and therefore precision.
	// LogN = 11 & LogQP = ~54 -> 128-bit secure.
	var paramsN11 ckks.Parameters
	if paramsN11, err = ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
		ParametersLiteral: rlwe.ParametersLiteral{
			LogN:     LogN - 1,
			Q:        Q[:1],
			P:        []uint64{0x42001},
			LogBase2: 12,
			H:        0,
			Sigma:    rlwe.DefaultSigma,
			RingType: ring.Standard,
		},
		LogSlots:     0,
		DefaultScale: 0,
	}); err != nil {
		panic(err)
	}

	// LUT interval
	a, b := -8.0, 8.0

	// Rescale inputs during Homomorphic Decoding by the normalization of the
	// LUT inputs and change of scale to ensure that upperbound on the homomorphic
	// decryption of LWE during the LUT evaluation X^{dec(lwe)} is smaller than N
	// to avoid negacyclic wrapping of X^{dec(lwe)}.
	diffScale := paramsN11.QiFloat64(0) / (4.0 * paramsN12.DefaultScale())
	normalization := 2.0 / (b - a) // all inputs are normalized before the LUT evaluation.

	// SlotsToCoeffsParameters homomorphic encoding parameters
	var SlotsToCoeffsParameters = ckksAdvanced.EncodingMatrixLiteral{
		LogN:                paramsN12.LogN(),
		LogSlots:            paramsN12.LogSlots(),
		Scaling:             normalization * diffScale,
		LinearTransformType: ckksAdvanced.SlotsToCoeffs,
		RepackImag2Real:     false,
		LevelStart:          1,     // starting level
		BSGSRatio:           4.0,   // ratio between n1/n2 for n1*n2 = slots
		BitReversed:         false, // bit-reversed input
		ScalingFactor: [][]float64{ // Decomposition level of the encoding matrix
			{paramsN12.QiFloat64(1)}, // Scale of the decoding matrix
		},
	}

	// CoeffsToSlotsParameters homomorphic decoding parameters
	var CoeffsToSlotsParameters = ckksAdvanced.EncodingMatrixLiteral{
		LinearTransformType: ckksAdvanced.CoeffsToSlots,
		RepackImag2Real:     false,
		LogN:                paramsN12.LogN(),
		LogSlots:            paramsN12.LogSlots(),
		Scaling:             1 / float64(paramsN12.Slots()),
		LevelStart:          1,     // starting level
		BSGSRatio:           4.0,   // ratio between n1/n2 for n1*n2 = slots
		BitReversed:         false, // bit-reversed input
		ScalingFactor: [][]float64{ // Decomposition level of the encoding matrix
			{paramsN12.QiFloat64(1)}, // Scale of the encoding matrix
		},
	}

	fmt.Printf("Generating LUT... ")
	now := time.Now()
	// Generate LUT, provide function, outputscale, ring and interval.
	LUTPoly := lut.InitLUT(sign, paramsN12.DefaultScale(), paramsN12.RingQ(), a, b)
	fmt.Printf("Done (%s)\n", time.Since(now))

	// Index of the LUT poly and repacking after evaluating the LUT.
	lutPolyMap := make(map[int]*ring.Poly) // Which slot to evaluate on the LUT
	repackIndex := make(map[int]int)       // Where to repack slots after the LUT
	gapN11 := paramsN11.N() / (2 * paramsN12.Slots())
	gapN12 := paramsN12.N() / (2 * paramsN12.Slots())

	for i := 0; i < paramsN12.Slots(); i++ {
		lutPolyMap[i*gapN11] = LUTPoly
		repackIndex[i*gapN11] = i * gapN12
	}

	kgenN12 := ckks.NewKeyGenerator(paramsN12)
	skN12 := kgenN12.GenSecretKey()
	encoderN12 := ckks.NewEncoder(paramsN12)
	encryptorN12 := ckks.NewEncryptor(paramsN12, skN12)
	decryptorN12 := ckks.NewDecryptor(paramsN12, skN12)

	kgenN11 := ckks.NewKeyGenerator(paramsN11)
	skN11 := kgenN11.GenSecretKey()
	//decryptorN11 := ckks.NewDecryptor(paramsN11, skN11)
	//encoderN11 := ckks.NewEncoder(paramsN11)

	// Switchingkey RLWEN12 -> RLWEN11
	swkN12ToN11 := ckks.NewKeyGenerator(paramsN12ToN11).GenSwitchingKey(skN12, skN11)

	fmt.Printf("Gen SlotsToCoeffs Matrices... ")
	now = time.Now()
	SlotsToCoeffsMatrix := ckksAdvanced.NewHomomorphicEncodingMatrixFromLiteral(SlotsToCoeffsParameters, encoderN12)
	CoeffsToSlotsMatrix := ckksAdvanced.NewHomomorphicEncodingMatrixFromLiteral(CoeffsToSlotsParameters, encoderN12)
	fmt.Printf("Done (%s)\n", time.Since(now))

	// Rotation Keys
	rotations := []int{}
	for i := 1; i < paramsN12.N(); i <<= 1 {
		rotations = append(rotations, i)
	}

	rotations = append(rotations, SlotsToCoeffsParameters.Rotations()...)
	rotations = append(rotations, CoeffsToSlotsParameters.Rotations()...)

	rotKey := kgenN12.GenRotationKeysForRotations(rotations, true, skN12)

	// LUT Evaluator
	evalLUT := lut.NewEvaluator(paramsN12.Parameters, paramsN11.Parameters, rotKey)

	// CKKS Evaluator
	evalCKKS := ckksAdvanced.NewEvaluator(paramsN12, rlwe.EvaluationKey{Rlk: nil, Rtks: rotKey})
	evalCKKSN12ToN11 := ckks.NewEvaluator(paramsN12ToN11, rlwe.EvaluationKey{})

	fmt.Printf("Encrypting bits of skLWE in RGSW... ")
	now = time.Now()
	LUTKEY := lut.GenKey(paramsN12.Parameters, skN12, paramsN11.Parameters, skN11) // Generate RGSW(sk_i) for all coefficients of sk
	fmt.Printf("Done (%s)\n", time.Since(now))

	// Generates the starting plaintext values.
	interval := (b - a) / float64(paramsN12.Slots())
	values := make([]float64, paramsN12.Slots())
	for i := 0; i < paramsN12.Slots(); i++ {
		values[i] = a + float64(i)*interval
	}

	pt := ckks.NewPlaintext(paramsN12, paramsN12.MaxLevel(), paramsN12.DefaultScale())
	encoderN12.EncodeSlots(values, pt, paramsN12.LogSlots())
	ctN12 := encryptorN12.EncryptNew(pt)

	fmt.Printf("Homomorphic Decoding... ")
	now = time.Now()
	// Homomorphic Decoding: [(a+bi), (c+di)] -> [a, c, b, d]
	ctN12 = evalCKKS.SlotsToCoeffsNew(ctN12, nil, SlotsToCoeffsMatrix)
	ctN12.Scale = paramsN11.QiFloat64(0) / 4.0

	// Key-Switch from LogN = 12 to LogN = 10
	evalCKKS.DropLevel(ctN12, ctN12.Level())                    // drop to LUT level
	ctTmp := evalCKKSN12ToN11.SwitchKeysNew(ctN12, swkN12ToN11) // key-switch to LWE degree
	ctN11 := ckks.NewCiphertext(paramsN11, 1, paramsN11.MaxLevel(), ctTmp.Scale)
	rlwe.SwitchCiphertextRingDegreeNTT(ctTmp.Ciphertext, paramsN11.RingQ(), paramsN12.RingQ(), ctN11.Ciphertext)
	fmt.Printf("Done (%s)\n", time.Since(now))

	//for i, v := range encoderN11.DecodeCoeffs(decryptorN11.DecryptNew(ctN11)){
	//	fmt.Printf("%3d: %7.4f\n", i, v)
	//}

	fmt.Printf("Evaluating LUT... ")
	now = time.Now()
	// Extracts & EvalLUT(LWEs, indexLUT) on the fly -> Repack(LWEs, indexRepack) -> RLWE
	ctN12.Ciphertext = evalLUT.EvaluateAndRepack(ctN11.Ciphertext, lutPolyMap, repackIndex, LUTKEY)
	ctN12.Scale = paramsN12.DefaultScale()
	fmt.Printf("Done (%s)\n", time.Since(now))

	//for i, v := range encoderN12.DecodeCoeffs(decryptorN12.DecryptNew(ctN12)){
	//	fmt.Printf("%3d: %7.4f\n", i, v)
	//}

	fmt.Printf("Homomorphic Encoding... ")
	now = time.Now()
	// Homomorphic Encoding: [LUT(a), LUT(c), LUT(b), LUT(d)] -> [(LUT(a)+LUT(b)i), (LUT(c)+LUT(d)i)]
	ctN12, _ = evalCKKS.CoeffsToSlotsNew(ctN12, CoeffsToSlotsMatrix)
	fmt.Printf("Done (%s)\n", time.Since(now))

	for i, v := range encoderN12.Decode(decryptorN12.DecryptNew(ctN12), paramsN12.LogSlots()) {
		fmt.Printf("%7.4f -> %7.4f\n", values[i], v)
	}
}
