package main

import (
	"fmt"
	"github.com/tuneinsight/lattigo/v3/ckks"
	ckksAdvanced "github.com/tuneinsight/lattigo/v3/ckks/advanced"
	"github.com/tuneinsight/lattigo/v3/lwe"
	"github.com/tuneinsight/lattigo/v3/ring"
	"github.com/tuneinsight/lattigo/v3/rlwe"
	"math"
	"time"
)

func main() {
	LUT()
}

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

func sigmoid(x float64) (y float64) {
	return 1.0 / (math.Exp(-x) + 1)
}

func identity(x float64) (y float64) {
	return x
}

func relu(x float64) (y float64) {
	if x < 0 {
		return 0
	}

	return x
}

// Q modulus Q
var Q = []uint64{0x80000000080001, 0x2000000e0001, 0x1fffffc20001}

// P modulus P
var P = []uint64{0x4000000008a0001}

// Starting RLWE params, size of these params
// determine the complexity of the LUT:
// each LUT takes N RGSW ciphertext-ciphetext mul.
var ckksParamsN12 = ckks.ParametersLiteral{
	LogN:         6,
	LogSlots:     4,
	Q:            Q,
	P:            P,
	DefaultScale: 1 << 40,
	Sigma:        rlwe.DefaultSigma,
	RingType:     ring.Standard,
}

// LUT RLWE params, N of these params determine
// the LUT poly and therefore precision.
var ckksParamsN10 = ckks.ParametersLiteral{
	LogN:     5,
	Q:        Q[:1],
	P:        P[:1],
	Sigma:    rlwe.DefaultSigma,
	RingType: ring.Standard,
}

// LUT example
func LUT() {
	var err error
	var paramsN12, paramsN10 ckks.Parameters
	if paramsN12, err = ckks.NewParametersFromLiteral(ckksParamsN12); err != nil {
		panic(err)
	}

	if paramsN10, err = ckks.NewParametersFromLiteral(ckksParamsN10); err != nil {
		panic(err)
	}

	// LUT interval
	a, b := -8.0, 8.0

	fmt.Printf("Generating LUT... ")
	now := time.Now()
	// Generate LUT, provide function, outputscale, ring and interval.
	LUTPoly := lwe.InitLUT(sign, paramsN12.DefaultScale(), paramsN12.RingQ(), a, b)
	fmt.Printf("Done (%s)\n", time.Since(now))

	lutPolyMap := make(map[int]*ring.Poly) // Which slot to evaluate on the LUT
	repackIndex := make(map[int]int)       // Where to repack slots after the LUT
	gapN10 := paramsN10.N() / (2 * paramsN12.Slots())
	gapN12 := paramsN12.N() / (2 * paramsN12.Slots())
	for i := 0; i < paramsN12.Slots(); i++ {
		lutPolyMap[i*gapN10] = LUTPoly
		repackIndex[i] = i * gapN12
	}

	kgenN12 := ckks.NewKeyGenerator(paramsN12)
	skN12 := kgenN12.GenSecretKey()
	encoderN12 := ckks.NewEncoder(paramsN12)
	encryptorN12 := ckks.NewEncryptor(paramsN12, skN12)
	decryptorN12 := ckks.NewDecryptor(paramsN12, skN12)

	kgenN10 := ckks.NewKeyGenerator(paramsN10)
	skN10 := kgenN10.GenSecretKey()
	//decryptorN10 := ckks.NewDecryptor(paramsN10, skN10)
	//encoderN10 := ckks.NewEncoder(paramsN10)

	// Switchingkey RLWEN12 -> RLWEN10
	swkN12ToN10 := kgenN12.GenSwitchingKey(skN12, skN10)

	fmt.Printf("Gen SlotsToCoeffs Matrices... ")
	now = time.Now()

	// Rescale inputs during Homomorphic Decoding by the normalization of the
	// LUT inputs and change of scale to ensure that upperbound on the homomorphic
	// decryption of LWE during the LUT evaluation X^{dec(lwe)} is smaller than N
	// to avoid negacyclic wrapping of X^{dec(lwe)}.
	diffScale := paramsN10.QiFloat64(0) / (4.0 * paramsN12.DefaultScale())
	normalization := 2.0 / (b - a)

	// SlotsToCoeffsParameters homomorphic encoding parameters
	var SlotsToCoeffsParameters = ckksAdvanced.EncodingMatrixLiteral{
		LogN:                paramsN12.LogN(),
		LogSlots:            paramsN12.LogSlots(),
		Scaling:             normalization * diffScale,
		LinearTransformType: ckksAdvanced.SlotsToCoeffs,
		RepackImag2Real:     false,
		LevelStart:          2,     // starting level
		BSGSRatio:           4.0,   // ratio between n1/n2 for n1*n2 = slots
		BitReversed:         false, // bit-reversed input
		ScalingFactor: [][]float64{ // Decomposition level of the encoding matrix
			{0x2000000e0001}, // Scale of the second matriox
			{0x1fffffc20001}, // Scale of the first matrix
		},
	}

	// CoeffsToSlotsParameters homomorphic decoding parameters
	var CoeffsToSlotsParameters = ckksAdvanced.EncodingMatrixLiteral{
		LinearTransformType: ckksAdvanced.CoeffsToSlots,
		RepackImag2Real:     false,
		LogN:                paramsN12.LogN(),
		LogSlots:            paramsN12.LogSlots(),
		Scaling:             1 / float64(paramsN12.Slots()),
		LevelStart:          2,     // starting level
		BSGSRatio:           4.0,   // ratio between n1/n2 for n1*n2 = slots
		BitReversed:         false, // bit-reversed input
		ScalingFactor: [][]float64{ // Decomposition level of the encoding matrix
			{0x2000000e0001}, // Scale of the second matriox
			{0x1fffffc20001}, // Scale of the first matrix
		},
	}

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

	for i := 0; i < 32; i++ {
		rotations = append(rotations, i)
	}

	rotKey := kgenN12.GenRotationKeysForRotations(rotations, true, skN12)

	// LUT handler
	handler := lwe.NewHandler(paramsN12.Parameters, paramsN10.Parameters, rotKey)

	eval := ckksAdvanced.NewEvaluator(paramsN12, rlwe.EvaluationKey{Rlk: nil, Rtks: rotKey})

	fmt.Printf("Encrypting bits of skLWE in RGSW... ")
	now = time.Now()
	LUTKEY := handler.GenLUTKey(skN12, skN10) // Generate RGSW(sk_i) for all coefficients of sk
	fmt.Printf("Done (%s)\n", time.Since(now))

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
	// Homomorphic decoding: [(a+bi), (c+di)] -> [a, c, b, d]
	ctN12 = eval.SlotsToCoeffsNew(ctN12, nil, SlotsToCoeffsMatrix)
	ctN12.Scale = paramsN10.QiFloat64(0) / 4.0
	eval.DropLevel(ctN12, ctN12.Level())            // drop to LUT level
	ctTmp := eval.SwitchKeysNew(ctN12, swkN12ToN10) // key-switch to LWE degree
	ctN10 := ckks.NewCiphertext(paramsN10, 1, paramsN10.MaxLevel(), ctTmp.Scale)
	rlwe.SwitchCiphertextRingDegreeNTT(ctTmp.Ciphertext, paramsN10.RingQ(), paramsN12.RingQ(), ctN10.Ciphertext)
	fmt.Printf("Done (%s)\n", time.Since(now))

	fmt.Printf("Evaluating LUT... ")
	now = time.Now()
	// Extracts & EvalLUT(LWEs, indexLUT) on the fly -> Repack(LWEs, indexRepack) -> RLWE
	ctN12.Ciphertext = handler.ExtractAndEvaluateLUTAndRepack(ctN10.Ciphertext, lutPolyMap, repackIndex, LUTKEY)
	ctN12.Scale = paramsN12.DefaultScale()
	fmt.Printf("Done (%s)\n", time.Since(now))

	fmt.Println("Homomorphic Encoding... ")
	now = time.Now()
	// [LUT(a), LUT(c), LUT(b), LUT(d)] -> [(LUT(a)+LUT(b)i), (LUT(c)+LUT(d)i)]
	ctN12, _ = eval.CoeffsToSlotsNew(ctN12, CoeffsToSlotsMatrix)
	fmt.Printf("Done (%s)\n", time.Since(now))

	v := encoderN12.Decode(decryptorN12.DecryptNew(ctN12), paramsN12.LogSlots())

	for i := range v {
		fmt.Printf("%7.4f -> %7.4f\n", values[i], v[i])
	}
}

// PrintPoly prints poly
func PrintPoly(pol *ring.Poly, scale float64, Q uint64) {
	fmt.Printf("[")
	for _, c := range pol.Coeffs[0][:1] {
		if c > Q>>1 {
			fmt.Printf("%8.4f, ", float64(int(c)-int(Q))/scale)
		} else {
			fmt.Printf("%8.4f, ", float64(int(c))/scale)
		}
	}
	fmt.Printf("]\n")
}

// DecryptAndPrint decrypts and prints the first N values.
func DecryptAndPrint(decryptor ckks.Decryptor, LogSlots int, ringQ *ring.Ring, ciphertext *ckks.Ciphertext, scale float64) {
	plaintext := decryptor.DecryptNew(ciphertext)
	ringQ.InvNTTLvl(ciphertext.Level(), plaintext.Value, plaintext.Value)

	v := make([]float64, 1<<LogSlots)

	gap := ringQ.N / (1 << LogSlots)
	for i, j := 0, 0; i < 1<<LogSlots; i, j = i+1, j+gap {
		if plaintext.Value.Coeffs[0][j] >= ringQ.Modulus[0]>>1 {
			v[i] = -float64(ringQ.Modulus[0] - plaintext.Value.Coeffs[0][j])
		} else {
			v[i] = float64(plaintext.Value.Coeffs[0][i])
		}

		v[i] /= scale
	}

	for i := 0; i < 1<<LogSlots; i++ {
		if i&15 == 0 {
			fmt.Printf("\n")
		}
		fmt.Printf("%7.4f ", v[i])

	}
	fmt.Printf("\n")
}

// DecryptAndCenter decrypts and prints
func DecryptAndCenter(n int, b, a, sk *ring.Poly, ringQ *ring.Ring, mForm bool, scale float64, slots int) {

	pt := ringQ.NewPolyLvl(0)
	ringQ.MulCoeffsMontgomeryLvl(0, a, sk, pt)
	ringQ.AddLvl(0, pt, b, pt)
	ringQ.InvNTTLvl(0, pt, pt)
	if mForm {
		ringQ.InvMFormLvl(0, pt, pt)
	}

	Q := ringQ.Modulus[0]
	fmt.Printf("[")
	for i, c := range pt.Coeffs[0][:n] {
		if i%slots == 0 {
			if c > Q>>1 {
				fmt.Printf("%10.6f, ", (float64(c)-float64(Q))/scale)
			} else {
				fmt.Printf("%10.6f, ", float64(c)/scale)
			}
		}
	}
	fmt.Printf("]\n")
}
