package main

import (
	"fmt"
	"github.com/tuneinsight/lattigo/v3/ckks"
	"github.com/tuneinsight/lattigo/v3/ring"
	"github.com/tuneinsight/lattigo/v3/rlwe"
	ckksAdvanced "github.com/tuneinsight/lattigo/v3/ckks/advanced"
)

func main() {
	LUT()
}

// Q modulus Q
var Q = []uint64{0x80000000080001, 0x2000000e0001, 0x1fffffc20001}

// P modulus P
var P = []uint64{0x4000000008a0001}

var ckksParamsN12 = ckks.ParametersLiteral{
	LogN:         5,
	LogSlots:     4,
	Q:            Q,
	P:            P,
	DefaultScale: 1 << 40,
	Sigma:        rlwe.DefaultSigma,
	RingType:     ring.Standard,
}

// LUT example
func LUT() {
	var err error
	var paramsN12 ckks.Parameters

	if paramsN12, err = ckks.NewParametersFromLiteral(ckksParamsN12); err != nil {
		panic(err)
	}

	kgenN12 := ckks.NewKeyGenerator(paramsN12)
	skN12 := kgenN12.GenSecretKey()
	encoderN12 := ckks.NewEncoder(paramsN12)
	encryptorN12 := ckks.NewEncryptor(paramsN12, skN12)
	decryptorN12 := ckks.NewDecryptor(paramsN12, skN12)

	fmt.Printf("Gen SlotsToCoeffs Matrices... ")

	// SlotsToCoeffsParameters homomorphic encoding parameters
	var SlotsToCoeffsParameters = ckksAdvanced.EncodingMatrixLiteral{
		LogN: paramsN12.LogN(),
		LogSlots: paramsN12.LogSlots(),
		Scaling: 1.0,
		LinearTransformType: ckksAdvanced.SlotsToCoeffs,
		RepackImag2Real: false,
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
		RepackImag2Real: false,
		LogN: paramsN12.LogN(),
		LogSlots: paramsN12.LogSlots(),
		Scaling: 1/16.0,
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

	// Rotation Keys
	rotations := []int{}
	for i := 1; i < paramsN12.N(); i <<= 1 {
		rotations = append(rotations, i)
	}

	rotations = append(rotations, SlotsToCoeffsParameters.Rotations()...)
	rotations = append(rotations, CoeffsToSlotsParameters.Rotations()...)

	for i := 0; i < 32; i++{
		rotations = append(rotations, i)
	}

	rotKey := kgenN12.GenRotationKeysForRotations(rotations, true, skN12)

	eval := ckksAdvanced.NewEvaluator(paramsN12, rlwe.EvaluationKey{Rlk: nil, Rtks: rotKey})

	values := make([]complex128, paramsN12.Slots())

	for i := 0; i < paramsN12.Slots(); i++ {
		values[i] = complex(float64(i+1), float64(i+1+paramsN12.Slots()))
	}

	fmt.Println("Before")
	ckks.SliceBitReverseInPlaceComplex128(values, paramsN12.Slots())
	for i := range values{
		fmt.Printf("%2d: %7.4f\n", i, values[i])
	}
	ckks.SliceBitReverseInPlaceComplex128(values, paramsN12.Slots())

	pt := ckks.NewPlaintext(paramsN12, paramsN12.MaxLevel(), paramsN12.DefaultScale())
	encoderN12.EncodeSlots(values, pt, paramsN12.LogSlots())
	ctN12 := encryptorN12.EncryptNew(pt)


	ctN12R := eval.SlotsToCoeffsNew(ctN12, nil, SlotsToCoeffsMatrix)

	valuesF := encoderN12.DecodeCoeffs(decryptorN12.DecryptNew(ctN12R))

	fmt.Println("After")
	for i, v := range valuesF{
		fmt.Printf("%2d: %7.4f\n", i, v)
	}

	encoderN12.EncodeCoeffs(valuesF, pt)
	ctN12 = encryptorN12.EncryptNew(pt)

	var ctN12I *ckks.Ciphertext
	ctN12R, ctN12I = eval.CoeffsToSlotsNew(ctN12, CoeffsToSlotsMatrix)

	if ctN12I != nil{
		eval.MultByi(ctN12I, ctN12I)
		eval.Add(ctN12R, ctN12I, ctN12R)
	}

	values = encoderN12.DecodeSlots(decryptorN12.DecryptNew(ctN12R), paramsN12.LogSlots())

	fmt.Println("After")
	for i, v := range values{
		fmt.Printf("%2d: %7.4f\n", i, v)
	}

}
