package main

import (
	"fmt"
	"github.com/ldsec/lattigo/v2/ckks"
	"github.com/ldsec/lattigo/v2/lwe"
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/rlwe"
	//"github.com/ldsec/lattigo/v2/utils"
	ckksAdvanced "github.com/ldsec/lattigo/v2/ckks/advanced"
	"math"
	"time"
)

func main() {
	LUT()
}

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
	} else {
		return x
	}
}

var Q = []uint64{0x80000000080001, 0x2000000e0001, 0x1fffffc20001}
var P = []uint64{0x4000000008a0001}

var ckksParamsN12 = ckks.ParametersLiteral{
	LogN:     10,
	LogSlots: 5,
	Q:        Q,
	P:        P,
	Scale: 1<<40,
	Sigma:    rlwe.DefaultSigma,
	RingType: rlwe.RingStandard,
}

var ckksParamsN10 = ckks.ParametersLiteral{
	LogN:     9,
	Q:        Q[:1],
	P:        P[:1],
	Sigma:    rlwe.DefaultSigma,
	RingType: rlwe.RingStandard,
}

var SlotsToCoeffsParameters = ckksAdvanced.EncodingMatrixLiteral{
	LinearTransformType: ckksAdvanced.SlotsToCoeffs,
	LevelStart:          2,     // starting level
	BSGSRatio:           4.0,  // ratio between n1/n2 for n1*n2 = slots
	BitReversed:         false, // bit-reversed input
	ScalingFactor: [][]float64{ // Decomposition level of the encoding matrix
		{0x2000000e0001}, // Scale of the second matriox
		{0x1fffffc20001}, // Scale of the first matrix
	},
}

var CoeffsToSlotsParameters = ckksAdvanced.EncodingMatrixLiteral{
	LinearTransformType: ckksAdvanced.CoeffsToSlots,
	LevelStart:          2,     // starting level
	BSGSRatio:           4.0,  // ratio between n1/n2 for n1*n2 = slots
	BitReversed:         false, // bit-reversed input
	ScalingFactor: [][]float64{ // Decomposition level of the encoding matrix
		{0x2000000e0001}, // Scale of the second matriox
		{0x1fffffc20001}, // Scale of the first matrix
	},
}

func LUT() {
	var err error
	var paramsN12, paramsN10 ckks.Parameters
	
	if paramsN12, err = ckks.NewParametersFromLiteral(ckksParamsN12); err != nil {
		panic(err)
	}

	if paramsN10, err = ckks.NewParametersFromLiteral(ckksParamsN10); err != nil {
		panic(err)
	}

	a, b := -1.0, 1.0

	fmt.Printf("Generating LUT... ")
	now := time.Now()
	LUTPoly := lwe.InitLUT(identity, paramsN12.Scale(), paramsN12.RingQ(), a, b)
	fmt.Printf("Done (%s)\n", time.Since(now))

	lutPolyMap := make(map[int]*ring.Poly)
	repackIndex := make(map[int]int)
	gap := paramsN10.N() / (2*paramsN12.Slots())
	for i := 0; i < paramsN10.N(); i += gap {
		lutPolyMap[i] = LUTPoly
		repackIndex[i] = i/gap
	}

	kgenN12 := ckks.NewKeyGenerator(paramsN12)
	skN12 := kgenN12.GenSecretKey()
	encoderN12 := ckks.NewEncoder(paramsN12)
	encryptorN12 := ckks.NewEncryptor(paramsN12, skN12)
	decryptorN12 := ckks.NewDecryptor(paramsN12, skN12)

	kgenN10 := ckks.NewKeyGenerator(paramsN10)
	skN10 := kgenN10.GenSecretKey()

	swkN12ToN10 := kgenN12.GenSwitchingKey(skN12, skN10)

	fmt.Printf("Gen SlotsToCoeffs Matrices... ")
	now = time.Now()
	diffScale := paramsN10.QiFloat64(0) / (4.0 * paramsN12.Scale())
	normalization := 2.0/(b-a)
	sf := math.Pow(normalization*diffScale, 0.5)
	SlotsToCoeffsMatrix := ckksAdvanced.NewHomomorphicEncodingMatrixFromLiteral(SlotsToCoeffsParameters, encoderN12, paramsN12.LogN(), paramsN12.LogSlots()+1, complex(sf, 0))
	CoeffsToSlotsMatrix := ckksAdvanced.NewHomomorphicEncodingMatrixFromLiteral(CoeffsToSlotsParameters, encoderN12, paramsN12.LogN(), paramsN12.LogSlots(), complex(math.Pow(1/float64(2*paramsN12.Slots()), 0.5), 0))
	fmt.Printf("Done (%s)\n", time.Since(now))

	// Rotation Keys
	rotations := []int{}
	for i := 1; i < paramsN12.N(); i <<= 1 {
		rotations = append(rotations, i)
	}

	rotations = append(rotations, SlotsToCoeffsParameters.Rotations(paramsN12.LogN(), paramsN12.LogSlots()+1)...)
	rotations = append(rotations, CoeffsToSlotsParameters.Rotations(paramsN12.LogN(), paramsN12.LogSlots())...)

	rotKey := kgenN12.GenRotationKeysForRotations(rotations, true, skN12)

	handler := lwe.NewHandler(paramsN12.Parameters, paramsN10.Parameters, rotKey)

	eval := ckksAdvanced.NewEvaluator(paramsN12, rlwe.EvaluationKey{Rlk: nil, Rtks: rotKey})

	fmt.Printf("Encrypting bits of skLWE in RGSW... ")
	now = time.Now()
	LUTKEY := handler.GenLUTKey(skN12, skN10)
	fmt.Printf("Done (%s)\n", time.Since(now))

	fmt.Printf("Generating Plaintext & Encrypting... ")
	now = time.Now()
	values := make([]float64, paramsN12.Slots())
	for i := 0; i < paramsN12.Slots(); i++ {
		values[i] = float64(i+1)/float64(paramsN12.Slots())
	}

	/*
	fmt.Println()
	ckks.SliceBitReverseInPlaceFloat64(values, paramsN12.Slots())
	for i := range values{
		fmt.Printf("%7.4f\n", values[i])
	}
	ckks.SliceBitReverseInPlaceFloat64(values, paramsN12.Slots())
	*/

	pt := ckks.NewPlaintext(paramsN12, paramsN12.MaxLevel(), paramsN12.Scale())
	encoderN12.Encode(pt, values, paramsN12.LogSlots())
	ctN12 := encryptorN12.EncryptNew(pt)
	fmt.Printf("Done (%s)\n", time.Since(now))


	fmt.Printf("Homomorphic Decoding... ")
	now = time.Now()
	//eval.AddConst(ctN12, 2*(-a - b), ctN12)
	ctN12 = eval.SlotsToCoeffsNew(ctN12, nil, SlotsToCoeffsMatrix)
	ctN12.Scale = paramsN10.QiFloat64(0) / 4.0
	eval.DropLevel(ctN12, ctN12.Level())
	ctTmp := eval.SwitchKeysNew(ctN12, swkN12ToN10)
	ctN10 := ckks.NewCiphertext(paramsN10, 1, paramsN10.MaxLevel(), ctTmp.Scale)
	paramsN12.RingQ().InvNTTLvl(ctTmp.Level(), ctTmp.Value[0], ctTmp.Value[0])
	paramsN12.RingQ().InvNTTLvl(ctTmp.Level(), ctTmp.Value[1], ctTmp.Value[1])
	rlwe.SwitchCiphertextRingDegree(ctTmp.El(), ctN10.El())
	paramsN10.RingQ().NTT(ctN10.Value[0], ctN10.Value[0])
	paramsN10.RingQ().NTT(ctN10.Value[1], ctN10.Value[1])
	fmt.Printf("Done (%s)\n", time.Since(now))

	fmt.Printf("Evaluating LUT... ")
	now = time.Now()
	ctN12.Ciphertext = handler.ExtractAndEvaluateLUTAndRepack(ctN10.Ciphertext, lutPolyMap, repackIndex, LUTKEY)
	ctN12.Scale = paramsN12.Scale()
	fmt.Printf("Done (%s)\n", time.Since(now))

	fmt.Println("Homomorphic Encoding... ")
	now = time.Now()
	ctN12, _ = eval.CoeffsToSlotsNew(ctN12, CoeffsToSlotsMatrix)
	fmt.Printf("Done (%s)\n", time.Since(now))

	//eval.Rotate(ctN12, -16, ctN12)

	v := encoderN12.Decode(decryptorN12.DecryptNew(ctN12), paramsN12.LogSlots())

	for i := range v{
		fmt.Printf("%7.4f -> %7.4f\n", values[i], v[i])
	}
}

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
		if i%slots==0{
			if c > Q>>1 {
				fmt.Printf("%8.4f, ", (float64(c)-float64(Q))/scale)
			} else {
				fmt.Printf("%8.4f, ", float64(c)/scale)
			}
		}
	}
	fmt.Printf("]\n")
}
