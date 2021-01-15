package main

import (
	"fmt"
	"math"
	"time"

	"github.com/ldsec/lattigo/v2/ckks"
	//"github.com/ldsec/lattigo/v2/utils"
)

func main() {

	var err error

	// Scheme params
	LogN := uint64(14)
	LogSlots := uint64(13)

	LogModuli := ckks.LogModuli{
		LogQi: []uint64{55, 40, 40, 40},
		LogPi: []uint64{50, 50},
	}

	Scale := float64(1 << 40)

	params, err := ckks.NewParametersFromLogModuli(LogN, &LogModuli)
	if err != nil {
		panic(err)
	}
	params.SetScale(Scale)
	params.SetLogSlots(LogSlots)

	encoder := ckks.NewEncoder(params)

	// Keys
	kgen := ckks.NewKeyGenerator(params)
	sk, pk := kgen.GenKeyPair()

	// Relinearization key
	rlk := kgen.GenRelinKey(sk)

	// Encryptor
	encryptor := ckks.NewEncryptorFromPk(params, pk)

	// Decryptor
	decryptor := ckks.NewDecryptor(params, sk)

	// Evaluator
	eval := ckks.NewEvaluator(params)

	// Values to encrypt
	// We encrypt the 32x128 Matrix as 4 32x32 matrices in one ciphertext of 16384 slots.
	// We do the same for the 128x32 matrix.
	//  __________ 32x128 ___________
	// |                             |
	//  _____   _____   _____   _____     _____  _      	_____
	// |	 | |	 | |	 | |	 |   |	   |  |    	   |	 |
	// |32x32| |32x32| |32x32| |32x32| x |32x32|  |    =   |32x32|
	// |_____| |_____| |_____| |_____|   |_____|  |   	   |_____|
	//  	                              _____   |
	// 									 |	   |  |
	// 	                                 |32x32|  |
	// 	                                 |_____|  |
	//  	                              _____   128x32
	// 									 |	   |  |
	// 	                                 |32x32|  |
	// 	                                 |_____|  |
	//  	                              _____   |
	// 									 |	   |  |
	// 	                                 |32x32|  |
	// 	                                 |_____| _|
	//
	// We then perform the 32x128 x 128x32 matrix multiplication in two steps :
	//
	//
	// 1) The parallel multiplication between each 4 matrices (slot wise)
	//
	//        _____   _____   _____   _____
	//  	 |     | |     | |     | |     |
	// ct0 = |32x32| |32x32| |32x32| |32x32|
	//       |_____| |_____| |_____| |_____|
	//
	//  x       x       x       x       x
	//        _____   _____   _____   _____
	//  	 |     | |     | |     | |     |
	// ct1 = |32x32| |32x32| |32x32| |32x32|
	//       |_____| |_____| |_____| |_____|
	//
	//  =       =       =       =       =
	//        _____   _____   _____   _____
	// 		 |     | |     | |     | |     |
	// res = |32x32| |32x32| |32x32| |32x32|
	//       |_____| |_____| |_____| |_____|
	//
	//
	//
	// 2) The inner sum between the result (slot wise)
	//
	//        _____     _____     _____     _____
	// 		 |     |   |     |   |     |   |     |
	// res = |32x32| + |32x32| + |32x32| + |32x32|
	//       |_____|   |_____|   |_____|   |_____|

	// Matrices parameters
	var rows uint64 = 32
	var cols uint64 = 32
	var nb uint64 = 4

	dxd := int(rows * cols)

	// Rotation-keys generation
	mmpt := ckks.GenMatMulLinTrans(params, params.MaxLevel(), rows, encoder)
	rotKeys := ckks.NewRotationKeys()
	kgen.GenMatMulRotKeys(mmpt, sk, rotKeys)

	for i := uint64(1); i < nb; i <<= 1 {
		kgen.GenRotationKey(ckks.RotationLeft, sk, uint64(dxd)*i, rotKeys)
	}

	// Matrices generation
	//
	// 32x128 and 128x32 matrix generation (each split into 4 32x32 square matrices)

	m32x128 := ckks.GenRandomComplexMatrices(rows, cols, nb)
	m128x32 := ckks.GenRandomComplexMatrices(rows, cols, nb)

	values := make([]complex128, params.Slots())

	for i, matrix := range m32x128 {
		for j, c := range matrix.M {
			values[i*dxd+j] = c
		}
	}

	plaintextm32x128 := encoder.EncodeNew(values, params.LogSlots())

	for i, matrix := range m128x32 {
		for j, c := range matrix.M {
			values[i*dxd+j] = c
		}
	}

	plaintextm128x32 := encoder.EncodeNew(values, params.LogSlots())

	// Encryption process

	ciphertextm32x128 := encryptor.EncryptNew(plaintextm32x128)
	ciphertextm128x32 := encryptor.EncryptNew(plaintextm128x32)

	start := time.Now()
	ct32x32 := eval.MulMatrix(ciphertextm32x128, ciphertextm128x32, mmpt, rlk, rotKeys)

	// Inner-sum using a baby-step giant-step approach.
	tmp := ckks.NewCiphertext(params, ct32x32.Degree(), ct32x32.Level(), ct32x32.Scale())
	for i := uint64(1); i < nb; i <<= 1 {
		eval.Rotate(ct32x32, uint64(dxd)*i, rotKeys, tmp)
		eval.Add(ct32x32, tmp, ct32x32)
	}

	fmt.Println("Done :", time.Since(start))

	// Computes the reference plaintext by doing the computation in the plaintext domain.
	for j := range m32x128 {
		m32x128[j].MulMat(m32x128[j], m128x32[j])

		if j > 0 {
			m32x128[0].Add(m32x128[0], m32x128[j])
		}
	}

	fmt.Println()
	fmt.Printf("Level: %d (logQ = %d)\n", ct32x32.Level(), params.LogQLvl(ct32x32.Level()))
	fmt.Printf("Scale: 2^%f\n", math.Log2(ct32x32.Scale()))

	// We only care about the first 32x32 matrix, which is the result of the computation
	valuesWant := m32x128[0].M
	valuesHave := encoder.Decode(decryptor.DecryptNew(ct32x32), params.LogSlots())[:dxd]

	// Print results and comparison
	printDebug(rows, cols, valuesWant, valuesHave, params, encoder)

}

func printDebug(rows, cols uint64, valuesWant, valuesHave []complex128, params *ckks.Parameters, encoder ckks.Encoder) {

	precStats := ckks.GetPrecisionStats(params, encoder, nil, valuesWant, valuesHave, math.Exp2(53))

	fmt.Println(precStats.String())

	fmt.Println("Have")
	for i := uint64(0); i < rows; i++ {
		for j := uint64(0); j < cols; j++ {
			fmt.Printf("%8.4f ", real(valuesHave[i*rows+j]))
		}
		fmt.Printf("\n")
	}

	fmt.Println("Want")
	for i := uint64(0); i < rows; i++ {
		for j := uint64(0); j < cols; j++ {
			fmt.Printf("%8.4f ", real(valuesWant[i*rows+j]))
		}
		fmt.Printf("\n")
	}
	fmt.Printf("\n")

	fmt.Println()

	return
}
