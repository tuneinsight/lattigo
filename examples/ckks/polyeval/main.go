package main

import (
	"fmt"
	"math"

	"github.com/tuneinsight/lattigo/v3/ckks"
	"github.com/tuneinsight/lattigo/v3/rlwe"
	"github.com/tuneinsight/lattigo/v3/utils"
)

func chebyshevinterpolation() {

	var err error

	// This example packs random 8192 float64 values in the range [-8, 8]
	// and approximates the function f = 1/(exp(-x) + 1) over the range [-8, 8]
	// for the even slots and the function g = f * (1-f) over the range [-8, 8]
	// for the odd slots.
	// The result is then parsed and compared to the expected result.

	// Scheme params are taken directly from the proposed defaults
	params, err := ckks.NewParametersFromLiteral(ckks.PN14QP438)
	if err != nil {
		panic(err)
	}

	encoder := ckks.NewEncoder(params)

	// Keys
	kgen := ckks.NewKeyGenerator(params)
	sk, pk := kgen.GenKeyPair()

	// Relinearization key
	rlk := kgen.GenRelinearizationKey(sk, 1)

	// Encryptor
	encryptor := ckks.NewEncryptor(params, pk)

	// Decryptor
	decryptor := ckks.NewDecryptor(params, sk)

	// Evaluator
	evaluator := ckks.NewEvaluator(params, rlwe.EvaluationKey{Rlk: rlk})

	// Values to encrypt
	values := make([]float64, params.Slots())
	for i := range values {
		values[i] = utils.RandFloat64(-8, 8)
	}

	fmt.Printf("CKKS parameters: logN = %d, logQ = %d, levels = %d, scale= %f, sigma = %f \n",
		params.LogN(), params.LogQP(), params.MaxLevel()+1, params.DefaultScale(), params.Sigma())

	fmt.Println()
	fmt.Printf("Values     : %6f %6f %6f %6f...\n",
		round(values[0]), round(values[1]), round(values[2]), round(values[3]))
	fmt.Println()

	// Plaintext creation and encoding process
	plaintext := encoder.EncodeNew(values, params.MaxLevel(), params.DefaultScale(), params.LogSlots())

	// Encryption process
	var ciphertext *ckks.Ciphertext
	ciphertext = encryptor.EncryptNew(plaintext)

	a, b := -8.0, 8.0
	deg := 63

	fmt.Printf("Evaluation of the function f(x) for even slots and g(x) for odd slots in the range [%0.2f, %0.2f] (degree of approximation: %d)\n", a, b, deg)

	// Evaluation process
	// We approximate f(x) in the range [-8, 8] with a Chebyshev interpolant of 33 coefficients (degree 32).
	approxF := ckks.Approximate(f, a, b, deg)
	approxG := ckks.Approximate(g, a, b, deg)

	// Map storing which polynomial has to be applied to which slot.
	slotsIndex := make(map[int][]int)

	idxF := make([]int, params.Slots()>>1)
	idxG := make([]int, params.Slots()>>1)
	for i := 0; i < params.Slots()>>1; i++ {
		idxF[i] = i * 2   // Index with all even slots
		idxG[i] = i*2 + 1 // Index with all odd slots
	}

	slotsIndex[0] = idxF // Assigns index of all even slots to poly[0] = f(x)
	slotsIndex[1] = idxG // Assigns index of all odd slots to poly[1] = g(x)

	// Change of variable
	evaluator.MultByConst(ciphertext, 2/(b-a), ciphertext)
	evaluator.AddConst(ciphertext, (-a-b)/(b-a), ciphertext)
	if err := evaluator.Rescale(ciphertext, params.DefaultScale(), ciphertext); err != nil {
		panic(err)
	}

	polyVec := ckks.PolynomialVector{
		Value:      []ckks.Polynomial{approxF, approxG},
		Encoder:    encoder,
		SlotsIndex: slotsIndex,
	}

	// We evaluate the interpolated Chebyshev interpolant on the ciphertext
	if ciphertext, err = evaluator.EvaluatePoly(ciphertext, polyVec, ciphertext.Scale); err != nil {
		panic(err)
	}

	fmt.Println("Done... Consumed levels:", params.MaxLevel()-ciphertext.Level())

	// Computation of the reference values
	for i := 0; i < params.Slots()>>1; i++ {
		values[i*2] = f(values[i*2])
		values[i*2+1] = g(values[i*2+1])
	}

	// Print results and comparison
	printDebug(params, ciphertext, values, decryptor, encoder)

}

func f(x float64) float64 {
	return 1 / (math.Exp(-x) + 1)
}

func g(x float64) float64 {
	return f(x) * (1 - f(x))
}

func round(x float64) float64 {
	return math.Round(x*100000000) / 100000000
}

func printDebug(params ckks.Parameters, ciphertext *ckks.Ciphertext, valuesWant []float64, decryptor ckks.Decryptor, encoder ckks.Encoder) (valuesTest []float64) {

	tmp := encoder.Decode(decryptor.DecryptNew(ciphertext), params.LogSlots())

	valuesTest = make([]float64, len(tmp))
	for i := range tmp {
		valuesTest[i] = real(tmp[i])
	}

	fmt.Println()
	fmt.Printf("Level: %d (logQ = %d)\n", ciphertext.Level(), params.LogQLvl(ciphertext.Level()))
	fmt.Printf("Scale: 2^%f\n", math.Log2(ciphertext.Scale()))
	fmt.Printf("ValuesTest: %6.10f %6.10f %6.10f %6.10f...\n", valuesTest[0], valuesTest[1], valuesTest[2], valuesTest[3])
	fmt.Printf("ValuesWant: %6.10f %6.10f %6.10f %6.10f...\n", valuesWant[0], valuesWant[1], valuesWant[2], valuesWant[3])
	fmt.Println()

	precStats := ckks.GetPrecisionStats(params, encoder, nil, valuesWant, valuesTest, params.LogSlots(), 0)

	fmt.Println(precStats.String())

	return
}

func main() {
	chebyshevinterpolation()
}
