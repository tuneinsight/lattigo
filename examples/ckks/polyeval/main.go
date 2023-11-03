package main

import (
	"fmt"
	"math"
	"math/big"

	"github.com/tuneinsight/lattigo/v4/ckks"
	"github.com/tuneinsight/lattigo/v4/he/float"
	"github.com/tuneinsight/lattigo/v4/rlwe"
	"github.com/tuneinsight/lattigo/v4/utils/bignum"
	"github.com/tuneinsight/lattigo/v4/utils/sampling"
)

func chebyshevinterpolation() {

	var err error

	// This example packs random 8192 float64 values in the range [-8, 8]
	// and approximates the function f = 1/(exp(-x) + 1) over the range [-8, 8]
	// for the even slots and the function g = f * (1-f) over the range [-8, 8]
	// for the odd slots.
	// The result is then parsed and compared to the expected result.

	// Scheme params are taken directly from the proposed defaults
	params, err := ckks.NewParametersFromLiteral(
		ckks.ParametersLiteral{
			LogN:            14,
			LogQ:            []int{55, 40, 40, 40, 40, 40, 40, 40},
			LogP:            []int{45, 45},
			LogDefaultScale: 40,
		})
	if err != nil {
		panic(err)
	}

	encoder := ckks.NewEncoder(params)

	// Keys
	kgen := ckks.NewKeyGenerator(params)
	sk, pk := kgen.GenKeyPairNew()

	// Encryptor
	encryptor := ckks.NewEncryptor(params, pk)

	// Decryptor
	decryptor := ckks.NewDecryptor(params, sk)

	// Evaluator with relinearization key
	evaluator := ckks.NewEvaluator(params, rlwe.NewMemEvaluationKeySet(kgen.GenRelinearizationKeyNew(sk)))

	// Values to encrypt
	slots := params.MaxSlots()
	values := make([]float64, slots)
	for i := range values {
		values[i] = sampling.RandFloat64(-8, 8)
	}

	fmt.Printf("CKKS parameters: logN = %d, logQ = %f, levels = %d, scale= %f, noise = %T %v \n",
		params.LogN(), params.LogQP(), params.MaxLevel()+1, params.DefaultScale().Float64(), params.Xe(), params.Xe())

	fmt.Println()
	fmt.Printf("Values     : %6f %6f %6f %6f...\n",
		round(values[0]), round(values[1]), round(values[2]), round(values[3]))
	fmt.Println()

	// Plaintext creation and encoding process
	plaintext := ckks.NewPlaintext(params, params.MaxLevel())
	if err := encoder.Encode(values, plaintext); err != nil {
		panic(err)
	}

	// Encryption process
	var ciphertext *rlwe.Ciphertext
	ciphertext, err = encryptor.EncryptNew(plaintext)
	if err != nil {
		panic(err)
	}

	a, b := -8.0, 8.0
	deg := 63

	fmt.Printf("Evaluation of the function f(x) for even slots and g(x) for odd slots in the range [%0.2f, %0.2f] (degree of approximation: %d)\n", a, b, deg)

	// Evaluation process
	// We approximate f(x) in the range [-8, 8] with a Chebyshev interpolant of 33 coefficients (degree 32).

	interval := bignum.Interval{
		Nodes: deg,
		A:     *new(big.Float).SetFloat64(a),
		B:     *new(big.Float).SetFloat64(b),
	}

	approxF := bignum.ChebyshevApproximation(f, interval)
	approxG := bignum.ChebyshevApproximation(g, interval)

	// Map storing which polynomial has to be applied to which slot.
	mapping := make(map[int][]int)

	idxF := make([]int, slots>>1)
	idxG := make([]int, slots>>1)
	for i := 0; i < slots>>1; i++ {
		idxF[i] = i * 2   // Index with all even slots
		idxG[i] = i*2 + 1 // Index with all odd slots
	}

	mapping[0] = idxF // Assigns index of all even slots to poly[0] = f(x)
	mapping[1] = idxG // Assigns index of all odd slots to poly[1] = g(x)

	// Change of variable
	if err := evaluator.Mul(ciphertext, 2/(b-a), ciphertext); err != nil {
		panic(err)
	}

	if err := evaluator.Add(ciphertext, (-a-b)/(b-a), ciphertext); err != nil {
		panic(err)
	}

	if err := evaluator.Rescale(ciphertext, ciphertext); err != nil {
		panic(err)
	}

	polyVec, err := float.NewPolynomialVector([]bignum.Polynomial{approxF, approxG}, mapping)
	if err != nil {
		panic(err)
	}

	polyEval := float.NewPolynomialEvaluator(params, evaluator)

	// We evaluate the interpolated Chebyshev interpolant on the ciphertext
	if ciphertext, err = polyEval.Evaluate(ciphertext, polyVec, ciphertext.Scale); err != nil {
		panic(err)
	}

	fmt.Println("Done... Consumed levels:", params.MaxLevel()-ciphertext.Level())

	// Computation of the reference values
	for i := 0; i < slots>>1; i++ {
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

func printDebug(params ckks.Parameters, ciphertext *rlwe.Ciphertext, valuesWant []float64, decryptor *rlwe.Decryptor, encoder *ckks.Encoder) (valuesTest []float64) {

	valuesTest = make([]float64, 1<<ciphertext.LogDimensions.Cols)

	if err := encoder.Decode(decryptor.DecryptNew(ciphertext), valuesTest); err != nil {
		panic(err)
	}

	fmt.Println()
	fmt.Printf("Level: %d (logQ = %d)\n", ciphertext.Level(), params.LogQLvl(ciphertext.Level()))
	fmt.Printf("Scale: 2^%f\n", math.Log2(ciphertext.Scale.Float64()))
	fmt.Printf("ValuesTest: %6.10f %6.10f %6.10f %6.10f...\n", valuesTest[0], valuesTest[1], valuesTest[2], valuesTest[3])
	fmt.Printf("ValuesWant: %6.10f %6.10f %6.10f %6.10f...\n", valuesWant[0], valuesWant[1], valuesWant[2], valuesWant[3])
	fmt.Println()

	precStats := ckks.GetPrecisionStats(params, encoder, nil, valuesWant, valuesTest, 0, false)

	fmt.Println(precStats.String())

	return
}

func main() {
	chebyshevinterpolation()
}
