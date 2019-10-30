package main

import (
	"fmt"
	"log"
	"math"
	"math/cmplx"
	"math/rand"
	"sort"
	"time"

	"github.com/ldsec/lattigo/ckks"
)

// Random float in [min, max)
func randomFloat(min, max float64) float64 {
	return min + rand.Float64()*(max-min)
}

func chebyshevInterpolation() {
	// This example will pack random 8192 float64 values in the range [-8, 8] and approximate the function 1/(exp(-x) + 1) over the range [-8, 8].
	// The result is then parsed and compared to the expected result.

	rand.Seed(time.Now().UnixNano())

	// Scheme params
	params := ckks.DefaultParams[14]

	// Context
	ckksContext, err := ckks.NewCkksContext(params)
	if err != nil {
		log.Fatal(err)
	}

	encoder := ckksContext.NewEncoder()

	// Keys
	kgen := ckksContext.NewKeyGenerator()
	sk, pk := kgen.NewKeyPair()

	// Relinearization key
	var rlk *ckks.EvaluationKey
	if rlk, err = kgen.NewRelinKey(sk, ckksContext.Scale()); err != nil {
		log.Fatal(err)
	}

	// Encryptor
	encryptor, err := ckksContext.NewEncryptorFromPk(pk)
	if err != nil {
		log.Fatal(err)
	}

	// Decryptor
	decryptor, err := ckksContext.NewDecryptor(sk)
	if err != nil {
		log.Fatal(err)
	}

	// Evaluator
	evaluator := ckksContext.NewEvaluator() // Cannot return error?

	// Values to encrypt
	values := make([]complex128, ckksContext.Slots())
	for i := range values {
		values[i] = complex(randomFloat(-8, 8), 0)
	}

	fmt.Printf("HEAAN parameters : logN = %d, logQ = %d, levels = %d, logScale = %d, sigma = %f\n",
		ckksContext.LogN(), ckksContext.LogQ(), ckksContext.Levels(), ckksContext.Scale(), ckksContext.Sigma())

	fmt.Printf("\nValues     : %6f %6f %6f %6f...\n\n",
		round(values[0]), round(values[1]), round(values[2]), round(values[3]))

	// Plaintext creation
	plaintext := ckksContext.NewPlaintext(ckksContext.Levels()-1, ckksContext.Scale())

	// Encoding process
	err = encoder.EncodeComplex(plaintext, values)
	if err != nil {
		log.Fatal(err)
	}

	// Encryption process
	ciphertext, err := encryptor.EncryptNew(plaintext)
	if err != nil {
		log.Fatal(err)
	}

	fmt.Println("Evaluation of the function 1/(exp(-x)+1) in the range [-8, 8] (degree of approximation : 32)")

	// Evaluation process
	// We ask to approximate f(x) in the range [-8, 8] with a chebyshev polynomial of 33 coefficients (degree 32).
	chebyApproximation := ckks.Approximate(f, -8, 8, 33)

	// We evaluate the interpolated chebyshev polynomial on the ciphertext
	resultCiphertext, err := evaluator.EvaluateCheby(ciphertext, chebyApproximation, rlk)
	if err != nil {
		log.Fatal(err)
	}

	fmt.Println("Done... Consumed levels :", ckksContext.Levels()-1-resultCiphertext.Level())

	// Decryption process
	resultPlaintext := decryptor.DecryptNew(resultCiphertext)

	// Decoding process
	computedValues := encoder.DecodeComplex(resultPlaintext)

	// Compute the reference values
	expectedValues := make([]complex128, len(values))
	for i := range values {
		expectedValues[i] = f(values[i])
	}

	// Printing results and comparison
	fmt.Println()
	fmt.Printf("Computed values : %6f %6f %6f %6f...\n",
		round(computedValues[0]), round(computedValues[1]), round(computedValues[2]), round(computedValues[3]))
	fmt.Printf("Expected values : %6f %6f %6f %6f...\n",
		round(expectedValues[0]), round(expectedValues[1]), round(expectedValues[2]), round(expectedValues[3]))

	verifyVector(expectedValues, computedValues)
}

func f(x complex128) complex128 {
	return 1 / (cmplx.Exp(-x) + 1)
}

func round(x complex128) complex128 {
	var factor float64
	factor = 100000000
	a := math.Round(real(x)*factor) / factor
	b := math.Round(imag(x)*factor) / factor

	return complex(a, b)
}

func verifyVector(expectedValues, computedValues []complex128) {
	var deltaReal, deltaImag float64
	var minPrec, maxPrec, meanPrec, medianPrec complex128

	diff := make([]complex128, len(expectedValues))

	minPrec = complex(0, 0)
	maxPrec = complex(1, 1)

	meanPrec = complex(0, 0)

	for i := range expectedValues {
		// Test the ratio for big values (> 1) and difference for small values (< 1)

		deltaReal = math.Abs(real(computedValues[i]) - real(expectedValues[i]))
		deltaImag = math.Abs(imag(computedValues[i]) - imag(expectedValues[i]))

		diff[i] += complex(deltaReal, 0)
		diff[i] += complex(0, deltaImag)

		meanPrec += diff[i]

		if real(diff[i]) > real(minPrec) {
			minPrec = complex(real(diff[i]), 0)
		}

		if imag(diff[i]) > imag(minPrec) {
			minPrec = complex(real(minPrec), imag(diff[i]))
		}

		if real(diff[i]) < real(maxPrec) {
			maxPrec = complex(real(diff[i]), 0)
		}

		if imag(diff[i]) < imag(maxPrec) {
			maxPrec = complex(real(maxPrec), imag(diff[i]))
		}
	}

	meanPrec /= complex(float64(len(expectedValues)), 0)
	medianPrec = computeMedian(diff)

	fmt.Println()
	fmt.Printf("Minimum precision : (%.2f, %.2f) bits\n",
		math.Log2(1/real(minPrec)), math.Log2(1/imag(minPrec)))
	fmt.Printf("Maximum precision : (%.2f, %.2f) bits\n",
		math.Log2(1/real(maxPrec)), math.Log2(1/imag(maxPrec)))
	fmt.Printf("Mean    precision : (%.2f, %.2f) bits\n",
		math.Log2(1/real(meanPrec)), math.Log2(1/imag(meanPrec)))
	fmt.Printf("Median  precision : (%.2f, %.2f) bits\n",
		math.Log2(1/real(medianPrec)), math.Log2(1/imag(medianPrec)))
	fmt.Println()
}

func computeMedian(values []complex128) (median complex128) {
	tmp := make([]float64, len(values))

	for i := range values {
		tmp[i] = real(values[i])
	}

	sort.Float64s(tmp)

	for i := range values {
		values[i] = complex(tmp[i], imag(values[i]))
	}

	for i := range values {
		tmp[i] = imag(values[i])
	}

	sort.Float64s(tmp)

	for i := range values {
		values[i] = complex(real(values[i]), tmp[i])
	}

	index := len(values) / 2

	if len(values)&1 == 1 {
		return values[index]
	}

	return (values[index] + values[index+1]) / 2
}

func main() {
	chebyshevInterpolation()
}
