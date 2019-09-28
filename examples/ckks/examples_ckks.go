package main

import (
	"fmt"
	"github.com/ldsec/lattigo/ckks"
	"log"
	"math"
	"math/cmplx"
	"math/rand"
)

func randomFloat(min, max float64) float64 {
	return min + rand.Float64()*(max-min)
}

func randomComplex(min, max float64) complex128 {
	return complex(randomFloat(min, max), randomFloat(min, max))
}

func chebyshevinterpolation() {
	var err error

	// Scheme params
	params := ckks.DefaultParams[14]

	// Context
	var ckkscontext *ckks.CkksContext
	if ckkscontext, err = ckks.NewCkksContext(params); err != nil {
		log.Fatal(err)
	}

	encoder := ckkscontext.NewEncoder()

	kgen := ckkscontext.NewKeyGenerator()

	// Keys
	var sk *ckks.SecretKey
	var pk *ckks.PublicKey
	sk, pk = kgen.NewKeyPair()

	// Relinearization key
	var rlk *ckks.EvaluationKey
	if rlk, err = kgen.NewRelinKey(sk, ckkscontext.Scale()); err != nil {
		log.Fatal(err)
	}

	// Encryptor
	var encryptor *ckks.Encryptor
	if encryptor, err = ckkscontext.NewEncryptorFromPk(pk); err != nil {
		log.Fatal(err)
	}

	// Decryptor
	var decryptor *ckks.Decryptor
	if decryptor, err = ckkscontext.NewDecryptor(sk); err != nil {
		log.Fatal(err)
	}

	//Evaluator
	var evaluator *ckks.Evaluator
	if evaluator = ckkscontext.NewEvaluator(); err != nil {
		log.Fatal(err)
	}

	// Values to encrypt
	values := make([]complex128, ckkscontext.Slots())
	for i := range values {
		values[i] = complex(randomFloat(-1, 1), randomFloat(-0.1, 0.1))
	}

	fmt.Printf("HEAAN parameters : logN = %d, logQ = %d, levels = %d, logScale = %d, sigma = %f \n", ckkscontext.LogN(), ckkscontext.LogQ(), ckkscontext.Levels(), ckkscontext.Scale(), ckkscontext.Sigma())

	fmt.Println()
	fmt.Printf("Values     : %6f %6f %6f %6f...\n", round(values[0]), round(values[1]), round(values[2]), round(values[3]))
	fmt.Println()

	// Plaintext creation and encoding process
	plaintext := ckkscontext.NewPlaintext(ckkscontext.Levels()-1, ckkscontext.Scale())
	if err = encoder.EncodeComplex(plaintext, values); err != nil {
		log.Fatal(err)
	}

	// Encryption process
	var ciphertext *ckks.Ciphertext
	if ciphertext, err = encryptor.EncryptNew(plaintext); err != nil {
		log.Fatal(err)
	}

	fmt.Println("Evaluation of the function cos(exp(2*pi*i*x)) in the range [-1, 1] (degree of approximation : 65)")
	// Evaluation process
	chebyapproximation := ckks.Approximate(f, -1, 1, 65)
	if ciphertext, err = evaluator.EvaluateCheby(ciphertext, chebyapproximation, rlk); err != nil {
		log.Fatal(err)
	}
	fmt.Println("Done... Consumed levels :", ckkscontext.Levels()-1-ciphertext.Level())

	// Decryption process
	if plaintext, err = decryptor.DecryptNew(ciphertext); err != nil {
		log.Fatal(err)
	}

	// Decoding process
	valuesTest := encoder.DecodeComplex(plaintext)

	// Computation of the reference values
	for i := range values {
		values[i] = f(values[i])
	}

	// Prints results to compare
	fmt.Println()
	fmt.Printf("ValuesTest : %6f %6f %6f %6f...\n", round(valuesTest[0]), round(valuesTest[1]), round(valuesTest[2]), round(valuesTest[3]))
	fmt.Printf("ValuesWant : %6f %6f %6f %6f...\n", round(values[0]), round(values[1]), round(values[2]), round(values[3]))

}

func f(x complex128) complex128 {
	return cmplx.Cos(cmplx.Exp(2 * 3.141592653589793 * complex(0, 1) * x)) // cos(exp(2*pi*i*x))
	//return cmplx.Exp(x)/(cmplx.Exp(x) + 1) // sigmoid function
}

func round(x complex128) complex128 {
	var factor float64
	factor = 100000000
	a := math.Round(real(x)*factor) / factor
	b := math.Round(imag(x)*factor) / factor
	return complex(a, b)
}

func main() {
	chebyshevinterpolation()
}
