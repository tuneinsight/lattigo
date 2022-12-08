package main

import (
	"fmt"
	"math"
	"math/cmplx"
	"time"

	"github.com/tuneinsight/lattigo/v4/ckks"
	"github.com/tuneinsight/lattigo/v4/rlwe"
	"github.com/tuneinsight/lattigo/v4/utils/bignum/polynomial"
)

func example() {

	var start time.Time
	var err error

	// Schemes parameters are created from scratch
	params, err := ckks.NewParametersFromLiteral(
		ckks.ParametersLiteral{
			LogN:     14,
			LogQ:     []int{55, 40, 40, 40, 40, 40, 40, 40},
			LogP:     []int{45, 45},
			LogSlots: 13,
			LogScale: 40,
		})
	if err != nil {
		panic(err)
	}

	fmt.Println()
	fmt.Println("=========================================")
	fmt.Println("         INSTANTIATING SCHEME            ")
	fmt.Println("=========================================")
	fmt.Println()

	start = time.Now()

	kgen := ckks.NewKeyGenerator(params)

	sk := kgen.GenSecretKeyNew()

	encryptor := ckks.NewEncryptor(params, sk)

	decryptor := ckks.NewDecryptor(params, sk)

	encoder := ckks.NewEncoder(params)

	evk := rlwe.NewEvaluationKeySet()
	evk.RelinearizationKey = kgen.GenRelinearizationKeyNew(sk)

	evaluator := ckks.NewEvaluator(params, evk)

	fmt.Printf("Done in %s \n", time.Since(start))

	fmt.Println()
	fmt.Printf("CKKS parameters: logN = %d, logSlots = %d, logQP = %f, levels = %d, scale= %f, noise = %T %v \n", params.LogN(), params.LogSlots(), params.LogQP(), params.MaxLevel()+1, params.DefaultScale().Float64(), params.Xe(), params.Xe())

	fmt.Println()
	fmt.Println("=========================================")
	fmt.Println("           PLAINTEXT CREATION            ")
	fmt.Println("=========================================")
	fmt.Println()

	start = time.Now()

	r := float64(16)

	pi := 3.141592653589793

	slots := params.Slots()

	values := make([]complex128, slots)
	for i := range values {
		values[i] = complex(2*pi, 0)
	}

	plaintext := ckks.NewPlaintext(params, params.MaxLevel())
	plaintext.Scale = plaintext.Scale.Div(rlwe.NewScale(r))
	encoder.Encode(values, plaintext, params.LogSlots())

	fmt.Printf("Done in %s \n", time.Since(start))

	fmt.Println()
	fmt.Println("=========================================")
	fmt.Println("              ENCRYPTION                 ")
	fmt.Println("=========================================")
	fmt.Println()

	start = time.Now()

	ciphertext := encryptor.EncryptNew(plaintext)

	fmt.Printf("Done in %s \n", time.Since(start))

	printDebug(params, ciphertext, values, decryptor, encoder)

	fmt.Println()
	fmt.Println("===============================================")
	fmt.Printf("        EVALUATION OF i*x on %d values\n", slots)
	fmt.Println("===============================================")
	fmt.Println()

	start = time.Now()

	evaluator.MultByConst(ciphertext, 1i, ciphertext)

	fmt.Printf("Done in %s \n", time.Since(start))

	for i := range values {
		values[i] *= 1i
	}

	printDebug(params, ciphertext, values, decryptor, encoder)

	fmt.Println()
	fmt.Println("===============================================")
	fmt.Printf("       EVALUATION of x/r on %d values\n", slots)
	fmt.Println("===============================================")
	fmt.Println()

	start = time.Now()

	ciphertext.Scale = ciphertext.Mul(rlwe.NewScale(r))

	fmt.Printf("Done in %s \n", time.Since(start))

	for i := range values {
		values[i] /= complex(r, 0)
	}

	printDebug(params, ciphertext, values, decryptor, encoder)

	fmt.Println()
	fmt.Println("===============================================")
	fmt.Printf("       EVALUATION of e^x on %d values\n", slots)
	fmt.Println("===============================================")
	fmt.Println()

	start = time.Now()

	coeffs := []complex128{
		1.0,
		1.0,
		1.0 / 2,
		1.0 / 6,
		1.0 / 24,
		1.0 / 120,
		1.0 / 720,
		1.0 / 5040,
	}

	poly := ckks.NewPoly(coeffs)

	if ciphertext, err = evaluator.EvaluatePoly(ciphertext, poly, ciphertext.Scale); err != nil {
		panic(err)
	}

	fmt.Printf("Done in %s \n", time.Since(start))

	for i := range values {
		values[i] = cmplx.Exp(values[i])
	}

	printDebug(params, ciphertext, values, decryptor, encoder)

	fmt.Println()
	fmt.Println("===============================================")
	fmt.Printf("       EVALUATION of x^r on %d values\n", slots)
	fmt.Println("===============================================")
	fmt.Println()

	start = time.Now()

	monomialBasis := ckks.NewPowerBasis(ciphertext, polynomial.Monomial)
	if err = monomialBasis.GenPower(int(r), false, params.DefaultScale(), evaluator); err != nil {
		panic(err)
	}
	ciphertext = monomialBasis.Value[int(r)]

	fmt.Printf("Done in %s \n", time.Since(start))

	for i := range values {
		values[i] = cmplx.Pow(values[i], complex(r, 0))
	}

	printDebug(params, ciphertext, values, decryptor, encoder)

	fmt.Println()
	fmt.Println("=========================================")
	fmt.Println("         DECRYPTION & DECODING           ")
	fmt.Println("=========================================")
	fmt.Println()

	start = time.Now()

	encoder.Decode(decryptor.DecryptNew(ciphertext), params.LogSlots())

	fmt.Printf("Done in %s \n", time.Since(start))

	printDebug(params, ciphertext, values, decryptor, encoder)

}

func printDebug(params ckks.Parameters, ciphertext *rlwe.Ciphertext, valuesWant []complex128, decryptor rlwe.Decryptor, encoder ckks.Encoder) (valuesTest []complex128) {

	valuesTest = encoder.Decode(decryptor.DecryptNew(ciphertext), params.LogSlots())

	fmt.Println()
	fmt.Printf("Level: %d (logQ = %d)\n", ciphertext.Level(), params.LogQLvl(ciphertext.Level()))
	fmt.Printf("Scale: 2^%f\n", math.Log2(ciphertext.Scale.Float64()))
	fmt.Printf("ValuesTest: %6.10f %6.10f %6.10f %6.10f...\n", valuesTest[0], valuesTest[1], valuesTest[2], valuesTest[3])
	fmt.Printf("ValuesWant: %6.10f %6.10f %6.10f %6.10f...\n", valuesWant[0], valuesWant[1], valuesWant[2], valuesWant[3])
	fmt.Println()

	precStats := ckks.GetPrecisionStats(params, encoder, nil, valuesWant, valuesTest, params.LogSlots(), nil)

	fmt.Println(precStats.String())

	return
}

func main() {
	example()
}
