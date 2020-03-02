package main

import (
	"fmt"
	"math"
	"math/cmplx"
	"sort"
	"time"

	"github.com/ldsec/lattigo/ckks"
)

func example() {

	var start time.Time

	params := &ckks.Parameters{
		LogN:     14,
		LogSlots: 13,
		LogModuli: ckks.LogModuli{
			LogQi: []uint64{55, 40, 40, 40, 40, 40, 40, 40},
			LogPi: []uint64{45, 45},
		},
		Scale: 1 << 40,
		Sigma: 3.2,
	}

	fmt.Println()
	fmt.Println("=========================================")
	fmt.Println("         INSTANTIATING SCHEME            ")
	fmt.Println("=========================================")
	fmt.Println()

	start = time.Now()

	params.GenFromLogModuli()

	kgen := ckks.NewKeyGenerator(params)

	sk := kgen.GenSecretKey()

	rlk := kgen.GenRelinKey(sk)

	encryptor := ckks.NewEncryptorFromSk(params, sk)

	decryptor := ckks.NewDecryptor(params, sk)

	encoder := ckks.NewEncoder(params)

	evaluator := ckks.NewEvaluator(params)

	fmt.Printf("Done in %s \n", time.Since(start))

	fmt.Println()
	fmt.Printf("CKKS parameters : logN = %d, logSlots = %d, logQP = %d, levels = %d, scale= %f, sigma = %f \n", params.LogN, params.LogSlots, params.LogQP(), params.MaxLevel()+1, params.Scale, params.Sigma)

	fmt.Println()
	fmt.Println("=========================================")
	fmt.Println("           PLAINTEXT CREATION            ")
	fmt.Println("=========================================")
	fmt.Println()

	start = time.Now()

	r := float64(16)

	pi := 3.141592653589793

	slots := uint64(1 << params.LogSlots)

	values := make([]complex128, slots)
	for i := range values {
		values[i] = complex(2*pi, 0)
	}

	plaintext := ckks.NewPlaintext(params, params.MaxLevel(), params.Scale/r)
	encoder.Encode(plaintext, values, slots)

	fmt.Printf("Done in %s \n", time.Since(start))

	fmt.Println()
	fmt.Println("=========================================")
	fmt.Println("              ENCRYPTION                 ")
	fmt.Println("=========================================")
	fmt.Println()

	start = time.Now()

	ciphertext := encryptor.EncryptNew(plaintext)

	fmt.Printf("Done in %s \n", time.Since(start))

	printDebug(ciphertext, values, decryptor, encoder)

	fmt.Println()
	fmt.Println("=========================================")
	fmt.Println("             EVALUATION OF i*x           ")
	fmt.Println("=========================================")
	fmt.Println()

	start = time.Now()

	evaluator.MultByi(ciphertext, ciphertext)

	fmt.Printf("Done in %s \n", time.Since(start))

	for i := range values {
		values[i] *= complex(0, 1)
	}

	printDebug(ciphertext, values, decryptor, encoder)

	fmt.Println()
	fmt.Println("=========================================")
	fmt.Println("            EVALUATION of x/r            ")
	fmt.Println("=========================================")
	fmt.Println()

	start = time.Now()

	ciphertext.MulScale(r)

	fmt.Printf("Done in %s \n", time.Since(start))

	for i := range values {
		values[i] /= complex(r, 0)
	}

	printDebug(ciphertext, values, decryptor, encoder)

	fmt.Println()
	fmt.Println("=========================================")
	fmt.Println("            EVALUATION of e^x            ")
	fmt.Println("=========================================")
	fmt.Println()

	start = time.Now()

	coeffs := []float64{1.0, 1.0, 1.0 / 2, 1.0 / 6, 1.0 / 24, 1.0 / 120, 1.0 / 720, 1.0 / 5040}

	ciphertext = evaluator.EvaluatePolyFast(ciphertext, coeffs, rlk)

	fmt.Printf("Done in %s \n", time.Since(start))

	for i := range values {
		values[i] = cmplx.Exp(values[i])
	}

	printDebug(ciphertext, values, decryptor, encoder)

	fmt.Println()
	fmt.Println("=========================================")
	fmt.Println("            EVALUATION of x^r            ")
	fmt.Println("=========================================")
	fmt.Println()

	start = time.Now()

	evaluator.Power(ciphertext, uint64(r), rlk, ciphertext)

	fmt.Printf("Done in %s \n", time.Since(start))

	for i := range values {
		values[i] = cmplx.Pow(values[i], complex(r, 0))
	}

	printDebug(ciphertext, values, decryptor, encoder)

	fmt.Println()
	fmt.Println("=========================================")
	fmt.Println("         DECRYPTION & DECODING           ")
	fmt.Println("=========================================")
	fmt.Println()

	start = time.Now()

	encoder.Decode(decryptor.DecryptNew(ciphertext), slots)

	fmt.Printf("Done in %s \n", time.Since(start))

	printDebug(ciphertext, values, decryptor, encoder)

}

func printDebug(ciphertext *ckks.Ciphertext, valuesWant []complex128, decryptor ckks.Decryptor, encoder ckks.Encoder) {

	slots := uint64(len(valuesWant))

	valuesTest := encoder.Decode(decryptor.DecryptNew(ciphertext), slots)

	fmt.Println()
	fmt.Printf("ValuesTest : %6.10f %6.10f %6.10f %6.10f...\n", valuesTest[0], valuesTest[1], valuesTest[2], valuesTest[3])
	fmt.Printf("ValuesWant : %6.10f %6.10f %6.10f %6.10f...\n", valuesWant[0], valuesWant[1], valuesWant[2], valuesWant[3])
	verifyVector(valuesWant, valuesTest)
}

func verifyVector(valuesWant, valuesTest []complex128) (err error) {

	var deltaReal, deltaImag float64
	var delta, minprec, maxprec, meanprec, medianprec complex128

	diff := make([]complex128, len(valuesWant))

	minprec = complex(0, 0)
	maxprec = complex(1, 1)

	meanprec = complex(0, 0)

	for i := range valuesWant {

		delta = valuesTest[i] - valuesWant[i]
		deltaReal = math.Abs(real(delta))
		deltaImag = math.Abs(imag(delta))

		diff[i] += complex(deltaReal, deltaImag)

		meanprec += diff[i]

		if deltaReal > real(minprec) {
			minprec = complex(deltaReal, imag(minprec))
		}

		if deltaImag > imag(minprec) {
			minprec = complex(real(minprec), deltaImag)
		}

		if deltaReal < real(maxprec) {
			maxprec = complex(deltaReal, imag(maxprec))
		}

		if deltaImag < imag(maxprec) {
			maxprec = complex(real(maxprec), deltaImag)
		}
	}

	meanprec /= complex(float64(len(valuesWant)), 0)
	medianprec = calcmedian(diff)

	fmt.Println()
	fmt.Printf("Minimum precision : (%.2f, %.2f) bits \n", math.Log2(1/real(minprec)), math.Log2(1/imag(minprec)))
	fmt.Printf("Maximum precision : (%.2f, %.2f) bits \n", math.Log2(1/real(maxprec)), math.Log2(1/imag(maxprec)))
	fmt.Printf("Mean    precision : (%.2f, %.2f) bits \n", math.Log2(1/real(meanprec)), math.Log2(1/imag(meanprec)))
	fmt.Printf("Median  precision : (%.2f, %.2f) bits \n", math.Log2(1/real(medianprec)), math.Log2(1/imag(medianprec)))
	fmt.Println()

	return nil
}

func calcmedian(values []complex128) (median complex128) {

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
	example()
}
