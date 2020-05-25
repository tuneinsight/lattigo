package ckks

import (
	"fmt"
	"math"
	"sort"
)

type PrecisionStats struct {
	Min, Max, Mean, Median complex128
}

func (prec PrecisionStats) String() string {
	return fmt.Sprintf("Minimum precision : (%.2f, %.2f) bits \n", math.Log2(1/real(prec.Min)), math.Log2(1/imag(prec.Min)))+
		fmt.Sprintf("Maximum precision : (%.2f, %.2f) bits \n", math.Log2(1/real(prec.Max)), math.Log2(1/imag(prec.Max)))+
		fmt.Sprintf("Mean    precision : (%.2f, %.2f) bits \n", math.Log2(1/real(prec.Mean)), math.Log2(1/imag(prec.Mean)))+
		fmt.Sprintf("Median  precision : (%.2f, %.2f) bits \n", math.Log2(1/real(prec.Max)), math.Log2(1/imag(prec.Median)))
}

func GetPrecisionStats(params *Parameters, encoder Encoder, decryptor Decryptor, valuesWant []complex128, element interface{}) (prec PrecisionStats) {
	var plaintextTest *Plaintext
	var valuesTest []complex128

	switch element.(type) {
	case *Ciphertext:
		plaintextTest = decryptor.DecryptNew(element.(*Ciphertext))
	case *Plaintext:
		plaintextTest = element.(*Plaintext)
	}

	valuesTest = encoder.Decode(plaintextTest, params.Slots)

	//fmt.Println(valuesTest[:4])
	//fmt.Println(valuesWant[:4])

	var deltaReal, deltaImag float64

	var delta complex128

	diff := make([]complex128, params.Slots)

	prec.Min = complex(0, 0)
	prec.Max = complex(1, 1)

	prec.Mean = complex(0, 0)

	distribReal := make(map[uint64]uint64)
	distribImag := make(map[uint64]uint64)

	distribPrec := float64(25)

	for i := range valuesWant {

		delta = valuesTest[i] - valuesWant[i]
		deltaReal = math.Abs(real(delta))
		deltaImag = math.Abs(imag(delta))

		diff[i] += complex(deltaReal, deltaImag)

		prec.Mean += diff[i]

		if deltaReal > real(prec.Min) {
			prec.Min = complex(deltaReal, imag(prec.Min))
		}

		if deltaImag > imag(prec.Min) {
			prec.Min = complex(real(prec.Min), deltaImag)
		}

		if deltaReal < real(prec.Max) {
			prec.Max = complex(deltaReal, imag(prec.Max))
		}

		if deltaImag < imag(prec.Max) {
			prec.Max = complex(real(prec.Max), deltaImag)
		}

		distribReal[uint64(math.Floor(distribPrec*math.Log2(1/deltaReal)))]++
		distribImag[uint64(math.Floor(distribPrec*math.Log2(1/deltaImag)))]++
	}

	prec.Mean /= complex(float64(params.Slots), 0)
	prec.Median = calcmedian(diff)
	return prec
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

	if index+1 == len(values) {
		return values[index]
	}

	return (values[index] + values[index+1]) / 2
}