package ckks

import (
	"fmt"
	"math"
	"sort"
)

type PrecisionStats struct {
	MinDelta, MaxDelta, Mean, Median complex128
	RealDist, ImagDist               []struct {
		Prec  float64
		Count int
	}

	cdfResol int
}

func (prec PrecisionStats) String() string {
	return fmt.Sprintf("\nMinimum precision : (%.2f, %.2f) bits \n", math.Log2(1/real(prec.MinDelta)), math.Log2(1/imag(prec.MinDelta))) +
		fmt.Sprintf("Maximum precision : (%.2f, %.2f) bits \n", math.Log2(1/real(prec.MaxDelta)), math.Log2(1/imag(prec.MaxDelta))) +
		fmt.Sprintf("Mean    precision : (%.2f, %.2f) bits \n", math.Log2(1/real(prec.Mean)), math.Log2(1/imag(prec.Mean))) +
		fmt.Sprintf("Median  precision : (%.2f, %.2f) bits \n", math.Log2(1/real(prec.Median)), math.Log2(1/imag(prec.Median)))
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

	prec.MinDelta = complex(0, 0)
	prec.MaxDelta = complex(1, 1)

	prec.Mean = complex(0, 0)

	prec.cdfResol = 500

	prec.RealDist = make([]struct {
		Prec  float64
		Count int
	}, prec.cdfResol)
	prec.ImagDist = make([]struct {
		Prec  float64
		Count int
	}, prec.cdfResol)

	precReal := make([]float64, len(valuesWant))
	precImag := make([]float64, len(valuesWant))

	//distribPrec := float64(25)

	for i := range valuesWant {

		delta = valuesTest[i] - valuesWant[i]
		deltaReal = math.Abs(real(delta))
		deltaImag = math.Abs(imag(delta))
		precReal[i] = math.Log2(1 / deltaReal)
		precImag[i] = math.Log2(1 / deltaImag)

		diff[i] += complex(deltaReal, deltaImag)

		prec.Mean += diff[i]

		if deltaReal > real(prec.MinDelta) { // TODO: Min is in fact the Min precision (so the MaxDelta...) to be clarified
			prec.MinDelta = complex(deltaReal, imag(prec.MinDelta))
		}

		if deltaImag > imag(prec.MinDelta) {
			prec.MinDelta = complex(real(prec.MinDelta), deltaImag)
		}

		if deltaReal < real(prec.MaxDelta) {
			prec.MaxDelta = complex(deltaReal, imag(prec.MaxDelta))
		}

		if deltaImag < imag(prec.MaxDelta) {
			prec.MaxDelta = complex(real(prec.MaxDelta), deltaImag)
		}
	}

	prec.calcCDF(precReal, prec.RealDist)
	prec.calcCDF(precImag, prec.ImagDist)

	prec.Mean /= complex(float64(params.Slots), 0)
	prec.Median = calcmedian(diff)
	return prec
}

func (prec *PrecisionStats) calcCDF(precs []float64, res []struct {
	Prec  float64
	Count int
}) {
	sortedPrecs := make([]float64, len(precs))
	copy(sortedPrecs, precs)
	sort.Float64s(sortedPrecs)
	minPrec := sortedPrecs[0]                  //math.Log2(1/real(prec.MinDelta))
	maxPrec := sortedPrecs[len(sortedPrecs)-1] //math.Log2(1/real(prec.MaxDelta))
	for i := 0; i < prec.cdfResol; i += 1 {
		curPrec := minPrec + float64(i)*(maxPrec-minPrec)/float64(prec.cdfResol)
		for countSmaller, p := range sortedPrecs {
			if p >= curPrec {
				res[i].Prec = curPrec
				res[i].Count = countSmaller
				break
			}
		}
	}
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
