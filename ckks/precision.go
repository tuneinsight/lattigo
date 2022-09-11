package ckks

import (
	"fmt"
	"math"
	"sort"
)

// PrecisionStats is a struct storing statistic about the precision of a CKKS plaintext
type PrecisionStats struct {
	MaxDelta        Stats
	MinDelta        Stats
	MaxPrecision    Stats
	MinPrecision    Stats
	MeanDelta       Stats
	MeanPrecision   Stats
	MedianDelta     Stats
	MedianPrecision Stats
	STDFreq         float64
	STDTime         float64

	RealDist, ImagDist, L2Dist []struct {
		Prec  float64
		Count int
	}

	cdfResol int
}

// Stats is a struct storing the real, imaginary and L2 norm (modulus)
// about the precision of a complex value.
type Stats struct {
	Real, Imag, L2 float64
}

func (prec PrecisionStats) String() string {
	return fmt.Sprintf(`
┌─────────┬───────┬───────┬───────┐
│    Log2 │ REAL  │ IMAG  │ L2    │
├─────────┼───────┼───────┼───────┤
│MIN Prec │ %5.2f │ %5.2f │ %5.2f │
│MAX Prec │ %5.2f │ %5.2f │ %5.2f │
│AVG Prec │ %5.2f │ %5.2f │ %5.2f │
│MED Prec │ %5.2f │ %5.2f │ %5.2f │
└─────────┴───────┴───────┴───────┘
Err STD Slots  : %5.2f Log2
Err STD Coeffs : %5.2f Log2
`,
		prec.MinPrecision.Real, prec.MinPrecision.Imag, prec.MinPrecision.L2,
		prec.MaxPrecision.Real, prec.MaxPrecision.Imag, prec.MaxPrecision.L2,
		prec.MeanPrecision.Real, prec.MeanPrecision.Imag, prec.MeanPrecision.L2,
		prec.MedianPrecision.Real, prec.MedianPrecision.Imag, prec.MedianPrecision.L2,
		math.Log2(prec.STDFreq),
		math.Log2(prec.STDTime))
}

// GetPrecisionStats generates a PrecisionStats struct from the reference values and the decrypted values
// vWant.(type) must be either []complex128 or []float64
// element.(type) must be either *Plaintext, *Ciphertext, []complex128 or []float64. If not *Ciphertext, then decryptor can be nil.
func GetPrecisionStats(params Parameters, encoder Encoder, decryptor Decryptor, vWant, element interface{}, logSlots int, sigma float64) (prec PrecisionStats) {

	var valuesTest []complex128

	switch element := element.(type) {
	case *Ciphertext:
		valuesTest = encoder.DecodePublic(decryptor.DecryptNew(element), logSlots, sigma)
	case *Plaintext:
		valuesTest = encoder.DecodePublic(element, logSlots, sigma)
	case []complex128:
		valuesTest = element
	case []float64:
		valuesTest = make([]complex128, len(element))
		for i := range element {
			valuesTest[i] = complex(element[i], 0)
		}
	}

	var valuesWant []complex128
	switch element := vWant.(type) {
	case []complex128:
		valuesWant = element
	case []float64:
		valuesWant = make([]complex128, len(element))
		for i := range element {
			valuesWant[i] = complex(element[i], 0)
		}
	}

	var deltaReal, deltaImag, deltaL2 float64

	slots := len(valuesWant)

	diff := make([]Stats, slots)

	prec.MaxDelta = Stats{0, 0, 0}
	prec.MinDelta = Stats{1, 1, 1}
	prec.MeanDelta = Stats{0, 0, 0}

	prec.cdfResol = 500

	prec.RealDist = make([]struct {
		Prec  float64
		Count int
	}, prec.cdfResol)
	prec.ImagDist = make([]struct {
		Prec  float64
		Count int
	}, prec.cdfResol)
	prec.L2Dist = make([]struct {
		Prec  float64
		Count int
	}, prec.cdfResol)

	precReal := make([]float64, len(valuesWant))
	precImag := make([]float64, len(valuesWant))
	precL2 := make([]float64, len(valuesWant))

	for i := range valuesWant {

		deltaReal = math.Abs(real(valuesTest[i]) - real(valuesWant[i]))
		deltaImag = math.Abs(imag(valuesTest[i]) - imag(valuesWant[i]))
		deltaL2 = math.Sqrt(deltaReal*deltaReal + deltaImag*deltaImag)
		precReal[i] = math.Log2(1 / deltaReal)
		precImag[i] = math.Log2(1 / deltaImag)
		precL2[i] = math.Log2(1 / deltaL2)

		diff[i].Real = deltaReal
		diff[i].Imag = deltaImag
		diff[i].L2 = deltaL2

		prec.MeanDelta.Real += deltaReal
		prec.MeanDelta.Imag += deltaImag
		prec.MeanDelta.L2 += deltaL2

		if deltaReal > prec.MaxDelta.Real {
			prec.MaxDelta.Real = deltaReal
		}

		if deltaImag > prec.MaxDelta.Imag {
			prec.MaxDelta.Imag = deltaImag
		}

		if deltaL2 > prec.MaxDelta.L2 {
			prec.MaxDelta.L2 = deltaL2
		}

		if deltaReal < prec.MinDelta.Real {
			prec.MinDelta.Real = deltaReal
		}

		if deltaImag < prec.MinDelta.Imag {
			prec.MinDelta.Imag = deltaImag
		}

		if deltaL2 < prec.MinDelta.L2 {
			prec.MinDelta.L2 = deltaL2
		}
	}

	prec.calcCDF(precReal, prec.RealDist)
	prec.calcCDF(precImag, prec.ImagDist)
	prec.calcCDF(precL2, prec.L2Dist)

	prec.MinPrecision = deltaToPrecision(prec.MaxDelta)
	prec.MaxPrecision = deltaToPrecision(prec.MinDelta)
	prec.MeanDelta.Real /= float64(slots)
	prec.MeanDelta.Imag /= float64(slots)
	prec.MeanDelta.L2 /= float64(slots)
	prec.MeanPrecision = deltaToPrecision(prec.MeanDelta)
	prec.MedianDelta = calcmedian(diff)
	prec.MedianPrecision = deltaToPrecision(prec.MedianDelta)
	prec.STDFreq = encoder.GetErrSTDSlotDomain(valuesWant[:], valuesTest[:], params.DefaultScale())
	prec.STDTime = encoder.GetErrSTDCoeffDomain(valuesWant, valuesTest, params.DefaultScale())
	return prec
}

func deltaToPrecision(c Stats) Stats {
	return Stats{math.Log2(1 / c.Real), math.Log2(1 / c.Imag), math.Log2(1 / c.L2)}
}

func (prec *PrecisionStats) calcCDF(precs []float64, res []struct {
	Prec  float64
	Count int
}) {
	sortedPrecs := make([]float64, len(precs))
	copy(sortedPrecs, precs)
	sort.Float64s(sortedPrecs)
	minPrec := sortedPrecs[0]
	maxPrec := sortedPrecs[len(sortedPrecs)-1]
	for i := 0; i < prec.cdfResol; i++ {
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

func calcmedian(values []Stats) (median Stats) {

	tmp := make([]float64, len(values))

	for i := range values {
		tmp[i] = values[i].Real
	}

	sort.Float64s(tmp)

	for i := range values {
		values[i].Real = tmp[i]
	}

	for i := range values {
		tmp[i] = values[i].Imag
	}

	sort.Float64s(tmp)

	for i := range values {
		values[i].Imag = tmp[i]
	}

	for i := range values {
		tmp[i] = values[i].L2
	}

	sort.Float64s(tmp)

	for i := range values {
		values[i].L2 = tmp[i]
	}

	index := len(values) / 2

	if len(values)&1 == 1 || index+1 == len(values) {
		return Stats{values[index].Real, values[index].Imag, values[index].L2}
	}

	return Stats{(values[index].Real + values[index+1].Real) / 2,
		(values[index].Imag + values[index+1].Imag) / 2,
		(values[index].L2 + values[index+1].L2) / 2}
}
