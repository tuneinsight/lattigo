package ckks

import (
	"fmt"
	"math"
	"sort"
)

// NoiseEstimator is a struct storing the necessary pre-computed
// to compute metrics on the plaintexts precision and error.
type NoiseEstimator struct {
	params *Parameters

	values      []complex128
	valuesfloat []float64

	m        uint64
	roots    []complex128
	rotGroup []uint64
}

// NewNoiseEstimator creates a new NoiseEstimator.
func NewNoiseEstimator(params *Parameters) (estimator *NoiseEstimator) {

	estimator = new(NoiseEstimator)

	m := 2 * params.N()

	rotGroup := make([]uint64, m>>1)
	fivePows := uint64(1)
	for i := uint64(0); i < m>>2; i++ {
		rotGroup[i] = fivePows
		fivePows *= GaloisGen
		fivePows &= (m - 1)
	}

	var angle float64
	roots := make([]complex128, m+1)
	for i := uint64(0); i < m; i++ {
		angle = 2 * 3.141592653589793 * float64(i) / float64(m)

		roots[i] = complex(math.Cos(angle), math.Sin(angle))
	}

	roots[m] = roots[0]

	return &NoiseEstimator{
		params:      params.Copy(),
		values:      make([]complex128, m>>2),
		valuesfloat: make([]float64, m>>1),
		m:           m,
		roots:       roots,
		rotGroup:    rotGroup,
	}
}

// StandardDeviationSlotDomain returns the scaled standard deviation of the difference between two complex vectors in the slot domains.
func (ne *NoiseEstimator) StandardDeviationSlotDomain(valuesWant, valuesHave []complex128, scale float64) (std float64) {

	var err complex128
	for i := range valuesWant {
		err = valuesWant[i] - valuesHave[i]
		ne.valuesfloat[2*i] = real(err)
		ne.valuesfloat[2*i+1] = imag(err)
	}

	return StandardDeviation(ne.valuesfloat[:len(valuesWant)*2], scale)
}

// StandardDeviationCoefDomain returns the scaled standard deviation of the [coefficient domain] of the difference between two complex vectors in the [slot domains].
func (ne *NoiseEstimator) StandardDeviationCoefDomain(valuesWant, valuesHave []complex128, scale float64) (std float64) {

	for i := range valuesHave {
		ne.values[i] = (valuesWant[i] - valuesHave[i])
	}

	invfft(ne.values, uint64(len(valuesWant)), ne.m, ne.rotGroup, ne.roots)

	for i := range valuesWant {
		ne.valuesfloat[2*i] = real(ne.values[i])
		ne.valuesfloat[2*i+1] = imag(ne.values[i])
	}

	return StandardDeviation(ne.valuesfloat[:len(valuesWant)*2], scale)
}

// PrecisionStats is a struct storing statistic about the precision of a CKKS plaintext
type PrecisionStats struct {
	MaxDelta        complex128
	MinDelta        complex128
	MaxPrecision    complex128
	MinPrecision    complex128
	MeanDelta       complex128
	MeanPrecision   complex128
	MedianDelta     complex128
	MedianPrecision complex128
	STDFreq         float64
	STDTime         float64

	RealDist, ImagDist []struct {
		Prec  float64
		Count int
	}

	cdfResol int
}

func (prec PrecisionStats) String() string {
	return fmt.Sprintf("\nMIN Prec : (%.2f, %.2f) Log2 \n", real(prec.MinPrecision), imag(prec.MinPrecision)) +
		fmt.Sprintf("MAX Prec : (%.2f, %.2f) Log2 \n", real(prec.MaxPrecision), imag(prec.MaxPrecision)) +
		fmt.Sprintf("AVG Prec : (%.2f, %.2f) Log2 \n", real(prec.MeanPrecision), imag(prec.MeanPrecision)) +
		fmt.Sprintf("MED Prec : (%.2f, %.2f) Log2 \n", real(prec.MedianPrecision), imag(prec.MedianPrecision)) +
		fmt.Sprintf("Err stdF : %5.2f Log2 \n", math.Log2(prec.STDFreq)) +
		fmt.Sprintf("Err stdT : %5.2f Log2 \n", math.Log2(prec.STDTime))

}

// PrecisionStats generates a PrecisionStats struct for a given vector based on a reference vector.
func (ne *NoiseEstimator) PrecisionStats(valuesWant, valuesHave []complex128) (prec PrecisionStats) {

	if len(valuesWant) != len(valuesHave) {
		panic("len(valuesWant) != len(valuesHave)")
	}

	var deltaReal, deltaImag float64
	var delta complex128

	diff := make([]complex128, len(valuesWant))

	prec.MaxDelta = complex(0, 0)
	prec.MinDelta = complex(1, 1)

	prec.MeanDelta = complex(0, 0)

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

	for i := range valuesWant {

		delta = valuesHave[i] - valuesWant[i]
		deltaReal = math.Abs(real(delta))
		deltaImag = math.Abs(imag(delta))
		precReal[i] = math.Log2(1 / deltaReal)
		precImag[i] = math.Log2(1 / deltaImag)

		diff[i] += complex(deltaReal, deltaImag)

		prec.MeanDelta += diff[i]

		if deltaReal > real(prec.MaxDelta) {
			prec.MaxDelta = complex(deltaReal, imag(prec.MaxDelta))
		}

		if deltaImag > imag(prec.MaxDelta) {
			prec.MaxDelta = complex(real(prec.MaxDelta), deltaImag)
		}

		if deltaReal < real(prec.MinDelta) {
			prec.MinDelta = complex(deltaReal, imag(prec.MinDelta))
		}

		if deltaImag < imag(prec.MinDelta) {
			prec.MinDelta = complex(real(prec.MinDelta), deltaImag)
		}
	}

	prec.calcCDF(precReal, prec.RealDist)
	prec.calcCDF(precImag, prec.ImagDist)

	prec.MinPrecision = deltaToPrecision(prec.MaxDelta)
	prec.MaxPrecision = deltaToPrecision(prec.MinDelta)
	prec.MeanDelta /= complex(float64(len(valuesWant)), 0)
	prec.MeanPrecision = deltaToPrecision(prec.MeanDelta)
	prec.MedianDelta = calcmedian(diff)
	prec.MedianPrecision = deltaToPrecision(prec.MedianDelta)
	prec.STDFreq = ne.StandardDeviationSlotDomain(valuesWant, valuesHave, ne.params.Scale())
	prec.STDTime = ne.StandardDeviationCoefDomain(valuesWant, valuesHave, ne.params.Scale())

	return
}

func deltaToPrecision(c complex128) complex128 {
	return complex(math.Log2(1/real(c)), math.Log2(1/imag(c)))
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
