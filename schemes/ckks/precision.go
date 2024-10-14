package ckks

import (
	"fmt"
	"math"
	"math/big"
	"sort"
	"testing"

	"github.com/stretchr/testify/require"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/utils"
	"github.com/tuneinsight/lattigo/v6/utils/bignum"
)

// PrecisionStats is a struct storing statistic about the precision of a CKKS plaintext
type PrecisionStats struct {
	MINLog2Prec Stats
	MAXLog2Prec Stats
	AVGLog2Prec Stats
	MEDLog2Prec Stats
	STDLog2Prec Stats

	MINLog2Err Stats
	MAXLog2Err Stats
	AVGLog2Err Stats
	MEDLog2Err Stats
	STDLog2Err Stats

	Log2Scale float64

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
│STD Prec │ %5.2f │ %5.2f │ %5.2f │
├─────────┼───────┼───────┼───────┤
│MIN Err  │ %5.2f │ %5.2f │ %5.2f │
│MAX Err  │ %5.2f │ %5.2f │ %5.2f │
│AVG Err  │ %5.2f │ %5.2f │ %5.2f │
│MED Err  │ %5.2f │ %5.2f │ %5.2f │
│STD Err  │ %5.2f │ %5.2f │ %5.2f │
└─────────┴───────┴───────┴───────┘
`,
		prec.MINLog2Prec.Real, prec.MINLog2Prec.Imag, prec.MINLog2Prec.L2,
		prec.MAXLog2Prec.Real, prec.MAXLog2Prec.Imag, prec.MAXLog2Prec.L2,
		prec.AVGLog2Prec.Real, prec.AVGLog2Prec.Imag, prec.AVGLog2Prec.L2,
		prec.MEDLog2Prec.Real, prec.MEDLog2Prec.Imag, prec.MEDLog2Prec.L2,
		prec.STDLog2Prec.Real, prec.STDLog2Prec.Imag, prec.STDLog2Prec.L2,
		prec.MINLog2Err.Real, prec.MINLog2Err.Imag, prec.MINLog2Err.L2,
		prec.MAXLog2Err.Real, prec.MAXLog2Err.Imag, prec.MAXLog2Err.L2,
		prec.AVGLog2Err.Real, prec.AVGLog2Err.Imag, prec.AVGLog2Err.L2,
		prec.MEDLog2Err.Real, prec.MEDLog2Err.Imag, prec.MEDLog2Err.L2,
		prec.STDLog2Err.Real, prec.STDLog2Err.Imag, prec.STDLog2Err.L2)
}

// GetPrecisionStats generates a [PrecisionStats] struct from the reference values and the decrypted values
// vWant.(type) must be either []complex128 or []float64
// element.(type) must be either *Plaintext, *Ciphertext, []complex128 or []float64. If not *Ciphertext, then decryptor can be nil.
func GetPrecisionStats(params Parameters, encoder *Encoder, decryptor *rlwe.Decryptor, want, have interface{}, logprec float64, computeDCF bool) (prec PrecisionStats) {
	return getPrecisionStats(params, encoder, decryptor, want, have, logprec, computeDCF)
}

func VerifyTestVectors(params Parameters, encoder *Encoder, decryptor *rlwe.Decryptor, valuesWant, valuesHave interface{}, log2MinPrec int, logprec float64, printPrecisionStats bool, t *testing.T) {

	precStats := GetPrecisionStats(params, encoder, decryptor, valuesWant, valuesHave, logprec, false)

	if printPrecisionStats {
		t.Log(precStats.String())
	}

	switch params.RingType() {
	case ring.Standard:
		log2MinPrec -= params.LogN() + 2 // Z[X]/(X^{N} + 1)
	case ring.ConjugateInvariant:
		log2MinPrec -= params.LogN() + 3 // Z[X + X^1]/(X^{2N} + 1)
	}
	if log2MinPrec < 0 {
		log2MinPrec = 0
	}

	require.GreaterOrEqual(t, precStats.AVGLog2Prec.Real, float64(log2MinPrec))
	require.GreaterOrEqual(t, precStats.AVGLog2Prec.Imag, float64(log2MinPrec))
}

func getPrecisionStats(params Parameters, encoder *Encoder, decryptor *rlwe.Decryptor, want, have interface{}, logprec float64, computeDCF bool) (prec PrecisionStats) {

	valuesWant, valuesHave := getRawVectors(params, encoder, decryptor, want, have, logprec)

	Log2Scale := params.DefaultScale().Log2()

	slots := float64(len(valuesWant))

	Log2PrecReal := make([]float64, len(valuesWant))
	Log2PrecImag := make([]float64, len(valuesWant))
	Log2PrecL2 := make([]float64, len(valuesWant))

	var log2PrecReal, log2PrecImag, log2PrecL2 float64
	var log2AVGPrecReal, log2AVGPrecImag, log2AVGPrecL2 float64
	var log2MAXPrecReal, log2MAXPrecImag, log2MAXPrecL2 float64
	var log2MINPrecReal, log2MINPrecImag, log2MINPrecL2 float64 = 1e10, 1e10, 1e10
	var log2STDPrecReal, log2STDPrecImag, log2STDPrecL2 float64

	errReal := new(big.Float)
	errImag := new(big.Float)
	errL2 := new(big.Float)
	tmp := new(big.Float)

	for i := range valuesWant {

		errReal.Sub(valuesHave[i][0], valuesWant[i][0])
		errImag.Sub(valuesHave[i][1], valuesWant[i][1])
		errL2.Mul(errReal, errReal)
		errL2.Add(errL2, tmp.Mul(errReal, errReal))

		log2PrecReal, _ = errReal.Float64()
		log2PrecImag, _ = errImag.Float64()
		log2PrecL2, _ = errL2.Float64()

		log2PrecL2 = math.Sqrt(log2PrecL2)

		log2PrecReal = -math.Log2(math.Abs(log2PrecReal))
		log2PrecImag = -math.Log2(math.Abs(log2PrecImag))
		log2PrecL2 = -math.Log2(log2PrecL2)

		if math.IsInf(log2PrecReal, 0) {
			log2PrecReal = Log2Scale
		}

		if math.IsInf(log2PrecImag, 0) {
			log2PrecImag = Log2Scale
		}

		if math.IsInf(log2PrecL2, 0) {
			log2PrecL2 = Log2Scale
		}

		Log2PrecReal[i] = log2PrecReal
		Log2PrecImag[i] = log2PrecImag
		Log2PrecL2[i] = log2PrecL2

		log2AVGPrecReal += log2PrecReal
		log2AVGPrecImag += log2PrecImag
		log2AVGPrecL2 += log2PrecL2

		log2MAXPrecReal = utils.Max(log2MAXPrecReal, log2PrecReal)
		log2MAXPrecImag = utils.Max(log2MAXPrecImag, log2PrecImag)
		log2MAXPrecL2 = utils.Max(log2MAXPrecL2, log2PrecL2)

		log2MINPrecReal = utils.Min(log2MINPrecReal, log2PrecReal)
		log2MINPrecImag = utils.Min(log2MINPrecImag, log2PrecImag)
		log2MINPrecL2 = utils.Min(log2MINPrecL2, log2PrecL2)
	}

	log2AVGPrecReal /= slots
	log2AVGPrecImag /= slots
	log2AVGPrecL2 /= slots

	for _, x := range Log2PrecReal {
		x -= log2AVGPrecReal
		log2STDPrecReal += x * x
	}

	for _, x := range Log2PrecImag {
		x -= log2AVGPrecImag
		log2STDPrecImag += x * x
	}

	for _, x := range Log2PrecL2 {
		x -= log2AVGPrecL2
		log2STDPrecL2 += x * x
	}

	prec.MINLog2Prec.Real = log2MINPrecReal
	prec.MINLog2Prec.Imag = log2MINPrecImag
	prec.MINLog2Prec.L2 = log2MINPrecL2

	prec.MAXLog2Prec.Real = utils.Min(log2MAXPrecReal, Log2Scale)
	prec.MAXLog2Prec.Imag = utils.Min(log2MAXPrecImag, Log2Scale)
	prec.MAXLog2Prec.L2 = utils.Min(log2MAXPrecL2, Log2Scale)

	prec.AVGLog2Prec.Real = log2AVGPrecReal
	prec.AVGLog2Prec.Imag = log2AVGPrecImag
	prec.AVGLog2Prec.L2 = log2AVGPrecL2

	prec.MEDLog2Prec.Real = calcmedian(Log2PrecReal)
	prec.MEDLog2Prec.Imag = calcmedian(Log2PrecImag)
	prec.MEDLog2Prec.L2 = calcmedian(Log2PrecL2)

	prec.STDLog2Prec.Real = math.Sqrt(log2STDPrecReal / (slots - 1))
	prec.STDLog2Prec.Imag = math.Sqrt(log2STDPrecImag / (slots - 1))
	prec.STDLog2Prec.L2 = math.Sqrt(log2STDPrecL2 / (slots - 1))

	prec.MAXLog2Err.Real = Log2Scale - prec.MINLog2Prec.Real
	prec.MAXLog2Err.Imag = Log2Scale - prec.MINLog2Prec.Imag
	prec.MAXLog2Err.L2 = Log2Scale - prec.MINLog2Prec.L2

	prec.MINLog2Err.Real = Log2Scale - prec.MAXLog2Prec.Real
	prec.MINLog2Err.Imag = Log2Scale - prec.MAXLog2Prec.Imag
	prec.MINLog2Err.L2 = Log2Scale - prec.MAXLog2Prec.L2

	prec.AVGLog2Err.Real = Log2Scale - prec.AVGLog2Prec.Real
	prec.AVGLog2Err.Imag = Log2Scale - prec.AVGLog2Prec.Imag
	prec.AVGLog2Err.L2 = Log2Scale - prec.AVGLog2Prec.L2

	prec.MEDLog2Err.Real = Log2Scale - prec.MEDLog2Prec.Real
	prec.MEDLog2Err.Imag = Log2Scale - prec.MEDLog2Prec.Imag
	prec.MEDLog2Err.L2 = Log2Scale - prec.MEDLog2Prec.L2

	prec.STDLog2Err.Real = prec.STDLog2Prec.Real
	prec.STDLog2Err.Imag = prec.STDLog2Prec.Imag
	prec.STDLog2Err.L2 = prec.STDLog2Prec.L2

	prec.Log2Scale = Log2Scale

	if computeDCF {

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

		prec.calcCDF(Log2PrecReal, prec.RealDist)
		prec.calcCDF(Log2PrecImag, prec.ImagDist)
		prec.calcCDF(Log2PrecL2, prec.L2Dist)
	}

	return prec
}

func getRawVectors(params Parameters, encoder *Encoder, decryptor *rlwe.Decryptor, want, have interface{}, logprec float64) (valuesWant, valuesHave []*bignum.Complex) {

	precision := encoder.Prec()

	switch want := want.(type) {
	case []complex128:
		valuesWant = make([]*bignum.Complex, len(want))
		for i := range want {
			valuesWant[i] = &bignum.Complex{
				new(big.Float).SetPrec(precision).SetFloat64(real(want[i])),
				new(big.Float).SetPrec(precision).SetFloat64(imag(want[i])),
			}
		}
	case []float64:
		valuesWant = make([]*bignum.Complex, len(want))
		for i := range want {
			valuesWant[i] = &bignum.Complex{
				new(big.Float).SetPrec(precision).SetFloat64(want[i]),
				new(big.Float).SetPrec(precision),
			}
		}
	case []*big.Float:
		valuesWant = make([]*bignum.Complex, len(want))
		for i := range want {
			valuesWant[i] = &bignum.Complex{
				want[i],
				new(big.Float).SetPrec(precision),
			}
		}
	case []*bignum.Complex:
		valuesWant = want

		for i := range valuesWant {
			if valuesWant[i] == nil {
				valuesWant[i] = &bignum.Complex{new(big.Float), new(big.Float)}
			}
		}
	}

	switch have := have.(type) {
	case *rlwe.Ciphertext:
		valuesHave = make([]*bignum.Complex, len(valuesWant))
		if err := encoder.DecodePublic(decryptor.DecryptNew(have), valuesHave, logprec); err != nil {
			// Sanity check, this error should never happen.
			panic(err)
		}
	case *rlwe.Plaintext:
		valuesHave = make([]*bignum.Complex, len(valuesWant))
		if err := encoder.DecodePublic(have, valuesHave, logprec); err != nil {
			// Sanity check, this error should never happen.
			panic(err)
		}
	case []complex128:
		valuesHave = make([]*bignum.Complex, len(have))
		for i := range have {
			valuesHave[i] = &bignum.Complex{
				new(big.Float).SetPrec(precision).SetFloat64(real(have[i])),
				new(big.Float).SetPrec(precision).SetFloat64(imag(have[i])),
			}
		}
	case []float64:
		valuesHave = make([]*bignum.Complex, len(have))
		for i := range have {
			valuesHave[i] = &bignum.Complex{
				new(big.Float).SetPrec(precision).SetFloat64(have[i]),
				new(big.Float).SetPrec(precision),
			}
		}
	case []*big.Float:
		valuesHave = make([]*bignum.Complex, len(have))
		for i := range have {
			valuesHave[i] = &bignum.Complex{
				have[i],
				new(big.Float).SetPrec(precision),
			}
		}
	case []*bignum.Complex:
		valuesHave = have
		for i := range valuesHave {
			if valuesHave[i] == nil {
				valuesHave[i] = &bignum.Complex{new(big.Float), new(big.Float)}
			}
		}
	}

	return
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

func calcmedian(values []float64) (median float64) {
	sort.Float64s(values)
	index := len(values) / 2
	if len(values)&1 == 1 || index+1 == len(values) {
		return values[index]
	}
	return (values[index-1] + values[index]) / 2
}
