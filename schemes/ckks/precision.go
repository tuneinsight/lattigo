package ckks

import (
	"fmt"
	"math"
	"math/big"
	"sort"
	"testing"

	"github.com/stretchr/testify/require"
	"github.com/tuneinsight/lattigo/v5/core/rlwe"
	"github.com/tuneinsight/lattigo/v5/ring"
	"github.com/tuneinsight/lattigo/v5/utils/bignum"
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

	RealDist, ImagDist, L2Dist []struct {
		Prec  *big.Float
		Count int
	}

	cdfResol int
}

// Stats is a struct storing the real, imaginary and L2 norm (modulus)
// about the precision of a complex value.
type Stats struct {
	Real, Imag, L2 *big.Float
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
`,
		prec.MinPrecision.Real, prec.MinPrecision.Imag, prec.MinPrecision.L2,
		prec.MaxPrecision.Real, prec.MaxPrecision.Imag, prec.MaxPrecision.L2,
		prec.MeanPrecision.Real, prec.MeanPrecision.Imag, prec.MeanPrecision.L2,
		prec.MedianPrecision.Real, prec.MedianPrecision.Imag, prec.MedianPrecision.L2)
}

// GetPrecisionStats generates a PrecisionStats struct from the reference values and the decrypted values
// vWant.(type) must be either []complex128 or []float64
// element.(type) must be either *Plaintext, *Ciphertext, []complex128 or []float64. If not *Ciphertext, then decryptor can be nil.
func GetPrecisionStats(params Parameters, encoder *Encoder, decryptor *rlwe.Decryptor, want, have interface{}, logprec float64, computeDCF bool) (prec PrecisionStats) {

	if encoder.Prec() <= 53 {
		return getPrecisionStatsF64(params, encoder, decryptor, want, have, logprec, computeDCF)
	}

	return getPrecisionStatsF128(params, encoder, decryptor, want, have, logprec, computeDCF)
}

func VerifyTestVectors(params Parameters, encoder *Encoder, decryptor *rlwe.Decryptor, valuesWant, valuesHave interface{}, log2MinPrec int, logprec float64, printPrecisionStats bool, t *testing.T) {

	precStats := GetPrecisionStats(params, encoder, decryptor, valuesWant, valuesHave, logprec, false)

	if printPrecisionStats {
		t.Log(precStats.String())
	}

	rf64, _ := precStats.MeanPrecision.Real.Float64()
	if64, _ := precStats.MeanPrecision.Imag.Float64()

	switch params.RingType() {
	case ring.Standard:
		log2MinPrec -= params.LogN() + 2 // Z[X]/(X^{N} + 1)
	case ring.ConjugateInvariant:
		log2MinPrec -= params.LogN() + 3 // Z[X + X^1]/(X^{2N} + 1)
	}
	if log2MinPrec < 0 {
		log2MinPrec = 0
	}

	require.GreaterOrEqual(t, rf64, float64(log2MinPrec))
	require.GreaterOrEqual(t, if64, float64(log2MinPrec))
}

func getPrecisionStatsF64(params Parameters, encoder *Encoder, decryptor *rlwe.Decryptor, want, have interface{}, logprec float64, computeDCF bool) (prec PrecisionStats) {

	precision := encoder.Prec()

	var valuesWant []complex128
	switch want := want.(type) {
	case []complex128:
		valuesWant = make([]complex128, len(want))
		copy(valuesWant, want)
	case []float64:
		valuesWant = make([]complex128, len(want))
		for i := range want {
			valuesWant[i] = complex(want[i], 0)
		}
	case []*big.Float:
		valuesWant = make([]complex128, len(want))
		for i := range want {
			if want[i] != nil {
				f64, _ := want[i].Float64()
				valuesWant[i] = complex(f64, 0)
			}
		}
	case []*bignum.Complex:
		valuesWant = make([]complex128, len(want))
		for i := range want {
			if want[i] != nil {
				valuesWant[i] = want[i].Complex128()
			}

		}
	}

	var valuesHave = make([]complex128, len(valuesWant))

	switch have := have.(type) {
	case *rlwe.Ciphertext:
		if err := encoder.DecodePublic(decryptor.DecryptNew(have), valuesHave, logprec); err != nil {
			// Sanity check, this error should never happen.
			panic(err)
		}
	case *rlwe.Plaintext:
		if err := encoder.DecodePublic(have, valuesHave, logprec); err != nil {
			// Sanity check, this error should never happen.
			panic(err)
		}
	case []complex128:
		copy(valuesHave, have)
	case []float64:
		for i := range have {
			valuesHave[i] = complex(have[i], 0)
		}
	case []*big.Float:
		for i := range have {
			if have[i] != nil {
				f64, _ := have[i].Float64()
				valuesHave[i] = complex(f64, 0)
			}
		}
	case []*bignum.Complex:
		for i := range have {
			if have[i] != nil {
				valuesHave[i] = have[i].Complex128()
			}
		}
	}

	slots := len(valuesWant)

	diff := make([]struct{ Real, Imag, L2 float64 }, slots)

	precReal := make([]float64, len(valuesWant))
	precImag := make([]float64, len(valuesWant))
	precL2 := make([]float64, len(valuesWant))

	var deltaReal, deltaImag, deltaL2 float64
	var MeanDeltaReal, MeanDeltaImag, MeanDeltaL2 float64
	var MaxDeltaReal, MaxDeltaImag, MaxDeltaL2 float64
	var MinDeltaReal, MinDeltaImag, MinDeltaL2 float64 = 1, 1, 1

	for i := range valuesWant {

		deltaReal = math.Abs(real(valuesHave[i]) - real(valuesWant[i]))
		deltaImag = math.Abs(imag(valuesHave[i]) - imag(valuesWant[i]))
		deltaL2 = math.Sqrt(deltaReal*deltaReal + deltaReal*deltaReal)

		precReal[i] = -math.Log2(deltaReal)
		precImag[i] = -math.Log2(deltaImag)
		precL2[i] = -math.Log2(deltaL2)

		diff[i].Real = deltaReal
		diff[i].Imag = deltaImag
		diff[i].L2 = deltaL2

		MeanDeltaReal += deltaReal
		MeanDeltaImag += deltaImag
		MeanDeltaL2 += deltaL2

		if deltaReal > MaxDeltaReal {
			MaxDeltaReal = deltaReal
		}

		if deltaImag < MaxDeltaImag {
			MaxDeltaImag = deltaImag
		}

		if deltaL2 < MaxDeltaL2 {
			MaxDeltaL2 = deltaL2
		}

		if deltaReal < MinDeltaReal {
			MinDeltaReal = deltaReal
		}

		if deltaImag < MinDeltaImag {
			MinDeltaImag = deltaImag
		}

		if deltaL2 < MinDeltaL2 {
			MinDeltaL2 = deltaL2
		}
	}

	if computeDCF {

		prec.cdfResol = 500

		prec.RealDist = make([]struct {
			Prec  *big.Float
			Count int
		}, prec.cdfResol)
		prec.ImagDist = make([]struct {
			Prec  *big.Float
			Count int
		}, prec.cdfResol)
		prec.L2Dist = make([]struct {
			Prec  *big.Float
			Count int
		}, prec.cdfResol)

		prec.calcCDFF64(precReal, prec.RealDist)
		prec.calcCDFF64(precImag, prec.ImagDist)
		prec.calcCDFF64(precL2, prec.L2Dist)
	}

	prec.MinPrecision = deltaToPrecisionF64(struct{ Real, Imag, L2 float64 }{Real: MaxDeltaReal, Imag: MaxDeltaImag, L2: MaxDeltaL2})
	prec.MaxPrecision = deltaToPrecisionF64(struct{ Real, Imag, L2 float64 }{Real: MinDeltaReal, Imag: MinDeltaImag, L2: MinDeltaL2})
	prec.MeanDelta.Real = new(big.Float).SetFloat64(MeanDeltaReal / float64(slots))
	prec.MeanDelta.Imag = new(big.Float).SetFloat64(MeanDeltaImag / float64(slots))
	prec.MeanDelta.L2 = new(big.Float).SetFloat64(MeanDeltaL2 / float64(slots))
	prec.MeanPrecision = deltaToPrecisionF64(struct{ Real, Imag, L2 float64 }{Real: MeanDeltaReal / float64(slots), Imag: MeanDeltaImag / float64(slots), L2: MeanDeltaL2 / float64(slots)})
	prec.MedianDelta = calcmedianF64(diff)
	prec.MedianPrecision = deltaToPrecisionF128(prec.MedianDelta, bignum.Log(new(big.Float).SetPrec(precision).SetInt64(2)))
	return prec
}

func deltaToPrecisionF64(c struct{ Real, Imag, L2 float64 }) Stats {

	return Stats{
		new(big.Float).SetFloat64(-math.Log2(c.Real)),
		new(big.Float).SetFloat64(-math.Log2(c.Imag)),
		new(big.Float).SetFloat64(-math.Log2(c.L2)),
	}
}

func (prec *PrecisionStats) calcCDFF64(precs []float64, res []struct {
	Prec  *big.Float
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
				res[i].Prec = new(big.Float).SetFloat64(curPrec)
				res[i].Count = countSmaller
				break
			}
		}
	}
}

func calcmedianF64(values []struct{ Real, Imag, L2 float64 }) (median Stats) {

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
		return Stats{
			new(big.Float).SetFloat64(values[index].Real),
			new(big.Float).SetFloat64(values[index].Imag),
			new(big.Float).SetFloat64(values[index].L2),
		}
	}

	return Stats{
		new(big.Float).SetFloat64((values[index-1].Real + values[index].Real) / 2),
		new(big.Float).SetFloat64((values[index-1].Imag + values[index].Imag) / 2),
		new(big.Float).SetFloat64((values[index-1].L2 + values[index].L2) / 2),
	}
}

func getPrecisionStatsF128(params Parameters, encoder *Encoder, decryptor *rlwe.Decryptor, want, have interface{}, logprec float64, computeDCF bool) (prec PrecisionStats) {
	precision := encoder.Prec()

	var valuesWant []*bignum.Complex
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

	var valuesHave []*bignum.Complex

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

	slots := len(valuesWant)

	diff := make([]Stats, slots)

	prec.MaxDelta = Stats{
		new(big.Float).SetPrec(precision),
		new(big.Float).SetPrec(precision),
		new(big.Float).SetPrec(precision),
	}
	prec.MinDelta = Stats{
		new(big.Float).SetPrec(precision).SetInt64(1),
		new(big.Float).SetPrec(precision).SetInt64(1),
		new(big.Float).SetPrec(precision).SetInt64(1),
	}
	prec.MeanDelta = Stats{
		new(big.Float).SetPrec(precision),
		new(big.Float).SetPrec(precision),
		new(big.Float).SetPrec(precision),
	}

	precReal := make([]*big.Float, len(valuesWant))
	precImag := make([]*big.Float, len(valuesWant))
	precL2 := make([]*big.Float, len(valuesWant))

	deltaReal := new(big.Float)
	deltaImag := new(big.Float)
	deltaL2 := new(big.Float)

	tmp := new(big.Float)

	ln2 := bignum.Log(new(big.Float).SetPrec(precision).SetInt64(2))

	for i := range valuesWant {

		deltaReal.Sub(valuesHave[i][0], valuesWant[i][0])
		deltaReal.Abs(deltaReal)

		deltaImag.Sub(valuesHave[i][1], valuesWant[i][1])
		deltaImag.Abs(deltaImag)

		deltaL2.Mul(deltaReal, deltaReal)
		deltaL2.Add(deltaL2, tmp.Mul(deltaImag, deltaImag))
		deltaL2.Sqrt(deltaL2)

		precReal[i] = bignum.Log(deltaReal)
		precReal[i].Quo(precReal[i], ln2)
		precReal[i].Neg(precReal[i])

		precImag[i] = bignum.Log(deltaImag)
		precImag[i].Quo(precImag[i], ln2)
		precImag[i].Neg(precImag[i])

		precL2[i] = bignum.Log(deltaL2)
		precL2[i].Quo(precL2[i], ln2)
		precL2[i].Neg(precL2[i])

		diff[i].Real = new(big.Float).Set(deltaReal)
		diff[i].Imag = new(big.Float).Set(deltaImag)
		diff[i].L2 = new(big.Float).Set(deltaL2)

		prec.MeanDelta.Real.Add(prec.MeanDelta.Real, deltaReal)
		prec.MeanDelta.Imag.Add(prec.MeanDelta.Imag, deltaImag)
		prec.MeanDelta.L2.Add(prec.MeanDelta.L2, deltaL2)

		if deltaReal.Cmp(prec.MaxDelta.Real) == 1 {
			prec.MaxDelta.Real.Set(deltaReal)
		}

		if deltaImag.Cmp(prec.MaxDelta.Imag) == 1 {
			prec.MaxDelta.Imag.Set(deltaImag)
		}

		if deltaL2.Cmp(prec.MaxDelta.L2) == 1 {
			prec.MaxDelta.L2.Set(deltaL2)
		}

		if deltaReal.Cmp(prec.MinDelta.Real) == -1 {
			prec.MinDelta.Real.Set(deltaReal)
		}

		if deltaImag.Cmp(prec.MinDelta.Imag) == -1 {
			prec.MinDelta.Imag.Set(deltaImag)
		}

		if deltaL2.Cmp(prec.MinDelta.L2) == -1 {
			prec.MinDelta.L2.Set(deltaL2)
		}
	}

	if computeDCF {

		prec.cdfResol = 500

		prec.RealDist = make([]struct {
			Prec  *big.Float
			Count int
		}, prec.cdfResol)
		prec.ImagDist = make([]struct {
			Prec  *big.Float
			Count int
		}, prec.cdfResol)
		prec.L2Dist = make([]struct {
			Prec  *big.Float
			Count int
		}, prec.cdfResol)

		prec.calcCDFF128(precReal, prec.RealDist)
		prec.calcCDFF128(precImag, prec.ImagDist)
		prec.calcCDFF128(precL2, prec.L2Dist)
	}

	prec.MinPrecision = deltaToPrecisionF128(prec.MaxDelta, ln2)
	prec.MaxPrecision = deltaToPrecisionF128(prec.MinDelta, ln2)
	prec.MeanDelta.Real.Quo(prec.MeanDelta.Real, new(big.Float).SetPrec(precision).SetInt64(int64(slots)))
	prec.MeanDelta.Imag.Quo(prec.MeanDelta.Imag, new(big.Float).SetPrec(precision).SetInt64(int64(slots)))
	prec.MeanDelta.L2.Quo(prec.MeanDelta.L2, new(big.Float).SetPrec(precision).SetInt64(int64(slots)))
	prec.MeanPrecision = deltaToPrecisionF128(prec.MeanDelta, ln2)
	prec.MedianDelta = calcmedianF128(diff)
	prec.MedianPrecision = deltaToPrecisionF128(prec.MedianDelta, ln2)
	return prec
}

func deltaToPrecisionF128(c Stats, ln2 *big.Float) Stats {

	real := bignum.Log(c.Real)
	real.Quo(real, ln2)
	real.Neg(real)

	imag := bignum.Log(c.Imag)
	imag.Quo(imag, ln2)
	imag.Neg(imag)

	l2 := bignum.Log(c.L2)
	l2.Quo(l2, ln2)
	l2.Neg(l2)

	return Stats{
		real,
		imag,
		l2,
	}
}

func (prec *PrecisionStats) calcCDFF128(precs []*big.Float, res []struct {
	Prec  *big.Float
	Count int
}) {
	sortedPrecs := make([]*big.Float, len(precs))
	copy(sortedPrecs, precs)

	sort.Slice(sortedPrecs, func(i, j int) bool {
		return sortedPrecs[i].Cmp(sortedPrecs[j]) > 0
	})

	minPrec := sortedPrecs[0]
	maxPrec := sortedPrecs[len(sortedPrecs)-1]

	curPrec := new(big.Float)

	a := new(big.Float).Sub(maxPrec, minPrec)
	a.Quo(a, new(big.Float).SetInt64(int64(prec.cdfResol)))

	b := new(big.Float).Quo(minPrec, new(big.Float).SetInt64(int64(prec.cdfResol)))

	for i := 0; i < prec.cdfResol; i++ {

		curPrec.Mul(new(big.Float).SetInt64(int64(i)), a)
		curPrec.Add(curPrec, b)

		for countSmaller, p := range sortedPrecs {
			if p.Cmp(curPrec) >= 0 {
				res[i].Prec = new(big.Float).Set(curPrec)
				res[i].Count = countSmaller
				break
			}
		}
	}
}

func calcmedianF128(values []Stats) (median Stats) {

	tmp := make([]*big.Float, len(values))

	for i := range values {
		tmp[i] = values[i].Real
	}

	sort.Slice(tmp, func(i, j int) bool {
		return tmp[i].Cmp(tmp[j]) > 0
	})

	for i := range values {
		values[i].Real.Set(tmp[i])
	}

	for i := range values {
		tmp[i] = values[i].Imag
	}

	sort.Slice(tmp, func(i, j int) bool {
		return tmp[i].Cmp(tmp[j]) > 0
	})

	for i := range values {
		values[i].Imag = tmp[i]
	}

	for i := range values {
		tmp[i] = values[i].L2
	}

	sort.Slice(tmp, func(i, j int) bool {
		return tmp[i].Cmp(tmp[j]) > 0
	})

	for i := range values {
		values[i].L2 = tmp[i]
	}

	index := len(values) / 2

	if len(values)&1 == 1 || index+1 == len(values) {
		return Stats{values[index].Real, values[index].Imag, values[index].L2}
	}

	real := new(big.Float).Add(values[index-1].Real, values[index].Real)
	real.Quo(real, new(big.Float).SetInt64(2))

	imag := new(big.Float).Add(values[index-1].Imag, values[index].Imag)
	imag.Quo(imag, new(big.Float).SetInt64(2))

	l2 := new(big.Float).Add(values[index-1].L2, values[index].L2)
	l2.Quo(l2, new(big.Float).SetInt64(2))

	return Stats{
		real,
		imag,
		l2,
	}
}
