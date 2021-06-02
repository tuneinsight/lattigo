package main

import (
	"math"
	"sort"

	"github.com/ldsec/lattigo/v2/ckks"
	"github.com/ldsec/lattigo/v2/mkckks"
	"github.com/ldsec/lattigo/v2/mkrlwe"
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/utils"
)

// example of computation of a simple circuit with 2 parties
func main() {

	paramLit := ckks.DefaultParams[0]
	params, err := ckks.NewParametersFromLiteral(paramLit)

	if err != nil {
		panic("Couldn't retrieve default ckks parameters")
	}

	// generation of crs
	prng, err := utils.NewKeyedPRNG([]byte{'l', 'a', 't', 't', 'i', 'g', 'o'})
	if err != nil {
		panic(err)
	}

	crs := mkrlwe.GenCommonPublicParam(&params.Parameters, prng)

	// Participant 1
	keys1 := mkrlwe.KeyGen(&params.Parameters, crs)
	encryptor1 := mkckks.NewMKEncryptor(keys1.PublicKey, &params)
	encoder1 := ckks.NewEncoder(params)
	decryptor1 := mkrlwe.NewMKDecryptor(&params.Parameters, 0.6)

	value1 := newTestValue(&params, complex(-1, -1), complex(1, 1))
	plaintext1 := encoder1.EncodeNTTAtLvlNew(params.MaxLevel(), value1, params.LogSlots())

	cipher1 := encryptor1.Encrypt(plaintext1)
	evk1 := keys1.EvalKey
	pk1 := keys1.PublicKey

	// Participant 2
	keys2 := mkrlwe.KeyGen(&params.Parameters, crs)
	encryptor2 := mkckks.NewMKEncryptor(keys2.PublicKey, &params)
	encoder2 := ckks.NewEncoder(params)
	decryptor2 := mkrlwe.NewMKDecryptor(&params.Parameters, 0.6)

	value2 := newTestValue(&params, complex(-1, -1), complex(1, 1))
	plaintext2 := encoder2.EncodeNTTAtLvlNew(params.MaxLevel(), value2, params.LogSlots())

	cipher2 := encryptor2.Encrypt(plaintext2)
	evk2 := keys2.EvalKey
	pk2 := keys2.PublicKey

	// Evaluator: evaluates (c1 - c2) * (c1 + c2)
	evaluator := mkckks.NewMKEvaluator(&params)

	// decide on an indexing method for the participants and their public material and ciphertexts
	ids := []uint64{1, 2}
	evk1.PeerID = 1
	evk2.PeerID = 2
	pk1.PeerID = 1
	pk2.PeerID = 2
	evalKeys := []*mkrlwe.MKEvaluationKey{evk1, evk2}
	pubKeys := []*mkrlwe.MKPublicKey{pk1, pk2}

	// convert the ckks ciphertexts into multi key ciphertexts
	ciphers := evaluator.ConvertToMKCiphertext([]*ckks.Ciphertext{cipher1, cipher2}, ids)

	// evaluate circuit
	res1 := evaluator.Sub(ciphers[0], ciphers[1])
	res2 := evaluator.Add(ciphers[0], ciphers[1])
	res := evaluator.Mul(res1, res2)
	evaluator.RelinInPlace(res, evalKeys, pubKeys)

	// convert the multi key result into ckks ciphertexts for all participants
	resCKKS := evaluator.ConvertToCKKSCiphertext(res)

	// Partial Decryption done by each participants
	ckksCipher1 := resCKKS[0]
	ckksCipher2 := resCKKS[1]

	part1 := decryptor1.PartDec(&ckksCipher1.El().Element, ckksCipher1.Level(), keys1.SecretKey, 0.6)
	part2 := decryptor2.PartDec(&ckksCipher2.El().Element, ckksCipher2.Level(), keys2.SecretKey, 0.6)

	// Final decryption using the partial shares
	decrypted := decryptor1.MergeDec(&ckksCipher1.El().Element, ckksCipher1.Level(), []*ring.Poly{part1, part2})

	// decode
	pt := ckks.NewPlaintext(params, ckksCipher1.Level(), ckksCipher1.Scale())
	pt.SetValue(decrypted)

	finalValues := encoder1.Decode(pt, params.LogSlots())

	// perform the operation in plaintext space and check correctness
	expected := make([]complex128, len(value1))

	for i := range value1 {
		expected[i] = (value1[i] - value2[i]) * (value1[i] + value2[i])
	}

	verifyTestVectors(&params, expected, finalValues)

	println("Computation Successful")

}

// returs ramdom slices of complex numbers
func newTestValue(params *ckks.Parameters, a, b complex128) []complex128 {

	logSlots := params.LogSlots()

	values := make([]complex128, 1<<logSlots)

	for i := uint64(0); i < 1<<logSlots; i++ {
		values[i] = complex(utils.RandFloat64(real(a), real(b)), utils.RandFloat64(imag(a), imag(b)))
	}

	values[0] = complex(0.607538, 0)

	return values
}

func verifyTestVectors(params *ckks.Parameters, valuesWant, valuesTest []complex128) {

	precStats := GetPrecisionStats(params, valuesWant, valuesTest)

	if real(precStats.MeanPrecision) < 13.0 || imag(precStats.MeanPrecision) < 13.0 {
		panic("Circuit Error")
	}
}

// GetPrecisionStats generates a PrecisionStats struct from the reference values and the decrypted values
func GetPrecisionStats(params *ckks.Parameters, valuesWant, valuesTest []complex128) (prec ckks.PrecisionStats) {

	logSlots := params.LogSlots()
	slots := uint64(1 << logSlots)

	var deltaReal, deltaImag float64

	var delta complex128

	diff := make([]complex128, slots)

	prec.MaxDelta = complex(0, 0)
	prec.MinDelta = complex(1, 1)

	prec.MeanDelta = complex(0, 0)

	cdfResol := 500

	prec.RealDist = make([]struct {
		Prec  float64
		Count int
	}, cdfResol)
	prec.ImagDist = make([]struct {
		Prec  float64
		Count int
	}, cdfResol)

	precReal := make([]float64, len(valuesWant))
	precImag := make([]float64, len(valuesWant))

	for i := range valuesWant {

		delta = valuesTest[i] - valuesWant[i]
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

	calcCDF(precReal, cdfResol, prec.RealDist)
	calcCDF(precImag, cdfResol, prec.ImagDist)

	prec.MinPrecision = deltaToPrecision(prec.MaxDelta)
	prec.MaxPrecision = deltaToPrecision(prec.MinDelta)
	prec.MeanDelta /= complex(float64(slots), 0)
	prec.MeanPrecision = deltaToPrecision(prec.MeanDelta)
	prec.MedianDelta = calcmedian(diff)
	prec.MedianPrecision = deltaToPrecision(prec.MedianDelta)
	return prec
}

func deltaToPrecision(c complex128) complex128 {
	return complex(math.Log2(1/real(c)), math.Log2(1/imag(c)))
}

func calcCDF(precs []float64, cdfResol int, res []struct {
	Prec  float64
	Count int
}) {
	sortedPrecs := make([]float64, len(precs))
	copy(sortedPrecs, precs)
	sort.Float64s(sortedPrecs)
	minPrec := sortedPrecs[0]
	maxPrec := sortedPrecs[len(sortedPrecs)-1]
	for i := 0; i < cdfResol; i++ {
		curPrec := minPrec + float64(i)*(maxPrec-minPrec)/float64(cdfResol)
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
