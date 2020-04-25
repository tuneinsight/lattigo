package ckks

import (
	"fmt"
	"math"
	"math/cmplx"
	"math/rand"
	"sort"
	"testing"
	"time"

	"github.com/ldsec/lattigo/ring"
	"github.com/ldsec/lattigo/utils"
	"github.com/stretchr/testify/require"
)

func testString(opname string, params *Parameters) string {
	return fmt.Sprintf("%slogN=%d/LogSlots=%d/logQP=%d/levels=%d/a=%d/b=%d",
		opname,
		params.LogN,
		params.LogSlots,
		params.LogQP,
		params.MaxLevel+1,
		params.Alpha,
		params.Beta)
}

type ckksParams struct {
	params      *Parameters
	ckkscontext *Context
	encoder     Encoder
	kgen        KeyGenerator
	sk          *SecretKey
	pk          *PublicKey
	encryptorPk Encryptor
	encryptorSk Encryptor
	decryptor   Decryptor
	evaluator   Evaluator
}

type ckksTestParameters struct {
	verbose    bool
	medianprec float64
	slots      uint64

	ckksParameters []*Parameters
}

var err error
var testParams = new(ckksTestParameters)

func init() {
	rand.Seed(time.Now().UnixNano())

	testParams.medianprec = 15
	testParams.verbose = true

	testParams.ckksParameters = DefaultParams[PN12QP109 : PN12QP109+4]
}

func TestCKKS(t *testing.T) {
	t.Run("Encoder", testEncoder)
	t.Run("Encryptor", testEncryptor)
	t.Run("Evaluator/Add", testEvaluatorAdd)
	t.Run("Evaluator/Sub", testEvaluatorSub)
	t.Run("Evaluator/Rescale", testEvaluatorRescale)
	t.Run("Evaluator/AddConst", testEvaluatorAddConst)
	t.Run("Evaluator/MultByConst", testEvaluatorMultByConst)
	t.Run("Evaluator/MultByConstAndAdd", testEvaluatorMultByConstAndAdd)
	t.Run("Evaluator/Mul", testEvaluatorMul)
	t.Run("Evaluator/Functions", testFunctions)
	t.Run("Evaluator/EvaluatePoly", testEvaluatePoly)
	t.Run("Evaluator/ChebyshevInterpolator", testChebyshevInterpolator)
	t.Run("Evaluator/SwitchKeys", testSwitchKeys)
	t.Run("Evaluator/Conjugate", testConjugate)
	t.Run("Evaluator/RotateColumns", testRotateColumns)
	t.Run("Marshalling", testMarshaller)
}

func genCkksParams(contextParameters *Parameters) (params *ckksParams) {

	params = new(ckksParams)

	params.params = contextParameters.Copy()

	params.ckkscontext = newContext(contextParameters)

	params.kgen = NewKeyGenerator(contextParameters)

	params.sk, params.pk = params.kgen.GenKeyPairSparse(128)

	params.encoder = NewEncoder(contextParameters)

	params.encryptorPk = NewEncryptorFromPk(contextParameters, params.pk)
	params.encryptorSk = NewEncryptorFromSk(contextParameters, params.sk)
	params.decryptor = NewDecryptor(contextParameters, params.sk)

	params.evaluator = NewEvaluator(contextParameters)

	return

}

func newTestVectors(contextParams *ckksParams, encryptor Encryptor, a float64, t *testing.T) (values []complex128, plaintext *Plaintext, ciphertext *Ciphertext) {

	slots := uint64(1 << contextParams.params.LogSlots)

	values = make([]complex128, slots)

	for i := uint64(0); i < slots; i++ {
		values[i] = randomComplex(-a, a)
	}

	values[0] = complex(0.607538, 0.555668)

	plaintext = NewPlaintext(contextParams.params, contextParams.params.MaxLevel, contextParams.params.Scale)

	contextParams.encoder.EncodeNTT(plaintext, values, slots)

	if encryptor != nil {
		ciphertext = encryptor.EncryptNew(plaintext)
	}

	return values, plaintext, ciphertext
}

func newTestVectorsReals(contextParams *ckksParams, encryptor Encryptor, a, b float64, t *testing.T) (values []complex128, plaintext *Plaintext, ciphertext *Ciphertext) {

	slots := uint64(1 << contextParams.params.LogSlots)

	values = make([]complex128, slots)

	for i := uint64(0); i < slots; i++ {
		values[i] = complex(randomFloat(a, b), 0)
	}

	values[0] = complex(0.607538, 0)

	plaintext = NewPlaintext(contextParams.params, contextParams.params.MaxLevel, contextParams.params.Scale)

	contextParams.encoder.EncodeNTT(plaintext, values, slots)

	if encryptor != nil {
		ciphertext = encryptor.EncryptNew(plaintext)
	}

	return values, plaintext, ciphertext
}

func newTestVectorsSineBoot(contextParams *ckksParams, encryptor Encryptor, a, b float64, t *testing.T) (values []complex128, plaintext *Plaintext, ciphertext *Ciphertext) {

	slots := uint64(1 << contextParams.params.LogSlots)

	values = make([]complex128, slots)

	for i := uint64(0); i < slots; i++ {
		values[i] = complex(math.Round(randomFloat(a, b))+randomFloat(-1, 1)/1000, 0)
	}

	plaintext = NewPlaintext(contextParams.params, contextParams.params.MaxLevel, contextParams.params.Scale)

	contextParams.encoder.EncodeNTT(plaintext, values, slots)

	if encryptor != nil {
		ciphertext = encryptor.EncryptNew(plaintext)
	}

	return values, plaintext, ciphertext
}

func verifyTestVectors(contextParams *ckksParams, decryptor Decryptor, valuesWant []complex128, element interface{}, t *testing.T) {

	var plaintextTest *Plaintext
	var valuesTest []complex128

	switch element.(type) {
	case *Ciphertext:
		plaintextTest = decryptor.DecryptNew(element.(*Ciphertext))
	case *Plaintext:
		plaintextTest = element.(*Plaintext)
	}

	slots := uint64(1 << contextParams.params.LogSlots)

	valuesTest = contextParams.encoder.Decode(plaintextTest, slots)

	//fmt.Println(valuesTest[:4])
	//fmt.Println(valuesWant[:4])

	var deltaReal, deltaImag float64

	var delta, minprec, maxprec, meanprec, medianprec complex128

	diff := make([]complex128, slots)

	minprec = complex(0, 0)
	maxprec = complex(1, 1)

	meanprec = complex(0, 0)

	distribReal := make(map[uint64]uint64)
	distribImag := make(map[uint64]uint64)

	distribPrec := float64(25)

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

		distribReal[uint64(math.Floor(distribPrec*math.Log2(1/deltaReal)))]++
		distribImag[uint64(math.Floor(distribPrec*math.Log2(1/deltaImag)))]++
	}

	meanprec /= complex(float64(slots), 0)
	medianprec = calcmedian(diff)

	if testParams.verbose {
		t.Logf("Minimum precision : (%.2f, %.2f) bits \n", math.Log2(1/real(minprec)), math.Log2(1/imag(minprec)))
		t.Logf("Maximum precision : (%.2f, %.2f) bits \n", math.Log2(1/real(maxprec)), math.Log2(1/imag(maxprec)))
		t.Logf("Mean    precision : (%.2f, %.2f) bits \n", math.Log2(1/real(meanprec)), math.Log2(1/imag(meanprec)))
		t.Logf("Median  precision : (%.2f, %.2f) bits \n", math.Log2(1/real(medianprec)), math.Log2(1/imag(medianprec)))
		t.Log()
	}

	require.GreaterOrEqual(t, math.Log2(1/real(medianprec)), testParams.medianprec)
	require.GreaterOrEqual(t, math.Log2(1/imag(medianprec)), testParams.medianprec)
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

func testEncoder(t *testing.T) {

	for _, parameters := range testParams.ckksParameters {

		params := genCkksParams(parameters)

		t.Run(testString("Encode/", parameters), func(t *testing.T) {

			values, plaintext, _ := newTestVectors(params, nil, 1, t)

			verifyTestVectors(params, params.decryptor, values, plaintext, t)
		})

		t.Run(testString("EncodeCoeffs/", parameters), func(t *testing.T) {

			slots := parameters.N

			valuesWant := make([]float64, slots)

			for i := uint64(0); i < slots; i++ {
				valuesWant[i] = randomFloat(-1, 1)
			}

			valuesWant[0] = 0.607538

			plaintext := NewPlaintext(parameters, parameters.MaxLevel, parameters.Scale)

			params.encoder.EncodeCoeffs(valuesWant, plaintext)

			valuesTest := params.encoder.DecodeCoeffs(plaintext)

			var meanprec float64

			for i := range valuesWant {
				meanprec += math.Abs(valuesTest[i] - valuesWant[i])
			}

			meanprec /= float64(slots)

			require.GreaterOrEqual(t, math.Log2(1/meanprec), testParams.medianprec)
		})

		t.Run(testString("EncodeBigComplex/", parameters), func(t *testing.T) {

			prec := uint64(parameters.LogQP >> 1)

			encoder := NewEncoderBigComplex(parameters, prec)

			scale := math.Pow(2, float64(prec))

			values := make([]*ring.Complex, parameters.Slots)

			for i := range values {
				values[i] = ring.NewComplex(ring.NewFloat(float64(i+1), prec), ring.NewFloat(0, prec))
			}

			plaintext := NewPlaintext(parameters, parameters.MaxLevel, scale)

			encoder.Encode(plaintext, values, parameters.Slots)

			valuesTest := encoder.Decode(plaintext, parameters.Slots)

			error := ring.NewFloat(math.Pow(2, float64(prec-1)), prec)

			for i := range values {
				values[i].Sub(values[i], valuesTest[i])
				values[i].Real().Abs(values[i].Real())
				values[i].Imag().Abs(values[i].Imag())

				/*
					if i == 0 {
						fmt.Println(prec)
						fmt.Println(values[i].Float64())
					}
				*/

				require.Greater(t, 1, values[i].Real().Cmp(error))
				require.Greater(t, 1, values[i].Imag().Cmp(error))
			}
		})
	}
}

func testEncryptor(t *testing.T) {

	for _, parameters := range testParams.ckksParameters {

		params := genCkksParams(parameters)

		t.Run(testString("EncryptFromPk/", parameters), func(t *testing.T) {

			values, _, ciphertext := newTestVectors(params, params.encryptorPk, 1, t)

			verifyTestVectors(params, params.decryptor, values, ciphertext, t)
		})

		t.Run(testString("EncryptFromPkFast/", parameters), func(t *testing.T) {

			slots := uint64(1 << params.params.LogSlots)

			values := make([]complex128, slots)

			for i := uint64(0); i < slots; i++ {
				values[i] = randomComplex(-1, 1)
			}

			values[0] = complex(0.607538, 0.555668)

			plaintext := NewPlaintext(params.params, params.params.MaxLevel, params.params.Scale)

			params.encoder.Encode(plaintext, values, slots)

			verifyTestVectors(params, params.decryptor, values, params.encryptorPk.EncryptFastNew(plaintext), t)
		})

		t.Run(testString("EncryptFromSk/", parameters), func(t *testing.T) {

			values, _, ciphertext := newTestVectors(params, params.encryptorSk, 1, t)

			verifyTestVectors(params, params.decryptor, values, ciphertext, t)
		})

		t.Run(testString("EncryptFromSkFast/", parameters), func(t *testing.T) {

			slots := uint64(1 << params.params.LogSlots)

			values := make([]complex128, slots)

			for i := uint64(0); i < slots; i++ {
				values[i] = randomComplex(-1, 1)
			}

			values[0] = complex(0.607538, 0.555668)

			plaintext := NewPlaintext(params.params, params.params.MaxLevel, params.params.Scale)

			params.encoder.Encode(plaintext, values, slots)

			verifyTestVectors(params, params.decryptor, values, params.encryptorSk.EncryptFastNew(plaintext), t)
		})
	}
}

func testEvaluatorAdd(t *testing.T) {

	for _, parameters := range testParams.ckksParameters {

		params := genCkksParams(parameters)

		t.Run(testString("CtCtInPlace/", parameters), func(t *testing.T) {

			values1, _, ciphertext1 := newTestVectors(params, params.encryptorSk, 1, t)
			values2, _, ciphertext2 := newTestVectors(params, params.encryptorSk, 1, t)

			for i := range values1 {
				values1[i] += values2[i]
			}

			params.evaluator.Add(ciphertext1, ciphertext2, ciphertext1)

			verifyTestVectors(params, params.decryptor, values1, ciphertext1, t)
		})

		t.Run(testString("CtCtNew/", parameters), func(t *testing.T) {

			values1, _, ciphertext1 := newTestVectors(params, params.encryptorSk, 1, t)
			values2, _, ciphertext2 := newTestVectors(params, params.encryptorSk, 1, t)

			for i := range values1 {
				values1[i] += values2[i]
			}

			ciphertext3 := params.evaluator.AddNew(ciphertext1, ciphertext2)

			verifyTestVectors(params, params.decryptor, values1, ciphertext3, t)
		})

		t.Run(testString("CtPlainInPlace/", parameters), func(t *testing.T) {

			values1, _, ciphertext1 := newTestVectors(params, params.encryptorSk, 1, t)
			values2, plaintext2, _ := newTestVectors(params, params.encryptorSk, 1, t)

			for i := range values1 {
				values1[i] += values2[i]
			}

			params.evaluator.Add(ciphertext1, plaintext2, ciphertext1)

			verifyTestVectors(params, params.decryptor, values1, ciphertext1, t)

			for i := range values1 {
				values1[i] += values2[i]
			}

			params.evaluator.Add(plaintext2, ciphertext1, ciphertext1)

			verifyTestVectors(params, params.decryptor, values1, ciphertext1, t)
		})

		t.Run(testString("CtPlainInPlaceNew/", parameters), func(t *testing.T) {

			values1, _, ciphertext1 := newTestVectors(params, params.encryptorSk, 1, t)
			values2, plaintext2, _ := newTestVectors(params, params.encryptorSk, 1, t)

			for i := range values1 {
				values1[i] += values2[i]
			}

			ciphertext3 := params.evaluator.AddNew(ciphertext1, plaintext2)

			verifyTestVectors(params, params.decryptor, values1, ciphertext3, t)

			ciphertext3 = params.evaluator.AddNew(plaintext2, ciphertext1)

			verifyTestVectors(params, params.decryptor, values1, ciphertext3, t)
		})
	}
}

func testEvaluatorSub(t *testing.T) {

	for _, parameters := range testParams.ckksParameters {

		params := genCkksParams(parameters)

		t.Run(testString("CtCtInPlace/", parameters), func(t *testing.T) {

			values1, _, ciphertext1 := newTestVectors(params, params.encryptorSk, 1, t)
			values2, _, ciphertext2 := newTestVectors(params, params.encryptorSk, 1, t)

			for i := range values1 {
				values1[i] -= values2[i]
			}

			params.evaluator.Sub(ciphertext1, ciphertext2, ciphertext1)

			verifyTestVectors(params, params.decryptor, values1, ciphertext1, t)
		})

		t.Run(testString("CtCtNew/", parameters), func(t *testing.T) {

			values1, _, ciphertext1 := newTestVectors(params, params.encryptorSk, 1, t)
			values2, _, ciphertext2 := newTestVectors(params, params.encryptorSk, 1, t)

			for i := range values1 {
				values1[i] -= values2[i]
			}

			ciphertext3 := params.evaluator.SubNew(ciphertext1, ciphertext2)

			verifyTestVectors(params, params.decryptor, values1, ciphertext3, t)
		})

		t.Run(testString("CtPlainInPlace/", parameters), func(t *testing.T) {

			values1, _, ciphertext1 := newTestVectors(params, params.encryptorSk, 1, t)
			values2, plaintext2, ciphertext2 := newTestVectors(params, params.encryptorSk, 1, t)

			valuesTest := make([]complex128, len(values1))
			for i := range values1 {
				valuesTest[i] = values1[i] - values2[i]
			}

			params.evaluator.Sub(ciphertext1, plaintext2, ciphertext2)

			verifyTestVectors(params, params.decryptor, valuesTest, ciphertext2, t)

			for i := range values1 {
				valuesTest[i] = values2[i] - values1[i]
			}

			params.evaluator.Sub(plaintext2, ciphertext1, ciphertext2)

			verifyTestVectors(params, params.decryptor, valuesTest, ciphertext2, t)
		})

		t.Run(testString("CtPlainNew/", parameters), func(t *testing.T) {

			values1, _, ciphertext1 := newTestVectors(params, params.encryptorSk, 1, t)
			values2, plaintext2, _ := newTestVectors(params, params.encryptorSk, 1, t)

			valuesTest := make([]complex128, len(values1))
			for i := range values1 {
				valuesTest[i] = values1[i] - values2[i]
			}

			ciphertext3 := params.evaluator.SubNew(ciphertext1, plaintext2)

			verifyTestVectors(params, params.decryptor, valuesTest, ciphertext3, t)

			for i := range values1 {
				valuesTest[i] = values2[i] - values1[i]
			}

			ciphertext3 = params.evaluator.SubNew(plaintext2, ciphertext1)

			verifyTestVectors(params, params.decryptor, valuesTest, ciphertext3, t)
		})
	}
}

func testEvaluatorRescale(t *testing.T) {

	for _, parameters := range testParams.ckksParameters {

		params := genCkksParams(parameters)

		t.Run(testString("Single/", parameters), func(t *testing.T) {

			values, _, ciphertext := newTestVectors(params, params.encryptorSk, 1, t)

			constant := params.ckkscontext.contextQ.Modulus[ciphertext.Level()]

			params.evaluator.MultByConst(ciphertext, constant, ciphertext)

			ciphertext.MulScale(float64(constant))

			params.evaluator.Rescale(ciphertext, parameters.Scale, ciphertext)

			verifyTestVectors(params, params.decryptor, values, ciphertext, t)
		})

		t.Run(testString("Many/", parameters), func(t *testing.T) {

			values, _, ciphertext := newTestVectors(params, params.encryptorSk, 1, t)

			nbRescales := parameters.MaxLevel

			for i := uint64(0); i < nbRescales; i++ {
				constant := params.ckkscontext.contextQ.Modulus[ciphertext.Level()]
				params.evaluator.MultByConst(ciphertext, constant, ciphertext)
				ciphertext.MulScale(float64(constant))
			}

			params.evaluator.RescaleMany(ciphertext, nbRescales, ciphertext)

			verifyTestVectors(params, params.decryptor, values, ciphertext, t)
		})

	}
}

func testEvaluatorAddConst(t *testing.T) {

	for _, parameters := range testParams.ckksParameters {

		params := genCkksParams(parameters)

		t.Run(testString("", parameters), func(t *testing.T) {

			values, _, ciphertext := newTestVectors(params, params.encryptorSk, 1, t)

			constant := complex(3.1415, -1.4142)

			for i := range values {
				values[i] += constant
			}

			params.evaluator.AddConst(ciphertext, constant, ciphertext)

			verifyTestVectors(params, params.decryptor, values, ciphertext, t)
		})
	}
}

func testEvaluatorMultByConst(t *testing.T) {

	for _, parameters := range testParams.ckksParameters {

		params := genCkksParams(parameters)

		t.Run(testString("", parameters), func(t *testing.T) {

			values, _, ciphertext := newTestVectors(params, params.encryptorSk, 1, t)

			constant := 1.0 / complex(3.1415, -1.4142)

			for i := range values {
				values[i] *= constant
			}

			params.evaluator.MultByConst(ciphertext, constant, ciphertext)

			verifyTestVectors(params, params.decryptor, values, ciphertext, t)
		})
	}
}

func testEvaluatorMultByConstAndAdd(t *testing.T) {

	for _, parameters := range testParams.ckksParameters {

		params := genCkksParams(parameters)

		t.Run(testString("", parameters), func(t *testing.T) {

			values1, _, ciphertext1 := newTestVectors(params, params.encryptorSk, 1, t)
			values2, _, ciphertext2 := newTestVectors(params, params.encryptorSk, 1, t)

			constant := 1.0 / complex(3.1415, -1.4142)

			for i := range values1 {
				values2[i] += (constant * values1[i])
			}

			params.evaluator.MultByConstAndAdd(ciphertext1, constant, ciphertext2)

			verifyTestVectors(params, params.decryptor, values2, ciphertext2, t)
		})
	}
}

func testEvaluatorMul(t *testing.T) {

	for _, parameters := range testParams.ckksParameters {

		params := genCkksParams(parameters)

		t.Run(testString("CtCtInPlace/", parameters), func(t *testing.T) {

			values1, _, ciphertext1 := newTestVectors(params, params.encryptorSk, 1, t)
			values2, _, ciphertext2 := newTestVectors(params, params.encryptorSk, 1, t)

			for i := range values1 {
				values2[i] *= values1[i]
			}

			params.evaluator.MulRelin(ciphertext1, ciphertext2, nil, ciphertext2)

			verifyTestVectors(params, params.decryptor, values2, ciphertext2, t)
		})

		t.Run(testString("CtCtNew/", parameters), func(t *testing.T) {

			values1, _, ciphertext1 := newTestVectors(params, params.encryptorSk, 1, t)
			values2, _, ciphertext2 := newTestVectors(params, params.encryptorSk, 1, t)

			for i := range values1 {
				values2[i] *= values1[i]
			}

			var ciphertext3 *Ciphertext
			ciphertext3 = params.evaluator.MulRelinNew(ciphertext1, ciphertext2, nil)

			verifyTestVectors(params, params.decryptor, values2, ciphertext3, t)
		})

		t.Run(testString("CtPlain/", parameters), func(t *testing.T) {

			values1, plaintext1, ciphertext1 := newTestVectors(params, params.encryptorSk, 1, t)
			values2, plaintext2, ciphertext2 := newTestVectors(params, params.encryptorSk, 1, t)

			for i := range values1 {
				values1[i] *= values1[i]
				values2[i] *= values2[i]
			}

			params.evaluator.MulRelin(ciphertext1, plaintext1, nil, ciphertext1)

			verifyTestVectors(params, params.decryptor, values1, ciphertext1, t)

			params.evaluator.MulRelin(plaintext2, ciphertext2, nil, ciphertext2)

			verifyTestVectors(params, params.decryptor, values2, ciphertext2, t)
		})

		t.Run(testString("Relinearize/", parameters), func(t *testing.T) {

			rlk := params.kgen.GenRelinKey(params.sk)

			values1, _, ciphertext1 := newTestVectors(params, params.encryptorSk, 1, t)
			values2, _, ciphertext2 := newTestVectors(params, params.encryptorSk, 1, t)

			for i := range values1 {
				values2[i] *= values1[i]
			}

			params.evaluator.MulRelin(ciphertext1, ciphertext2, nil, ciphertext2)

			params.evaluator.Relinearize(ciphertext2, rlk, ciphertext2)

			require.Equal(t, ciphertext2.Degree(), uint64(1))

			verifyTestVectors(params, params.decryptor, values2, ciphertext2, t)
		})
	}
}

func testFunctions(t *testing.T) {

	for _, parameters := range testParams.ckksParameters {

		params := genCkksParams(parameters)

		rlk := params.kgen.GenRelinKey(params.sk)

		if parameters.MaxLevel > 2 {

			t.Run(testString("PowerOf2/", parameters), func(t *testing.T) {

				values, _, ciphertext := newTestVectors(params, params.encryptorSk, 1, t)

				n := uint64(2)

				valuesWant := make([]complex128, len(values))
				for i := 0; i < len(valuesWant); i++ {
					valuesWant[i] = values[i]
				}

				for i := uint64(0); i < n; i++ {
					for j := 0; j < len(valuesWant); j++ {
						valuesWant[j] *= valuesWant[j]
					}
				}

				params.evaluator.PowerOf2(ciphertext, n, rlk, ciphertext)

				verifyTestVectors(params, params.decryptor, valuesWant, ciphertext, t)
			})
		}

		if parameters.MaxLevel > 3 {

			t.Run(testString("Power/", parameters), func(t *testing.T) {

				values, _, ciphertext := newTestVectors(params, params.encryptorSk, 1, t)

				n := uint64(3)

				for i := range values {
					values[i] = cmplx.Pow(values[i], complex(float64(n), 0))
				}

				params.evaluator.Power(ciphertext, n, rlk, ciphertext)

				verifyTestVectors(params, params.decryptor, values, ciphertext, t)
			})
		}

		if parameters.MaxLevel > 6 {
			t.Run(testString("Inverse/", parameters), func(t *testing.T) {

				values, _, ciphertext := newTestVectorsReals(params, params.encryptorSk, 0.1, 1, t)

				n := uint64(7)

				for i := range values {
					values[i] = 1.0 / values[i]
				}

				ciphertext = params.evaluator.InverseNew(ciphertext, n, rlk)

				verifyTestVectors(params, params.decryptor, values, ciphertext, t)
			})
		}
	}
}

func testEvaluatePoly(t *testing.T) {

	for _, parameters := range testParams.ckksParameters {

		params := genCkksParams(parameters)

		rlk := params.kgen.GenRelinKey(params.sk)

		if parameters.MaxLevel > 4 {

			t.Run(testString("Fast/Exp/", parameters), func(t *testing.T) {

				values, _, ciphertext := newTestVectorsReals(params, params.encryptorSk, -1, 1, t)

				coeffs := []float64{1.0, 1.0, 1.0 / 2, 1.0 / 6, 1.0 / 24, 1.0 / 120, 1.0 / 720, 1.0 / 5040}

				for i := range values {
					values[i] = cmplx.Exp(values[i])
				}

				ciphertext = params.evaluator.EvaluatePolyFast(ciphertext, coeffs, rlk)

				verifyTestVectors(params, params.decryptor, values, ciphertext, t)
			})
		}

		if parameters.MaxLevel > 3 {

			t.Run(testString("Eco/Exp/", parameters), func(t *testing.T) {

				values, _, ciphertext := newTestVectorsReals(params, params.encryptorSk, -1, 1, t)

				coeffs := []float64{1.0, 1.0, 1.0 / 2, 1.0 / 6, 1.0 / 24, 1.0 / 120, 1.0 / 720, 1.0 / 5040}

				for i := range values {
					values[i] = cmplx.Exp(values[i])
				}

				ciphertext = params.evaluator.EvaluatePolyEco(ciphertext, coeffs, rlk)

				verifyTestVectors(params, params.decryptor, values, ciphertext, t)
			})
		}
	}
}

func testChebyshevInterpolator(t *testing.T) {

	for _, parameters := range testParams.ckksParameters {

		params := genCkksParams(parameters)

		rlk := params.kgen.GenRelinKey(params.sk)

		if parameters.MaxLevel > 5 {

			t.Run(testString("Fast/Sin/", parameters), func(t *testing.T) {

				values, _, ciphertext := newTestVectorsReals(params, params.encryptorSk, -1, 1, t)

				cheby := Approximate(cmplx.Sin, complex(-3, 0), complex(3, 0), 15)

				for i := range values {
					values[i] = cmplx.Sin(values[i])
				}

				ciphertext = params.evaluator.EvaluateChebyFast(ciphertext, cheby, rlk)

				verifyTestVectors(params, params.decryptor, values, ciphertext, t)
			})
		}

		if parameters.MaxLevel > 4 {

			t.Run(testString("Eco/Sin/", parameters), func(t *testing.T) {

				values, _, ciphertext := newTestVectorsReals(params, params.encryptorSk, -3, 3, t)

				cheby := Approximate(cmplx.Sin, complex(-3, 0), complex(3, 0), 15)

				for i := range values {
					values[i] = cmplx.Sin(values[i])
				}

				ciphertext = params.evaluator.EvaluateChebyEco(ciphertext, cheby, rlk)

				verifyTestVectors(params, params.decryptor, values, ciphertext, t)
			})
		}
	}
}

func testSwitchKeys(t *testing.T) {

	for _, parameters := range testParams.ckksParameters {

		params := genCkksParams(parameters)

		sk2 := params.kgen.GenSecretKey()
		decryptorSk2 := NewDecryptor(parameters, sk2)
		switchingKey := params.kgen.GenSwitchingKey(params.sk, sk2)

		t.Run(testString("InPlace/", parameters), func(t *testing.T) {

			values, _, ciphertext := newTestVectorsReals(params, params.encryptorSk, -1, 1, t)

			params.evaluator.SwitchKeys(ciphertext, switchingKey, ciphertext)

			verifyTestVectors(params, decryptorSk2, values, ciphertext, t)
		})

		t.Run(testString("New/", parameters), func(t *testing.T) {

			values, _, ciphertext := newTestVectorsReals(params, params.encryptorSk, -1, 1, t)

			ciphertext = params.evaluator.SwitchKeysNew(ciphertext, switchingKey)

			verifyTestVectors(params, decryptorSk2, values, ciphertext, t)
		})
	}
}

func testConjugate(t *testing.T) {

	for _, parameters := range testParams.ckksParameters {

		params := genCkksParams(parameters)

		rotKey := NewRotationKeys()
		params.kgen.GenRot(Conjugate, params.sk, 0, rotKey)

		t.Run(testString("InPlace/", parameters), func(t *testing.T) {

			values, _, ciphertext := newTestVectorsReals(params, params.encryptorSk, -1, 1, t)

			for i := range values {
				values[i] = complex(real(values[i]), -imag(values[i]))
			}

			params.evaluator.Conjugate(ciphertext, rotKey, ciphertext)

			verifyTestVectors(params, params.decryptor, values, ciphertext, t)
		})

		t.Run(testString("New/", parameters), func(t *testing.T) {

			values, _, ciphertext := newTestVectorsReals(params, params.encryptorSk, -1, 1, t)

			for i := range values {
				values[i] = complex(real(values[i]), -imag(values[i]))
			}

			ciphertext = params.evaluator.ConjugateNew(ciphertext, rotKey)

			verifyTestVectors(params, params.decryptor, values, ciphertext, t)
		})
	}
}

func testRotateColumns(t *testing.T) {

	for _, parameters := range testParams.ckksParameters {

		params := genCkksParams(parameters)

		rotKey := params.kgen.GenRotationKeysPow2(params.sk)

		t.Run(testString("InPlace/", parameters), func(t *testing.T) {

			values1, _, ciphertext1 := newTestVectorsReals(params, params.encryptorSk, -1, 1, t)

			values2 := make([]complex128, len(values1))
			ciphertext2 := NewCiphertext(parameters, ciphertext1.Degree(), ciphertext1.Level(), ciphertext1.Scale())

			for n := 1; n < len(values1); n <<= 1 {

				// Applies the column rotation to the values
				for i := range values1 {
					values2[i] = values1[(i+n)%len(values1)]
				}

				params.evaluator.RotateColumns(ciphertext1, uint64(n), rotKey, ciphertext2)

				verifyTestVectors(params, params.decryptor, values2, ciphertext2, t)
			}

		})

		t.Run(testString("New/", parameters), func(t *testing.T) {

			values1, _, ciphertext1 := newTestVectorsReals(params, params.encryptorSk, -1, 1, t)

			values2 := make([]complex128, len(values1))
			ciphertext2 := NewCiphertext(parameters, ciphertext1.Degree(), ciphertext1.Level(), ciphertext1.Scale())

			for n := 1; n < len(values1); n <<= 1 {

				// Applies the column rotation to the values
				for i := range values1 {
					values2[i] = values1[(i+n)%len(values1)]
				}

				ciphertext2 = params.evaluator.RotateColumnsNew(ciphertext1, uint64(n), rotKey)

				verifyTestVectors(params, params.decryptor, values2, ciphertext2, t)
			}

		})

		t.Run(testString("Random/", parameters), func(t *testing.T) {

			values1, _, ciphertext1 := newTestVectorsReals(params, params.encryptorSk, -1, 1, t)

			values2 := make([]complex128, len(values1))
			ciphertext2 := NewCiphertext(parameters, ciphertext1.Degree(), ciphertext1.Level(), ciphertext1.Scale())

			for n := 1; n < 4; n++ {

				rand := rand.Uint64() % uint64(len(values1))

				// Applies the column rotation to the values
				for i := range values1 {
					values2[i] = values1[(i+int(rand))%len(values1)]
				}

				params.evaluator.RotateColumns(ciphertext1, rand, rotKey, ciphertext2)

				verifyTestVectors(params, params.decryptor, values2, ciphertext2, t)
			}

		})

		t.Run(testString("Hoisted/", parameters), func(t *testing.T) {

			values1, _, ciphertext1 := newTestVectorsReals(params, params.encryptorSk, -1, 1, t)

			values2 := make([]complex128, len(values1))
			rotations := []uint64{1, 2, 3, 4, 5, 17, 31, 97, 127, 511}

			for _, n := range rotations {
				params.kgen.GenRot(RotationLeft, params.sk, n, rotKey)
			}

			ciphertexts := params.evaluator.RotateHoisted(ciphertext1, rotations, rotKey)

			for _, n := range rotations {

				for i := range values1 {
					values2[i] = values1[(i+int(n))%len(values1)]
				}

				verifyTestVectors(params, params.decryptor, values2, ciphertexts[n], t)
			}
		})
	}
}

func testMarshaller(t *testing.T) {

	for _, parameters := range testParams.ckksParameters {

		params := genCkksParams(parameters)

		contextQP := params.ckkscontext.contextQP

		t.Run(testString("Ciphertext/", parameters), func(t *testing.T) {
			t.Run(testString("Ciphertext/EndToEnd", parameters), func(t *testing.T) {
				t.Parallel()

				ciphertextWant := NewCiphertextRandom(parameters, 2, parameters.MaxLevel, parameters.Scale)

				marshalledCiphertext, err := ciphertextWant.MarshalBinary()
				require.NoError(t, err)

				ciphertextTest := new(Ciphertext)
				require.NoError(t, ciphertextTest.UnmarshalBinary(marshalledCiphertext))

				require.Equal(t, ciphertextWant.Degree(), ciphertextTest.Degree())
				require.Equal(t, ciphertextWant.Level(), ciphertextTest.Level())
				require.Equal(t, ciphertextWant.Scale(), ciphertextTest.Scale())

				for i := range ciphertextWant.value {
					require.True(t, params.ckkscontext.contextQ.EqualLvl(ciphertextWant.Level(), ciphertextWant.Value()[i], ciphertextTest.Value()[i]))
				}
			})

			t.Run(testString("Ciphertext/Minimal", parameters), func(t *testing.T) {
				t.Parallel()

				ciphertext := NewCiphertextRandom(parameters, 0, parameters.MaxLevel, parameters.Scale)

				marshalledCiphertext, err := ciphertext.MarshalBinary()
				require.NoError(t, err)

				ciphertextTest := new(Ciphertext)
				require.Error(t, ciphertextTest.UnmarshalBinary(nil))
				require.NoError(t, ciphertextTest.UnmarshalBinary(marshalledCiphertext))

				require.Equal(t, ciphertext.Degree(), uint64(0))
				require.Equal(t, ciphertext.Level(), parameters.MaxLevel)
				require.Equal(t, ciphertext.Scale(), parameters.Scale)
				require.Equal(t, len(ciphertext.Value()), 1)
			})
		})

		t.Run(testString("Sk", parameters), func(t *testing.T) {

			marshalledSk, err := params.sk.MarshalBinary()
			require.NoError(t, err)

			sk := new(SecretKey)
			err = sk.UnmarshalBinary(marshalledSk)
			require.NoError(t, err)

			require.True(t, contextQP.Equal(sk.sk, params.sk.sk))

		})

		t.Run(testString("Pk", parameters), func(t *testing.T) {

			marshalledPk, err := params.pk.MarshalBinary()
			require.NoError(t, err)

			pk := new(PublicKey)
			err = pk.UnmarshalBinary(marshalledPk)
			require.NoError(t, err)

			for k := range params.pk.pk {
				require.Truef(t, contextQP.Equal(pk.pk[k], params.pk.pk[k]), "Marshal PublicKey element [%d]", k)
			}
		})

		t.Run(testString("EvaluationKey", parameters), func(t *testing.T) {

			evalKey := params.kgen.GenRelinKey(params.sk)
			data, err := evalKey.MarshalBinary()
			require.NoError(t, err)

			resEvalKey := new(EvaluationKey)
			err = resEvalKey.UnmarshalBinary(data)
			require.NoError(t, err)

			evakeyWant := evalKey.evakey.evakey
			evakeyTest := resEvalKey.evakey.evakey

			for j := range evakeyWant {
				for k := range evakeyWant[j] {
					require.Truef(t, contextQP.Equal(evakeyWant[j][k], evakeyTest[j][k]), "Marshal EvaluationKey element [%d][%d]", j, k)
				}
			}
		})

		t.Run(testString("SwitchingKey", parameters), func(t *testing.T) {

			skOut := params.kgen.GenSecretKey()

			switchingKey := params.kgen.GenSwitchingKey(params.sk, skOut)
			data, err := switchingKey.MarshalBinary()
			require.NoError(t, err)

			resSwitchingKey := new(SwitchingKey)
			err = resSwitchingKey.UnmarshalBinary(data)
			require.NoError(t, err)

			evakeyWant := switchingKey.evakey
			evakeyTest := resSwitchingKey.evakey

			for j := range evakeyWant {
				for k := range evakeyWant[j] {
					require.True(t, contextQP.Equal(evakeyWant[j][k], evakeyTest[j][k]))
				}
			}
		})

		t.Run(testString("RotationKey", parameters), func(t *testing.T) {

			rotationKey := NewRotationKeys()

			params.kgen.GenRot(Conjugate, params.sk, 0, rotationKey)
			params.kgen.GenRot(RotationLeft, params.sk, 1, rotationKey)
			params.kgen.GenRot(RotationLeft, params.sk, 2, rotationKey)
			params.kgen.GenRot(RotationRight, params.sk, 3, rotationKey)
			params.kgen.GenRot(RotationRight, params.sk, 5, rotationKey)

			data, err := rotationKey.MarshalBinary()
			require.NoError(t, err)

			resRotationKey := new(RotationKeys)
			err = resRotationKey.UnmarshalBinary(data)
			require.NoError(t, err)

			for i := uint64(1); i < params.ckkscontext.n>>1; i++ {

				if rotationKey.evakeyRotColLeft[i] != nil {

					evakeyWant := rotationKey.evakeyRotColLeft[i].evakey
					evakeyTest := resRotationKey.evakeyRotColLeft[i].evakey

					evakeyNTTIndexWant := rotationKey.permuteNTTLeftIndex[i]
					evakeyNTTIndexTest := resRotationKey.permuteNTTLeftIndex[i]

					require.True(t, utils.EqualSliceUint64(evakeyNTTIndexWant, evakeyNTTIndexTest))

					for j := range evakeyWant {
						for k := range evakeyWant[j] {
							require.Truef(t, contextQP.Equal(evakeyWant[j][k], evakeyTest[j][k]), "Marshal RotationKey RotateLeft %d element [%d][%d]", i, j, k)
						}
					}
				}

				if rotationKey.evakeyRotColRight[i] != nil {

					evakeyWant := rotationKey.evakeyRotColRight[i].evakey
					evakeyTest := resRotationKey.evakeyRotColRight[i].evakey

					evakeyNTTIndexWant := rotationKey.permuteNTTRightIndex[i]
					evakeyNTTIndexTest := resRotationKey.permuteNTTRightIndex[i]

					require.True(t, utils.EqualSliceUint64(evakeyNTTIndexWant, evakeyNTTIndexTest))

					for j := range evakeyWant {
						for k := range evakeyWant[j] {
							require.Truef(t, contextQP.Equal(evakeyWant[j][k], evakeyTest[j][k]), "Marshal RotationKey RotateRight %d element [%d][%d]", i, j, k)
						}
					}
				}
			}

			if rotationKey.evakeyConjugate != nil {

				evakeyWant := rotationKey.evakeyConjugate.evakey
				evakeyTest := resRotationKey.evakeyConjugate.evakey

				evakeyNTTIndexWant := rotationKey.permuteNTTConjugateIndex
				evakeyNTTIndexTest := resRotationKey.permuteNTTConjugateIndex

				require.True(t, utils.EqualSliceUint64(evakeyNTTIndexWant, evakeyNTTIndexTest))

				for j := range evakeyWant {
					for k := range evakeyWant[j] {
						require.Truef(t, contextQP.Equal(evakeyWant[j][k], evakeyTest[j][k]), "Marshal RotationKey RotateRow element [%d][%d]", j, k)
					}
				}
			}
		})
	}
}
