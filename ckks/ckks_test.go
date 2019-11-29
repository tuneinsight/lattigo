package ckks

import (
	"fmt"
	"log"
	"math"
	"math/cmplx"
	"math/rand"
	"sort"
	"testing"
	"time"
)

func check(t *testing.T, err error) {
	if err != nil {
		t.Error(err)
		log.Fatal(err)
	}
}

func testString(opname string, params *ckksParams) string {
	return fmt.Sprintf("%slogN=%d/logSlots=%d/logQ=%d/levels=%d/a=%d/b=%d",
		opname,
		params.ckkscontext.logN,
		params.ckkscontext.logN-1,
		params.ckkscontext.logQ,
		params.ckkscontext.levels,
		params.ckkscontext.alpha,
		params.ckkscontext.beta)
}

type ckksParams struct {
	params      *Parameters
	ckkscontext *Context
	encoder     *Encoder
	kgen        *KeyGenerator
	sk          *SecretKey
	pk          *PublicKey
	encryptorPk *Encryptor
	encryptorSk *Encryptor
	decryptor   *Decryptor
	evaluator   *Evaluator
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
	testParams.verbose = false

	testParams.ckksParameters = []*Parameters{
		DefaultParams[13],
		DefaultParams[14],
		//DefaultParams[15],
		//DefaultParams[16],
	}
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

	params.ckkscontext = NewContext(contextParameters)

	params.kgen = NewKeyGenerator(contextParameters)

	params.sk, params.pk = params.kgen.NewKeyPairSparse(128)

	params.encoder = NewEncoder(contextParameters)

	params.encryptorPk = NewEncryptorFromPk(contextParameters, params.pk)
	params.encryptorSk = NewEncryptorFromSk(contextParameters, params.sk)
	params.decryptor = NewDecryptor(contextParameters, params.sk)

	params.evaluator = NewEvaluator(contextParameters)

	return

}

func newTestVectors(contextParams *ckksParams, encryptor *Encryptor, a float64, t *testing.T) (values []complex128, plaintext *Plaintext, ciphertext *Ciphertext) {

	slots := contextParams.ckkscontext.Slots()

	values = make([]complex128, slots)

	for i := uint64(0); i < slots; i++ {
		values[i] = randomComplex(-a, a)
	}

	values[0] = complex(0.607538, 0.555668)

	plaintext = NewPlaintext(contextParams.params, contextParams.ckkscontext.Levels()-1, contextParams.ckkscontext.Scale())

	contextParams.encoder.Encode(plaintext, values, slots)

	if encryptor != nil {
		ciphertext = encryptor.EncryptNew(plaintext)
	}

	return values, plaintext, ciphertext
}

func newTestVectorsReals(contextParams *ckksParams, encryptor *Encryptor, a, b float64, t *testing.T) (values []complex128, plaintext *Plaintext, ciphertext *Ciphertext) {

	slots := contextParams.ckkscontext.Slots()

	values = make([]complex128, slots)

	for i := uint64(0); i < slots; i++ {
		values[i] = complex(randomFloat(a, b), 0)
	}

	values[0] = complex(0.607538, 0)

	plaintext = NewPlaintext(contextParams.params, contextParams.ckkscontext.Levels()-1, contextParams.ckkscontext.Scale())

	contextParams.encoder.Encode(plaintext, values, slots)

	if encryptor != nil {
		ciphertext = encryptor.EncryptNew(plaintext)
	}

	return values, plaintext, ciphertext
}

func verifyTestVectors(contextParams *ckksParams, decryptor *Decryptor, valuesWant []complex128, element interface{}, t *testing.T) {

	var plaintextTest *Plaintext
	var valuesTest []complex128

	switch element.(type) {
	case *Ciphertext:
		plaintextTest = decryptor.DecryptNew(element.(*Ciphertext))
	case *Plaintext:
		plaintextTest = element.(*Plaintext)
	}

	valuesTest = contextParams.encoder.Decode(plaintextTest, contextParams.ckkscontext.Slots())

	var deltaReal, deltaImag float64

	var minprec, maxprec, meanprec, medianprec complex128

	diff := make([]complex128, contextParams.ckkscontext.Slots())

	minprec = complex(0, 0)
	maxprec = complex(1, 1)

	meanprec = complex(0, 0)

	distribReal := make(map[uint64]uint64)
	distribImag := make(map[uint64]uint64)

	for i := range valuesWant {

		deltaReal = math.Abs(real(valuesTest[i]) - real(valuesWant[i]))
		deltaImag = math.Abs(imag(valuesTest[i]) - imag(valuesWant[i]))

		diff[i] += complex(deltaReal, 0)
		diff[i] += complex(0, deltaImag)

		meanprec += diff[i]

		if real(diff[i]) > real(minprec) {
			minprec = complex(real(diff[i]), 0)
		}

		if imag(diff[i]) > imag(minprec) {
			minprec = complex(real(minprec), imag(diff[i]))
		}

		if real(diff[i]) < real(maxprec) {
			maxprec = complex(real(diff[i]), 0)
		}

		if imag(diff[i]) < imag(maxprec) {
			maxprec = complex(real(maxprec), imag(diff[i]))
		}

		distribReal[uint64(math.Floor(math.Log2(1/real(diff[i]))))]++
		distribImag[uint64(math.Floor(math.Log2(1/imag(diff[i]))))]++
	}

	meanprec /= complex(float64(contextParams.ckkscontext.Slots()), 0)
	medianprec = calcmedian(diff)

	if testParams.verbose {
		t.Logf("Minimum precision : (%.2f, %.2f) bits \n", math.Log2(1/real(minprec)), math.Log2(1/imag(minprec)))
		t.Logf("Maximum precision : (%.2f, %.2f) bits \n", math.Log2(1/real(maxprec)), math.Log2(1/imag(maxprec)))
		t.Logf("Mean    precision : (%.2f, %.2f) bits \n", math.Log2(1/real(meanprec)), math.Log2(1/imag(meanprec)))
		t.Logf("Median  precision : (%.2f, %.2f) bits \n", math.Log2(1/real(medianprec)), math.Log2(1/imag(medianprec)))
		t.Log()
	}

	if math.Log2(1/real(medianprec)) < testParams.medianprec || math.Log2(1/imag(medianprec)) < testParams.medianprec {
		t.Errorf("Mean precision error : target (%.2f, %.2f) > result (%.2f, %.2f)", testParams.medianprec, testParams.medianprec, math.Log2(1/real(medianprec)), math.Log2(1/imag(medianprec)))
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

func testEncoder(t *testing.T) {

	for _, parameters := range testParams.ckksParameters {

		params := genCkksParams(parameters)

		t.Run(testString("Encode/", params), func(t *testing.T) {

			values, plaintext, _ := newTestVectors(params, nil, 1, t)

			verifyTestVectors(params, params.decryptor, values, plaintext, t)
		})
	}
}

func testEncryptor(t *testing.T) {

	for _, parameters := range testParams.ckksParameters {

		params := genCkksParams(parameters)

		t.Run(testString("EncryptFromPk/", params), func(t *testing.T) {

			values, _, ciphertext := newTestVectors(params, params.encryptorPk, 1, t)

			verifyTestVectors(params, params.decryptor, values, ciphertext, t)
		})
	}

	for _, parameters := range testParams.ckksParameters {

		params := genCkksParams(parameters)

		t.Run(testString("EncryptFromSk/", params), func(t *testing.T) {

			values, _, ciphertext := newTestVectors(params, params.encryptorSk, 1, t)

			verifyTestVectors(params, params.decryptor, values, ciphertext, t)
		})
	}
}

func testEvaluatorAdd(t *testing.T) {

	for _, parameters := range testParams.ckksParameters {

		params := genCkksParams(parameters)

		t.Run(testString("CtCtInPlace/", params), func(t *testing.T) {

			values1, _, ciphertext1 := newTestVectors(params, params.encryptorSk, 1, t)
			values2, _, ciphertext2 := newTestVectors(params, params.encryptorSk, 1, t)

			for i := range values1 {
				values1[i] += values2[i]
			}

			params.evaluator.Add(ciphertext1, ciphertext2, ciphertext1)

			verifyTestVectors(params, params.decryptor, values1, ciphertext1, t)
		})

		t.Run(testString("CtCtNew/", params), func(t *testing.T) {

			values1, _, ciphertext1 := newTestVectors(params, params.encryptorSk, 1, t)
			values2, _, ciphertext2 := newTestVectors(params, params.encryptorSk, 1, t)

			for i := range values1 {
				values1[i] += values2[i]
			}

			ciphertext3 := params.evaluator.AddNew(ciphertext1, ciphertext2)

			verifyTestVectors(params, params.decryptor, values1, ciphertext3, t)
		})

		t.Run(testString("CtPlainInPlace/", params), func(t *testing.T) {

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

		t.Run(testString("CtPlainInPlaceNew/", params), func(t *testing.T) {

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

		t.Run(testString("CtCtInPlace/", params), func(t *testing.T) {

			values1, _, ciphertext1 := newTestVectors(params, params.encryptorSk, 1, t)
			values2, _, ciphertext2 := newTestVectors(params, params.encryptorSk, 1, t)

			for i := range values1 {
				values1[i] -= values2[i]
			}

			params.evaluator.Sub(ciphertext1, ciphertext2, ciphertext1)

			verifyTestVectors(params, params.decryptor, values1, ciphertext1, t)
		})

		t.Run(testString("CtCtNew/", params), func(t *testing.T) {

			values1, _, ciphertext1 := newTestVectors(params, params.encryptorSk, 1, t)
			values2, _, ciphertext2 := newTestVectors(params, params.encryptorSk, 1, t)

			for i := range values1 {
				values1[i] -= values2[i]
			}

			ciphertext3 := params.evaluator.SubNew(ciphertext1, ciphertext2)

			verifyTestVectors(params, params.decryptor, values1, ciphertext3, t)
		})

		t.Run(testString("CtPlainInPlace/", params), func(t *testing.T) {

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

		t.Run(testString("CtPlainNew/", params), func(t *testing.T) {

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

		t.Run(testString("Single/", params), func(t *testing.T) {

			values, _, ciphertext := newTestVectors(params, params.encryptorSk, 1, t)

			constant := params.ckkscontext.contextQ.Modulus[ciphertext.Level()]

			params.evaluator.MultByConst(ciphertext, constant, ciphertext)

			ciphertext.MulScale(float64(constant))

			params.evaluator.Rescale(ciphertext, params.ckkscontext.scale, ciphertext)

			verifyTestVectors(params, params.decryptor, values, ciphertext, t)
		})

		t.Run(testString("Many/", params), func(t *testing.T) {

			values, _, ciphertext := newTestVectors(params, params.encryptorSk, 1, t)

			nbRescales := uint64(2)

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

		t.Run(testString("", params), func(t *testing.T) {

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

		t.Run(testString("", params), func(t *testing.T) {

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

		t.Run(testString("", params), func(t *testing.T) {

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

		t.Run(testString("CtCtInPlace/", params), func(t *testing.T) {

			values1, _, ciphertext1 := newTestVectors(params, params.encryptorSk, 1, t)
			values2, _, ciphertext2 := newTestVectors(params, params.encryptorSk, 1, t)

			for i := range values1 {
				values2[i] *= values1[i]
			}

			params.evaluator.MulRelin(ciphertext1, ciphertext2, nil, ciphertext2)

			verifyTestVectors(params, params.decryptor, values2, ciphertext2, t)
		})

		t.Run(testString("CtCtNew/", params), func(t *testing.T) {

			values1, _, ciphertext1 := newTestVectors(params, params.encryptorSk, 1, t)
			values2, _, ciphertext2 := newTestVectors(params, params.encryptorSk, 1, t)

			for i := range values1 {
				values2[i] *= values1[i]
			}

			var ciphertext3 *Ciphertext
			ciphertext3 = params.evaluator.MulRelinNew(ciphertext1, ciphertext2, nil)

			verifyTestVectors(params, params.decryptor, values2, ciphertext3, t)
		})

		t.Run(testString("CtPlain/", params), func(t *testing.T) {

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

		t.Run(testString("Relinearize/", params), func(t *testing.T) {

			rlk := params.kgen.NewRelinKey(params.sk)

			values1, _, ciphertext1 := newTestVectors(params, params.encryptorSk, 1, t)
			values2, _, ciphertext2 := newTestVectors(params, params.encryptorSk, 1, t)

			for i := range values1 {
				values2[i] *= values1[i]
			}

			params.evaluator.MulRelin(ciphertext1, ciphertext2, nil, ciphertext2)

			params.evaluator.Relinearize(ciphertext2, rlk, ciphertext2)

			if ciphertext2.Degree() != 1 {
				t.Errorf("Relinearize error")
			}

			verifyTestVectors(params, params.decryptor, values2, ciphertext2, t)
		})
	}
}

func testFunctions(t *testing.T) {

	for _, parameters := range testParams.ckksParameters {

		params := genCkksParams(parameters)

		rlk := params.kgen.NewRelinKey(params.sk)

		t.Run(testString("PowerOf2/", params), func(t *testing.T) {

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

		t.Run(testString("Power/", params), func(t *testing.T) {

			values, _, ciphertext := newTestVectors(params, params.encryptorSk, 1, t)

			n := uint64(7)

			for i := range values {
				values[i] = cmplx.Pow(values[i], complex(float64(n), 0))
			}

			params.evaluator.Power(ciphertext, n, rlk, ciphertext)

			verifyTestVectors(params, params.decryptor, values, ciphertext, t)
		})

		if params.ckkscontext.levels > 7 {
			t.Run(testString("Inverse/", params), func(t *testing.T) {

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

		rlk := params.kgen.NewRelinKey(params.sk)

		t.Run(testString("Fast/Exp/", params), func(t *testing.T) {

			values, _, ciphertext := newTestVectorsReals(params, params.encryptorSk, -1, 1, t)

			coeffs := []float64{1.0, 1.0, 1.0 / 2, 1.0 / 6, 1.0 / 24, 1.0 / 120, 1.0 / 720, 1.0 / 5040}

			for i := range values {
				values[i] = cmplx.Exp(values[i])
			}

			ciphertext = params.evaluator.EvaluatePolyFast(ciphertext, coeffs, rlk)

			verifyTestVectors(params, params.decryptor, values, ciphertext, t)
		})

		t.Run(testString("Eco/Exp/", params), func(t *testing.T) {

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

func testChebyshevInterpolator(t *testing.T) {

	for _, parameters := range testParams.ckksParameters {

		params := genCkksParams(parameters)

		rlk := params.kgen.NewRelinKey(params.sk)

		t.Run(testString("Fast/Sin/", params), func(t *testing.T) {

			values, _, ciphertext := newTestVectorsReals(params, params.encryptorSk, -1, 1, t)

			cheby := Approximate(cmplx.Sin, complex(-1, 0), complex(1, 0), 16)

			for i := range values {
				values[i] = cmplx.Sin(values[i])
			}

			ciphertext = params.evaluator.EvaluateChebyFast(ciphertext, cheby, rlk)

			verifyTestVectors(params, params.decryptor, values, ciphertext, t)
		})

		t.Run(testString("Eco/Sin/", params), func(t *testing.T) {

			values, _, ciphertext := newTestVectorsReals(params, params.encryptorSk, -1, 1, t)

			cheby := Approximate(cmplx.Sin, complex(-1, 0), complex(1, 0), 16)

			for i := range values {
				values[i] = cmplx.Sin(values[i])
			}

			ciphertext = params.evaluator.EvaluateChebyEco(ciphertext, cheby, rlk)

			verifyTestVectors(params, params.decryptor, values, ciphertext, t)
		})
	}
}

func testSwitchKeys(t *testing.T) {

	for _, parameters := range testParams.ckksParameters {

		params := genCkksParams(parameters)

		sk2 := params.kgen.NewSecretKey()
		decryptorSk2 := NewDecryptor(parameters, sk2)
		switchingKey := params.kgen.NewSwitchingKey(params.sk, sk2)

		t.Run(testString("InPlace/", params), func(t *testing.T) {

			values, _, ciphertext := newTestVectorsReals(params, params.encryptorSk, -1, 1, t)

			params.evaluator.SwitchKeys(ciphertext, switchingKey, ciphertext)

			verifyTestVectors(params, decryptorSk2, values, ciphertext, t)
		})

		t.Run(testString("New/", params), func(t *testing.T) {

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

		t.Run(testString("InPlace/", params), func(t *testing.T) {

			values, _, ciphertext := newTestVectorsReals(params, params.encryptorSk, -1, 1, t)

			for i := range values {
				values[i] = complex(real(values[i]), -imag(values[i]))
			}

			params.evaluator.Conjugate(ciphertext, rotKey, ciphertext)

			verifyTestVectors(params, params.decryptor, values, ciphertext, t)
		})

		t.Run(testString("New/", params), func(t *testing.T) {

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

		rotKey := params.kgen.NewRotationKeysPow2(params.sk)

		t.Run(testString("InPlace/", params), func(t *testing.T) {

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

		t.Run(testString("New/", params), func(t *testing.T) {

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

		t.Run(testString("Random/", params), func(t *testing.T) {

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
	}
}

func testMarshaller(t *testing.T) {

	for _, parameters := range testParams.ckksParameters {

		params := genCkksParams(parameters)

		contextQP := params.ckkscontext.contextQP

		t.Run(testString("Ciphertext/", params), func(t *testing.T) {

			ciphertextWant := NewCiphertextRandom(parameters, 2, params.ckkscontext.levels-1, params.ckkscontext.scale)

			marshalledCiphertext, err := ciphertextWant.MarshalBinary()
			check(t, err)

			ciphertextTest := NewCkksElement().Ciphertext()
			err = ciphertextTest.UnmarshalBinary(marshalledCiphertext)
			check(t, err)

			if ciphertextWant.Degree() != ciphertextTest.Degree() {
				t.Errorf("Marshal Cipehrtext Degree")
			}

			if ciphertextWant.Level() != ciphertextTest.Level() {
				t.Errorf("Marshal Cipehrtext Level")
			}

			if ciphertextWant.Scale() != ciphertextTest.Scale() {
				t.Errorf("Marshal Cipehrtext Scale")
			}

			for i := range ciphertextWant.value {
				if !params.ckkscontext.contextQ.EqualLvl(ciphertextWant.Level(), ciphertextWant.Value()[i], ciphertextTest.Value()[i]) {
					t.Errorf("Marshal Ciphertext Coefficients")
				}
			}
		})

		t.Run(testString("Sk", params), func(t *testing.T) {

			marshalledSk, err := params.sk.MarshalBinary()
			check(t, err)

			sk := new(SecretKey)
			err = sk.UnmarshalBinary(marshalledSk)
			check(t, err)

			if !contextQP.Equal(sk.sk, params.sk.sk) {
				t.Errorf("Marshal SecretKey")
			}

		})

		t.Run(testString("Pk", params), func(t *testing.T) {

			marshalledPk, err := params.pk.MarshalBinary()
			check(t, err)

			pk := new(PublicKey)
			err = pk.UnmarshalBinary(marshalledPk)
			check(t, err)

			for k := range params.pk.pk {
				if !contextQP.Equal(pk.pk[k], params.pk.pk[k]) {
					t.Errorf("Marshal PublicKey element [%d]", k)
				}
			}
		})

		t.Run(testString("EvaluationKey", params), func(t *testing.T) {

			evalKey := params.kgen.NewRelinKey(params.sk)
			data, err := evalKey.MarshalBinary()
			check(t, err)

			resEvalKey := new(EvaluationKey)
			err = resEvalKey.UnmarshalBinary(data)
			check(t, err)

			evakeyWant := evalKey.evakey.evakey
			evakeyTest := resEvalKey.evakey.evakey

			for j := range evakeyWant {

				for k := range evakeyWant[j] {
					if !contextQP.Equal(evakeyWant[j][k], evakeyTest[j][k]) {
						t.Errorf("Marshal EvaluationKey element [%d][%d]", j, k)
					}
				}
			}
		})

		t.Run(testString("SwitchingKey", params), func(t *testing.T) {

			skOut := params.kgen.NewSecretKey()

			switchingKey := params.kgen.NewSwitchingKey(params.sk, skOut)
			data, err := switchingKey.MarshalBinary()
			check(t, err)

			resSwitchingKey := new(SwitchingKey)
			err = resSwitchingKey.UnmarshalBinary(data)
			check(t, err)

			evakeyWant := switchingKey.evakey
			evakeyTest := resSwitchingKey.evakey

			for j := range evakeyWant {

				for k := range evakeyWant[j] {
					if !contextQP.Equal(evakeyWant[j][k], evakeyTest[j][k]) {
						t.Errorf("Marshal Switchkey element [%d][%d]", j, k)
					}
				}
			}
		})

		t.Run(testString("RotationKey", params), func(t *testing.T) {

			rotationKey := NewRotationKeys()

			params.kgen.GenRot(Conjugate, params.sk, 0, rotationKey)
			params.kgen.GenRot(RotationLeft, params.sk, 1, rotationKey)
			params.kgen.GenRot(RotationLeft, params.sk, 2, rotationKey)
			params.kgen.GenRot(RotationRight, params.sk, 3, rotationKey)
			params.kgen.GenRot(RotationRight, params.sk, 5, rotationKey)

			data, err := rotationKey.MarshalBinary()
			check(t, err)

			resRotationKey := new(RotationKeys)
			err = resRotationKey.UnmarshalBinary(data)
			check(t, err)

			for i := uint64(1); i < params.ckkscontext.n>>1; i++ {

				if rotationKey.evakeyRotColLeft[i] != nil {

					evakeyWant := rotationKey.evakeyRotColLeft[i].evakey
					evakeyTest := resRotationKey.evakeyRotColLeft[i].evakey

					for j := range evakeyWant {

						for k := range evakeyWant[j] {
							if !contextQP.Equal(evakeyWant[j][k], evakeyTest[j][k]) {
								t.Errorf("Marshal RotKey RotateLeft %d element [%d][%d]", i, j, k)
							}
						}
					}
				}

				if rotationKey.evakeyRotColRight[i] != nil {

					evakeyWant := rotationKey.evakeyRotColRight[i].evakey
					evakeyTest := resRotationKey.evakeyRotColRight[i].evakey

					for j := range evakeyWant {

						for k := range evakeyWant[j] {
							if !contextQP.Equal(evakeyWant[j][k], evakeyTest[j][k]) {
								t.Errorf("Marshal RotKey RotateRight %d element [%d][%d]", i, j, k)
							}
						}
					}
				}
			}

			if rotationKey.evakeyConjugate != nil {

				evakeyWant := rotationKey.evakeyConjugate.evakey
				evakeyTest := resRotationKey.evakeyConjugate.evakey

				for j := range evakeyWant {

					for k := range evakeyWant[j] {
						if !contextQP.Equal(evakeyWant[j][k], evakeyTest[j][k]) {
							t.Errorf("Marshal RotKey RotateRow element [%d][%d]", j, k)
						}
					}
				}
			}
		})
	}
}
