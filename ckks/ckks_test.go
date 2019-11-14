package ckks

import (
	"fmt"
	//"github.com/ldsec/lattigo/ring"
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
	return fmt.Sprintf("%slogN=%d/logQ=%d/levels=%d/a=%d/b=%d", opname, params.ckksContext.logN, params.ckksContext.logQ, params.ckksContext.levels, params.ckksContext.alpha, params.ckksContext.beta)
}

type ckksParams struct {
	ckksContext *CkksContext
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

func Test_CKKS(t *testing.T) {
	t.Run("Encoder", testEncoder)
	t.Run("Encryptor", testEncryptor)
	t.Run("Evaluator/Add", testEvaluatorAdd)
	t.Run("Evaluator/Sub", testEvaluatorSub)
	t.Run("Evaluator/Rescale", testEvaluatorRescale)
	t.Run("Evaluator/AddConst", testEvaluatorAddConst)
	t.Run("Evaluator/MultConst", testEvaluatorMultConst)
	t.Run("Evaluator/MultConstAndAdd", testEvaluatorMultConstAndAdd)
	t.Run("Evaluator/Mul", testEvaluatorMul)
	t.Run("Evaluator/Functions", testFunctions)
	t.Run("Evaluator/ChebyshevInterpolator", testChebyshevInterpolator)
	t.Run("Evaluator/SwitchKeys", testSwitchKeys)
	t.Run("Evaluator/Conjugate", testConjugate)
	t.Run("Evaluator/RotateColumns", testRotateColumns)
}

func genCkksParams(contextParameters *Parameters) (params *ckksParams) {

	params = new(ckksParams)

	if params.ckksContext, err = NewCkksContext(contextParameters); err != nil {
		log.Fatal(err)
	}

	params.kgen = params.ckksContext.NewKeyGenerator()

	params.sk, params.pk = params.kgen.NewKeyPairSparse(128)

	params.encoder = params.ckksContext.NewEncoder()

	if params.encryptorPk, err = params.ckksContext.NewEncryptorFromPk(params.pk); err != nil {
		log.Fatal(err)
	}

	if params.encryptorSk, err = params.ckksContext.NewEncryptorFromSk(params.sk); err != nil {
		log.Fatal(err)
	}

	if params.decryptor, err = params.ckksContext.NewDecryptor(params.sk); err != nil {
		log.Fatal(err)
	}

	params.evaluator = params.ckksContext.NewEvaluator()

	return

}

func new_test_vectors(contextParams *ckksParams, encryptor *Encryptor, a float64, t *testing.T) (values []complex128, plaintext *Plaintext, ciphertext *Ciphertext) {

	slots := contextParams.ckksContext.Slots()

	values = make([]complex128, slots)

	for i := uint64(0); i < slots; i++ {
		values[i] = randomComplex(-a, a)
	}

	values[0] = complex(0.607538, 0.555668)

	plaintext = contextParams.ckksContext.NewPlaintext(contextParams.ckksContext.Levels()-1, contextParams.ckksContext.Scale())

	check(t, contextParams.encoder.Encode(plaintext, values, slots))

	if encryptor != nil {
		ciphertext, err = encryptor.EncryptNew(plaintext)
		check(t, err)
	}

	return values, plaintext, ciphertext
}

func new_test_vectors_reals(contextParams *ckksParams, encryptor *Encryptor, a, b float64, t *testing.T) (values []complex128, plaintext *Plaintext, ciphertext *Ciphertext) {

	slots := contextParams.ckksContext.Slots()

	values = make([]complex128, slots)

	for i := uint64(0); i < slots; i++ {
		values[i] = complex(randomFloat(a, b), 0)
	}

	values[0] = complex(0.607538, 0)

	plaintext = contextParams.ckksContext.NewPlaintext(contextParams.ckksContext.Levels()-1, contextParams.ckksContext.Scale())

	check(t, contextParams.encoder.Encode(plaintext, values, slots))

	if encryptor != nil {
		ciphertext, err = encryptor.EncryptNew(plaintext)
		check(t, err)
	}

	return values, plaintext, ciphertext
}

func verify_test_vectors(contextParams *ckksParams, decryptor *Decryptor, valuesWant []complex128, element interface{}, t *testing.T) {

	var plaintextTest *Plaintext
	var valuesTest []complex128

	switch element.(type) {
	case *Ciphertext:
		plaintextTest = decryptor.DecryptNew(element.(*Ciphertext))
	case *Plaintext:
		plaintextTest = element.(*Plaintext)
	}

	valuesTest = contextParams.encoder.Decode(plaintextTest, contextParams.ckksContext.Slots())

	var deltaReal, deltaImag float64

	var minprec, maxprec, meanprec, medianprec complex128

	diff := make([]complex128, contextParams.ckksContext.Slots())

	minprec = complex(0, 0)
	maxprec = complex(1, 1)

	meanprec = complex(0, 0)

	distrib_real := make(map[uint64]uint64)
	distrib_imag := make(map[uint64]uint64)

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

		distrib_real[uint64(math.Floor(math.Log2(1/real(diff[i]))))] += 1
		distrib_imag[uint64(math.Floor(math.Log2(1/imag(diff[i]))))] += 1
	}

	meanprec /= complex(float64(contextParams.ckksContext.Slots()), 0)
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

			values, plaintext, _ := new_test_vectors(params, nil, 1, t)

			verify_test_vectors(params, params.decryptor, values, plaintext, t)
		})
	}
}

func testEncryptor(t *testing.T) {

	for _, parameters := range testParams.ckksParameters {

		params := genCkksParams(parameters)

		t.Run(testString("EncryptFromPk/", params), func(t *testing.T) {

			values, _, ciphertext := new_test_vectors(params, params.encryptorPk, 1, t)

			verify_test_vectors(params, params.decryptor, values, ciphertext, t)
		})
	}

	for _, parameters := range testParams.ckksParameters {

		params := genCkksParams(parameters)

		t.Run(testString("EncryptFromSk/", params), func(t *testing.T) {

			values, _, ciphertext := new_test_vectors(params, params.encryptorSk, 1, t)

			verify_test_vectors(params, params.decryptor, values, ciphertext, t)
		})
	}
}

func testEvaluatorAdd(t *testing.T) {

	for _, parameters := range testParams.ckksParameters {

		params := genCkksParams(parameters)

		t.Run(testString("CtCtInPlace/", params), func(t *testing.T) {

			values1, _, ciphertext1 := new_test_vectors(params, params.encryptorSk, 1, t)
			values2, _, ciphertext2 := new_test_vectors(params, params.encryptorSk, 1, t)

			for i := range values1 {
				values1[i] += values2[i]
			}

			check(t, params.evaluator.Add(ciphertext1, ciphertext2, ciphertext1))

			verify_test_vectors(params, params.decryptor, values1, ciphertext1, t)
		})

		t.Run(testString("CtCtNew/", params), func(t *testing.T) {

			values1, _, ciphertext1 := new_test_vectors(params, params.encryptorSk, 1, t)
			values2, _, ciphertext2 := new_test_vectors(params, params.encryptorSk, 1, t)

			for i := range values1 {
				values1[i] += values2[i]
			}

			ciphertext3, err := params.evaluator.AddNew(ciphertext1, ciphertext2)

			check(t, err)

			verify_test_vectors(params, params.decryptor, values1, ciphertext3, t)
		})

		t.Run(testString("CtPlainInPlace/", params), func(t *testing.T) {

			values1, _, ciphertext1 := new_test_vectors(params, params.encryptorSk, 1, t)
			values2, plaintext2, _ := new_test_vectors(params, params.encryptorSk, 1, t)

			for i := range values1 {
				values1[i] += values2[i]
			}

			check(t, params.evaluator.Add(ciphertext1, plaintext2, ciphertext1))

			verify_test_vectors(params, params.decryptor, values1, ciphertext1, t)

			for i := range values1 {
				values1[i] += values2[i]
			}

			check(t, params.evaluator.Add(plaintext2, ciphertext1, ciphertext1))

			verify_test_vectors(params, params.decryptor, values1, ciphertext1, t)
		})

		t.Run(testString("CtPlainInPlaceNew/", params), func(t *testing.T) {

			values1, _, ciphertext1 := new_test_vectors(params, params.encryptorSk, 1, t)
			values2, plaintext2, _ := new_test_vectors(params, params.encryptorSk, 1, t)

			for i := range values1 {
				values1[i] += values2[i]
			}

			ciphertext3, err := params.evaluator.AddNew(ciphertext1, plaintext2)

			check(t, err)

			verify_test_vectors(params, params.decryptor, values1, ciphertext3, t)

			ciphertext3, err = params.evaluator.AddNew(plaintext2, ciphertext1)

			check(t, err)

			verify_test_vectors(params, params.decryptor, values1, ciphertext3, t)
		})
	}
}

func testEvaluatorSub(t *testing.T) {

	for _, parameters := range testParams.ckksParameters {

		params := genCkksParams(parameters)

		t.Run(testString("CtCtInPlace/", params), func(t *testing.T) {

			values1, _, ciphertext1 := new_test_vectors(params, params.encryptorSk, 1, t)
			values2, _, ciphertext2 := new_test_vectors(params, params.encryptorSk, 1, t)

			for i := range values1 {
				values1[i] -= values2[i]
			}

			check(t, params.evaluator.Sub(ciphertext1, ciphertext2, ciphertext1))

			verify_test_vectors(params, params.decryptor, values1, ciphertext1, t)
		})

		t.Run(testString("CtCtNew/", params), func(t *testing.T) {

			values1, _, ciphertext1 := new_test_vectors(params, params.encryptorSk, 1, t)
			values2, _, ciphertext2 := new_test_vectors(params, params.encryptorSk, 1, t)

			for i := range values1 {
				values1[i] -= values2[i]
			}

			ciphertext3, err := params.evaluator.SubNew(ciphertext1, ciphertext2)

			check(t, err)

			verify_test_vectors(params, params.decryptor, values1, ciphertext3, t)
		})

		t.Run(testString("CtPlainInPlace/", params), func(t *testing.T) {

			values1, _, ciphertext1 := new_test_vectors(params, params.encryptorSk, 1, t)
			values2, plaintext2, ciphertext2 := new_test_vectors(params, params.encryptorSk, 1, t)

			valuesTest := make([]complex128, len(values1))
			for i := range values1 {
				valuesTest[i] = values1[i] - values2[i]
			}

			check(t, params.evaluator.Sub(ciphertext1, plaintext2, ciphertext2))

			verify_test_vectors(params, params.decryptor, valuesTest, ciphertext2, t)

			for i := range values1 {
				valuesTest[i] = values2[i] - values1[i]
			}

			check(t, params.evaluator.Sub(plaintext2, ciphertext1, ciphertext2))

			verify_test_vectors(params, params.decryptor, valuesTest, ciphertext2, t)
		})

		t.Run(testString("CtPlainNew/", params), func(t *testing.T) {

			values1, _, ciphertext1 := new_test_vectors(params, params.encryptorSk, 1, t)
			values2, plaintext2, _ := new_test_vectors(params, params.encryptorSk, 1, t)

			valuesTest := make([]complex128, len(values1))
			for i := range values1 {
				valuesTest[i] = values1[i] - values2[i]
			}

			ciphertext3, err := params.evaluator.SubNew(ciphertext1, plaintext2)

			check(t, err)

			verify_test_vectors(params, params.decryptor, valuesTest, ciphertext3, t)

			for i := range values1 {
				valuesTest[i] = values2[i] - values1[i]
			}

			ciphertext3, err = params.evaluator.SubNew(plaintext2, ciphertext1)

			check(t, err)

			verify_test_vectors(params, params.decryptor, valuesTest, ciphertext3, t)
		})
	}
}

func testEvaluatorRescale(t *testing.T) {

	for _, parameters := range testParams.ckksParameters {

		params := genCkksParams(parameters)

		t.Run(testString("Single/", params), func(t *testing.T) {

			values, _, ciphertext := new_test_vectors(params, params.encryptorSk, 1, t)

			constant := params.ckksContext.moduli[ciphertext.Level()]

			check(t, params.evaluator.MultConst(ciphertext, constant, ciphertext))

			ciphertext.MulScale(float64(constant))

			check(t, params.evaluator.Rescale(ciphertext, params.ckksContext.scale, ciphertext))

			verify_test_vectors(params, params.decryptor, values, ciphertext, t)
		})

		t.Run(testString("Many/", params), func(t *testing.T) {

			values, _, ciphertext := new_test_vectors(params, params.encryptorSk, 1, t)

			nbRescales := uint64(2)

			for i := uint64(0); i < nbRescales; i++ {
				constant := params.ckksContext.moduli[ciphertext.Level()-i]
				check(t, params.evaluator.MultConst(ciphertext, constant, ciphertext))
				ciphertext.MulScale(float64(constant))
			}

			check(t, params.evaluator.RescaleMany(ciphertext, nbRescales, ciphertext))

			verify_test_vectors(params, params.decryptor, values, ciphertext, t)
		})
	}
}

func testEvaluatorAddConst(t *testing.T) {

	for _, parameters := range testParams.ckksParameters {

		params := genCkksParams(parameters)

		t.Run(testString("", params), func(t *testing.T) {

			values, _, ciphertext := new_test_vectors(params, params.encryptorSk, 1, t)

			constant := complex(3.1415, -1.4142)

			for i := range values {
				values[i] += constant
			}

			check(t, params.evaluator.AddConst(ciphertext, constant, ciphertext))

			verify_test_vectors(params, params.decryptor, values, ciphertext, t)
		})
	}
}

func testEvaluatorMultConst(t *testing.T) {

	for _, parameters := range testParams.ckksParameters {

		params := genCkksParams(parameters)

		t.Run(testString("", params), func(t *testing.T) {

			values, _, ciphertext := new_test_vectors(params, params.encryptorSk, 1, t)

			constant := 1.0 / complex(3.1415, -1.4142)

			for i := range values {
				values[i] *= constant
			}

			check(t, params.evaluator.MultConst(ciphertext, constant, ciphertext))

			verify_test_vectors(params, params.decryptor, values, ciphertext, t)
		})
	}
}

func testEvaluatorMultConstAndAdd(t *testing.T) {

	for _, parameters := range testParams.ckksParameters {

		params := genCkksParams(parameters)

		t.Run(testString("", params), func(t *testing.T) {

			values1, _, ciphertext1 := new_test_vectors(params, params.encryptorSk, 1, t)
			values2, _, ciphertext2 := new_test_vectors(params, params.encryptorSk, 1, t)

			constant := 1.0 / complex(3.1415, -1.4142)

			for i := range values1 {
				values2[i] += (constant * values1[i])
			}

			check(t, params.evaluator.MultByConstAndAdd(ciphertext1, constant, ciphertext2))

			verify_test_vectors(params, params.decryptor, values2, ciphertext2, t)
		})
	}
}

func testEvaluatorMul(t *testing.T) {

	for _, parameters := range testParams.ckksParameters {

		params := genCkksParams(parameters)

		t.Run(testString("CtCtInPlace/", params), func(t *testing.T) {

			values1, _, ciphertext1 := new_test_vectors(params, params.encryptorSk, 1, t)
			values2, _, ciphertext2 := new_test_vectors(params, params.encryptorSk, 1, t)

			for i := range values1 {
				values2[i] *= values1[i]
			}

			check(t, params.evaluator.MulRelin(ciphertext1, ciphertext2, nil, ciphertext2))

			verify_test_vectors(params, params.decryptor, values2, ciphertext2, t)
		})

		t.Run(testString("CtCtNew/", params), func(t *testing.T) {

			values1, _, ciphertext1 := new_test_vectors(params, params.encryptorSk, 1, t)
			values2, _, ciphertext2 := new_test_vectors(params, params.encryptorSk, 1, t)

			for i := range values1 {
				values2[i] *= values1[i]
			}

			var ciphertext3 *Ciphertext
			ciphertext3, err = params.evaluator.MulRelinNew(ciphertext1, ciphertext2, nil)

			verify_test_vectors(params, params.decryptor, values2, ciphertext3, t)
		})

		t.Run(testString("CtPlain/", params), func(t *testing.T) {

			values1, plaintext1, ciphertext1 := new_test_vectors(params, params.encryptorSk, 1, t)
			values2, plaintext2, ciphertext2 := new_test_vectors(params, params.encryptorSk, 1, t)

			for i := range values1 {
				values1[i] *= values1[i]
				values2[i] *= values2[i]
			}

			check(t, params.evaluator.MulRelin(ciphertext1, plaintext1, nil, ciphertext1))

			verify_test_vectors(params, params.decryptor, values1, ciphertext1, t)

			check(t, params.evaluator.MulRelin(plaintext2, ciphertext2, nil, ciphertext2))

			verify_test_vectors(params, params.decryptor, values2, ciphertext2, t)
		})

		t.Run(testString("Relinearize/", params), func(t *testing.T) {

			rlk := params.kgen.NewRelinKey(params.sk)

			values1, _, ciphertext1 := new_test_vectors(params, params.encryptorSk, 1, t)
			values2, _, ciphertext2 := new_test_vectors(params, params.encryptorSk, 1, t)

			for i := range values1 {
				values2[i] *= values1[i]
			}

			check(t, params.evaluator.MulRelin(ciphertext1, ciphertext2, nil, ciphertext2))

			check(t, params.evaluator.Relinearize(ciphertext2, rlk, ciphertext2))

			if ciphertext2.Degree() != 1 {
				t.Errorf("Relinearize error")
			}

			verify_test_vectors(params, params.decryptor, values2, ciphertext2, t)
		})
	}
}

func testFunctions(t *testing.T) {

	for _, parameters := range testParams.ckksParameters {

		params := genCkksParams(parameters)

		rlk := params.kgen.NewRelinKey(params.sk)

		t.Run(testString("PowerOf2/", params), func(t *testing.T) {

			values, _, ciphertext := new_test_vectors(params, params.encryptorSk, 1, t)

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

			check(t, params.evaluator.PowerOf2(ciphertext, n, rlk, ciphertext))

			verify_test_vectors(params, params.decryptor, valuesWant, ciphertext, t)
		})

		t.Run(testString("Power/", params), func(t *testing.T) {

			values, _, ciphertext := new_test_vectors(params, params.encryptorSk, 1, t)

			n := uint64(7)

			for i := range values {
				values[i] = cmplx.Pow(values[i], complex(float64(n), 0))
			}

			check(t, params.evaluator.Power(ciphertext, n, rlk, ciphertext))

			verify_test_vectors(params, params.decryptor, values, ciphertext, t)
		})

		if params.ckksContext.levels > 7 {
			t.Run(testString("Inverse/", params), func(t *testing.T) {

				values, _, ciphertext := new_test_vectors_reals(params, params.encryptorSk, 0.1, 1, t)

				n := uint64(7)

				for i := range values {
					values[i] = 1.0 / values[i]
				}

				ciphertext, err = params.evaluator.InverseNew(ciphertext, n, rlk)

				check(t, err)

				verify_test_vectors(params, params.decryptor, values, ciphertext, t)
			})
		}
	}
}

func testChebyshevInterpolator(t *testing.T) {

	for _, parameters := range testParams.ckksParameters {

		params := genCkksParams(parameters)

		rlk := params.kgen.NewRelinKey(params.sk)

		t.Run(testString("Sin/", params), func(t *testing.T) {

			values, _, ciphertext := new_test_vectors_reals(params, params.encryptorSk, -1, 1, t)

			cheby := Approximate(cmplx.Sin, complex(-1, 0), complex(1, 0), 16)

			for i := range values {
				values[i] = cmplx.Sin(values[i])
			}

			ciphertext, err = params.evaluator.EvaluateCheby(ciphertext, cheby, rlk)
			check(t, err)

			verify_test_vectors(params, params.decryptor, values, ciphertext, t)
		})
	}
}

func testSwitchKeys(t *testing.T) {

	for _, parameters := range testParams.ckksParameters {

		params := genCkksParams(parameters)

		sk2 := params.kgen.NewSecretKey()
		decryptorSk2, err := params.ckksContext.NewDecryptor(sk2)
		check(t, err)
		switchingKey, err := params.kgen.NewSwitchingKey(params.sk, sk2)
		check(t, err)

		t.Run(testString("InPlace/", params), func(t *testing.T) {

			values, _, ciphertext := new_test_vectors_reals(params, params.encryptorSk, -1, 1, t)

			check(t, params.evaluator.SwitchKeys(ciphertext, switchingKey, ciphertext))

			verify_test_vectors(params, decryptorSk2, values, ciphertext, t)
		})

		t.Run(testString("New/", params), func(t *testing.T) {

			values, _, ciphertext := new_test_vectors_reals(params, params.encryptorSk, -1, 1, t)

			ciphertext, err = params.evaluator.SwitchKeysNew(ciphertext, switchingKey)
			check(t, err)

			verify_test_vectors(params, decryptorSk2, values, ciphertext, t)
		})
	}
}

func testConjugate(t *testing.T) {

	for _, parameters := range testParams.ckksParameters {

		params := genCkksParams(parameters)

		rotKey := params.ckksContext.NewRotationKeys()
		params.kgen.GenRot(Conjugate, params.sk, 0, rotKey)

		t.Run(testString("InPlace/", params), func(t *testing.T) {

			values, _, ciphertext := new_test_vectors_reals(params, params.encryptorSk, -1, 1, t)

			for i := range values {
				values[i] = complex(real(values[i]), -imag(values[i]))
			}

			check(t, params.evaluator.Conjugate(ciphertext, rotKey, ciphertext))

			verify_test_vectors(params, params.decryptor, values, ciphertext, t)
		})

		t.Run(testString("New/", params), func(t *testing.T) {

			values, _, ciphertext := new_test_vectors_reals(params, params.encryptorSk, -1, 1, t)

			for i := range values {
				values[i] = complex(real(values[i]), -imag(values[i]))
			}

			ciphertext, err = params.evaluator.ConjugateNew(ciphertext, rotKey)
			check(t, err)

			verify_test_vectors(params, params.decryptor, values, ciphertext, t)
		})
	}
}

func testRotateColumns(t *testing.T) {

	for _, parameters := range testParams.ckksParameters {

		params := genCkksParams(parameters)

		rotKey := params.kgen.NewRotationKeysPow2(params.sk)

		t.Run(testString("InPlace/", params), func(t *testing.T) {

			values1, _, ciphertext1 := new_test_vectors_reals(params, params.encryptorSk, -1, 1, t)

			values2 := make([]complex128, len(values1))
			ciphertext2 := params.ckksContext.NewCiphertext(ciphertext1.Degree(), ciphertext1.Level(), ciphertext1.Scale())

			for n := 1; n < len(values1); n <<= 1 {

				// Applies the column rotation to the values
				for i := range values1 {
					values2[i] = values1[(i+n)%len(values1)]
				}

				check(t, params.evaluator.RotateColumns(ciphertext1, uint64(n), rotKey, ciphertext2))

				verify_test_vectors(params, params.decryptor, values2, ciphertext2, t)
			}

		})

		t.Run(testString("New/", params), func(t *testing.T) {

			values1, _, ciphertext1 := new_test_vectors_reals(params, params.encryptorSk, -1, 1, t)

			values2 := make([]complex128, len(values1))
			ciphertext2 := params.ckksContext.NewCiphertext(ciphertext1.Degree(), ciphertext1.Level(), ciphertext1.Scale())

			for n := 1; n < len(values1); n <<= 1 {

				// Applies the column rotation to the values
				for i := range values1 {
					values2[i] = values1[(i+n)%len(values1)]
				}

				ciphertext2, err = params.evaluator.RotateColumnsNew(ciphertext1, uint64(n), rotKey)

				verify_test_vectors(params, params.decryptor, values2, ciphertext2, t)
			}

		})

		t.Run(testString("Random/", params), func(t *testing.T) {

			values1, _, ciphertext1 := new_test_vectors_reals(params, params.encryptorSk, -1, 1, t)

			values2 := make([]complex128, len(values1))
			ciphertext2 := params.ckksContext.NewCiphertext(ciphertext1.Degree(), ciphertext1.Level(), ciphertext1.Scale())

			for n := 1; n < 4; n++ {

				rand := rand.Uint64() % uint64(len(values1))

				// Applies the column rotation to the values
				for i := range values1 {
					values2[i] = values1[(i+int(rand))%len(values1)]
				}

				check(t, params.evaluator.RotateColumns(ciphertext1, rand, rotKey, ciphertext2))

				verify_test_vectors(params, params.decryptor, values2, ciphertext2, t)
			}

		})
	}
}
