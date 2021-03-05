package ckks

import (
	"flag"
	"fmt"
	"math"
	"math/cmplx"
	"runtime"
	"testing"

	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/utils"
	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

var flagLongTest = flag.Bool("long", false, "run the long test suite (all parameters + secure bootstrapping). Overrides -short and requires -timeout=0.")
var printPrecisionStats = flag.Bool("print-precision", false, "print precision stats")
var testBootstrapping = flag.Bool("test-bootstrapping", false, "run the bootstrapping tests (memory intensive)")

var minPrec float64 = 15.0

func testString(testContext *testParams, opname string) string {
	return fmt.Sprintf("%slogN=%d/LogSlots=%d/logQP=%d/levels=%d/a=%d/b=%d",
		opname,
		testContext.params.LogN(),
		testContext.params.LogSlots(),
		testContext.params.LogQP(),
		testContext.params.MaxLevel()+1,
		testContext.params.Alpha(),
		testContext.params.Beta())
}

type testParams struct {
	params      *Parameters
	ringQ       *ring.Ring
	ringP       *ring.Ring
	ringQP      *ring.Ring
	prng        utils.PRNG
	encoder     Encoder
	kgen        KeyGenerator
	sk          *SecretKey
	pk          *PublicKey
	rlk         *RelinearizationKey
	encryptorPk Encryptor
	encryptorSk Encryptor
	decryptor   Decryptor
	evaluator   Evaluator
}

func TestCKKS(t *testing.T) {

	var err error

	var defaultParams = DefaultParams[PN12QP109 : PN12QP109+4] // the default test runs for ring degree N=2^12, 2^13, 2^14, 2^15

	if testing.Short() {
		defaultParams = DefaultParams[PN12QP109 : PN12QP109+1] // the short test suite runs for ring degree N=2^12, 2^13
	}

	if *flagLongTest {
		defaultParams = DefaultParams // the long test suite runs for all default parameters
	}

	for _, defaultParam := range defaultParams {
		var testContext *testParams
		if testContext, err = genTestParams(defaultParam, 0); err != nil {
			panic(err)
		}

		for _, testSet := range []func(testContext *testParams, t *testing.T){
			testParameters,
			testEncoder,
			testEncryptor,
			testEvaluatorAdd,
			testEvaluatorSub,
			testEvaluatorRescale,
			testEvaluatorAddConst,
			testEvaluatorMultByConst,
			testEvaluatorMultByConstAndAdd,
			testEvaluatorMul,
			testFunctions,
			testDecryptPublic,
			testEvaluatePoly,
			testChebyshevInterpolator,
			testSwitchKeys,
			testAutomorphisms,
			testMarshaller,
		} {
			testSet(testContext, t)
			runtime.GC()
		}
	}
}

func genTestParams(defaultParam *Parameters, hw uint64) (testContext *testParams, err error) {

	testContext = new(testParams)

	testContext.params = defaultParam.Copy()

	testContext.kgen = NewKeyGenerator(testContext.params)

	if hw == 0 {
		testContext.sk, testContext.pk = testContext.kgen.GenKeyPair()
	} else {
		testContext.sk, testContext.pk = testContext.kgen.GenKeyPairSparse(hw)
	}

	if testContext.ringQ, err = ring.NewRing(testContext.params.N(), testContext.params.qi); err != nil {
		return nil, err
	}

	if testContext.ringQP, err = ring.NewRing(testContext.params.N(), append(testContext.params.qi, testContext.params.pi...)); err != nil {
		return nil, err
	}

	if testContext.params.PiCount() != 0 {
		if testContext.ringP, err = ring.NewRing(testContext.params.N(), testContext.params.pi); err != nil {
			return nil, err
		}

		testContext.rlk = testContext.kgen.GenRelinearizationKey(testContext.sk)
	}

	if testContext.prng, err = utils.NewPRNG(); err != nil {
		return nil, err
	}

	testContext.encoder = NewEncoder(testContext.params)

	testContext.encryptorPk = NewEncryptorFromPk(testContext.params, testContext.pk)
	testContext.encryptorSk = NewEncryptorFromSk(testContext.params, testContext.sk)
	testContext.decryptor = NewDecryptor(testContext.params, testContext.sk)

	testContext.evaluator = NewEvaluator(testContext.params, EvaluationKey{testContext.rlk, nil})

	return testContext, nil

}

func newTestVectors(testContext *testParams, encryptor Encryptor, a, b complex128, t *testing.T) (values []complex128, plaintext *Plaintext, ciphertext *Ciphertext) {

	logSlots := testContext.params.LogSlots()

	values = make([]complex128, 1<<logSlots)

	for i := uint64(0); i < 1<<logSlots; i++ {
		values[i] = complex(utils.RandFloat64(real(a), real(b)), utils.RandFloat64(imag(a), imag(b)))
	}

	values[0] = complex(0.607538, 0)

	plaintext = testContext.encoder.EncodeNTTAtLvlNew(testContext.params.MaxLevel(), values, logSlots)

	if encryptor != nil {

		switch encryptor := encryptor.(type) {
		case *pkEncryptor:
			if testContext.params.PiCount() != 0 {
				ciphertext = encryptor.EncryptNew(plaintext)
			} else {
				ciphertext = encryptor.EncryptFastNew(plaintext)
			}
		case *skEncryptor:
			ciphertext = encryptor.EncryptNew(plaintext)
		}
	}

	return values, plaintext, ciphertext
}

func verifyTestVectors(testContext *testParams, decryptor Decryptor, valuesWant []complex128, element interface{}, t *testing.T, bound float64) {

	precStats := GetPrecisionStats(testContext.params, testContext.encoder, decryptor, valuesWant, element, bound)

	if *printPrecisionStats {
		t.Log(precStats.String())
	}

	require.GreaterOrEqual(t, real(precStats.MeanPrecision), minPrec)
	require.GreaterOrEqual(t, imag(precStats.MeanPrecision), minPrec)
}

func testParameters(testContext *testParams, t *testing.T) {

	t.Run("Parameters/NewParametersFromModuli/", func(t *testing.T) {
		p, err := NewParametersFromModuli(testContext.params.LogN(), testContext.params.Moduli())
		p.SetLogSlots(testContext.params.LogSlots())
		p.SetScale(testContext.params.Scale())
		assert.NoError(t, err)
		assert.True(t, p.Equals(testContext.params))
	})

	t.Run("Parameters/NewParametersFromLogModuli/", func(t *testing.T) {
		p, err := NewParametersFromLogModuli(testContext.params.LogN(), testContext.params.LogModuli())
		p.SetLogSlots(testContext.params.LogSlots())
		p.SetScale(testContext.params.Scale())
		assert.NoError(t, err)
		assert.True(t, p.Equals(testContext.params))
	})
}

func testEncoder(testContext *testParams, t *testing.T) {

	t.Run(testString(testContext, "Encoder/Slots/"), func(t *testing.T) {

		values, plaintext, _ := newTestVectors(testContext, nil, complex(-1, -1), complex(1, 1), t)

		verifyTestVectors(testContext, testContext.decryptor, values, plaintext, t, 0)
	})

	t.Run(testString(testContext, "Encoder/Coeffs/"), func(t *testing.T) {

		slots := testContext.params.N()

		valuesWant := make([]float64, slots)

		for i := uint64(0); i < slots; i++ {
			valuesWant[i] = utils.RandFloat64(-1, 1)
		}

		valuesWant[0] = 0.607538

		plaintext := NewPlaintext(testContext.params, testContext.params.MaxLevel(), testContext.params.Scale())

		testContext.encoder.EncodeCoeffs(valuesWant, plaintext)

		valuesTest := testContext.encoder.DecodeCoeffs(plaintext)

		var meanprec float64

		for i := range valuesWant {
			meanprec += math.Abs(valuesTest[i] - valuesWant[i])
		}

		meanprec /= float64(slots)

		if *printPrecisionStats {
			t.Log(fmt.Sprintf("\nMean    precision : %.2f \n", math.Log2(1/meanprec)))
		}

		require.GreaterOrEqual(t, math.Log2(1/meanprec), minPrec)
	})

}

func testEncryptor(testContext *testParams, t *testing.T) {

	t.Run(testString(testContext, "Encryptor/EncryptFromPk/Lvl=Max/"), func(t *testing.T) {

		if testContext.params.PiCount() == 0 {
			t.Skip("#Pi is empty")
		}

		values, _, ciphertext := newTestVectors(testContext, testContext.encryptorPk, complex(-1, -1), complex(1, 1), t)

		verifyTestVectors(testContext, testContext.decryptor, values, ciphertext, t, 0)
	})

	t.Run(testString(testContext, "Encryptor/EncryptFromPkFast/Lvl=Max/"), func(t *testing.T) {

		logSlots := testContext.params.LogSlots()

		values := make([]complex128, 1<<logSlots)

		for i := uint64(0); i < 1<<logSlots; i++ {
			values[i] = utils.RandComplex128(-1, 1)
		}

		values[0] = complex(0.607538, 0.555668)

		plaintext := testContext.encoder.EncodeNew(values, logSlots)

		verifyTestVectors(testContext, testContext.decryptor, values, testContext.encryptorPk.EncryptFastNew(plaintext), t, 0)
	})

	t.Run(testString(testContext, "Encryptor/EncryptFromSk/Lvl=Max/"), func(t *testing.T) {

		values, _, ciphertext := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)

		verifyTestVectors(testContext, testContext.decryptor, values, ciphertext, t, 0)
	})

	t.Run(testString(testContext, "Encryptor/EncryptFromPk/Lvl=1/"), func(t *testing.T) {

		if testContext.params.PiCount() == 0 {
			t.Skip("#Pi is empty")
		}

		if testContext.params.MaxLevel() < 1 {
			t.Skip("skipping test for params max level < 1")
		}

		logSlots := testContext.params.LogSlots()

		values := make([]complex128, 1<<logSlots)

		for i := uint64(0); i < 1<<logSlots; i++ {
			values[i] = utils.RandComplex128(-1, 1)
		}

		values[0] = complex(0.607538, 0.555668)

		plaintext := testContext.encoder.EncodeAtLvlNew(1, values, logSlots)

		verifyTestVectors(testContext, testContext.decryptor, values, testContext.encryptorPk.EncryptNew(plaintext), t, 0)
	})

	t.Run(testString(testContext, "Encryptor/EncryptFromPkFast/Lvl=1/"), func(t *testing.T) {

		if testContext.params.MaxLevel() < 1 {
			t.Skip("skipping test for params max level < 1")
		}

		logSlots := testContext.params.LogSlots()

		values := make([]complex128, 1<<logSlots)

		for i := uint64(0); i < 1<<logSlots; i++ {
			values[i] = utils.RandComplex128(-1, 1)
		}

		values[0] = complex(0.607538, 0.555668)

		plaintext := testContext.encoder.EncodeAtLvlNew(1, values, logSlots)

		verifyTestVectors(testContext, testContext.decryptor, values, testContext.encryptorPk.EncryptFastNew(plaintext), t, 0)
	})

	t.Run(testString(testContext, "Encryptor/EncryptFromSk/Lvl=1/"), func(t *testing.T) {

		if testContext.params.MaxLevel() < 1 {
			t.Skip("skipping test for params max level < 1")
		}

		logSlots := testContext.params.LogSlots()

		values := make([]complex128, 1<<logSlots)

		for i := uint64(0); i < 1<<logSlots; i++ {
			values[i] = utils.RandComplex128(-1, 1)
		}

		values[0] = complex(0.607538, 0.555668)

		plaintext := testContext.encoder.EncodeAtLvlNew(1, values, logSlots)

		verifyTestVectors(testContext, testContext.decryptor, values, testContext.encryptorSk.EncryptNew(plaintext), t, 0)
	})

}

func testEvaluatorAdd(testContext *testParams, t *testing.T) {

	t.Run(testString(testContext, "EvaluatorAdd/CtCtInPlace/"), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)
		values2, _, ciphertext2 := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)

		for i := range values1 {
			values1[i] += values2[i]
		}

		testContext.evaluator.Add(ciphertext1, ciphertext2, ciphertext1)

		verifyTestVectors(testContext, testContext.decryptor, values1, ciphertext1, t, 0)
	})

	t.Run(testString(testContext, "EvaluatorAdd/CtCtNew/"), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)
		values2, _, ciphertext2 := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)

		for i := range values1 {
			values1[i] += values2[i]
		}

		ciphertext3 := testContext.evaluator.AddNew(ciphertext1, ciphertext2)

		verifyTestVectors(testContext, testContext.decryptor, values1, ciphertext3, t, 0)
	})

	t.Run(testString(testContext, "EvaluatorAdd/CtPlainInPlace/"), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)
		values2, plaintext2, _ := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)

		for i := range values1 {
			values1[i] += values2[i]
		}

		testContext.evaluator.Add(ciphertext1, plaintext2, ciphertext1)

		verifyTestVectors(testContext, testContext.decryptor, values1, ciphertext1, t, 0)

		for i := range values1 {
			values1[i] += values2[i]
		}

		testContext.evaluator.Add(plaintext2, ciphertext1, ciphertext1)

		verifyTestVectors(testContext, testContext.decryptor, values1, ciphertext1, t, 0)
	})

	t.Run(testString(testContext, "EvaluatorAdd/CtPlainInPlaceNew/"), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)
		values2, plaintext2, _ := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)

		for i := range values1 {
			values1[i] += values2[i]
		}

		ciphertext3 := testContext.evaluator.AddNew(ciphertext1, plaintext2)

		verifyTestVectors(testContext, testContext.decryptor, values1, ciphertext3, t, 0)

		ciphertext3 = testContext.evaluator.AddNew(plaintext2, ciphertext1)

		verifyTestVectors(testContext, testContext.decryptor, values1, ciphertext3, t, 0)
	})

}

func testEvaluatorSub(testContext *testParams, t *testing.T) {

	t.Run(testString(testContext, "EvaluatorSub/CtCtInPlace/"), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)
		values2, _, ciphertext2 := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)

		for i := range values1 {
			values1[i] -= values2[i]
		}

		testContext.evaluator.Sub(ciphertext1, ciphertext2, ciphertext1)

		verifyTestVectors(testContext, testContext.decryptor, values1, ciphertext1, t, 0)
	})

	t.Run(testString(testContext, "EvaluatorSub/CtCtNew/"), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)
		values2, _, ciphertext2 := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)

		for i := range values1 {
			values1[i] -= values2[i]
		}

		ciphertext3 := testContext.evaluator.SubNew(ciphertext1, ciphertext2)

		verifyTestVectors(testContext, testContext.decryptor, values1, ciphertext3, t, 0)
	})

	t.Run(testString(testContext, "EvaluatorSub/CtPlainInPlace/"), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)
		values2, plaintext2, ciphertext2 := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)

		valuesTest := make([]complex128, len(values1))
		for i := range values1 {
			valuesTest[i] = values1[i] - values2[i]
		}

		testContext.evaluator.Sub(ciphertext1, plaintext2, ciphertext2)

		verifyTestVectors(testContext, testContext.decryptor, valuesTest, ciphertext2, t, 0)

		for i := range values1 {
			valuesTest[i] = values2[i] - values1[i]
		}

		testContext.evaluator.Sub(plaintext2, ciphertext1, ciphertext2)

		verifyTestVectors(testContext, testContext.decryptor, valuesTest, ciphertext2, t, 0)
	})

	t.Run(testString(testContext, "EvaluatorSub/CtPlainNew/"), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)
		values2, plaintext2, _ := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)

		valuesTest := make([]complex128, len(values1))
		for i := range values1 {
			valuesTest[i] = values1[i] - values2[i]
		}

		ciphertext3 := testContext.evaluator.SubNew(ciphertext1, plaintext2)

		verifyTestVectors(testContext, testContext.decryptor, valuesTest, ciphertext3, t, 0)

		for i := range values1 {
			valuesTest[i] = values2[i] - values1[i]
		}

		ciphertext3 = testContext.evaluator.SubNew(plaintext2, ciphertext1)

		verifyTestVectors(testContext, testContext.decryptor, valuesTest, ciphertext3, t, 0)
	})

}

func testEvaluatorRescale(testContext *testParams, t *testing.T) {

	t.Run(testString(testContext, "Evaluator/Rescale/Single/"), func(t *testing.T) {

		if testContext.params.PiCount() == 0 {
			t.Skip("#Pi is empty")
		}

		values, _, ciphertext := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)

		constant := testContext.ringQ.Modulus[ciphertext.Level()]

		testContext.evaluator.MultByConst(ciphertext, constant, ciphertext)

		ciphertext.MulScale(float64(constant))

		testContext.evaluator.Rescale(ciphertext, testContext.params.Scale(), ciphertext)

		verifyTestVectors(testContext, testContext.decryptor, values, ciphertext, t, 0)
	})

	t.Run(testString(testContext, "Evaluator/Rescale/Many/"), func(t *testing.T) {

		if testContext.params.PiCount() == 0 {
			t.Skip("#Pi is empty")
		}

		values, _, ciphertext := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)

		nbRescales := testContext.params.MaxLevel()
		if nbRescales > 5 {
			nbRescales = 5
		}

		for i := uint64(0); i < nbRescales; i++ {
			constant := testContext.ringQ.Modulus[ciphertext.Level()-i]
			testContext.evaluator.MultByConst(ciphertext, constant, ciphertext)
			ciphertext.MulScale(float64(constant))
		}

		testContext.evaluator.Rescale(ciphertext, testContext.params.Scale(), ciphertext)

		verifyTestVectors(testContext, testContext.decryptor, values, ciphertext, t, 0)
	})
}

func testEvaluatorAddConst(testContext *testParams, t *testing.T) {

	t.Run(testString(testContext, "Evaluator/AddConst/"), func(t *testing.T) {

		values, _, ciphertext := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)

		constant := complex(3.1415, -1.4142)

		for i := range values {
			values[i] += constant
		}

		testContext.evaluator.AddConst(ciphertext, constant, ciphertext)

		verifyTestVectors(testContext, testContext.decryptor, values, ciphertext, t, 0)
	})

}

func testEvaluatorMultByConst(testContext *testParams, t *testing.T) {

	t.Run(testString(testContext, "Evaluator/MultByConst/"), func(t *testing.T) {

		values, _, ciphertext := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)

		constant := 1.0 / complex(3.1415, -1.4142)

		for i := range values {
			values[i] *= constant
		}

		testContext.evaluator.MultByConst(ciphertext, constant, ciphertext)

		verifyTestVectors(testContext, testContext.decryptor, values, ciphertext, t, 0)
	})

}

func testEvaluatorMultByConstAndAdd(testContext *testParams, t *testing.T) {

	t.Run(testString(testContext, "Evaluator/MultByConstAndAdd/"), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)
		values2, _, ciphertext2 := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)

		constant := 1.0 / complex(3.1415, -1.4142)

		for i := range values1 {
			values2[i] += (constant * values1[i])
		}

		testContext.evaluator.MultByConstAndAdd(ciphertext1, constant, ciphertext2)

		verifyTestVectors(testContext, testContext.decryptor, values2, ciphertext2, t, 0)
	})

}

func testEvaluatorMul(testContext *testParams, t *testing.T) {

	t.Run(testString(testContext, "Evaluator/Mul/ct0*pt->ct0/"), func(t *testing.T) {

		values1, plaintext1, ciphertext1 := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)

		for i := range values1 {
			values1[i] *= values1[i]
		}

		testContext.evaluator.MulRelin(ciphertext1, plaintext1, ciphertext1)

		verifyTestVectors(testContext, testContext.decryptor, values1, ciphertext1, t, 0)
	})

<<<<<<< HEAD
	t.Run(testString(testContext, "Evaluator/Mul/pt*ct0->ct0/"), func(t *testing.T) {
=======
	t.Run(testString(testContext, "Evaluato/rMul/pt*ct0->ct0/"), func(t *testing.T) {
>>>>>>> dev_rlwe_layer

		values1, plaintext1, ciphertext1 := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)

		for i := range values1 {
			values1[i] *= values1[i]
		}

		testContext.evaluator.MulRelin(ciphertext1, plaintext1, ciphertext1)

		verifyTestVectors(testContext, testContext.decryptor, values1, ciphertext1, t, 0)
	})

	t.Run(testString(testContext, "Evaluator/Mul/ct0*pt->ct1/"), func(t *testing.T) {

		values1, plaintext1, ciphertext1 := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)

		for i := range values1 {
			values1[i] *= values1[i]
		}

		ciphertext2 := testContext.evaluator.MulRelinNew(ciphertext1, plaintext1)

		verifyTestVectors(testContext, testContext.decryptor, values1, ciphertext2, t, 0)
	})

	t.Run(testString(testContext, "Evaluator/Mul/ct0*ct1->ct0/"), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)
		values2, _, ciphertext2 := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)

		for i := range values1 {
			values2[i] *= values1[i]
		}

		testContext.evaluator.MulRelin(ciphertext1, ciphertext2, ciphertext1)

		verifyTestVectors(testContext, testContext.decryptor, values2, ciphertext1, t, 0)
	})

	t.Run(testString(testContext, "Evaluator/Mul/ct0*ct1->ct1/"), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)
		values2, _, ciphertext2 := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)

		for i := range values1 {
			values2[i] *= values1[i]
		}

		testContext.evaluator.MulRelin(ciphertext1, ciphertext2, ciphertext2)

		verifyTestVectors(testContext, testContext.decryptor, values2, ciphertext2, t, 0)
	})

	t.Run(testString(testContext, "Evaluator/Mul/ct0*ct1->ct2/"), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)
		values2, _, ciphertext2 := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)

		for i := range values1 {
			values2[i] *= values1[i]
		}

		ciphertext3 := testContext.evaluator.MulRelinNew(ciphertext1, ciphertext2)

		verifyTestVectors(testContext, testContext.decryptor, values2, ciphertext3, t, 0)
	})

	t.Run(testString(testContext, "Evaluator/Mul/ct0*ct0->ct0/"), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)

		for i := range values1 {
			values1[i] *= values1[i]
		}

		testContext.evaluator.MulRelin(ciphertext1, ciphertext1, ciphertext1)

		verifyTestVectors(testContext, testContext.decryptor, values1, ciphertext1, t, 0)
	})

	t.Run(testString(testContext, "Evaluator/Mul/ct0*ct0->ct1/"), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)

		for i := range values1 {
			values1[i] *= values1[i]
		}

		ciphertext2 := testContext.evaluator.MulRelinNew(ciphertext1, ciphertext1)

		verifyTestVectors(testContext, testContext.decryptor, values1, ciphertext2, t, 0)
	})

	t.Run(testString(testContext, "Evaluator/Mul/Relinearize(ct0*ct1->ct0)/"), func(t *testing.T) {

		if testContext.params.PiCount() == 0 {
			t.Skip("#Pi is empty")
		}

		values1, _, ciphertext1 := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)
		values2, _, ciphertext2 := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)

		for i := range values1 {
			values1[i] *= values2[i]
		}

		testContext.evaluator.Mul(ciphertext1, ciphertext2, ciphertext1)
<<<<<<< HEAD

=======
		require.Equal(t, ciphertext1.Degree(), uint64(2))
>>>>>>> dev_rlwe_layer
		testContext.evaluator.Relinearize(ciphertext1, ciphertext1)
		require.Equal(t, ciphertext1.Degree(), uint64(1))

		verifyTestVectors(testContext, testContext.decryptor, values1, ciphertext1, t, 0)
	})

	t.Run(testString(testContext, "Evaluator/Mul/Relinearize(ct0*ct1->ct1)/"), func(t *testing.T) {

		if testContext.params.PiCount() == 0 {
			t.Skip("#Pi is empty")
		}

		values1, _, ciphertext1 := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)
		values2, _, ciphertext2 := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)

		for i := range values1 {
			values2[i] *= values1[i]
		}

		testContext.evaluator.Mul(ciphertext1, ciphertext2, ciphertext2)
<<<<<<< HEAD

=======
		require.Equal(t, ciphertext2.Degree(), uint64(2))
>>>>>>> dev_rlwe_layer
		testContext.evaluator.Relinearize(ciphertext2, ciphertext2)
		require.Equal(t, ciphertext2.Degree(), uint64(1))

		verifyTestVectors(testContext, testContext.decryptor, values2, ciphertext2, t, 0)
	})

}

func testFunctions(testContext *testParams, t *testing.T) {

	t.Run(testString(testContext, "Functions/PowerOf2/"), func(t *testing.T) {

		if testContext.params.PiCount() == 0 {
			t.Skip("#Pi is empty")
		}

		if testContext.params.MaxLevel() < 3 {
			t.Skip("skipping test for params max level < 3")
		}

		values, _, ciphertext := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)

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

		testContext.evaluator.PowerOf2(ciphertext, n, ciphertext)

		verifyTestVectors(testContext, testContext.decryptor, valuesWant, ciphertext, t, 0)
	})

	t.Run(testString(testContext, "Functions/Power/"), func(t *testing.T) {

		if testContext.params.PiCount() == 0 {
			t.Skip("#Pi is empty")
		}

		if testContext.params.MaxLevel() < 4 {
			t.Skip("skipping test for params max level < 4")
		}

		values, _, ciphertext := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)

		n := uint64(3)

		for i := range values {
			values[i] = cmplx.Pow(values[i], complex(float64(n), 0))
		}

		testContext.evaluator.Power(ciphertext, n, ciphertext)

		verifyTestVectors(testContext, testContext.decryptor, values, ciphertext, t, 0)
	})

	t.Run(testString(testContext, "Functions/Inverse/"), func(t *testing.T) {

		if testContext.params.PiCount() == 0 {
			t.Skip("#Pi is empty")
		}

		if testContext.params.MaxLevel() < 7 {
			t.Skip("skipping test for params max level < 7")
		}

		values, _, ciphertext := newTestVectors(testContext, testContext.encryptorSk, complex(0.1, 0), complex(1, 0), t)

		n := uint64(7)

		for i := range values {
			values[i] = 1.0 / values[i]
		}

		ciphertext = testContext.evaluator.InverseNew(ciphertext, n)

		verifyTestVectors(testContext, testContext.decryptor, values, ciphertext, t, 0)
	})
}

func testEvaluatePoly(testContext *testParams, t *testing.T) {

	var err error

	t.Run(testString(testContext, "EvaluatePoly/Exp/"), func(t *testing.T) {

		if testContext.params.PiCount() == 0 {
			t.Skip("#Pi is empty")
		}

		if testContext.params.MaxLevel() < 3 {
			t.Skip("skipping test for params max level < 3")
		}

		values, _, ciphertext := newTestVectors(testContext, testContext.encryptorSk, complex(-1, 0), complex(1, 0), t)

		coeffs := []complex128{
			complex(1.0, 0),
			complex(1.0, 0),
			complex(1.0/2, 0),
			complex(1.0/6, 0),
			complex(1.0/24, 0),
			complex(1.0/120, 0),
			complex(1.0/720, 0),
			complex(1.0/5040, 0),
		}

		poly := NewPoly(coeffs)

		for i := range values {
			values[i] = cmplx.Exp(values[i])
		}

		if ciphertext, err = testContext.evaluator.EvaluatePoly(ciphertext, poly, ciphertext.Scale()); err != nil {
			t.Error(err)
		}

		verifyTestVectors(testContext, testContext.decryptor, values, ciphertext, t, 0)
	})
}

func testChebyshevInterpolator(testContext *testParams, t *testing.T) {

	var err error

	t.Run(testString(testContext, "ChebyshevInterpolator/Sin/"), func(t *testing.T) {

		if testContext.params.PiCount() == 0 {
			t.Skip("#Pi is empty")
		}

		if testContext.params.MaxLevel() < 5 {
			t.Skip("skipping test for params max level < 5")
		}

		eval := testContext.evaluator

		values, _, ciphertext := newTestVectors(testContext, testContext.encryptorSk, complex(-1, 0), complex(1, 0), t)

		cheby := Approximate(cmplx.Sin, complex(-1.5, 0), complex(1.5, 0), 15)

		for i := range values {
			values[i] = cmplx.Sin(values[i])
		}

		eval.MultByConst(ciphertext, 2/(cheby.b-cheby.a), ciphertext)
		eval.AddConst(ciphertext, (-cheby.a-cheby.b)/(cheby.b-cheby.a), ciphertext)
		eval.Rescale(ciphertext, eval.(*evaluator).scale, ciphertext)

		if ciphertext, err = eval.EvaluateCheby(ciphertext, cheby, ciphertext.Scale()); err != nil {
			t.Error(err)
		}

		verifyTestVectors(testContext, testContext.decryptor, values, ciphertext, t, 0)
	})
}

func testDecryptPublic(testContext *testParams, t *testing.T) {

	var err error

	t.Run(testString(testContext, "DecryptPublic/Sin/"), func(t *testing.T) {

		if testContext.params.PiCount() == 0 {
			t.Skip("#Pi is empty")
		}

		if testContext.params.MaxLevel() < 5 {
			t.Skip("skipping test for params max level < 5")
		}

		eval := testContext.evaluator

		values, _, ciphertext := newTestVectors(testContext, testContext.encryptorSk, complex(-1, 0), complex(1, 0), t)

		cheby := Approximate(cmplx.Sin, complex(-1.5, 0), complex(1.5, 0), 15)

		for i := range values {
			values[i] = cmplx.Sin(values[i])
		}

		eval.MultByConst(ciphertext, 2/(cheby.b-cheby.a), ciphertext)
		eval.AddConst(ciphertext, (-cheby.a-cheby.b)/(cheby.b-cheby.a), ciphertext)
		eval.Rescale(ciphertext, eval.(*evaluator).scale, ciphertext)

		if ciphertext, err = eval.EvaluateCheby(ciphertext, cheby, ciphertext.Scale()); err != nil {
			t.Error(err)
		}

		plaintext := testContext.decryptor.DecryptNew(ciphertext)

		valuesHave := testContext.encoder.Decode(plaintext, testContext.params.LogSlots())

		verifyTestVectors(testContext, nil, values, valuesHave, t, 0)

		sigma := testContext.encoder.GetErrSTDTimeDom(values, valuesHave, plaintext.Scale())

		valuesHave = testContext.encoder.DecodePublic(plaintext, testContext.params.LogSlots(), sigma)

		verifyTestVectors(testContext, nil, values, valuesHave, t, 0)
	})
}

func testSwitchKeys(testContext *testParams, t *testing.T) {

	var sk2 *SecretKey
	var decryptorSk2 Decryptor
	var switchingKey *SwitchingKey

	if testContext.params.PiCount() != 0 {
		sk2 = testContext.kgen.GenSecretKey()
		decryptorSk2 = NewDecryptor(testContext.params, sk2)
		switchingKey = testContext.kgen.GenSwitchingKey(testContext.sk, sk2)
	}

	t.Run(testString(testContext, "SwitchKeys/InPlace/"), func(t *testing.T) {

		if testContext.params.PiCount() == 0 {
			t.Skip("#Pi is empty")
		}

		values, _, ciphertext := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)

		testContext.evaluator.SwitchKeys(ciphertext, switchingKey, ciphertext)

		verifyTestVectors(testContext, decryptorSk2, values, ciphertext, t, 0)
	})

	t.Run(testString(testContext, "SwitchKeys/New/"), func(t *testing.T) {

		if testContext.params.PiCount() == 0 {
			t.Skip("#Pi is empty")
		}

		values, _, ciphertext := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)

		ciphertext = testContext.evaluator.SwitchKeysNew(ciphertext, switchingKey)

		verifyTestVectors(testContext, decryptorSk2, values, ciphertext, t, 0)
	})

}

func testAutomorphisms(testContext *testParams, t *testing.T) {

	if testContext.params.PiCount() == 0 {
		t.Skip("#Pi is empty")
	}
	rots := []int{0, 1, -1, 4, -4, 63, -63}
	rotKey := testContext.kgen.GenRotationKeysForRotations(rots, true, testContext.sk)
	evaluator := testContext.evaluator.ShallowCopyWithKey(EvaluationKey{testContext.rlk, rotKey})

	t.Run(testString(testContext, "Conjugate/InPlace/"), func(t *testing.T) {

		if testContext.params.PiCount() == 0 {
			t.Skip("#Pi is empty")
		}

		values, _, ciphertext := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)

		for i := range values {
			values[i] = complex(real(values[i]), -imag(values[i]))
		}

		evaluator.Conjugate(ciphertext, ciphertext)

		verifyTestVectors(testContext, testContext.decryptor, values, ciphertext, t, 0)
	})

	t.Run(testString(testContext, "Conjugate/New/"), func(t *testing.T) {

		if testContext.params.PiCount() == 0 {
			t.Skip("#Pi is empty")
		}

		values, _, ciphertext := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)

		for i := range values {
			values[i] = complex(real(values[i]), -imag(values[i]))
		}

		ciphertext = evaluator.ConjugateNew(ciphertext)

		verifyTestVectors(testContext, testContext.decryptor, values, ciphertext, t, 0)
	})

	t.Run(testString(testContext, "RotateColumns/InPlace/"), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)

		ciphertext2 := NewCiphertext(testContext.params, ciphertext1.Degree(), ciphertext1.Level(), ciphertext1.Scale())

		for _, n := range rots {
			evaluator.Rotate(ciphertext1, n, ciphertext2)
			verifyTestVectors(testContext, testContext.decryptor, utils.RotateComplex128Slice(values1, n), ciphertext2, t, 0)
		}
	})

	t.Run(testString(testContext, "RotateColumns/New/"), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)

		for _, n := range rots {
			verifyTestVectors(testContext, testContext.decryptor, utils.RotateComplex128Slice(values1, n), evaluator.RotateNew(ciphertext1, n), t, 0)
		}

	})

	t.Run(testString(testContext, "RotateHoisted/"), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)

		ciphertexts := evaluator.RotateHoisted(ciphertext1, rots)

		for _, n := range rots {
			verifyTestVectors(testContext, testContext.decryptor, utils.RotateComplex128Slice(values1, n), ciphertexts[n], t, 0)
		}
	})
}

func testMarshaller(testContext *testParams, t *testing.T) {

	ringQP := testContext.ringQP

	t.Run("Marshaller/Ciphertext/", func(t *testing.T) {
		t.Run(testString(testContext, "EndToEnd/"), func(t *testing.T) {

			ciphertextWant := NewCiphertextRandom(testContext.prng, testContext.params, 2, testContext.params.MaxLevel(), testContext.params.Scale())

			marshalledCiphertext, err := ciphertextWant.MarshalBinary()
			require.NoError(t, err)

			ciphertextTest := new(Ciphertext)
			require.NoError(t, ciphertextTest.UnmarshalBinary(marshalledCiphertext))

			require.Equal(t, ciphertextWant.Degree(), ciphertextTest.Degree())
			require.Equal(t, ciphertextWant.Level(), ciphertextTest.Level())
			require.Equal(t, ciphertextWant.Scale(), ciphertextTest.Scale())

			for i := range ciphertextWant.value {
				require.True(t, testContext.ringQ.EqualLvl(ciphertextWant.Level(), ciphertextWant.Value()[i], ciphertextTest.Value()[i]))
			}
		})

		t.Run(testString(testContext, "Minimal/"), func(t *testing.T) {

			ciphertext := NewCiphertextRandom(testContext.prng, testContext.params, 0, testContext.params.MaxLevel(), testContext.params.Scale())

			marshalledCiphertext, err := ciphertext.MarshalBinary()
			require.NoError(t, err)

			ciphertextTest := new(Ciphertext)
			require.Error(t, ciphertextTest.UnmarshalBinary(nil))
			require.NoError(t, ciphertextTest.UnmarshalBinary(marshalledCiphertext))

			require.Equal(t, ciphertext.Degree(), uint64(0))
			require.Equal(t, ciphertext.Level(), testContext.params.MaxLevel())
			require.Equal(t, ciphertext.Scale(), testContext.params.Scale())
			require.Equal(t, len(ciphertext.Value()), 1)
		})
	})

	t.Run(testString(testContext, "Marshaller/Sk/"), func(t *testing.T) {

		marshalledSk, err := testContext.sk.MarshalBinary()
		require.NoError(t, err)

		sk := new(SecretKey)
		err = sk.UnmarshalBinary(marshalledSk)
		require.NoError(t, err)

		require.True(t, ringQP.Equal(sk.Value, testContext.sk.Value))

	})

	t.Run(testString(testContext, "Marshaller/Pk/"), func(t *testing.T) {

		marshalledPk, err := testContext.pk.MarshalBinary()
		require.NoError(t, err)

		pk := new(PublicKey)
		err = pk.UnmarshalBinary(marshalledPk)
		require.NoError(t, err)

		for k := range testContext.pk.Value {
			require.Truef(t, ringQP.Equal(pk.Value[k], testContext.pk.Value[k]), "Marshal PublicKey element [%d]", k)
		}
	})

	t.Run(testString(testContext, "Marshaller/EvaluationKey/"), func(t *testing.T) {

		if testContext.params.PiCount() == 0 {
			t.Skip("#Pi is empty")
		}

		evalKey := testContext.kgen.GenRelinearizationKey(testContext.sk)
		data, err := evalKey.MarshalBinary()
		require.NoError(t, err)

		resEvalKey := new(RelinearizationKey)
		err = resEvalKey.UnmarshalBinary(data)
		require.NoError(t, err)

		evakeyWant := evalKey.Keys[0].Value
		evakeyTest := resEvalKey.Keys[0].Value

		for j := range evakeyWant {
			for k := range evakeyWant[j] {
				require.Truef(t, ringQP.Equal(evakeyWant[j][k], evakeyTest[j][k]), "Marshal EvaluationKey element [%d][%d]", j, k)
			}
		}
	})

	t.Run(testString(testContext, "Marshaller/SwitchingKey/"), func(t *testing.T) {

		if testContext.params.PiCount() == 0 {
			t.Skip("#Pi is empty")
		}

		skOut := testContext.kgen.GenSecretKey()

		switchingKey := testContext.kgen.GenSwitchingKey(testContext.sk, skOut)
		data, err := switchingKey.MarshalBinary()
		require.NoError(t, err)

		resSwitchingKey := new(SwitchingKey)
		err = resSwitchingKey.UnmarshalBinary(data)
		require.NoError(t, err)

		evakeyWant := switchingKey.Value
		evakeyTest := resSwitchingKey.Value

		for j := range evakeyWant {
			for k := range evakeyWant[j] {
				require.True(t, ringQP.Equal(evakeyWant[j][k], evakeyTest[j][k]))
			}
		}
	})

	t.Run(testString(testContext, "Marshaller/RotationKey/"), func(t *testing.T) {

		if testContext.params.PiCount() == 0 {
			t.Skip("#Pi is empty")
		}

		rots := []int{1, -1, 63, -63}
		galEls := []uint64{testContext.params.GaloisElementForRowRotation()}
		for _, n := range rots {
			galEls = append(galEls, testContext.params.GaloisElementForColumnRotationBy(n))
		}

		rotationKey := testContext.kgen.GenRotationKeys(galEls, testContext.sk)

		data, err := rotationKey.MarshalBinary()
		require.NoError(t, err)

		resRotationKey := new(RotationKeySet)
		err = resRotationKey.UnmarshalBinary(data)
		require.NoError(t, err)

		for _, galEl := range galEls {

			evakeyWant := rotationKey.Keys[galEl].Value
			evakeyTest := resRotationKey.Keys[galEl].Value

			for j := range evakeyWant {
				for k := range evakeyWant[j] {
					require.Truef(t, ringQP.Equal(evakeyWant[j][k], evakeyTest[j][k]), "Marshal RotationKey RotateLeft %d element [%d][%d]", galEl, j, k)
				}
			}
		}
	})
}
