package ckks

import (
	"encoding/json"
	"flag"
	"fmt"
	"math"
	"math/cmplx"
	"runtime"
	"testing"

	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/rlwe"
	"github.com/ldsec/lattigo/v2/utils"
	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

var flagLongTest = flag.Bool("long", false, "run the long test suite (all parameters + secure bootstrapping). Overrides -short and requires -timeout=0.")
var flagPostQuantum = flag.Bool("pq", false, "run post quantum test suite (does not run non-PQ parameters).")
var flagParamString = flag.String("params", "", "specify the test cryptographic parameters as a JSON string. Overrides -short and -long.")
var printPrecisionStats = flag.Bool("print-precision", false, "print precision stats")

var minPrec float64 = 15.0

func GetTestName(params Parameters, opname string) string {
	return fmt.Sprintf("RingType=%s/logN=%d/logQP=%d/LogSlots=%d/levels=%d/alpha=%d/beta=%d/%s",
		params.RingType(),
		params.LogN(),
		params.LogQP(),
		params.LogSlots(),
		params.MaxLevel()+1,
		params.PCount(),
		params.Beta(),
		opname)
}

type testContext struct {
	params      Parameters
	ringQ       *ring.Ring
	ringP       *ring.Ring
	prng        utils.PRNG
	encoder     Encoder
	kgen        rlwe.KeyGenerator
	sk          *rlwe.SecretKey
	pk          *rlwe.PublicKey
	rlk         *rlwe.RelinearizationKey
	encryptorPk Encryptor
	encryptorSk Encryptor
	decryptor   Decryptor
	evaluator   Evaluator
}

func TestCKKS(t *testing.T) {

	var testParams []ParametersLiteral
	switch {
	case *flagParamString != "": // the custom test suite reads the parameters from the -params flag
		testParams = append(testParams, ParametersLiteral{})
		json.Unmarshal([]byte(*flagParamString), &testParams[0])
	case *flagLongTest:
		for _, pls := range [][]ParametersLiteral{
			DefaultParams,
			DefaultConjugateInvariantParams,
			DefaultPostQuantumParams,
			DefaultPostQuantumConjugateInvariantParams} {
			testParams = append(testParams, pls...)
		}
	case *flagPostQuantum && testing.Short():
		testParams = append(DefaultPostQuantumParams[:2], DefaultPostQuantumConjugateInvariantParams[:2]...)
	case *flagPostQuantum:
		testParams = append(DefaultPostQuantumParams[:4], DefaultPostQuantumConjugateInvariantParams[:4]...)
	case testing.Short():
		testParams = append(DefaultParams[:2], DefaultConjugateInvariantParams[:2]...)
	default:
		testParams = append(DefaultParams[:4], DefaultConjugateInvariantParams[:4]...)
	}

	for _, paramsLiteral := range testParams[:] {

		params, err := NewParametersFromLiteral(paramsLiteral)

		if err != nil {
			panic(err)
		}
		var tc *testContext
		if tc, err = genTestParams(params, 0); err != nil {
			panic(err)
		}

		for _, testSet := range []func(tc *testContext, t *testing.T){
			testParameters,
			testEncoder,
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
			testBridge,
			testAutomorphisms,
			testInnerSum,
			testReplicate,
			testLinearTransform,
			testMarshaller,
		} {
			testSet(tc, t)
			runtime.GC()
		}
	}
}

func genTestParams(defaultParam Parameters, hw int) (tc *testContext, err error) {

	tc = new(testContext)

	tc.params = defaultParam

	tc.kgen = NewKeyGenerator(tc.params)

	if hw == 0 {
		tc.sk, tc.pk = tc.kgen.GenKeyPair()
	} else {
		tc.sk, tc.pk = tc.kgen.GenKeyPairSparse(hw)
	}

	tc.ringQ = defaultParam.RingQ()
	if tc.params.PCount() != 0 {
		tc.ringP = defaultParam.RingP()
		tc.rlk = tc.kgen.GenRelinearizationKey(tc.sk, 2)
	}

	if tc.prng, err = utils.NewPRNG(); err != nil {
		return nil, err
	}

	tc.encoder = NewEncoder(tc.params)

	tc.encryptorPk = NewEncryptor(tc.params, tc.pk)
	tc.encryptorSk = NewEncryptor(tc.params, tc.sk)
	tc.decryptor = NewDecryptor(tc.params, tc.sk)

	tc.evaluator = NewEvaluator(tc.params, rlwe.EvaluationKey{Rlk: tc.rlk})

	return tc, nil

}

func newTestVectors(tc *testContext, encryptor Encryptor, a, b complex128, t *testing.T) (values []complex128, plaintext *Plaintext, ciphertext *Ciphertext) {

	logSlots := tc.params.LogSlots()

	values = make([]complex128, 1<<logSlots)

	switch tc.params.RingType() {
	case ring.Standard:
		for i := 0; i < 1<<logSlots; i++ {
			values[i] = complex(utils.RandFloat64(real(a), real(b)), utils.RandFloat64(imag(a), imag(b)))
		}
	case ring.ConjugateInvariant:
		for i := 0; i < 1<<logSlots; i++ {
			values[i] = complex(utils.RandFloat64(real(a), real(b)), 0)
		}
	default:
		panic("invalid ring type")
	}

	values[0] = complex(0.607538, 0)

	plaintext = tc.encoder.EncodeNew(values, tc.params.MaxLevel(), tc.params.DefaultScale(), logSlots)

	if encryptor != nil {
		ciphertext = encryptor.EncryptNew(plaintext)
	}

	return values, plaintext, ciphertext
}

func randomConst(tp ring.Type, a, b complex128) (constant complex128) {
	switch tp {
	case ring.Standard:
		constant = complex(utils.RandFloat64(real(a), real(b)), utils.RandFloat64(imag(a), imag(b)))
	case ring.ConjugateInvariant:
		constant = complex(utils.RandFloat64(real(a), real(b)), 0)
	default:
		panic("invalid ring type")
	}
	return
}

func verifyTestVectors(params Parameters, encoder Encoder, decryptor Decryptor, valuesWant []complex128, element interface{}, logSlots int, bound float64, t *testing.T) {

	precStats := GetPrecisionStats(params, encoder, decryptor, valuesWant, element, logSlots, bound)

	if *printPrecisionStats {
		t.Log(precStats.String())
	}

	require.GreaterOrEqual(t, precStats.MeanPrecision.Real, minPrec)
	require.GreaterOrEqual(t, precStats.MeanPrecision.Imag, minPrec)
}

func testParameters(tc *testContext, t *testing.T) {

	t.Run(GetTestName(tc.params, "Parameters/NewParameters"), func(t *testing.T) {
		params, err := NewParametersFromLiteral(ParametersLiteral{LogN: 4, LogQ: []int{60, 60}, LogP: []int{60}})
		require.NoError(t, err)
		require.Equal(t, ring.Standard, params.RingType())   // Default ring type should be standard
		require.Equal(t, rlwe.DefaultSigma, params.Sigma())  // Default error std should be rlwe.DefaultSigma
		require.Equal(t, params.LogN()-1, params.LogSlots()) // Default number of slots should be N/2
	})

	t.Run(GetTestName(tc.params, "Parameters/CopyNew/"), func(t *testing.T) {
		params1, params2 := tc.params.CopyNew(), tc.params.CopyNew()
		assert.True(t, params1.Equals(tc.params) && params2.Equals(tc.params))
		params1.logSlots = 7
		params1.defaultScale = 3.14
		assert.False(t, params1.Equals(tc.params))
		assert.True(t, params2.Equals(tc.params))
	})

	t.Run(GetTestName(tc.params, "Parameters/StandardRing/"), func(t *testing.T) {
		params, err := tc.params.StandardParameters()
		switch tc.params.RingType() {
		case ring.Standard:
			require.True(t, params.Equals(tc.params))
			require.NoError(t, err)
		case ring.ConjugateInvariant:
			require.Equal(t, params.LogN(), tc.params.LogN()+1)
			require.Equal(t, params.LogSlots(), tc.params.LogSlots())
			require.NoError(t, err)
		default:
			t.Error("invalid RingType")
		}
	})
}

func testEncoder(tc *testContext, t *testing.T) {

	t.Run(GetTestName(tc.params, "Encoder/EncodeSlots/"), func(t *testing.T) {

		values, plaintext, _ := newTestVectors(tc, nil, complex(-1, -1), complex(1, 1), t)

		verifyTestVectors(tc.params, tc.encoder, tc.decryptor, values, plaintext, tc.params.LogSlots(), 0, t)
	})

	t.Run(GetTestName(tc.params, "Encoder/EncodeCoeffs/"), func(t *testing.T) {

		slots := tc.params.N()

		valuesWant := make([]float64, slots)

		for i := 0; i < slots; i++ {
			valuesWant[i] = utils.RandFloat64(-1, 1)
		}

		valuesWant[0] = 0.607538

		plaintext := tc.encoder.EncodeCoeffsNew(valuesWant, tc.params.MaxLevel(), tc.params.DefaultScale())

		valuesTest := tc.encoder.DecodeCoeffs(plaintext)

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

func testEvaluatorAdd(tc *testContext, t *testing.T) {

	t.Run(GetTestName(tc.params, "Evaluator/Add/CtCt/"), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectors(tc, tc.encryptorSk, complex(-1, -1), complex(1, 1), t)
		values2, _, ciphertext2 := newTestVectors(tc, tc.encryptorSk, complex(-1, -1), complex(1, 1), t)

		for i := range values1 {
			values1[i] += values2[i]
		}

		tc.evaluator.Add(ciphertext1, ciphertext2, ciphertext1)

		verifyTestVectors(tc.params, tc.encoder, tc.decryptor, values1, ciphertext1, tc.params.LogSlots(), 0, t)
	})

	t.Run(GetTestName(tc.params, "Evaluator/AddNew/CtCt/"), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectors(tc, tc.encryptorSk, complex(-1, -1), complex(1, 1), t)
		values2, _, ciphertext2 := newTestVectors(tc, tc.encryptorSk, complex(-1, -1), complex(1, 1), t)

		for i := range values1 {
			values1[i] += values2[i]
		}

		ciphertext3 := tc.evaluator.AddNew(ciphertext1, ciphertext2)

		verifyTestVectors(tc.params, tc.encoder, tc.decryptor, values1, ciphertext3, tc.params.LogSlots(), 0, t)
	})

	t.Run(GetTestName(tc.params, "Evaluator/Add/CtPlain/"), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectors(tc, tc.encryptorSk, complex(-1, -1), complex(1, 1), t)
		values2, plaintext2, _ := newTestVectors(tc, tc.encryptorSk, complex(-1, -1), complex(1, 1), t)

		for i := range values1 {
			values1[i] += values2[i]
		}

		tc.evaluator.Add(ciphertext1, plaintext2, ciphertext1)

		verifyTestVectors(tc.params, tc.encoder, tc.decryptor, values1, ciphertext1, tc.params.LogSlots(), 0, t)

		for i := range values1 {
			values1[i] += values2[i]
		}

		tc.evaluator.Add(plaintext2, ciphertext1, ciphertext1)

		verifyTestVectors(tc.params, tc.encoder, tc.decryptor, values1, ciphertext1, tc.params.LogSlots(), 0, t)
	})

	t.Run(GetTestName(tc.params, "Evaluator/AddNew/CtPlain/"), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectors(tc, tc.encryptorSk, complex(-1, -1), complex(1, 1), t)
		values2, plaintext2, _ := newTestVectors(tc, tc.encryptorSk, complex(-1, -1), complex(1, 1), t)

		for i := range values1 {
			values1[i] += values2[i]
		}

		ciphertext3 := tc.evaluator.AddNew(ciphertext1, plaintext2)

		verifyTestVectors(tc.params, tc.encoder, tc.decryptor, values1, ciphertext3, tc.params.LogSlots(), 0, t)

		ciphertext3 = tc.evaluator.AddNew(plaintext2, ciphertext1)

		verifyTestVectors(tc.params, tc.encoder, tc.decryptor, values1, ciphertext3, tc.params.LogSlots(), 0, t)
	})

}

func testEvaluatorSub(tc *testContext, t *testing.T) {

	t.Run(GetTestName(tc.params, "Evaluator/Sub/CtCt/"), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectors(tc, tc.encryptorSk, complex(-1, -1), complex(1, 1), t)
		values2, _, ciphertext2 := newTestVectors(tc, tc.encryptorSk, complex(-1, -1), complex(1, 1), t)

		for i := range values1 {
			values1[i] -= values2[i]
		}

		tc.evaluator.Sub(ciphertext1, ciphertext2, ciphertext1)

		verifyTestVectors(tc.params, tc.encoder, tc.decryptor, values1, ciphertext1, tc.params.LogSlots(), 0, t)
	})

	t.Run(GetTestName(tc.params, "Evaluator/SubNew/CtCt/"), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectors(tc, tc.encryptorSk, complex(-1, -1), complex(1, 1), t)
		values2, _, ciphertext2 := newTestVectors(tc, tc.encryptorSk, complex(-1, -1), complex(1, 1), t)

		for i := range values1 {
			values1[i] -= values2[i]
		}

		ciphertext3 := tc.evaluator.SubNew(ciphertext1, ciphertext2)

		verifyTestVectors(tc.params, tc.encoder, tc.decryptor, values1, ciphertext3, tc.params.LogSlots(), 0, t)
	})

	t.Run(GetTestName(tc.params, "Evaluator/Sub/CtPlain/"), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectors(tc, tc.encryptorSk, complex(-1, -1), complex(1, 1), t)
		values2, plaintext2, ciphertext2 := newTestVectors(tc, tc.encryptorSk, complex(-1, -1), complex(1, 1), t)

		valuesTest := make([]complex128, len(values1))
		for i := range values1 {
			valuesTest[i] = values1[i] - values2[i]
		}

		tc.evaluator.Sub(ciphertext1, plaintext2, ciphertext2)

		verifyTestVectors(tc.params, tc.encoder, tc.decryptor, valuesTest, ciphertext2, tc.params.LogSlots(), 0, t)

		for i := range values1 {
			valuesTest[i] = values2[i] - values1[i]
		}

		tc.evaluator.Sub(plaintext2, ciphertext1, ciphertext2)

		verifyTestVectors(tc.params, tc.encoder, tc.decryptor, valuesTest, ciphertext2, tc.params.LogSlots(), 0, t)
	})

	t.Run(GetTestName(tc.params, "Evaluator/SubNew/CtPlain/"), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectors(tc, tc.encryptorSk, complex(-1, -1), complex(1, 1), t)
		values2, plaintext2, _ := newTestVectors(tc, tc.encryptorSk, complex(-1, -1), complex(1, 1), t)

		valuesTest := make([]complex128, len(values1))
		for i := range values1 {
			valuesTest[i] = values1[i] - values2[i]
		}

		ciphertext3 := tc.evaluator.SubNew(ciphertext1, plaintext2)

		verifyTestVectors(tc.params, tc.encoder, tc.decryptor, valuesTest, ciphertext3, tc.params.LogSlots(), 0, t)

		for i := range values1 {
			valuesTest[i] = values2[i] - values1[i]
		}

		ciphertext3 = tc.evaluator.SubNew(plaintext2, ciphertext1)

		verifyTestVectors(tc.params, tc.encoder, tc.decryptor, valuesTest, ciphertext3, tc.params.LogSlots(), 0, t)
	})

}

func testEvaluatorRescale(tc *testContext, t *testing.T) {

	t.Run(GetTestName(tc.params, "Evaluator/Rescale/Single/"), func(t *testing.T) {

		if tc.params.PCount() == 0 {
			t.Skip("method is unsuported when params.PCount() == 0")
		}

		values, _, ciphertext := newTestVectors(tc, tc.encryptorSk, complex(-1, -1), complex(1, 1), t)

		constant := tc.ringQ.Modulus[ciphertext.Level()]

		tc.evaluator.MultByConst(ciphertext, constant, ciphertext)

		ciphertext.Scale *= float64(constant)

		tc.evaluator.Rescale(ciphertext, tc.params.DefaultScale(), ciphertext)

		verifyTestVectors(tc.params, tc.encoder, tc.decryptor, values, ciphertext, tc.params.LogSlots(), 0, t)
	})

	t.Run(GetTestName(tc.params, "Evaluator/Rescale/Many/"), func(t *testing.T) {

		if tc.params.PCount() == 0 {
			t.Skip("method is unsuported when params.PCount() == 0")
		}

		values, _, ciphertext := newTestVectors(tc, tc.encryptorSk, complex(-1, -1), complex(1, 1), t)

		nbRescales := tc.params.MaxLevel()
		if nbRescales > 5 {
			nbRescales = 5
		}

		for i := 0; i < nbRescales; i++ {
			constant := tc.ringQ.Modulus[ciphertext.Level()-i]
			tc.evaluator.MultByConst(ciphertext, constant, ciphertext)
			ciphertext.Scale *= float64(constant)
		}

		tc.evaluator.Rescale(ciphertext, tc.params.DefaultScale(), ciphertext)

		verifyTestVectors(tc.params, tc.encoder, tc.decryptor, values, ciphertext, tc.params.LogSlots(), 0, t)
	})
}

func testEvaluatorAddConst(tc *testContext, t *testing.T) {

	t.Run(GetTestName(tc.params, "Evaluator/AddConst/"), func(t *testing.T) {

		values, _, ciphertext := newTestVectors(tc, tc.encryptorSk, complex(-1, -1), complex(1, 1), t)

		constant := randomConst(tc.params.RingType(), complex(-1, 1), complex(-1, 1))

		for i := range values {
			values[i] += constant
		}

		tc.evaluator.AddConst(ciphertext, constant, ciphertext)

		verifyTestVectors(tc.params, tc.encoder, tc.decryptor, values, ciphertext, tc.params.LogSlots(), 0, t)
	})
}

func testEvaluatorMultByConst(tc *testContext, t *testing.T) {

	t.Run(GetTestName(tc.params, "Evaluator/MultByConst/"), func(t *testing.T) {

		values, _, ciphertext := newTestVectors(tc, tc.encryptorSk, complex(-1, -1), complex(1, 1), t)

		constant := randomConst(tc.params.RingType(), complex(-1, 1), complex(-1, 1))

		for i := range values {
			values[i] *= constant
		}

		tc.evaluator.MultByConst(ciphertext, constant, ciphertext)

		verifyTestVectors(tc.params, tc.encoder, tc.decryptor, values, ciphertext, tc.params.LogSlots(), 0, t)
	})

}

func testEvaluatorMultByConstAndAdd(tc *testContext, t *testing.T) {

	t.Run(GetTestName(tc.params, "Evaluator/MultByConstAndAdd/"), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectors(tc, tc.encryptorSk, complex(-1, -1), complex(1, 1), t)
		values2, _, ciphertext2 := newTestVectors(tc, tc.encryptorSk, complex(-1, -1), complex(1, 1), t)

		constant := randomConst(tc.params.RingType(), complex(-1, 1), complex(-1, 1))

		for i := range values1 {
			values2[i] += (constant * values1[i])
		}

		tc.evaluator.MultByConstAndAdd(ciphertext1, constant, ciphertext2)

		verifyTestVectors(tc.params, tc.encoder, tc.decryptor, values2, ciphertext2, tc.params.LogSlots(), 0, t)
	})

}

func testEvaluatorMul(tc *testContext, t *testing.T) {

	t.Run(GetTestName(tc.params, "Evaluator/Mul/ct0*pt->ct0/"), func(t *testing.T) {

		values1, plaintext1, ciphertext1 := newTestVectors(tc, tc.encryptorSk, complex(-1, -1), complex(1, 1), t)

		for i := range values1 {
			values1[i] *= values1[i]
		}

		tc.evaluator.MulRelin(ciphertext1, plaintext1, ciphertext1)

		verifyTestVectors(tc.params, tc.encoder, tc.decryptor, values1, ciphertext1, tc.params.LogSlots(), 0, t)
	})

	t.Run(GetTestName(tc.params, "Evaluator/Mul/pt*ct0->ct0/"), func(t *testing.T) {

		values1, plaintext1, ciphertext1 := newTestVectors(tc, tc.encryptorSk, complex(-1, -1), complex(1, 1), t)

		for i := range values1 {
			values1[i] *= values1[i]
		}

		tc.evaluator.MulRelin(ciphertext1, plaintext1, ciphertext1)

		verifyTestVectors(tc.params, tc.encoder, tc.decryptor, values1, ciphertext1, tc.params.LogSlots(), 0, t)
	})

	t.Run(GetTestName(tc.params, "Evaluator/Mul/ct0*pt->ct1/"), func(t *testing.T) {

		values1, plaintext1, ciphertext1 := newTestVectors(tc, tc.encryptorSk, complex(-1, -1), complex(1, 1), t)

		for i := range values1 {
			values1[i] *= values1[i]
		}

		ciphertext2 := tc.evaluator.MulRelinNew(ciphertext1, plaintext1)

		verifyTestVectors(tc.params, tc.encoder, tc.decryptor, values1, ciphertext2, tc.params.LogSlots(), 0, t)
	})

	t.Run(GetTestName(tc.params, "Evaluator/Mul/ct0*ct1->ct0/"), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectors(tc, tc.encryptorSk, complex(-1, -1), complex(1, 1), t)
		values2, _, ciphertext2 := newTestVectors(tc, tc.encryptorSk, complex(-1, -1), complex(1, 1), t)

		for i := range values1 {
			values2[i] *= values1[i]
		}

		tc.evaluator.MulRelin(ciphertext1, ciphertext2, ciphertext1)

		verifyTestVectors(tc.params, tc.encoder, tc.decryptor, values2, ciphertext1, tc.params.LogSlots(), 0, t)
	})

	t.Run(GetTestName(tc.params, "Evaluator/Mul/ct0*ct1->ct1/"), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectors(tc, tc.encryptorSk, complex(-1, -1), complex(1, 1), t)
		values2, _, ciphertext2 := newTestVectors(tc, tc.encryptorSk, complex(-1, -1), complex(1, 1), t)

		for i := range values1 {
			values2[i] *= values1[i]
		}

		tc.evaluator.MulRelin(ciphertext1, ciphertext2, ciphertext2)

		verifyTestVectors(tc.params, tc.encoder, tc.decryptor, values2, ciphertext2, tc.params.LogSlots(), 0, t)
	})

	t.Run(GetTestName(tc.params, "Evaluator/Mul/ct0*ct1->ct2/"), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectors(tc, tc.encryptorSk, complex(-1, -1), complex(1, 1), t)
		values2, _, ciphertext2 := newTestVectors(tc, tc.encryptorSk, complex(-1, -1), complex(1, 1), t)

		for i := range values1 {
			values2[i] *= values1[i]
		}

		ciphertext3 := tc.evaluator.MulRelinNew(ciphertext1, ciphertext2)

		verifyTestVectors(tc.params, tc.encoder, tc.decryptor, values2, ciphertext3, tc.params.LogSlots(), 0, t)
	})

	t.Run(GetTestName(tc.params, "Evaluator/Mul/ct0*ct0->ct0/"), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectors(tc, tc.encryptorSk, complex(-1, -1), complex(1, 1), t)

		for i := range values1 {
			values1[i] *= values1[i]
		}

		tc.evaluator.MulRelin(ciphertext1, ciphertext1, ciphertext1)

		verifyTestVectors(tc.params, tc.encoder, tc.decryptor, values1, ciphertext1, tc.params.LogSlots(), 0, t)
	})

	t.Run(GetTestName(tc.params, "Evaluator/Mul/ct0*ct0->ct1/"), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectors(tc, tc.encryptorSk, complex(-1, -1), complex(1, 1), t)

		for i := range values1 {
			values1[i] *= values1[i]
		}

		ciphertext2 := tc.evaluator.MulRelinNew(ciphertext1, ciphertext1)

		verifyTestVectors(tc.params, tc.encoder, tc.decryptor, values1, ciphertext2, tc.params.LogSlots(), 0, t)
	})

	t.Run(GetTestName(tc.params, "Evaluator/Mul/Relinearize(ct0*ct1->ct0)/"), func(t *testing.T) {

		if tc.params.PCount() == 0 {
			t.Skip("method is unsuported when params.PCount() == 0")
		}

		values1, _, ciphertext1 := newTestVectors(tc, tc.encryptorSk, complex(-1, -1), complex(1, 1), t)
		values2, _, ciphertext2 := newTestVectors(tc, tc.encryptorSk, complex(-1, -1), complex(1, 1), t)

		for i := range values1 {
			values1[i] *= values2[i]
		}

		tc.evaluator.Mul(ciphertext1, ciphertext2, ciphertext1)
		require.Equal(t, ciphertext1.Degree(), 2)
		tc.evaluator.Relinearize(ciphertext1, ciphertext1)
		require.Equal(t, ciphertext1.Degree(), 1)

		verifyTestVectors(tc.params, tc.encoder, tc.decryptor, values1, ciphertext1, tc.params.LogSlots(), 0, t)
	})

	t.Run(GetTestName(tc.params, "Evaluator/Mul/Relinearize(ct0*ct1->ct1)/"), func(t *testing.T) {

		if tc.params.PCount() == 0 {
			t.Skip("method is unsuported when params.PCount() == 0")
		}

		values1, _, ciphertext1 := newTestVectors(tc, tc.encryptorSk, complex(-1, -1), complex(1, 1), t)
		values2, _, ciphertext2 := newTestVectors(tc, tc.encryptorSk, complex(-1, -1), complex(1, 1), t)

		for i := range values1 {
			values2[i] *= values1[i]
		}

		tc.evaluator.Mul(ciphertext1, ciphertext2, ciphertext2)
		require.Equal(t, ciphertext2.Degree(), 2)
		tc.evaluator.Relinearize(ciphertext2, ciphertext2)
		require.Equal(t, ciphertext2.Degree(), 1)

		verifyTestVectors(tc.params, tc.encoder, tc.decryptor, values2, ciphertext2, tc.params.LogSlots(), 0, t)
	})

}

func testFunctions(tc *testContext, t *testing.T) {

	t.Run(GetTestName(tc.params, "Evaluator/PowerOf2/"), func(t *testing.T) {

		if tc.params.PCount() == 0 {
			t.Skip("method is unsuported when params.PCount() == 0")
		}

		if tc.params.MaxLevel() < 3 {
			t.Skip("skipping test for params max level < 3")
		}

		values, _, ciphertext := newTestVectors(tc, tc.encryptorSk, complex(-1, -1), complex(1, 1), t)

		n := 2

		valuesWant := make([]complex128, len(values))
		for i := 0; i < len(valuesWant); i++ {
			valuesWant[i] = values[i]
		}

		for i := 0; i < n; i++ {
			for j := 0; j < len(valuesWant); j++ {
				valuesWant[j] *= valuesWant[j]
			}
		}

		tc.evaluator.PowerOf2(ciphertext, n, ciphertext)

		verifyTestVectors(tc.params, tc.encoder, tc.decryptor, valuesWant, ciphertext, tc.params.LogSlots(), 0, t)
	})

	t.Run(GetTestName(tc.params, "Evaluator/Power/"), func(t *testing.T) {

		if tc.params.PCount() == 0 {
			t.Skip("method is unsuported when params.PCount() == 0")
		}

		if tc.params.MaxLevel() < 4 {
			t.Skip("skipping test for params max level < 4")
		}

		values, _, ciphertext := newTestVectors(tc, tc.encryptorSk, complex(-1, -1), complex(1, 1), t)

		n := 3

		for i := range values {
			values[i] = cmplx.Pow(values[i], complex(float64(n), 0))
		}

		tc.evaluator.Power(ciphertext, n, ciphertext)

		verifyTestVectors(tc.params, tc.encoder, tc.decryptor, values, ciphertext, tc.params.LogSlots(), 0, t)
	})

	t.Run(GetTestName(tc.params, "Evaluator/Inverse/"), func(t *testing.T) {

		if tc.params.PCount() == 0 {
			t.Skip("method is unsuported when params.PCount() == 0")
		}

		if tc.params.MaxLevel() < 7 {
			t.Skip("skipping test for params max level < 7")
		}

		values, _, ciphertext := newTestVectors(tc, tc.encryptorSk, complex(0.1, 0), complex(1, 0), t)

		n := 7

		for i := range values {
			values[i] = 1.0 / values[i]
		}

		ciphertext = tc.evaluator.InverseNew(ciphertext, n)

		verifyTestVectors(tc.params, tc.encoder, tc.decryptor, values, ciphertext, tc.params.LogSlots(), 0, t)
	})
}

func testEvaluatePoly(tc *testContext, t *testing.T) {

	var err error

	t.Run(GetTestName(tc.params, "EvaluatePoly/Exp/"), func(t *testing.T) {

		if tc.params.PCount() == 0 {
			t.Skip("method is unsuported when params.PCount() == 0")
		}

		if tc.params.MaxLevel() < 3 {
			t.Skip("skipping test for params max level < 3")
		}

		values, _, ciphertext := newTestVectors(tc, tc.encryptorSk, complex(-1, 0), complex(1, 0), t)

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

		if ciphertext, err = tc.evaluator.EvaluatePoly(ciphertext, poly, ciphertext.Scale); err != nil {
			t.Error(err)
		}

		verifyTestVectors(tc.params, tc.encoder, tc.decryptor, values, ciphertext, tc.params.LogSlots(), 0, t)
	})
}

func testChebyshevInterpolator(tc *testContext, t *testing.T) {

	var err error

	t.Run(GetTestName(tc.params, "ChebyshevInterpolator/Sin/"), func(t *testing.T) {

		if tc.params.PCount() == 0 {
			t.Skip("method is unsuported when params.PCount() == 0")
		}

		if tc.params.MaxLevel() < 5 {
			t.Skip("skipping test for params max level < 5")
		}

		eval := tc.evaluator

		values, _, ciphertext := newTestVectors(tc, tc.encryptorSk, complex(-1, 0), complex(1, 0), t)

		poly := Approximate(cmplx.Sin, complex(-1.5, 0), complex(1.5, 0), 15)

		for i := range values {
			values[i] = cmplx.Sin(values[i])
		}

		eval.MultByConst(ciphertext, 2/(poly.B-poly.A), ciphertext)
		eval.AddConst(ciphertext, (-poly.A-poly.B)/(poly.B-poly.A), ciphertext)
		eval.Rescale(ciphertext, tc.params.DefaultScale(), ciphertext)

		if ciphertext, err = eval.EvaluatePoly(ciphertext, poly, ciphertext.Scale); err != nil {
			t.Error(err)
		}

		verifyTestVectors(tc.params, tc.encoder, tc.decryptor, values, ciphertext, tc.params.LogSlots(), 0, t)
	})
}

func testDecryptPublic(tc *testContext, t *testing.T) {

	var err error

	t.Run(GetTestName(tc.params, "DecryptPublic/Sin/"), func(t *testing.T) {

		if tc.params.PCount() == 0 {
			t.Skip("method is unsuported when params.PCount() == 0")
		}

		if tc.params.MaxLevel() < 5 {
			t.Skip("skipping test for params max level < 5")
		}

		eval := tc.evaluator

		values, _, ciphertext := newTestVectors(tc, tc.encryptorSk, complex(-1, 0), complex(1, 0), t)

		poly := Approximate(cmplx.Sin, complex(-1.5, 0), complex(1.5, 0), 15)

		for i := range values {
			values[i] = cmplx.Sin(values[i])
		}

		eval.MultByConst(ciphertext, 2/(poly.B-poly.A), ciphertext)
		eval.AddConst(ciphertext, (-poly.A-poly.B)/(poly.B-poly.A), ciphertext)
		eval.Rescale(ciphertext, tc.params.DefaultScale(), ciphertext)

		if ciphertext, err = eval.EvaluatePoly(ciphertext, poly, ciphertext.Scale); err != nil {
			t.Error(err)
		}

		plaintext := tc.decryptor.DecryptNew(ciphertext)

		valuesHave := tc.encoder.Decode(plaintext, tc.params.LogSlots())

		verifyTestVectors(tc.params, tc.encoder, nil, values, valuesHave, tc.params.LogSlots(), 0, t)

		sigma := tc.encoder.GetErrSTDCoeffDomain(values, valuesHave, plaintext.Scale)

		valuesHave = tc.encoder.DecodePublic(plaintext, tc.params.LogSlots(), sigma)

		verifyTestVectors(tc.params, tc.encoder, nil, values, valuesHave, tc.params.LogSlots(), 0, t)
	})
}

func testSwitchKeys(tc *testContext, t *testing.T) {

	var sk2 *rlwe.SecretKey
	var decryptorSk2 Decryptor
	var switchingKey *rlwe.SwitchingKey

	if tc.params.PCount() != 0 {
		sk2 = tc.kgen.GenSecretKey()
		decryptorSk2 = NewDecryptor(tc.params, sk2)
		switchingKey = tc.kgen.GenSwitchingKey(tc.sk, sk2)
	}

	t.Run(GetTestName(tc.params, "SwitchKeys/"), func(t *testing.T) {

		if tc.params.PCount() == 0 {
			t.Skip("method is unsuported when params.PCount() == 0")
		}

		values, _, ciphertext := newTestVectors(tc, tc.encryptorSk, complex(-1, -1), complex(1, 1), t)

		tc.evaluator.SwitchKeys(ciphertext, switchingKey, ciphertext)

		verifyTestVectors(tc.params, tc.encoder, decryptorSk2, values, ciphertext, tc.params.LogSlots(), 0, t)
	})

	t.Run(GetTestName(tc.params, "SwitchKeysNew/"), func(t *testing.T) {

		if tc.params.PCount() == 0 {
			t.Skip("method is unsuported when params.PCount() == 0")
		}

		values, _, ciphertext := newTestVectors(tc, tc.encryptorSk, complex(-1, -1), complex(1, 1), t)

		ciphertext = tc.evaluator.SwitchKeysNew(ciphertext, switchingKey)

		verifyTestVectors(tc.params, tc.encoder, decryptorSk2, values, ciphertext, tc.params.LogSlots(), 0, t)
	})
}

func testBridge(tc *testContext, t *testing.T) {

	t.Run(GetTestName(tc.params, "Bridge/"), func(t *testing.T) {

		if tc.params.RingType() != ring.ConjugateInvariant {
			t.Skip("only tested for params.RingType() == ring.ConjugateInvariant")
		}

		ciParams := tc.params
		var stdParams Parameters
		var err error
		if stdParams, err = ciParams.StandardParameters(); err != nil {
			t.Errorf("all Conjugate Invariant parameters should have a standard counterpart but got: %f", err)
		}

		stdKeyGen := NewKeyGenerator(stdParams)
		stdSK := stdKeyGen.GenSecretKey()
		stdDecryptor := NewDecryptor(stdParams, stdSK)
		stdEncoder := NewEncoder(stdParams)
		stdEvaluator := NewEvaluator(stdParams, rlwe.EvaluationKey{Rlk: nil, Rtks: nil})

		swkCtR, swkRtC := stdKeyGen.GenSwitchingKeysForBridge(stdSK, tc.sk)

		switcher, err := NewDomainSwitcher(stdParams, swkCtR, swkRtC)
		if err != nil {
			t.Error(err)
		}

		values, _, ctCI := newTestVectors(tc, tc.encryptorSk, complex(-1, -1), complex(1, 1), t)

		stdCTHave := NewCiphertext(stdParams, ctCI.Degree(), ctCI.Level(), ctCI.Scale)

		switcher.RealToComplex(ctCI, stdCTHave)

		verifyTestVectors(stdParams, stdEncoder, stdDecryptor, values, stdCTHave, stdParams.LogSlots(), 0, t)

		stdCTImag := stdEvaluator.MultByiNew(stdCTHave)
		stdEvaluator.Add(stdCTHave, stdCTImag, stdCTHave)

		ciCTHave := NewCiphertext(ciParams, 1, stdCTHave.Level(), stdCTHave.Scale)
		switcher.ComplexToReal(stdCTHave, ciCTHave)

		verifyTestVectors(tc.params, tc.encoder, tc.decryptor, values, ciCTHave, ciParams.LogSlots(), 0, t)
	})
}

func testAutomorphisms(tc *testContext, t *testing.T) {

	params := tc.params

	rots := []int{0, 1, -1, 4, -4, 63, -63}
	var rotKey *rlwe.RotationKeySet

	if tc.params.PCount() != 0 {
		rotKey = tc.kgen.GenRotationKeysForRotations(rots, params.RingType() == ring.Standard, tc.sk)
	}

	evaluator := tc.evaluator.WithKey(rlwe.EvaluationKey{Rlk: tc.rlk, Rtks: rotKey})

	t.Run(GetTestName(params, "Conjugate/"), func(t *testing.T) {

		if params.PCount() == 0 {
			t.Skip("method is unsuported when params.PCount() == 0")
		}

		if params.RingType() != ring.Standard {
			t.Skip("Conjugate not defined in real-CKKS")
		}

		values, _, ciphertext := newTestVectors(tc, tc.encryptorSk, complex(-1, -1), complex(1, 1), t)

		for i := range values {
			values[i] = complex(real(values[i]), -imag(values[i]))
		}

		evaluator.Conjugate(ciphertext, ciphertext)

		verifyTestVectors(params, tc.encoder, tc.decryptor, values, ciphertext, params.LogSlots(), 0, t)
	})

	t.Run(GetTestName(params, "ConjugateNew/"), func(t *testing.T) {

		if params.PCount() == 0 {
			t.Skip("method is unsuported when params.PCount() == 0")
		}

		if params.RingType() != ring.Standard {
			t.Skip("Conjugate not defined in real-CKKS")
		}

		values, _, ciphertext := newTestVectors(tc, tc.encryptorSk, complex(-1, -1), complex(1, 1), t)

		for i := range values {
			values[i] = complex(real(values[i]), -imag(values[i]))
		}

		ciphertext = evaluator.ConjugateNew(ciphertext)

		verifyTestVectors(params, tc.encoder, tc.decryptor, values, ciphertext, params.LogSlots(), 0, t)
	})

	t.Run(GetTestName(tc.params, "Rotate/"), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectors(tc, tc.encryptorSk, complex(-1, -1), complex(1, 1), t)

		ciphertext2 := NewCiphertext(tc.params, ciphertext1.Degree(), ciphertext1.Level(), ciphertext1.Scale)

		for _, n := range rots {
			evaluator.Rotate(ciphertext1, n, ciphertext2)
			verifyTestVectors(tc.params, tc.encoder, tc.decryptor, utils.RotateComplex128Slice(values1, n), ciphertext2, tc.params.LogSlots(), 0, t)
		}
	})

	t.Run(GetTestName(tc.params, "RotateNew/"), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectors(tc, tc.encryptorSk, complex(-1, -1), complex(1, 1), t)

		for _, n := range rots {
			verifyTestVectors(tc.params, tc.encoder, tc.decryptor, utils.RotateComplex128Slice(values1, n), evaluator.RotateNew(ciphertext1, n), tc.params.LogSlots(), 0, t)
		}

	})

	t.Run(GetTestName(tc.params, "RotateHoisted/"), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectors(tc, tc.encryptorSk, complex(-1, -1), complex(1, 1), t)

		ciphertexts := evaluator.RotateHoistedNew(ciphertext1, rots)

		for _, n := range rots {
			verifyTestVectors(tc.params, tc.encoder, tc.decryptor, utils.RotateComplex128Slice(values1, n), ciphertexts[n], tc.params.LogSlots(), 0, t)
		}
	})
}

func testInnerSum(tc *testContext, t *testing.T) {

	t.Run(GetTestName(tc.params, "InnerSum/"), func(t *testing.T) {

		if tc.params.PCount() == 0 {
			t.Skip("method is unsuported when params.PCount() == 0")
		}

		batch := 7
		n := 35

		rotKey := tc.kgen.GenRotationKeysForRotations(tc.params.RotationsForInnerSum(batch, n), false, tc.sk)
		eval := tc.evaluator.WithKey(rlwe.EvaluationKey{Rlk: tc.rlk, Rtks: rotKey})

		values1, _, ciphertext1 := newTestVectors(tc, tc.encryptorSk, complex(-1, -1), complex(1, 1), t)

		eval.InnerSum(ciphertext1, batch, n, ciphertext1)

		tmp0 := make([]complex128, len(values1))
		copy(tmp0, values1)

		for i := 1; i < n; i++ {

			tmp1 := utils.RotateComplex128Slice(tmp0, i*batch)

			for j := range values1 {
				values1[j] += tmp1[j]
			}
		}

		verifyTestVectors(tc.params, tc.encoder, tc.decryptor, values1, ciphertext1, tc.params.LogSlots(), 0, t)
	})

	t.Run(GetTestName(tc.params, "InnerSumLog/"), func(t *testing.T) {

		if tc.params.PCount() == 0 {
			t.Skip("method is unsuported when params.PCount() == 0")
		}

		batch := 512
		n := tc.params.Slots() / batch

		rotKey := tc.kgen.GenRotationKeysForRotations(tc.params.RotationsForInnerSumLog(batch, n), false, tc.sk)
		eval := tc.evaluator.WithKey(rlwe.EvaluationKey{Rlk: tc.rlk, Rtks: rotKey})

		values1, _, ciphertext1 := newTestVectors(tc, tc.encryptorSk, complex(-1, -1), complex(1, 1), t)

		eval.InnerSumLog(ciphertext1, batch, n, ciphertext1)

		tmp0 := make([]complex128, len(values1))
		copy(tmp0, values1)

		for i := 1; i < n; i++ {

			tmp1 := utils.RotateComplex128Slice(tmp0, i*batch)

			for j := range values1 {
				values1[j] += tmp1[j]
			}
		}

		verifyTestVectors(tc.params, tc.encoder, tc.decryptor, values1, ciphertext1, tc.params.LogSlots(), 0, t)

	})

	t.Run(GetTestName(tc.params, "Average/"), func(t *testing.T) {

		if tc.params.PCount() == 0 {
			t.Skip("method is unsuported when params.PCount() == 0")
		}

		logBatch := 9
		batch := 1 << logBatch
		n := tc.params.Slots() / batch

		rotKey := tc.kgen.GenRotationKeysForRotations(tc.params.RotationsForInnerSumLog(batch, n), false, tc.sk)
		eval := tc.evaluator.WithKey(rlwe.EvaluationKey{Rlk: tc.rlk, Rtks: rotKey})

		values1, _, ciphertext1 := newTestVectors(tc, tc.encryptorSk, complex(-1, -1), complex(1, 1), t)

		eval.Average(ciphertext1, logBatch, ciphertext1)

		tmp0 := make([]complex128, len(values1))
		copy(tmp0, values1)

		for i := 1; i < n; i++ {

			tmp1 := utils.RotateComplex128Slice(tmp0, i*batch)

			for j := range values1 {
				values1[j] += tmp1[j]
			}
		}

		for i := range values1 {
			values1[i] /= complex(float64(n), 0)
		}

		verifyTestVectors(tc.params, tc.encoder, tc.decryptor, values1, ciphertext1, tc.params.LogSlots(), 0, t)

	})
}

func testReplicate(tc *testContext, t *testing.T) {

	if tc.params.PCount() == 0 {
		t.Skip("method is unsuported when params.PCount() == 0")
	}

	t.Run(GetTestName(tc.params, "Replicate/"), func(t *testing.T) {
		batch := 2
		n := 35

		rotKey := tc.kgen.GenRotationKeysForRotations(tc.params.RotationsForReplicate(batch, n), false, tc.sk)
		eval := tc.evaluator.WithKey(rlwe.EvaluationKey{Rlk: tc.rlk, Rtks: rotKey})

		values1, _, ciphertext1 := newTestVectors(tc, tc.encryptorSk, complex(-1, -1), complex(1, 1), t)

		eval.Replicate(ciphertext1, batch, n, ciphertext1)

		tmp0 := make([]complex128, len(values1))
		copy(tmp0, values1)

		for i := 1; i < n; i++ {

			tmp1 := utils.RotateComplex128Slice(tmp0, i*-batch)

			for j := range values1 {
				values1[j] += tmp1[j]
			}
		}

		verifyTestVectors(tc.params, tc.encoder, tc.decryptor, values1, ciphertext1, tc.params.LogSlots(), 0, t)
	})

	t.Run(GetTestName(tc.params, "ReplicateLog/"), func(t *testing.T) {

		batch := 3
		n := 15

		rotKey := tc.kgen.GenRotationKeysForRotations(tc.params.RotationsForReplicateLog(batch, n), false, tc.sk)
		eval := tc.evaluator.WithKey(rlwe.EvaluationKey{Rlk: tc.rlk, Rtks: rotKey})

		values1, _, ciphertext1 := newTestVectors(tc, tc.encryptorSk, complex(-1, -1), complex(1, 1), t)

		eval.ReplicateLog(ciphertext1, batch, n, ciphertext1)

		tmp0 := make([]complex128, len(values1))
		copy(tmp0, values1)

		for i := 1; i < n; i++ {

			tmp1 := utils.RotateComplex128Slice(tmp0, i*-batch)

			for j := range values1 {
				values1[j] += tmp1[j]
			}
		}

		verifyTestVectors(tc.params, tc.encoder, tc.decryptor, values1, ciphertext1, tc.params.LogSlots(), 0, t)

	})
}

func testLinearTransform(tc *testContext, t *testing.T) {

	if tc.params.PCount() == 0 {
		t.Skip("method is unsuported when params.PCount() == 0")
	}

	t.Run(GetTestName(tc.params, "LinearTransform/BSGS/"), func(t *testing.T) {

		params := tc.params

		values1, _, ciphertext1 := newTestVectors(tc, tc.encryptorSk, complex(-1, -1), complex(1, 1), t)

		diagMatrix := make(map[int][]complex128)

		diagMatrix[-15] = make([]complex128, params.Slots())
		diagMatrix[-4] = make([]complex128, params.Slots())
		diagMatrix[-1] = make([]complex128, params.Slots())
		diagMatrix[0] = make([]complex128, params.Slots())
		diagMatrix[1] = make([]complex128, params.Slots())
		diagMatrix[2] = make([]complex128, params.Slots())
		diagMatrix[3] = make([]complex128, params.Slots())
		diagMatrix[4] = make([]complex128, params.Slots())
		diagMatrix[15] = make([]complex128, params.Slots())

		for i := 0; i < params.Slots(); i++ {
			diagMatrix[-15][i] = complex(1, 0)
			diagMatrix[-4][i] = complex(1, 0)
			diagMatrix[-1][i] = complex(1, 0)
			diagMatrix[0][i] = complex(1, 0)
			diagMatrix[1][i] = complex(1, 0)
			diagMatrix[2][i] = complex(1, 0)
			diagMatrix[3][i] = complex(1, 0)
			diagMatrix[4][i] = complex(1, 0)
			diagMatrix[15][i] = complex(1, 0)
		}

		linTransf := NewLinearTransformBSGS(tc.encoder, diagMatrix, params.MaxLevel(), params.DefaultScale(), 1.0, params.logSlots)

		rots := tc.params.RotationsForDiagMatrixMult(linTransf)

		rotKey := tc.kgen.GenRotationKeysForRotations(rots, false, tc.sk)

		eval := tc.evaluator.WithKey(rlwe.EvaluationKey{Rlk: tc.rlk, Rtks: rotKey})

		eval.LinearTransform(ciphertext1, linTransf, []*Ciphertext{ciphertext1})

		tmp := make([]complex128, params.Slots())
		copy(tmp, values1)

		for i := 0; i < params.Slots(); i++ {
			values1[i] += tmp[(i-15+params.Slots())%params.Slots()]
			values1[i] += tmp[(i-4+params.Slots())%params.Slots()]
			values1[i] += tmp[(i-1+params.Slots())%params.Slots()]
			values1[i] += tmp[(i+1)%params.Slots()]
			values1[i] += tmp[(i+2)%params.Slots()]
			values1[i] += tmp[(i+3)%params.Slots()]
			values1[i] += tmp[(i+4)%params.Slots()]
			values1[i] += tmp[(i+15)%params.Slots()]
		}

		verifyTestVectors(tc.params, tc.encoder, tc.decryptor, values1, ciphertext1, tc.params.LogSlots(), 0, t)
	})

	t.Run(GetTestName(tc.params, "LinearTransform/Naive/"), func(t *testing.T) {

		params := tc.params

		values1, _, ciphertext1 := newTestVectors(tc, tc.encryptorSk, complex(-1, -1), complex(1, 1), t)

		diagMatrix := make(map[int][]complex128)

		diagMatrix[-1] = make([]complex128, params.Slots())
		diagMatrix[0] = make([]complex128, params.Slots())

		for i := 0; i < params.Slots(); i++ {
			diagMatrix[-1][i] = complex(1, 0)
			diagMatrix[0][i] = complex(1, 0)
		}

		linTransf := NewLinearTransform(tc.encoder, diagMatrix, params.MaxLevel(), params.DefaultScale(), params.LogSlots())

		rots := tc.params.RotationsForDiagMatrixMult(linTransf)

		rotKey := tc.kgen.GenRotationKeysForRotations(rots, false, tc.sk)

		eval := tc.evaluator.WithKey(rlwe.EvaluationKey{Rlk: tc.rlk, Rtks: rotKey})

		eval.LinearTransform(ciphertext1, linTransf, []*Ciphertext{ciphertext1})

		tmp := make([]complex128, params.Slots())
		copy(tmp, values1)

		for i := 0; i < params.Slots(); i++ {
			values1[i] += tmp[(i-1+params.Slots())%params.Slots()]
		}

		verifyTestVectors(tc.params, tc.encoder, tc.decryptor, values1, ciphertext1, tc.params.LogSlots(), 0, t)
	})
}

func testMarshaller(testctx *testContext, t *testing.T) {

	t.Run(GetTestName(testctx.params, "Marshaller/Parameters/Binary"), func(t *testing.T) {
		bytes, err := testctx.params.MarshalBinary()
		assert.Nil(t, err)
		var p Parameters
		err = p.UnmarshalBinary(bytes)
		assert.Nil(t, err)
		assert.Equal(t, testctx.params, p)
		assert.Equal(t, testctx.params.RingQ(), p.RingQ())
	})

	t.Run(GetTestName(testctx.params, "Marshaller/Parameters/JSON"), func(t *testing.T) {
		// checks that parameters can be marshalled without error
		data, err := json.Marshal(testctx.params)
		assert.Nil(t, err)
		assert.NotNil(t, data)

		// checks that ckks.Parameters can be unmarshalled without error
		var paramsRec Parameters
		err = json.Unmarshal(data, &paramsRec)
		assert.Nil(t, err)
		assert.True(t, testctx.params.Equals(paramsRec))

		// checks that ckks.Paramters can be unmarshalled with log-moduli definition without error
		dataWithLogModuli := []byte(fmt.Sprintf(`{"LogN":%d,"LogQ":[50,50],"LogP":[60],"Sigma":3.2,"T":65537}`, testctx.params.LogN()))
		var paramsWithLogModuli Parameters
		err = json.Unmarshal(dataWithLogModuli, &paramsWithLogModuli)
		assert.Nil(t, err)
		assert.Equal(t, 2, paramsWithLogModuli.QCount())
		assert.Equal(t, 1, paramsWithLogModuli.PCount())
		assert.Equal(t, ring.Standard, paramsWithLogModuli.RingType()) // Omitting the RingType field should result in a standard instance

		// checks that ckks.Paramters can be unmarshalled with log-moduli definition with empty P without error
		dataWithLogModuliNoP := []byte(fmt.Sprintf(`{"LogN":%d,"LogQ":[50,50],"LogP":[],"Sigma":3.2,"T":65537, "RingType": "ConjugateInvariant"}`, testctx.params.LogN()))
		var paramsWithLogModuliNoP Parameters
		err = json.Unmarshal(dataWithLogModuliNoP, &paramsWithLogModuliNoP)
		assert.Nil(t, err)
		assert.Equal(t, 2, paramsWithLogModuliNoP.QCount())
		assert.Equal(t, 0, paramsWithLogModuliNoP.PCount())
		assert.Equal(t, ring.ConjugateInvariant, paramsWithLogModuliNoP.RingType())
	})

	t.Run("Marshaller/Ciphertext/", func(t *testing.T) {
		t.Run(GetTestName(testctx.params, "EndToEnd/"), func(t *testing.T) {

			ciphertextWant := NewCiphertextRandom(testctx.prng, testctx.params, 2, testctx.params.MaxLevel(), testctx.params.DefaultScale())

			marshalledCiphertext, err := ciphertextWant.MarshalBinary()
			require.NoError(t, err)

			ciphertextTest := new(Ciphertext)
			require.NoError(t, ciphertextTest.UnmarshalBinary(marshalledCiphertext))

			require.Equal(t, ciphertextWant.Degree(), ciphertextTest.Degree())
			require.Equal(t, ciphertextWant.Level(), ciphertextTest.Level())
			require.Equal(t, ciphertextWant.Scale, ciphertextTest.Scale)

			for i := range ciphertextWant.Value {
				require.True(t, testctx.ringQ.EqualLvl(ciphertextWant.Level(), ciphertextWant.Value[i], ciphertextTest.Value[i]))
			}
		})

		t.Run(GetTestName(testctx.params, "Minimal/"), func(t *testing.T) {

			ciphertext := NewCiphertextRandom(testctx.prng, testctx.params, 0, testctx.params.MaxLevel(), testctx.params.DefaultScale())

			marshalledCiphertext, err := ciphertext.MarshalBinary()
			require.NoError(t, err)

			ciphertextTest := new(Ciphertext)
			require.Error(t, ciphertextTest.UnmarshalBinary(nil))
			require.NoError(t, ciphertextTest.UnmarshalBinary(marshalledCiphertext))

			require.Equal(t, ciphertext.Degree(), 0)
			require.Equal(t, ciphertext.Level(), testctx.params.MaxLevel())
			require.Equal(t, ciphertext.Scale, testctx.params.DefaultScale())
			require.Equal(t, len(ciphertext.Value), 1)
		})
	})
}
