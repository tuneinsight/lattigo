package ckks

import (
	"encoding/json"
	"flag"
	"fmt"
	"math"
	"math/big"
	"math/bits"
	"runtime"
	"testing"

	"github.com/stretchr/testify/require"
	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/ring/distribution"
	"github.com/tuneinsight/lattigo/v4/rlwe"
	"github.com/tuneinsight/lattigo/v4/utils"
	"github.com/tuneinsight/lattigo/v4/utils/bignum"
	"github.com/tuneinsight/lattigo/v4/utils/bignum/approximation"
	"github.com/tuneinsight/lattigo/v4/utils/bignum/polynomial"
	"github.com/tuneinsight/lattigo/v4/utils/sampling"
)

var flagParamString = flag.String("params", "", "specify the test cryptographic parameters as a JSON string. Overrides -short and -long.")
var printPrecisionStats = flag.Bool("print-precision", false, "print precision stats")

func GetTestName(params Parameters, opname string) string {
	return fmt.Sprintf("%s/RingType=%s/logN=%d/logQP=%d/Qi=%d/Pi=%d/LogPlaintextScale=%d",
		opname,
		params.RingType(),
		params.LogN(),
		int(math.Round(params.LogQP())),
		params.QCount(),
		params.PCount(),
		int(math.Log2(params.PlaintextScale().Float64())))
}

type testContext struct {
	params      Parameters
	ringQ       *ring.Ring
	ringP       *ring.Ring
	prng        sampling.PRNG
	encoder     *Encoder
	kgen        *rlwe.KeyGenerator
	sk          *rlwe.SecretKey
	pk          *rlwe.PublicKey
	encryptorPk rlwe.EncryptorInterface
	encryptorSk rlwe.EncryptorInterface
	decryptor   *rlwe.Decryptor
	evaluator   *Evaluator
}

func TestCKKS(t *testing.T) {

	var err error

	var testParams []ParametersLiteral
	switch {
	case *flagParamString != "": // the custom test suite reads the parameters from the -params flag
		testParams = append(testParams, ParametersLiteral{})
		if err = json.Unmarshal([]byte(*flagParamString), &testParams[0]); err != nil {
			t.Fatal(err)
		}
	default:
		testParams = TestParamsLiteral
	}

	for _, ringType := range []ring.Type{ring.Standard, ring.ConjugateInvariant} {

		for _, paramsLiteral := range testParams {

			paramsLiteral.RingType = ringType

			var params Parameters
			if params, err = NewParametersFromLiteral(paramsLiteral); err != nil {
				t.Fatal(err)
			}

			var tc *testContext
			if tc, err = genTestParams(params); err != nil {
				t.Fatal(err)
			}

			for _, testSet := range []func(tc *testContext, t *testing.T){
				testParameters,
				testEncoder,
				testEvaluatorAdd,
				testEvaluatorSub,
				testEvaluatorRescale,
				testEvaluatorMul,
				testEvaluatorMulThenAdd,
				testFunctions,
				testDecryptPublic,
				testEvaluatePoly,
				testChebyshevInterpolator,
				testBridge,
				testLinearTransform,
				testMarshaller,
			} {
				testSet(tc, t)
				runtime.GC()
			}
		}
	}

}

func genTestParams(defaultParam Parameters) (tc *testContext, err error) {

	tc = new(testContext)

	tc.params = defaultParam

	tc.kgen = NewKeyGenerator(tc.params)

	tc.sk, tc.pk = tc.kgen.GenKeyPairNew()

	tc.ringQ = defaultParam.RingQ()
	if tc.params.PCount() != 0 {
		tc.ringP = defaultParam.RingP()
	}

	if tc.prng, err = sampling.NewPRNG(); err != nil {
		return nil, err
	}

	tc.encoder = NewEncoder(tc.params)

	tc.encryptorPk = NewEncryptor(tc.params, tc.pk)
	tc.encryptorSk = NewEncryptor(tc.params, tc.sk)
	tc.decryptor = NewDecryptor(tc.params, tc.sk)

	tc.evaluator = NewEvaluator(tc.params, &rlwe.EvaluationKeySet{RelinearizationKey: tc.kgen.GenRelinearizationKeyNew(tc.sk)})

	return tc, nil

}

func newTestVectors(tc *testContext, encryptor rlwe.EncryptorInterface, a, b complex128, t *testing.T) (values []*bignum.Complex, pt *rlwe.Plaintext, ct *rlwe.Ciphertext) {

	prec := tc.encoder.Prec()

	pt = NewPlaintext(tc.params, tc.params.MaxLevel())

	values = make([]*bignum.Complex, pt.PlaintextSlots())

	switch tc.params.RingType() {
	case ring.Standard:
		for i := range values {
			values[i] = &bignum.Complex{
				bignum.NewFloat(sampling.RandFloat64(real(a), real(b)), prec),
				bignum.NewFloat(sampling.RandFloat64(imag(a), imag(b)), prec),
			}
		}
	case ring.ConjugateInvariant:
		for i := range values {
			values[i] = &bignum.Complex{
				bignum.NewFloat(sampling.RandFloat64(real(a), real(b)), prec),
				new(big.Float),
			}
		}
	default:
		panic("invalid ring type")
	}

	tc.encoder.Encode(values, pt)

	if encryptor != nil {
		ct = encryptor.EncryptNew(pt)
	}

	return values, pt, ct
}

func randomConst(tp ring.Type, prec uint, a, b complex128) (constant *bignum.Complex) {
	switch tp {
	case ring.Standard:
		constant = &bignum.Complex{
			bignum.NewFloat(sampling.RandFloat64(real(a), real(b)), prec),
			bignum.NewFloat(sampling.RandFloat64(imag(a), imag(b)), prec),
		}
	case ring.ConjugateInvariant:
		constant = &bignum.Complex{
			bignum.NewFloat(sampling.RandFloat64(real(a), real(b)), prec),
			new(big.Float),
		}
	default:
		panic("invalid ring type")
	}
	return
}

func verifyTestVectors(params Parameters, encoder *Encoder, decryptor *rlwe.Decryptor, valuesWant, valuesHave interface{}, noise distribution.Distribution, t *testing.T) {

	precStats := GetPrecisionStats(params, encoder, decryptor, valuesWant, valuesHave, noise, false)

	if *printPrecisionStats {
		t.Log(precStats.String())
	}

	rf64, _ := precStats.MeanPrecision.Real.Float64()
	if64, _ := precStats.MeanPrecision.Imag.Float64()

	minPrec := math.Log2(params.PlaintextScale().Float64()) - float64(params.LogN()+2)
	if minPrec < 0 {
		minPrec = 0
	}

	require.GreaterOrEqual(t, rf64, minPrec)
	require.GreaterOrEqual(t, if64, minPrec)
}

func testParameters(tc *testContext, t *testing.T) {

	t.Run(GetTestName(tc.params, "Parameters/NewParameters"), func(t *testing.T) {
		params, err := NewParametersFromLiteral(ParametersLiteral{
			LogN:              4,
			LogQ:              []int{60, 60},
			LogP:              []int{60},
			LogPlaintextScale: 0,
		})
		require.NoError(t, err)
		require.Equal(t, ring.Standard, params.RingType()) // Default ring type should be standard
		require.Equal(t, &rlwe.DefaultXe, params.Xe())     // Default error std should be rlwe.DefaultSigma
	})

	t.Run(GetTestName(tc.params, "Parameters/StandardRing"), func(t *testing.T) {
		params, err := tc.params.StandardParameters()
		switch tc.params.RingType() {
		case ring.Standard:
			require.True(t, params.Equal(tc.params))
			require.NoError(t, err)
		case ring.ConjugateInvariant:
			require.Equal(t, params.LogN(), tc.params.LogN()+1)
			require.NoError(t, err)
		default:
			t.Fatal("invalid RingType")
		}
	})
}

func testEncoder(tc *testContext, t *testing.T) {

	t.Run(GetTestName(tc.params, "Encoder/FrequencyDomain"), func(t *testing.T) {

		values, plaintext, _ := newTestVectors(tc, nil, -1-1i, 1+1i, t)

		verifyTestVectors(tc.params, tc.encoder, nil, values, plaintext, nil, t)
	})

	t.Run(GetTestName(tc.params, "Encoder/CoefficientDomain"), func(t *testing.T) {

		slots := tc.params.N()

		valuesWant := make([]float64, slots)

		for i := 0; i < slots; i++ {
			valuesWant[i] = sampling.RandFloat64(-1, 1)
		}

		valuesWant[0] = 0.607538

		pt := NewPlaintext(tc.params, tc.params.MaxLevel())
		pt.EncodingDomain = rlwe.TimeDomain

		tc.encoder.Encode(valuesWant, pt)

		valuesTest := make([]float64, len(valuesWant))

		tc.encoder.Decode(pt, valuesTest)

		var meanprec float64

		for i := range valuesWant {
			meanprec += math.Abs(valuesTest[i] - valuesWant[i])
		}

		meanprec /= float64(slots)

		if *printPrecisionStats {
			t.Logf("\nMean    precision : %.2f \n", math.Log2(1/meanprec))
		}

		minPrec := math.Log2(tc.params.PlaintextScale().Float64()) - float64(tc.params.LogN()+2)
		if minPrec < 0 {
			minPrec = 0
		}

		require.GreaterOrEqual(t, math.Log2(1/meanprec), minPrec)
	})

}

func testEvaluatorAdd(tc *testContext, t *testing.T) {

	t.Run(GetTestName(tc.params, "Evaluator/AddNew/Ct"), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectors(tc, tc.encryptorSk, -1-1i, 1+1i, t)
		values2, _, ciphertext2 := newTestVectors(tc, tc.encryptorSk, -1-1i, 1+1i, t)

		for i := range values1 {
			values1[i].Add(values1[i], values2[i])
		}

		ciphertext3 := tc.evaluator.AddNew(ciphertext1, ciphertext2)

		verifyTestVectors(tc.params, tc.encoder, tc.decryptor, values1, ciphertext3, nil, t)
	})

	t.Run(GetTestName(tc.params, "Evaluator/Add/Ct"), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectors(tc, tc.encryptorSk, -1-1i, 1+1i, t)
		values2, _, ciphertext2 := newTestVectors(tc, tc.encryptorSk, -1-1i, 1+1i, t)

		for i := range values1 {
			values1[i].Add(values1[i], values2[i])
		}

		tc.evaluator.Add(ciphertext1, ciphertext2, ciphertext1)

		verifyTestVectors(tc.params, tc.encoder, tc.decryptor, values1, ciphertext1, nil, t)
	})

	t.Run(GetTestName(tc.params, "Evaluator/Add/Pt"), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectors(tc, tc.encryptorSk, -1-1i, 1+1i, t)
		values2, plaintext2, _ := newTestVectors(tc, tc.encryptorSk, -1-1i, 1+1i, t)

		for i := range values1 {
			values1[i].Add(values1[i], values2[i])
		}

		tc.evaluator.Add(ciphertext1, plaintext2, ciphertext1)

		verifyTestVectors(tc.params, tc.encoder, tc.decryptor, values1, ciphertext1, nil, t)
	})

	t.Run(GetTestName(tc.params, "Evaluator/Add/Scalar"), func(t *testing.T) {

		values, _, ciphertext := newTestVectors(tc, tc.encryptorSk, -1-1i, 1+1i, t)

		constant := randomConst(tc.params.RingType(), tc.encoder.Prec(), -1+1i, -1+1i)

		for i := range values {
			values[i].Add(values[i], constant)
		}

		tc.evaluator.Add(ciphertext, constant, ciphertext)

		verifyTestVectors(tc.params, tc.encoder, tc.decryptor, values, ciphertext, nil, t)
	})

	t.Run(GetTestName(tc.params, "Evaluator/Add/Vector"), func(t *testing.T) {

		values1, _, ciphertext := newTestVectors(tc, tc.encryptorSk, -1-1i, 1+1i, t)
		values2, _, _ := newTestVectors(tc, tc.encryptorSk, -1-1i, 1+1i, t)

		for i := range values1 {
			values1[i].Add(values1[i], values2[i])
		}

		tc.evaluator.Add(ciphertext, values2, ciphertext)

		verifyTestVectors(tc.params, tc.encoder, tc.decryptor, values1, ciphertext, nil, t)
	})
}

func testEvaluatorSub(tc *testContext, t *testing.T) {

	t.Run(GetTestName(tc.params, "Evaluator/SubNew/Ct"), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectors(tc, tc.encryptorSk, -1-1i, 1+1i, t)
		values2, _, ciphertext2 := newTestVectors(tc, tc.encryptorSk, -1-1i, 1+1i, t)

		for i := range values1 {
			values1[i].Sub(values1[i], values2[i])
		}

		ciphertext3 := tc.evaluator.SubNew(ciphertext1, ciphertext2)

		verifyTestVectors(tc.params, tc.encoder, tc.decryptor, values1, ciphertext3, nil, t)
	})

	t.Run(GetTestName(tc.params, "Evaluator/Sub/Ct"), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectors(tc, tc.encryptorSk, -1-1i, 1+1i, t)
		values2, _, ciphertext2 := newTestVectors(tc, tc.encryptorSk, -1-1i, 1+1i, t)

		for i := range values1 {
			values1[i].Sub(values1[i], values2[i])
		}

		tc.evaluator.Sub(ciphertext1, ciphertext2, ciphertext1)

		verifyTestVectors(tc.params, tc.encoder, tc.decryptor, values1, ciphertext1, nil, t)
	})

	t.Run(GetTestName(tc.params, "Evaluator/Sub/Pt"), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectors(tc, tc.encryptorSk, -1-1i, 1+1i, t)
		values2, plaintext2, ciphertext2 := newTestVectors(tc, tc.encryptorSk, -1-1i, 1+1i, t)

		valuesTest := make([]*bignum.Complex, len(values1))
		for i := range values1 {
			valuesTest[i] = bignum.NewComplex()
			valuesTest[i].Sub(values1[i], values2[i])
		}

		tc.evaluator.Sub(ciphertext1, plaintext2, ciphertext2)

		verifyTestVectors(tc.params, tc.encoder, tc.decryptor, valuesTest, ciphertext2, nil, t)
	})

	t.Run(GetTestName(tc.params, "Evaluator/Sub/Scalar"), func(t *testing.T) {

		values, _, ciphertext := newTestVectors(tc, tc.encryptorSk, -1-1i, 1+1i, t)

		constant := randomConst(tc.params.RingType(), tc.encoder.Prec(), -1+1i, -1+1i)

		for i := range values {
			values[i].Sub(values[i], constant)
		}

		ciphertext = tc.evaluator.SubNew(ciphertext, constant)

		verifyTestVectors(tc.params, tc.encoder, tc.decryptor, values, ciphertext, nil, t)
	})

	t.Run(GetTestName(tc.params, "Evaluator/Sub/Vector"), func(t *testing.T) {

		values1, _, ciphertext := newTestVectors(tc, tc.encryptorSk, -1-1i, 1+1i, t)
		values2, _, _ := newTestVectors(tc, tc.encryptorSk, -1-1i, 1+1i, t)

		for i := range values1 {
			values1[i].Sub(values1[i], values2[i])
		}

		tc.evaluator.Sub(ciphertext, values2, ciphertext)

		verifyTestVectors(tc.params, tc.encoder, tc.decryptor, values1, ciphertext, nil, t)
	})
}

func testEvaluatorRescale(tc *testContext, t *testing.T) {

	t.Run(GetTestName(tc.params, "Evaluator/Rescale/Single"), func(t *testing.T) {

		if tc.params.MaxLevel() < 2 {
			t.Skip("skipping test for params max level < 2")
		}

		values, _, ciphertext := newTestVectors(tc, tc.encryptorSk, -1-1i, 1+1i, t)

		constant := tc.ringQ.SubRings[ciphertext.Level()].Modulus

		tc.evaluator.Mul(ciphertext, constant, ciphertext)

		ciphertext.PlaintextScale = ciphertext.PlaintextScale.Mul(rlwe.NewScale(constant))

		if err := tc.evaluator.Rescale(ciphertext, tc.params.PlaintextScale(), ciphertext); err != nil {
			t.Fatal(err)
		}

		verifyTestVectors(tc.params, tc.encoder, tc.decryptor, values, ciphertext, nil, t)
	})

	t.Run(GetTestName(tc.params, "Evaluator/Rescale/Many"), func(t *testing.T) {

		if tc.params.MaxLevel() < 2 {
			t.Skip("skipping test for params max level < 2")
		}

		values, _, ciphertext := newTestVectors(tc, tc.encryptorSk, -1-1i, 1+1i, t)

		nbRescales := tc.params.MaxLevel()
		if nbRescales > 5 {
			nbRescales = 5
		}

		for i := 0; i < nbRescales; i++ {
			constant := tc.ringQ.SubRings[ciphertext.Level()-i].Modulus
			tc.evaluator.Mul(ciphertext, constant, ciphertext)
			ciphertext.PlaintextScale = ciphertext.PlaintextScale.Mul(rlwe.NewScale(constant))
		}

		if err := tc.evaluator.Rescale(ciphertext, tc.params.PlaintextScale(), ciphertext); err != nil {
			t.Fatal(err)
		}

		verifyTestVectors(tc.params, tc.encoder, tc.decryptor, values, ciphertext, nil, t)
	})
}

func testEvaluatorMul(tc *testContext, t *testing.T) {

	t.Run(GetTestName(tc.params, "Evaluator/MulNew/Ct/Pt"), func(t *testing.T) {

		values1, plaintext1, ciphertext1 := newTestVectors(tc, tc.encryptorSk, -1-1i, 1+1i, t)

		mul := bignum.NewComplexMultiplier()

		for i := range values1 {
			mul.Mul(values1[i], values1[i], values1[i])
		}

		ciphertext2 := tc.evaluator.MulNew(ciphertext1, plaintext1)

		verifyTestVectors(tc.params, tc.encoder, tc.decryptor, values1, ciphertext2, nil, t)
	})

	t.Run(GetTestName(tc.params, "Evaluator/Mul/Ct/Scalar"), func(t *testing.T) {

		values, _, ciphertext := newTestVectors(tc, tc.encryptorSk, -1-1i, 1+1i, t)

		constant := randomConst(tc.params.RingType(), tc.encoder.Prec(), -1+1i, -1+1i)

		mul := bignum.NewComplexMultiplier()

		for i := range values {
			mul.Mul(values[i], constant, values[i])
		}

		tc.evaluator.Mul(ciphertext, constant, ciphertext)

		verifyTestVectors(tc.params, tc.encoder, tc.decryptor, values, ciphertext, nil, t)
	})

	t.Run(GetTestName(tc.params, "Evaluator/Mul/Ct/Vector"), func(t *testing.T) {

		values1, _, ciphertext := newTestVectors(tc, tc.encryptorSk, -1-1i, 1+1i, t)
		values2, _, _ := newTestVectors(tc, tc.encryptorSk, -1-1i, 1+1i, t)

		mul := bignum.NewComplexMultiplier()

		for i := range values1 {
			mul.Mul(values1[i], values2[i], values1[i])
		}

		tc.evaluator.Mul(ciphertext, values2, ciphertext)

		verifyTestVectors(tc.params, tc.encoder, tc.decryptor, values1, ciphertext, nil, t)
	})

	t.Run(GetTestName(tc.params, "Evaluator/Mul/Ct/Pt"), func(t *testing.T) {

		values1, plaintext1, ciphertext1 := newTestVectors(tc, tc.encryptorSk, -1-1i, 1+1i, t)

		mul := bignum.NewComplexMultiplier()

		for i := range values1 {
			mul.Mul(values1[i], values1[i], values1[i])
		}

		tc.evaluator.MulRelin(ciphertext1, plaintext1, ciphertext1)

		verifyTestVectors(tc.params, tc.encoder, tc.decryptor, values1, ciphertext1, nil, t)
	})

	t.Run(GetTestName(tc.params, "Evaluator/Mul/Ct/Ct/Degree0"), func(t *testing.T) {

		values1, plaintext1, _ := newTestVectors(tc, tc.encryptorSk, -1-1i, 1+1i, t)
		values2, _, ciphertext2 := newTestVectors(tc, tc.encryptorSk, -1-1i, 1+1i, t)

		mul := bignum.NewComplexMultiplier()

		for i := range values1 {
			mul.Mul(values2[i], values1[i], values2[i])
		}

		ciphertext1 := &rlwe.Ciphertext{}
		ciphertext1.Value = []*ring.Poly{plaintext1.Value}
		ciphertext1.MetaData = plaintext1.MetaData

		tc.evaluator.MulRelin(ciphertext1, ciphertext2, ciphertext1)

		verifyTestVectors(tc.params, tc.encoder, tc.decryptor, values2, ciphertext1, nil, t)
	})

	t.Run(GetTestName(tc.params, "Evaluator/MulRelin/Ct/Ct"), func(t *testing.T) {

		// op0 <- op0 * op1
		values1, _, ciphertext1 := newTestVectors(tc, tc.encryptorSk, -1-1i, 1+1i, t)
		values2, _, ciphertext2 := newTestVectors(tc, tc.encryptorSk, -1-1i, 1+1i, t)

		mul := bignum.NewComplexMultiplier()

		for i := range values1 {
			mul.Mul(values1[i], values2[i], values1[i])
		}

		tc.evaluator.MulRelin(ciphertext1, ciphertext2, ciphertext1)
		require.Equal(t, ciphertext1.Degree(), 1)

		verifyTestVectors(tc.params, tc.encoder, tc.decryptor, values1, ciphertext1, nil, t)

		// op1 <- op0 * op1
		values1, _, ciphertext1 = newTestVectors(tc, tc.encryptorSk, -1-1i, 1+1i, t)
		values2, _, ciphertext2 = newTestVectors(tc, tc.encryptorSk, -1-1i, 1+1i, t)

		for i := range values1 {
			mul.Mul(values2[i], values1[i], values2[i])
		}

		tc.evaluator.MulRelin(ciphertext1, ciphertext2, ciphertext2)
		require.Equal(t, ciphertext2.Degree(), 1)

		verifyTestVectors(tc.params, tc.encoder, tc.decryptor, values2, ciphertext2, nil, t)

		// op0 <- op0 * op0
		for i := range values1 {
			mul.Mul(values1[i], values1[i], values1[i])
		}

		tc.evaluator.MulRelin(ciphertext1, ciphertext1, ciphertext1)
		require.Equal(t, ciphertext1.Degree(), 1)

		verifyTestVectors(tc.params, tc.encoder, tc.decryptor, values1, ciphertext1, nil, t)
	})
}

func testEvaluatorMulThenAdd(tc *testContext, t *testing.T) {

	t.Run(GetTestName(tc.params, "Evaluator/MulThenAdd/Scalar"), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectors(tc, tc.encryptorSk, -1-1i, 1+1i, t)
		values2, _, ciphertext2 := newTestVectors(tc, tc.encryptorSk, -1-1i, 1+1i, t)

		constant := randomConst(tc.params.RingType(), tc.encoder.Prec(), -1+1i, -1+1i)

		mul := bignum.NewComplexMultiplier()

		tmp := new(bignum.Complex)
		tmp[0] = new(big.Float)
		tmp[1] = new(big.Float)

		for i := range values1 {
			mul.Mul(values1[i], constant, tmp)
			values2[i].Add(values2[i], tmp)
		}

		tc.evaluator.MulThenAdd(ciphertext1, constant, ciphertext2)

		verifyTestVectors(tc.params, tc.encoder, tc.decryptor, values2, ciphertext2, nil, t)
	})

	t.Run(GetTestName(tc.params, "Evaluator/MulThenAdd/Vector"), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectors(tc, tc.encryptorSk, -1, 1, t)
		values2, _, ciphertext2 := newTestVectors(tc, tc.encryptorSk, -1, 1, t)

		tc.evaluator.MulThenAdd(ciphertext2, values1, ciphertext1)

		mul := bignum.NewComplexMultiplier()

		tmp := new(bignum.Complex)
		tmp[0] = new(big.Float)
		tmp[1] = new(big.Float)

		for i := range values1 {
			mul.Mul(values2[i], values1[i], tmp)
			values1[i].Add(values1[i], tmp)
		}

		require.Equal(t, ciphertext1.Degree(), 1)

		verifyTestVectors(tc.params, tc.encoder, tc.decryptor, values1, ciphertext1, nil, t)
	})

	t.Run(GetTestName(tc.params, "Evaluator/MulThenAdd/Pt"), func(t *testing.T) {

		values1, plaintext1, ciphertext1 := newTestVectors(tc, tc.encryptorSk, -1, 1, t)
		values2, _, ciphertext2 := newTestVectors(tc, tc.encryptorSk, -1, 1, t)

		mul := bignum.NewComplexMultiplier()

		tmp := new(bignum.Complex)
		tmp[0] = new(big.Float)
		tmp[1] = new(big.Float)

		for i := range values1 {
			mul.Mul(values2[i], values1[i], tmp)
			values1[i].Add(values1[i], tmp)
		}

		tc.evaluator.MulThenAdd(ciphertext2, plaintext1, ciphertext1)

		require.Equal(t, ciphertext1.Degree(), 1)

		verifyTestVectors(tc.params, tc.encoder, tc.decryptor, values1, ciphertext1, nil, t)
	})

	t.Run(GetTestName(tc.params, "Evaluator/MulRelinThenAdd/Ct"), func(t *testing.T) {

		// op2 = op2 + op1 * op0
		values1, _, ciphertext1 := newTestVectors(tc, tc.encryptorSk, -1-1i, 1+1i, t)
		values2, _, ciphertext2 := newTestVectors(tc, tc.encryptorSk, -1-1i, 1+1i, t)

		mul := bignum.NewComplexMultiplier()

		for i := range values1 {
			mul.Mul(values1[i], values2[i], values2[i])
		}

		ciphertext3 := NewCiphertext(tc.params, 2, ciphertext1.Level())

		ciphertext3.PlaintextScale = ciphertext1.PlaintextScale.Mul(ciphertext2.PlaintextScale)

		tc.evaluator.MulThenAdd(ciphertext1, ciphertext2, ciphertext3)

		require.Equal(t, ciphertext3.Degree(), 2)

		tc.evaluator.Relinearize(ciphertext3, ciphertext3)

		require.Equal(t, ciphertext3.Degree(), 1)

		verifyTestVectors(tc.params, tc.encoder, tc.decryptor, values2, ciphertext3, nil, t)

		// op1 = op1 + op0*op0
		values1, _, ciphertext1 = newTestVectors(tc, tc.encryptorSk, -1-1i, 1+1i, t)
		values2, _, ciphertext2 = newTestVectors(tc, tc.encryptorSk, -1-1i, 1+1i, t)

		tmp := bignum.NewComplex()
		for i := range values1 {
			mul.Mul(values2[i], values2[i], tmp)
			values1[i].Add(values1[i], tmp)
		}

		tc.evaluator.MulRelinThenAdd(ciphertext2, ciphertext2, ciphertext1)

		require.Equal(t, ciphertext1.Degree(), 1)

		verifyTestVectors(tc.params, tc.encoder, tc.decryptor, values1, ciphertext1, nil, t)
	})
}

func testFunctions(tc *testContext, t *testing.T) {

	t.Run(GetTestName(tc.params, "Evaluator/GoldschmidtDivisionNew"), func(t *testing.T) {

		min := 0.1

		values, _, ciphertext := newTestVectors(tc, tc.encryptorSk, complex(min, 0), complex(2-min, 0), t)

		one := new(big.Float).SetInt64(1)
		for i := range values {
			values[i][0].Quo(one, values[i][0])
		}

		logPrec := math.Log2(tc.params.PlaintextScale().Float64()) - float64(tc.params.LogN()-1)

		var err error
		if ciphertext, err = tc.evaluator.GoldschmidtDivisionNew(ciphertext, min, logPrec, NewSecretKeyBootstrapper(tc.params, tc.sk)); err != nil {
			t.Fatal(err)
		}

		verifyTestVectors(tc.params, tc.encoder, tc.decryptor, values, ciphertext, nil, t)
	})
}

func testEvaluatePoly(tc *testContext, t *testing.T) {

	var err error

	t.Run(GetTestName(tc.params, "EvaluatePoly/PolySingle/Exp"), func(t *testing.T) {

		if tc.params.MaxLevel() < 3 {
			t.Skip("skipping test for params max level < 3")
		}

		values, _, ciphertext := newTestVectors(tc, tc.encryptorSk, -1, 1, t)

		prec := tc.encoder.Prec()

		coeffs := []*big.Float{
			bignum.NewFloat(1, prec),
			bignum.NewFloat(1, prec),
			new(big.Float).Quo(bignum.NewFloat(1, prec), bignum.NewFloat(2, prec)),
			new(big.Float).Quo(bignum.NewFloat(1, prec), bignum.NewFloat(6, prec)),
			new(big.Float).Quo(bignum.NewFloat(1, prec), bignum.NewFloat(24, prec)),
			new(big.Float).Quo(bignum.NewFloat(1, prec), bignum.NewFloat(120, prec)),
			new(big.Float).Quo(bignum.NewFloat(1, prec), bignum.NewFloat(720, prec)),
			new(big.Float).Quo(bignum.NewFloat(1, prec), bignum.NewFloat(5040, prec)),
		}

		poly := polynomial.NewPolynomial(polynomial.Monomial, coeffs, nil)

		for i := range values {
			values[i] = poly.Evaluate(values[i])
		}

		if ciphertext, err = tc.evaluator.Polynomial(ciphertext, poly, ciphertext.PlaintextScale); err != nil {
			t.Fatal(err)
		}

		verifyTestVectors(tc.params, tc.encoder, tc.decryptor, values, ciphertext, nil, t)
	})

	t.Run(GetTestName(tc.params, "Polynomial/PolyVector/Exp"), func(t *testing.T) {

		if tc.params.MaxLevel() < 3 {
			t.Skip("skipping test for params max level < 3")
		}

		values, _, ciphertext := newTestVectors(tc, tc.encryptorSk, -1, 1, t)

		prec := tc.encoder.Prec()

		coeffs := []*big.Float{
			bignum.NewFloat(1, prec),
			bignum.NewFloat(1, prec),
			new(big.Float).Quo(bignum.NewFloat(1, prec), bignum.NewFloat(2, prec)),
			new(big.Float).Quo(bignum.NewFloat(1, prec), bignum.NewFloat(6, prec)),
			new(big.Float).Quo(bignum.NewFloat(1, prec), bignum.NewFloat(24, prec)),
			new(big.Float).Quo(bignum.NewFloat(1, prec), bignum.NewFloat(120, prec)),
			new(big.Float).Quo(bignum.NewFloat(1, prec), bignum.NewFloat(720, prec)),
			new(big.Float).Quo(bignum.NewFloat(1, prec), bignum.NewFloat(5040, prec)),
		}

		poly := polynomial.NewPolynomial(polynomial.Monomial, coeffs, nil)

		slots := ciphertext.PlaintextSlots()

		slotIndex := make(map[int][]int)
		idx := make([]int, slots>>1)
		for i := 0; i < slots>>1; i++ {
			idx[i] = 2 * i
		}

		slotIndex[0] = idx

		valuesWant := make([]*bignum.Complex, slots)
		for _, j := range idx {
			valuesWant[j] = poly.Evaluate(values[j])
		}

		polyVector := rlwe.NewPolynomialVector([]*rlwe.Polynomial{rlwe.NewPolynomial(poly)}, slotIndex)

		if ciphertext, err = tc.evaluator.Polynomial(ciphertext, polyVector, ciphertext.PlaintextScale); err != nil {
			t.Fatal(err)
		}

		verifyTestVectors(tc.params, tc.encoder, tc.decryptor, valuesWant, ciphertext, nil, t)
	})
}

func testChebyshevInterpolator(tc *testContext, t *testing.T) {

	var err error

	t.Run(GetTestName(tc.params, "ChebyshevInterpolator/Sin"), func(t *testing.T) {

		degree := 13

		if tc.params.MaxDepth() < bits.Len64(uint64(degree)) {
			t.Skip("skipping test: not enough levels")
		}

		eval := tc.evaluator

		values, _, ciphertext := newTestVectors(tc, tc.encryptorSk, -1, 1, t)

		prec := tc.params.PlaintextPrecision()

		sin := func(x *bignum.Complex) (y *bignum.Complex) {
			xf64, _ := x[0].Float64()
			y = bignum.NewComplex()
			y.SetPrec(prec)
			y[0].SetFloat64(math.Sin(xf64))
			return
		}

		interval := polynomial.Interval{
			A: *new(big.Float).SetPrec(prec).SetFloat64(-8),
			B: *new(big.Float).SetPrec(prec).SetFloat64(8),
		}

		poly := rlwe.NewPolynomial(approximation.Chebyshev(sin, interval, degree))

		scalar, constant := poly.ChangeOfBasis()
		eval.Mul(ciphertext, scalar, ciphertext)
		eval.Add(ciphertext, constant, ciphertext)
		if err = eval.Rescale(ciphertext, tc.params.PlaintextScale(), ciphertext); err != nil {
			t.Fatal(err)
		}

		if ciphertext, err = eval.Polynomial(ciphertext, poly, ciphertext.PlaintextScale); err != nil {
			t.Fatal(err)
		}

		for i := range values {
			values[i] = poly.Evaluate(values[i])
		}

		verifyTestVectors(tc.params, tc.encoder, tc.decryptor, values, ciphertext, nil, t)
	})
}

func testDecryptPublic(tc *testContext, t *testing.T) {

	var err error

	t.Run(GetTestName(tc.params, "DecryptPublic/Sin"), func(t *testing.T) {

		degree := 7
		a, b := -1.5, 1.5

		if tc.params.MaxDepth() < bits.Len64(uint64(degree)) {
			t.Skip("skipping test: not enough levels")
		}

		eval := tc.evaluator

		values, _, ciphertext := newTestVectors(tc, tc.encryptorSk, complex(a, 0), complex(b, 0), t)

		prec := tc.params.PlaintextPrecision()

		sin := func(x *bignum.Complex) (y *bignum.Complex) {
			xf64, _ := x[0].Float64()
			y = bignum.NewComplex()
			y.SetPrec(prec)
			y[0].SetFloat64(math.Sin(xf64))
			return
		}

		interval := polynomial.Interval{
			A: *new(big.Float).SetPrec(prec).SetFloat64(a),
			B: *new(big.Float).SetPrec(prec).SetFloat64(b),
		}

		poly := approximation.Chebyshev(sin, interval, degree)

		for i := range values {
			values[i] = poly.Evaluate(values[i])
		}

		scalar, constant := poly.ChangeOfBasis()

		eval.Mul(ciphertext, scalar, ciphertext)
		eval.Add(ciphertext, constant, ciphertext)
		if err := eval.Rescale(ciphertext, tc.params.PlaintextScale(), ciphertext); err != nil {
			t.Fatal(err)
		}

		if ciphertext, err = eval.Polynomial(ciphertext, poly, ciphertext.PlaintextScale); err != nil {
			t.Fatal(err)
		}

		plaintext := tc.decryptor.DecryptNew(ciphertext)

		valuesHave := make([]*big.Float, plaintext.PlaintextSlots())

		tc.encoder.Decode(plaintext, valuesHave)

		verifyTestVectors(tc.params, tc.encoder, nil, values, valuesHave, nil, t)

		for i := range valuesHave {
			valuesHave[i].Sub(valuesHave[i], values[i][0])
		}

		// This should make it lose at most ~0.5 bit or precision.
		sigma := StandardDeviation(valuesHave, rlwe.NewScale(plaintext.PlaintextScale.Float64()/math.Sqrt(float64(len(values)))))

		tc.encoder.DecodePublic(plaintext, valuesHave, &distribution.DiscreteGaussian{Sigma: sigma, Bound: 2.5066282746310002 * sigma})

		verifyTestVectors(tc.params, tc.encoder, nil, values, valuesHave, nil, t)
	})
}

func testBridge(tc *testContext, t *testing.T) {

	t.Run(GetTestName(tc.params, "Bridge"), func(t *testing.T) {

		if tc.params.RingType() != ring.ConjugateInvariant {
			t.Skip("only tested for params.RingType() == ring.ConjugateInvariant")
		}

		ciParams := tc.params
		var err error
		if _, err = ciParams.StandardParameters(); err != nil {
			t.Fatalf("all Conjugate Invariant parameters should have a standard counterpart but got: %f", err)
		}

		// Create equivalent parameters with RingStandard ring type and different auxiliary modulus P
		stdParamsLit := ciParams.ParametersLiteral()
		stdParamsLit.LogN = ciParams.LogN() + 1
		stdParamsLit.P = []uint64{0x1ffffffff6c80001, 0x1ffffffff6140001} // Assigns new P to ensure that independence from auxiliary P is tested
		stdParamsLit.RingType = ring.Standard
		stdParams, err := NewParametersFromLiteral(stdParamsLit)
		require.Nil(t, err)

		stdKeyGen := NewKeyGenerator(stdParams)
		stdSK := stdKeyGen.GenSecretKeyNew()
		stdDecryptor := NewDecryptor(stdParams, stdSK)
		stdEncoder := NewEncoder(stdParams)
		stdEvaluator := NewEvaluator(stdParams, nil)

		evkCtR, evkRtC := stdKeyGen.GenEvaluationKeysForRingSwapNew(stdSK, tc.sk)

		switcher, err := NewDomainSwitcher(stdParams, evkCtR, evkRtC)
		if err != nil {
			t.Fatal(err)
		}

		evalStandar := NewEvaluator(stdParams, nil)

		values, _, ctCI := newTestVectors(tc, tc.encryptorSk, -1-1i, 1+1i, t)

		stdCTHave := NewCiphertext(stdParams, ctCI.Degree(), ctCI.Level())

		switcher.RealToComplex(evalStandar, ctCI, stdCTHave)

		verifyTestVectors(stdParams, stdEncoder, stdDecryptor, values, stdCTHave, nil, t)

		stdCTImag := stdEvaluator.MulNew(stdCTHave, 1i)
		stdEvaluator.Add(stdCTHave, stdCTImag, stdCTHave)

		ciCTHave := NewCiphertext(ciParams, 1, stdCTHave.Level())
		switcher.ComplexToReal(evalStandar, stdCTHave, ciCTHave)

		verifyTestVectors(tc.params, tc.encoder, tc.decryptor, values, ciCTHave, nil, t)
	})
}

func testLinearTransform(tc *testContext, t *testing.T) {

	t.Run(GetTestName(tc.params, "Average"), func(t *testing.T) {

		values, _, ciphertext := newTestVectors(tc, tc.encryptorSk, -1-1i, 1+1i, t)

		slots := ciphertext.PlaintextSlots()

		logBatch := 9
		batch := 1 << logBatch
		n := slots / batch

		evk := rlwe.NewEvaluationKeySet()
		for _, galEl := range tc.params.GaloisElementsForInnerSum(batch, n) {
			evk.GaloisKeys[galEl] = tc.kgen.GenGaloisKeyNew(galEl, tc.sk)
		}

		eval := tc.evaluator.WithKey(evk)

		eval.Average(ciphertext, logBatch, ciphertext)

		tmp0 := make([]*bignum.Complex, len(values))
		for i := range tmp0 {
			tmp0[i] = values[i].Clone()
		}

		for i := 1; i < n; i++ {

			tmp1 := utils.RotateSlice(tmp0, i*batch)

			for j := range values {
				values[j].Add(values[j], tmp1[j])
			}
		}

		nB := new(big.Float).SetFloat64(float64(n))

		for i := range values {
			values[i][0].Quo(values[i][0], nB)
			values[i][1].Quo(values[i][1], nB)
		}

		verifyTestVectors(tc.params, tc.encoder, tc.decryptor, values, ciphertext, nil, t)
	})

	t.Run(GetTestName(tc.params, "LinearTransform/BSGS=True"), func(t *testing.T) {

		params := tc.params

		values, _, ciphertext := newTestVectors(tc, tc.encryptorSk, -1-1i, 1+1i, t)

		slots := ciphertext.PlaintextSlots()

		nonZeroDiags := []int{-15, -4, -1, 0, 1, 2, 3, 4, 15}

		one := new(big.Float).SetInt64(1)
		zero := new(big.Float)

		diagMatrix := make(map[int][]*bignum.Complex)
		for _, i := range nonZeroDiags {
			diagMatrix[i] = make([]*bignum.Complex, slots)

			for j := 0; j < slots; j++ {
				diagMatrix[i][j] = &bignum.Complex{one, zero}
			}
		}

		LogBSGSRatio := 1

		linTransf, err := GenLinearTransform(diagMatrix, tc.encoder, params.MaxLevel(), rlwe.NewScale(params.Q()[params.MaxLevel()]), ciphertext.PlaintextLogDimensions[1], LogBSGSRatio)
		require.NoError(t, err)

		galEls := params.GaloisElementsForLinearTransform(nonZeroDiags, ciphertext.PlaintextLogDimensions[1], LogBSGSRatio)

		evk := rlwe.NewEvaluationKeySet()
		for _, galEl := range galEls {
			evk.GaloisKeys[galEl] = tc.kgen.GenGaloisKeyNew(galEl, tc.sk)
		}

		eval := tc.evaluator.WithKey(evk)

		eval.LinearTransform(ciphertext, linTransf, []*rlwe.Ciphertext{ciphertext})

		tmp := make([]*bignum.Complex, len(values))
		for i := range tmp {
			tmp[i] = values[i].Clone()
		}

		for i := 0; i < slots; i++ {
			values[i].Add(values[i], tmp[(i-15+slots)%slots])
			values[i].Add(values[i], tmp[(i-4+slots)%slots])
			values[i].Add(values[i], tmp[(i-1+slots)%slots])
			values[i].Add(values[i], tmp[(i+1)%slots])
			values[i].Add(values[i], tmp[(i+2)%slots])
			values[i].Add(values[i], tmp[(i+3)%slots])
			values[i].Add(values[i], tmp[(i+4)%slots])
			values[i].Add(values[i], tmp[(i+15)%slots])
		}

		verifyTestVectors(tc.params, tc.encoder, tc.decryptor, values, ciphertext, nil, t)
	})

	t.Run(GetTestName(tc.params, "LinearTransform/BSGS=False"), func(t *testing.T) {

		params := tc.params

		values, _, ciphertext := newTestVectors(tc, tc.encryptorSk, -1-1i, 1+1i, t)

		slots := ciphertext.PlaintextSlots()

		diagMatrix := make(map[int][]*bignum.Complex)

		diagMatrix[-1] = make([]*bignum.Complex, slots)
		diagMatrix[0] = make([]*bignum.Complex, slots)

		one := new(big.Float).SetInt64(1)
		zero := new(big.Float)

		for i := 0; i < slots; i++ {
			diagMatrix[-1][i] = &bignum.Complex{one, zero}
			diagMatrix[0][i] = &bignum.Complex{one, zero}
		}

		linTransf, err := GenLinearTransform(diagMatrix, tc.encoder, params.MaxLevel(), rlwe.NewScale(params.Q()[params.MaxLevel()]), ciphertext.PlaintextLogDimensions[1], -1)
		require.NoError(t, err)

		galEls := params.GaloisElementsForLinearTransform([]int{-1, 0}, ciphertext.PlaintextLogDimensions[1], -1)

		evk := rlwe.NewEvaluationKeySet()
		for _, galEl := range galEls {
			evk.GaloisKeys[galEl] = tc.kgen.GenGaloisKeyNew(galEl, tc.sk)
		}

		eval := tc.evaluator.WithKey(evk)

		eval.LinearTransform(ciphertext, linTransf, []*rlwe.Ciphertext{ciphertext})

		tmp := make([]*bignum.Complex, slots)
		for i := range tmp {
			tmp[i] = values[i].Clone()
		}

		for i := 0; i < slots; i++ {
			values[i].Add(values[i], tmp[(i-1+slots)%slots])
		}

		verifyTestVectors(tc.params, tc.encoder, tc.decryptor, values, ciphertext, nil, t)
	})
}

func testMarshaller(tc *testContext, t *testing.T) {

	/*
		t.Run(GetTestName(tc.params, "Marshaller/Parameters/JSON"), func(t *testing.T) {
			// checks that parameters can be marshalled without error
			data, err := json.Marshal(tc.params)
			require.Nil(t, err)
			require.NotNil(t, data)

				// checks that ckks.Parameters can be unmarshalled without error
				var paramsRec Parameters
				err = json.Unmarshal(data, &paramsRec)
				require.Nil(t, err)
				require.True(t, tc.params.Equals(paramsRec))

				// checks that ckks.Parameters can be unmarshalled with log-moduli definition without error
				dataWithLogModuli := []byte(fmt.Sprintf(`{"LogN":%d,"LogQ":[50,50],"LogP":[60], "DefaultScale":1.0}`, tc.params.LogN()))
				var paramsWithLogModuli Parameters
				err = json.Unmarshal(dataWithLogModuli, &paramsWithLogModuli)
				require.Nil(t, err)
				require.Equal(t, 2, paramsWithLogModuli.QCount())
				require.Equal(t, 1, paramsWithLogModuli.PCount())
				require.Equal(t, ring.Standard, paramsWithLogModuli.RingType())  // Omitting the RingType field should result in a standard instance
				require.Equal(t, rlwe.DefaultSigma, paramsWithLogModuli.Sigma()) // Omitting sigma should result in Default being used

				// checks that ckks.Parameters can be unmarshalled with log-moduli definition with empty P without error
				dataWithLogModuliNoP := []byte(fmt.Sprintf(`{"LogN":%d,"LogQ":[50,50],"LogP":[],"DefaultScale":1.0,"RingType": "ConjugateInvariant"}`, tc.params.LogN()))
				var paramsWithLogModuliNoP Parameters
				err = json.Unmarshal(dataWithLogModuliNoP, &paramsWithLogModuliNoP)
				require.Nil(t, err)
				require.Equal(t, 2, paramsWithLogModuliNoP.QCount())
				require.Equal(t, 0, paramsWithLogModuliNoP.PCount())
				require.Equal(t, ring.ConjugateInvariant, paramsWithLogModuliNoP.RingType())

			// checks that one can provide custom parameters for the secret-key and error distributions
			dataWithCustomSecrets := []byte(fmt.Sprintf(`{"LogN":%d,"LogQ":[50,50],"LogP":[60],"DefaultScale":1.0,"H": 192, "Sigma": 6.6}`, tc.params.LogN()))
			var paramsWithCustomSecrets Parameters
			err = json.Unmarshal(dataWithCustomSecrets, &paramsWithCustomSecrets)
			require.Nil(t, err)
			require.Equal(t, 6.6, paramsWithCustomSecrets.Sigma())
			require.Equal(t, 192, paramsWithCustomSecrets.HammingWeight())
		})
	*/
}
