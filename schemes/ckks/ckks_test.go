package ckks

import (
	"encoding/json"
	"flag"
	"fmt"
	"math"
	"math/big"
	"runtime"
	"testing"

	"github.com/stretchr/testify/require"
	"github.com/tuneinsight/lattigo/v5/ring"

	"github.com/tuneinsight/lattigo/v5/core/rlwe"
	"github.com/tuneinsight/lattigo/v5/utils/bignum"
	"github.com/tuneinsight/lattigo/v5/utils/sampling"
)

var flagParamString = flag.String("params", "", "specify the test cryptographic parameters as a JSON string. Overrides -short and -long.")
var printPrecisionStats = flag.Bool("print-precision", false, "print precision stats")

func GetTestName(params Parameters, opname string) string {
	return fmt.Sprintf("%s/RingType=%s/logN=%d/logQP=%d/Qi=%d/Pi=%d/LogScale=%d",
		opname,
		params.RingType(),
		params.LogN(),
		int(math.Round(params.LogQP())),
		params.QCount(),
		params.PCount(),
		int(math.Log2(params.DefaultScale().Float64())))
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
	encryptorPk *rlwe.Encryptor
	encryptorSk *rlwe.Encryptor
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
		testParams = testParamsLiteral
	}

	for _, ringType := range []ring.Type{ring.Standard, ring.ConjugateInvariant} {

		for _, paramsLiteral := range testParams {

			paramsLiteral.RingType = ringType

			if testing.Short() {
				paramsLiteral.LogN = 10
			}

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
				testBridge,
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

	tc.evaluator = NewEvaluator(tc.params, rlwe.NewMemEvaluationKeySet(tc.kgen.GenRelinearizationKeyNew(tc.sk)))

	return tc, nil

}

func newTestVectors(tc *testContext, encryptor *rlwe.Encryptor, a, b complex128, t *testing.T) (values []*bignum.Complex, pt *rlwe.Plaintext, ct *rlwe.Ciphertext) {

	prec := tc.encoder.Prec()

	pt = NewPlaintext(tc.params, tc.params.MaxLevel())

	values = make([]*bignum.Complex, pt.Slots())

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
		t.Fatal("invalid ring type")
	}

	tc.encoder.Encode(values, pt)

	if encryptor != nil {
		var err error
		ct, err = encryptor.EncryptNew(pt)
		if err != nil {
			panic(err)
		}
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

func testParameters(tc *testContext, t *testing.T) {

	t.Run(GetTestName(tc.params, "Parameters/NewParameters"), func(t *testing.T) {
		params, err := NewParametersFromLiteral(ParametersLiteral{
			LogN:            4,
			LogQ:            []int{60, 60},
			LogP:            []int{60},
			LogDefaultScale: 0,
		})
		require.NoError(t, err)
		require.Equal(t, ring.Standard, params.RingType()) // Default ring type should be standard
		require.Equal(t, rlwe.DefaultXe, params.Xe())
		require.Equal(t, rlwe.DefaultXs, params.Xs())
	})

	t.Run(GetTestName(tc.params, "Parameters/StandardRing"), func(t *testing.T) {
		params, err := tc.params.StandardParameters()
		switch tc.params.RingType() {
		case ring.Standard:
			require.True(t, params.Equal(&tc.params))
			require.NoError(t, err)
		case ring.ConjugateInvariant:
			require.Equal(t, params.LogN(), tc.params.LogN()+1)
			require.NoError(t, err)
		default:
			t.Fatal("invalid RingType")
		}
	})

	t.Run(GetTestName(tc.params, "Parameters/Marshaller/Binary"), func(t *testing.T) {

		bytes, err := tc.params.MarshalBinary()
		require.Nil(t, err)
		var p Parameters
		require.Nil(t, p.UnmarshalBinary(bytes))
		require.True(t, tc.params.Equal(&p))
	})

	t.Run(GetTestName(tc.params, "Parameters/Marshaller/JSON"), func(t *testing.T) {
		// checks that parameters can be marshalled without error
		data, err := json.Marshal(tc.params)
		require.Nil(t, err)
		require.NotNil(t, data)

		// checks that ckks.Parameters can be unmarshalled without error
		var paramsRec Parameters
		err = json.Unmarshal(data, &paramsRec)
		require.Nil(t, err)
		require.True(t, tc.params.Equal(&paramsRec))

		// checks that ckks.Parameters can be unmarshalled with log-moduli definition without error
		dataWithLogModuli := []byte(fmt.Sprintf(`{"LogN":%d,"LogQ":[50,50],"LogP":[60], "LogDefaultScale":30}`, tc.params.LogN()))
		var paramsWithLogModuli Parameters
		err = json.Unmarshal(dataWithLogModuli, &paramsWithLogModuli)
		require.Nil(t, err)
		require.Equal(t, 2, paramsWithLogModuli.QCount())
		require.Equal(t, 1, paramsWithLogModuli.PCount())
		require.Equal(t, ring.Standard, paramsWithLogModuli.RingType()) // Omitting the RingType field should result in a standard instance
		require.Equal(t, rlwe.DefaultXe, paramsWithLogModuli.Xe())      // Omitting Xe should result in Default being used
		require.Equal(t, float64(1<<30), paramsWithLogModuli.DefaultScale().Float64())

		// checks that ckks.Parameters can be unmarshalled with log-moduli definition with empty P without error
		dataWithLogModuliNoP := []byte(fmt.Sprintf(`{"LogN":%d,"LogQ":[50,50],"LogP":[], "RingType": "ConjugateInvariant"}`, tc.params.LogN()))
		var paramsWithLogModuliNoP Parameters
		err = json.Unmarshal(dataWithLogModuliNoP, &paramsWithLogModuliNoP)
		require.Nil(t, err)
		require.Equal(t, 2, paramsWithLogModuliNoP.QCount())
		require.Equal(t, 0, paramsWithLogModuliNoP.PCount())
		require.Equal(t, ring.ConjugateInvariant, paramsWithLogModuliNoP.RingType())

		// checks that one can provide custom parameters for the secret-key and error distributions
		dataWithCustomSecrets := []byte(fmt.Sprintf(`{"LogN":%d,"LogQ":[50,50],"LogP":[60], "Xs": {"Type": "Ternary", "H": 192}, "Xe": {"Type": "DiscreteGaussian", "Sigma": 6.6, "Bound": 39.6}}`, tc.params.LogN()))
		var paramsWithCustomSecrets Parameters
		err = json.Unmarshal(dataWithCustomSecrets, &paramsWithCustomSecrets)
		require.Nil(t, err)
		require.Equal(t, ring.DiscreteGaussian{Sigma: 6.6, Bound: 39.6}, paramsWithCustomSecrets.Xe())
		require.Equal(t, ring.Ternary{H: 192}, paramsWithCustomSecrets.Xs())
	})
}

func testEncoder(tc *testContext, t *testing.T) {

	t.Run(GetTestName(tc.params, "Encoder/IsBatched=true"), func(t *testing.T) {

		values, plaintext, _ := newTestVectors(tc, nil, -1-1i, 1+1i, t)

		VerifyTestVectors(tc.params, tc.encoder, nil, values, plaintext, tc.params.LogDefaultScale(), 0, *printPrecisionStats, t)
	})

	logprec := float64(tc.params.LogDefaultScale()) / 2

	t.Run(GetTestName(tc.params, "Encoder/IsBatched=true/DecodePublic/[]float64"), func(t *testing.T) {

		values, plaintext, _ := newTestVectors(tc, nil, -1-1i, 1+1i, t)

		have := make([]float64, len(values))

		require.NoError(t, tc.encoder.DecodePublic(plaintext, have, logprec))

		want := make([]float64, len(values))
		for i := range want {
			want[i], _ = values[i][0].Float64()
			want[i] -= have[i]
		}

		// Allows for a 10% error over the expected standard deviation of the error
		require.GreaterOrEqual(t, StandardDeviation(want, rlwe.NewScale(1)), math.Exp2(-logprec)/math.Sqrt(12)*0.9)
	})

	t.Run(GetTestName(tc.params, "Encoder/IsBatched=true/DecodePublic/[]complex128"), func(t *testing.T) {

		if tc.params.RingType() == ring.ConjugateInvariant {
			t.Skip()
		}
		values, plaintext, _ := newTestVectors(tc, nil, -1-1i, 1+1i, t)

		have := make([]complex128, len(values))
		require.NoError(t, tc.encoder.DecodePublic(plaintext, have, logprec))

		wantReal := make([]float64, len(values))
		wantImag := make([]float64, len(values))

		for i := range have {
			wantReal[i], _ = values[i][0].Float64()
			wantImag[i], _ = values[i][1].Float64()

			wantReal[i] -= real(have[i])
			wantImag[i] -= imag(have[i])
		}

		// Allows for a 10% error over the expected standard deviation of the error
		require.GreaterOrEqual(t, StandardDeviation(wantReal, rlwe.NewScale(1)), math.Exp2(-logprec)/math.Sqrt(12)*0.9)
		require.GreaterOrEqual(t, StandardDeviation(wantImag, rlwe.NewScale(1)), math.Exp2(-logprec)/math.Sqrt(12)*0.9)
	})

	t.Run(GetTestName(tc.params, "Encoder/IsBatched=true/DecodePublic/[]big.Float"), func(t *testing.T) {
		values, plaintext, _ := newTestVectors(tc, nil, -1-1i, 1+1i, t)
		have := make([]*big.Float, len(values))
		require.NoError(t, tc.encoder.DecodePublic(plaintext, have, logprec))

		want := make([]*big.Float, len(values))
		for i := range want {
			want[i] = values[i][0].Sub(values[i][0], have[i])
		}

		// Allows for a 10% error over the expected standard deviation of the error
		require.GreaterOrEqual(t, StandardDeviation(want, rlwe.NewScale(1)), math.Exp2(-logprec)/math.Sqrt(12)*0.9)
	})

	t.Run(GetTestName(tc.params, "Encoder/IsBatched=true/DecodePublic/[]bignum.Complex"), func(t *testing.T) {
		if tc.params.RingType() == ring.ConjugateInvariant {
			t.Skip()
		}
		values, plaintext, _ := newTestVectors(tc, nil, -1-1i, 1+1i, t)
		have := make([]*bignum.Complex, len(values))
		require.NoError(t, tc.encoder.DecodePublic(plaintext, have, logprec))

		wantReal := make([]*big.Float, len(values))
		wantImag := make([]*big.Float, len(values))

		for i := range have {
			wantReal[i] = values[i][0].Sub(values[i][0], have[i][0])
			wantImag[i] = values[i][1].Sub(values[i][1], have[i][1])
		}

		// Allows for a 10% error over the expected standard deviation of the error
		require.GreaterOrEqual(t, StandardDeviation(wantReal, rlwe.NewScale(1)), math.Exp2(-logprec)/math.Sqrt(12)*0.9)
		require.GreaterOrEqual(t, StandardDeviation(wantImag, rlwe.NewScale(1)), math.Exp2(-logprec)/math.Sqrt(12)*0.9)
	})

	t.Run(GetTestName(tc.params, "Encoder/IsBatched=false"), func(t *testing.T) {

		slots := tc.params.N()

		valuesWant := make([]float64, slots)

		for i := 0; i < slots; i++ {
			valuesWant[i] = sampling.RandFloat64(-1, 1)
		}

		valuesWant[0] = 0.607538

		pt := NewPlaintext(tc.params, tc.params.MaxLevel())
		pt.IsBatched = false

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

		minPrec := math.Log2(tc.params.DefaultScale().Float64()) - float64(tc.params.LogN()+2)
		if minPrec < 0 {
			minPrec = 0
		}

		require.GreaterOrEqual(t, math.Log2(1/meanprec), minPrec)

		// Also tests at level 0
		pt = NewPlaintext(tc.params, tc.params.LevelsConsumedPerRescaling()-1)
		pt.IsBatched = false

		tc.encoder.Encode(valuesWant, pt)

		tc.encoder.Decode(pt, valuesTest)

		meanprec = 0
		for i := range valuesWant {
			meanprec += math.Abs(valuesTest[i] - valuesWant[i])
		}

		meanprec /= float64(slots)

		if *printPrecisionStats {
			t.Logf("\nMean    precision : %.2f \n", math.Log2(1/meanprec))
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

		ciphertext3, err := tc.evaluator.AddNew(ciphertext1, ciphertext2)
		require.NoError(t, err)

		VerifyTestVectors(tc.params, tc.encoder, tc.decryptor, values1, ciphertext3, tc.params.LogDefaultScale(), 0, *printPrecisionStats, t)
	})

	t.Run(GetTestName(tc.params, "Evaluator/Add/Ct"), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectors(tc, tc.encryptorSk, -1-1i, 1+1i, t)
		values2, _, ciphertext2 := newTestVectors(tc, tc.encryptorSk, -1-1i, 1+1i, t)

		for i := range values1 {
			values1[i].Add(values1[i], values2[i])
		}

		require.NoError(t, tc.evaluator.Add(ciphertext1, ciphertext2, ciphertext1))

		VerifyTestVectors(tc.params, tc.encoder, tc.decryptor, values1, ciphertext1, tc.params.LogDefaultScale(), 0, *printPrecisionStats, t)
	})

	t.Run(GetTestName(tc.params, "Evaluator/Add/Pt"), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectors(tc, tc.encryptorSk, -1-1i, 1+1i, t)
		values2, plaintext2, _ := newTestVectors(tc, tc.encryptorSk, -1-1i, 1+1i, t)

		for i := range values1 {
			values1[i].Add(values1[i], values2[i])
		}

		require.NoError(t, tc.evaluator.Add(ciphertext1, plaintext2, ciphertext1))

		VerifyTestVectors(tc.params, tc.encoder, tc.decryptor, values1, ciphertext1, tc.params.LogDefaultScale(), 0, *printPrecisionStats, t)
	})

	t.Run(GetTestName(tc.params, "Evaluator/Add/Scalar"), func(t *testing.T) {

		values, _, ciphertext := newTestVectors(tc, tc.encryptorSk, -1-1i, 1+1i, t)

		constant := randomConst(tc.params.RingType(), tc.encoder.Prec(), -1+1i, -1+1i)

		for i := range values {
			values[i].Add(values[i], constant)
		}

		require.NoError(t, tc.evaluator.Add(ciphertext, constant, ciphertext))

		VerifyTestVectors(tc.params, tc.encoder, tc.decryptor, values, ciphertext, tc.params.LogDefaultScale(), 0, *printPrecisionStats, t)
	})

	t.Run(GetTestName(tc.params, "Evaluator/Add/Vector"), func(t *testing.T) {

		values1, _, ciphertext := newTestVectors(tc, tc.encryptorSk, -1-1i, 1+1i, t)
		values2, _, _ := newTestVectors(tc, tc.encryptorSk, -1-1i, 1+1i, t)

		for i := range values1 {
			values1[i].Add(values1[i], values2[i])
		}

		require.NoError(t, tc.evaluator.Add(ciphertext, values2, ciphertext))

		VerifyTestVectors(tc.params, tc.encoder, tc.decryptor, values1, ciphertext, tc.params.LogDefaultScale(), 0, *printPrecisionStats, t)
	})
}

func testEvaluatorSub(tc *testContext, t *testing.T) {

	t.Run(GetTestName(tc.params, "Evaluator/SubNew/Ct"), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectors(tc, tc.encryptorSk, -1-1i, 1+1i, t)
		values2, _, ciphertext2 := newTestVectors(tc, tc.encryptorSk, -1-1i, 1+1i, t)

		for i := range values1 {
			values1[i].Sub(values1[i], values2[i])
		}

		ciphertext3, err := tc.evaluator.SubNew(ciphertext1, ciphertext2)
		require.NoError(t, err)

		VerifyTestVectors(tc.params, tc.encoder, tc.decryptor, values1, ciphertext3, tc.params.LogDefaultScale(), 0, *printPrecisionStats, t)
	})

	t.Run(GetTestName(tc.params, "Evaluator/Sub/Ct"), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectors(tc, tc.encryptorSk, -1-1i, 1+1i, t)
		values2, _, ciphertext2 := newTestVectors(tc, tc.encryptorSk, -1-1i, 1+1i, t)

		for i := range values1 {
			values1[i].Sub(values1[i], values2[i])
		}

		require.NoError(t, tc.evaluator.Sub(ciphertext1, ciphertext2, ciphertext1))

		VerifyTestVectors(tc.params, tc.encoder, tc.decryptor, values1, ciphertext1, tc.params.LogDefaultScale(), 0, *printPrecisionStats, t)
	})

	t.Run(GetTestName(tc.params, "Evaluator/Sub/Pt"), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectors(tc, tc.encryptorSk, -1-1i, 1+1i, t)
		values2, plaintext2, ciphertext2 := newTestVectors(tc, tc.encryptorSk, -1-1i, 1+1i, t)

		valuesTest := make([]*bignum.Complex, len(values1))
		for i := range values1 {
			valuesTest[i] = bignum.NewComplex()
			valuesTest[i].Sub(values1[i], values2[i])
		}

		require.NoError(t, tc.evaluator.Sub(ciphertext1, plaintext2, ciphertext2))

		VerifyTestVectors(tc.params, tc.encoder, tc.decryptor, valuesTest, ciphertext2, tc.params.LogDefaultScale(), 0, *printPrecisionStats, t)
	})

	t.Run(GetTestName(tc.params, "Evaluator/Sub/Scalar"), func(t *testing.T) {

		values, _, ciphertext := newTestVectors(tc, tc.encryptorSk, -1-1i, 1+1i, t)

		constant := randomConst(tc.params.RingType(), tc.encoder.Prec(), -1+1i, -1+1i)

		for i := range values {
			values[i].Sub(values[i], constant)
		}

		require.NoError(t, tc.evaluator.Sub(ciphertext, constant, ciphertext))

		VerifyTestVectors(tc.params, tc.encoder, tc.decryptor, values, ciphertext, tc.params.LogDefaultScale(), 0, *printPrecisionStats, t)
	})

	t.Run(GetTestName(tc.params, "Evaluator/Sub/Vector"), func(t *testing.T) {

		values1, _, ciphertext := newTestVectors(tc, tc.encryptorSk, -1-1i, 1+1i, t)
		values2, _, _ := newTestVectors(tc, tc.encryptorSk, -1-1i, 1+1i, t)

		for i := range values1 {
			values1[i].Sub(values1[i], values2[i])
		}

		require.NoError(t, tc.evaluator.Sub(ciphertext, values2, ciphertext))

		VerifyTestVectors(tc.params, tc.encoder, tc.decryptor, values1, ciphertext, tc.params.LogDefaultScale(), 0, *printPrecisionStats, t)
	})
}

func testEvaluatorRescale(tc *testContext, t *testing.T) {

	t.Run(GetTestName(tc.params, "Evaluator/RescaleTo/Single"), func(t *testing.T) {

		if tc.params.MaxLevel() < 2 {
			t.Skip("skipping test for params max level < 2")
		}

		values, _, ciphertext := newTestVectors(tc, tc.encryptorSk, -1-1i, 1+1i, t)

		constant := tc.ringQ.SubRings[ciphertext.Level()].Modulus

		require.NoError(t, tc.evaluator.Mul(ciphertext, constant, ciphertext))

		ciphertext.Scale = ciphertext.Scale.Mul(rlwe.NewScale(constant))

		if err := tc.evaluator.RescaleTo(ciphertext, tc.params.DefaultScale(), ciphertext); err != nil {
			t.Fatal(err)
		}

		VerifyTestVectors(tc.params, tc.encoder, tc.decryptor, values, ciphertext, tc.params.LogDefaultScale(), 0, *printPrecisionStats, t)
	})

	t.Run(GetTestName(tc.params, "Evaluator/RescaleTo/Many"), func(t *testing.T) {

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
			require.NoError(t, tc.evaluator.Mul(ciphertext, constant, ciphertext))
			ciphertext.Scale = ciphertext.Scale.Mul(rlwe.NewScale(constant))
		}

		if err := tc.evaluator.RescaleTo(ciphertext, tc.params.DefaultScale(), ciphertext); err != nil {
			t.Fatal(err)
		}

		VerifyTestVectors(tc.params, tc.encoder, tc.decryptor, values, ciphertext, tc.params.LogDefaultScale(), 0, *printPrecisionStats, t)
	})
}

func testEvaluatorMul(tc *testContext, t *testing.T) {

	t.Run(GetTestName(tc.params, "Evaluator/MulNew/Ct/Pt"), func(t *testing.T) {

		values1, plaintext1, ciphertext1 := newTestVectors(tc, tc.encryptorSk, -1-1i, 1+1i, t)

		mul := bignum.NewComplexMultiplier()

		for i := range values1 {
			mul.Mul(values1[i], values1[i], values1[i])
		}

		ciphertext2, err := tc.evaluator.MulNew(ciphertext1, plaintext1)
		require.NoError(t, err)

		VerifyTestVectors(tc.params, tc.encoder, tc.decryptor, values1, ciphertext2, tc.params.LogDefaultScale(), 0, *printPrecisionStats, t)
	})

	t.Run(GetTestName(tc.params, "Evaluator/Mul/Ct/Scalar"), func(t *testing.T) {

		values, _, ciphertext := newTestVectors(tc, tc.encryptorSk, -1-1i, 1+1i, t)

		constant := randomConst(tc.params.RingType(), tc.encoder.Prec(), -1+1i, -1+1i)

		mul := bignum.NewComplexMultiplier()

		for i := range values {
			mul.Mul(values[i], constant, values[i])
		}

		require.NoError(t, tc.evaluator.Mul(ciphertext, constant, ciphertext))

		VerifyTestVectors(tc.params, tc.encoder, tc.decryptor, values, ciphertext, tc.params.LogDefaultScale(), 0, *printPrecisionStats, t)
	})

	t.Run(GetTestName(tc.params, "Evaluator/Mul/Ct/Vector"), func(t *testing.T) {

		values1, _, ciphertext := newTestVectors(tc, tc.encryptorSk, -1-1i, 1+1i, t)
		values2, _, _ := newTestVectors(tc, tc.encryptorSk, -1-1i, 1+1i, t)

		mul := bignum.NewComplexMultiplier()

		for i := range values1 {
			mul.Mul(values1[i], values2[i], values1[i])
		}

		tc.evaluator.Mul(ciphertext, values2, ciphertext)

		VerifyTestVectors(tc.params, tc.encoder, tc.decryptor, values1, ciphertext, tc.params.LogDefaultScale(), 0, *printPrecisionStats, t)
	})

	t.Run(GetTestName(tc.params, "Evaluator/Mul/Ct/Pt"), func(t *testing.T) {

		values1, plaintext1, ciphertext1 := newTestVectors(tc, tc.encryptorSk, -1-1i, 1+1i, t)

		mul := bignum.NewComplexMultiplier()

		for i := range values1 {
			mul.Mul(values1[i], values1[i], values1[i])
		}

		require.NoError(t, tc.evaluator.MulRelin(ciphertext1, plaintext1, ciphertext1))

		VerifyTestVectors(tc.params, tc.encoder, tc.decryptor, values1, ciphertext1, tc.params.LogDefaultScale(), 0, *printPrecisionStats, t)
	})

	t.Run(GetTestName(tc.params, "Evaluator/Mul/Ct/Ct/Degree0"), func(t *testing.T) {

		values1, plaintext1, _ := newTestVectors(tc, tc.encryptorSk, -1-1i, 1+1i, t)
		values2, _, ciphertext2 := newTestVectors(tc, tc.encryptorSk, -1-1i, 1+1i, t)

		mul := bignum.NewComplexMultiplier()

		for i := range values1 {
			mul.Mul(values2[i], values1[i], values2[i])
		}

		ciphertext1 := &rlwe.Ciphertext{}
		ciphertext1.Value = []ring.Poly{plaintext1.Value}
		ciphertext1.MetaData = plaintext1.MetaData.CopyNew()

		require.NoError(t, tc.evaluator.MulRelin(ciphertext1, ciphertext2, ciphertext1))

		VerifyTestVectors(tc.params, tc.encoder, tc.decryptor, values2, ciphertext1, tc.params.LogDefaultScale(), 0, *printPrecisionStats, t)
	})

	t.Run(GetTestName(tc.params, "Evaluator/MulRelin/Ct/Ct"), func(t *testing.T) {

		// op0 <- op0 * op1
		values1, _, ciphertext1 := newTestVectors(tc, tc.encryptorSk, -1-1i, 1+1i, t)
		values2, _, ciphertext2 := newTestVectors(tc, tc.encryptorSk, -1-1i, 1+1i, t)

		mul := bignum.NewComplexMultiplier()

		for i := range values1 {
			mul.Mul(values1[i], values2[i], values1[i])
		}

		require.NoError(t, tc.evaluator.MulRelin(ciphertext1, ciphertext2, ciphertext1))
		require.Equal(t, ciphertext1.Degree(), 1)

		VerifyTestVectors(tc.params, tc.encoder, tc.decryptor, values1, ciphertext1, tc.params.LogDefaultScale(), 0, *printPrecisionStats, t)

		// op1 <- op0 * op1
		values1, _, ciphertext1 = newTestVectors(tc, tc.encryptorSk, -1-1i, 1+1i, t)
		values2, _, ciphertext2 = newTestVectors(tc, tc.encryptorSk, -1-1i, 1+1i, t)

		for i := range values1 {
			mul.Mul(values2[i], values1[i], values2[i])
		}

		require.NoError(t, tc.evaluator.MulRelin(ciphertext1, ciphertext2, ciphertext2))
		require.Equal(t, ciphertext2.Degree(), 1)

		VerifyTestVectors(tc.params, tc.encoder, tc.decryptor, values2, ciphertext2, tc.params.LogDefaultScale(), 0, *printPrecisionStats, t)

		// op0 <- op0 * op0
		for i := range values1 {
			mul.Mul(values1[i], values1[i], values1[i])
		}

		require.NoError(t, tc.evaluator.MulRelin(ciphertext1, ciphertext1, ciphertext1))
		require.Equal(t, ciphertext1.Degree(), 1)

		VerifyTestVectors(tc.params, tc.encoder, tc.decryptor, values1, ciphertext1, tc.params.LogDefaultScale(), 0, *printPrecisionStats, t)
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

		require.NoError(t, tc.evaluator.MulThenAdd(ciphertext1, constant, ciphertext2))

		VerifyTestVectors(tc.params, tc.encoder, tc.decryptor, values2, ciphertext2, tc.params.LogDefaultScale(), 0, *printPrecisionStats, t)
	})

	t.Run(GetTestName(tc.params, "Evaluator/MulThenAdd/Vector"), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectors(tc, tc.encryptorSk, -1, 1, t)
		values2, _, ciphertext2 := newTestVectors(tc, tc.encryptorSk, -1, 1, t)

		require.NoError(t, tc.evaluator.MulThenAdd(ciphertext2, values1, ciphertext1))

		mul := bignum.NewComplexMultiplier()

		tmp := new(bignum.Complex)
		tmp[0] = new(big.Float)
		tmp[1] = new(big.Float)

		for i := range values1 {
			mul.Mul(values2[i], values1[i], tmp)
			values1[i].Add(values1[i], tmp)
		}

		require.Equal(t, ciphertext1.Degree(), 1)

		VerifyTestVectors(tc.params, tc.encoder, tc.decryptor, values1, ciphertext1, tc.params.LogDefaultScale(), 0, *printPrecisionStats, t)
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

		require.NoError(t, tc.evaluator.MulThenAdd(ciphertext2, plaintext1, ciphertext1))

		require.Equal(t, ciphertext1.Degree(), 1)

		VerifyTestVectors(tc.params, tc.encoder, tc.decryptor, values1, ciphertext1, tc.params.LogDefaultScale(), 0, *printPrecisionStats, t)
	})

	t.Run(GetTestName(tc.params, "Evaluator/MulRelinThenAdd/Ct"), func(t *testing.T) {

		// opOut = opOut + op1 * op0
		values1, _, ciphertext1 := newTestVectors(tc, tc.encryptorSk, -1-1i, 1+1i, t)
		values2, _, ciphertext2 := newTestVectors(tc, tc.encryptorSk, -1-1i, 1+1i, t)

		mul := bignum.NewComplexMultiplier()

		for i := range values1 {
			mul.Mul(values1[i], values2[i], values2[i])
		}

		ciphertext3 := NewCiphertext(tc.params, 2, ciphertext1.Level())

		ciphertext3.Scale = ciphertext1.Scale.Mul(ciphertext2.Scale)

		require.NoError(t, tc.evaluator.MulThenAdd(ciphertext1, ciphertext2, ciphertext3))

		require.Equal(t, ciphertext3.Degree(), 2)

		require.NoError(t, tc.evaluator.Relinearize(ciphertext3, ciphertext3))

		require.Equal(t, ciphertext3.Degree(), 1)

		VerifyTestVectors(tc.params, tc.encoder, tc.decryptor, values2, ciphertext3, tc.params.LogDefaultScale(), 0, *printPrecisionStats, t)

		// op1 = op1 + op0*op0
		values1, _, ciphertext1 = newTestVectors(tc, tc.encryptorSk, -1-1i, 1+1i, t)
		values2, _, ciphertext2 = newTestVectors(tc, tc.encryptorSk, -1-1i, 1+1i, t)

		tmp := bignum.NewComplex()
		for i := range values1 {
			mul.Mul(values2[i], values2[i], tmp)
			values1[i].Add(values1[i], tmp)
		}

		require.NoError(t, tc.evaluator.MulRelinThenAdd(ciphertext2, ciphertext2, ciphertext1))

		require.Equal(t, ciphertext1.Degree(), 1)

		VerifyTestVectors(tc.params, tc.encoder, tc.decryptor, values1, ciphertext1, tc.params.LogDefaultScale(), 0, *printPrecisionStats, t)
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

		VerifyTestVectors(stdParams, stdEncoder, stdDecryptor, values, stdCTHave, tc.params.LogDefaultScale(), 0, *printPrecisionStats, t)

		stdCTImag, err := stdEvaluator.MulNew(stdCTHave, 1i)
		require.NoError(t, err)
		require.NoError(t, stdEvaluator.Add(stdCTHave, stdCTImag, stdCTHave))

		ciCTHave := NewCiphertext(ciParams, 1, stdCTHave.Level())
		switcher.ComplexToReal(evalStandar, stdCTHave, ciCTHave)

		VerifyTestVectors(tc.params, tc.encoder, tc.decryptor, values, ciCTHave, tc.params.LogDefaultScale(), 0, *printPrecisionStats, t)
	})
}
