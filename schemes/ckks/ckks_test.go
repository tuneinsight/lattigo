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

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/utils/bignum"
	"github.com/tuneinsight/lattigo/v6/utils/sampling"
)

var flagParamString = flag.String("params", "", "specify the test cryptographic parameters as a JSON string. Overrides -short and -long.")
var printPrecisionStats = flag.Bool("print-precision", false, "print precision stats")

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
		testParams = testParametersLiteral
	}

	for _, ringType := range []ring.Type{ring.Standard, ring.ConjugateInvariant} {

		for _, paramsLiteral := range testParams {

			paramsLiteral.RingType = ringType

			if testing.Short() {
				paramsLiteral.LogN = 10
			}

			tc := NewTestContext(paramsLiteral)

			for _, testSet := range []func(tc *TestContext, t *testing.T){
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

func testParameters(tc *TestContext, t *testing.T) {

	t.Run(name("Parameters/NewParameters", tc), func(t *testing.T) {
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

	t.Run(name("Parameters/StandardRing", tc), func(t *testing.T) {
		params, err := tc.Params.StandardParameters()
		switch tc.Params.RingType() {
		case ring.Standard:
			require.True(t, params.Equal(&tc.Params))
			require.NoError(t, err)
		case ring.ConjugateInvariant:
			require.Equal(t, params.LogN(), tc.Params.LogN()+1)
			require.NoError(t, err)
		default:
			t.Fatal("invalid RingType")
		}
	})

	t.Run(name("Parameters/Marshaller/Binary", tc), func(t *testing.T) {

		bytes, err := tc.Params.MarshalBinary()
		require.Nil(t, err)
		var p Parameters
		require.Nil(t, p.UnmarshalBinary(bytes))
		require.True(t, tc.Params.Equal(&p))
	})

	t.Run(name("Parameters/Marshaller/JSON", tc), func(t *testing.T) {
		// checks that parameters can be marshalled without error
		data, err := json.Marshal(tc.Params)
		require.Nil(t, err)
		require.NotNil(t, data)

		// checks that ckks.Parameters can be unmarshalled without error
		var paramsRec Parameters
		err = json.Unmarshal(data, &paramsRec)
		require.Nil(t, err)
		require.True(t, tc.Params.Equal(&paramsRec))

		// checks that ckks.Parameters can be unmarshalled with log-moduli definition without error
		dataWithLogModuli := []byte(fmt.Sprintf(`{"LogN":%d,"LogQ":[50,50],"LogP":[60], "LogDefaultScale":30}`, tc.Params.LogN()))
		var paramsWithLogModuli Parameters
		err = json.Unmarshal(dataWithLogModuli, &paramsWithLogModuli)
		require.Nil(t, err)
		require.Equal(t, 2, paramsWithLogModuli.QCount())
		require.Equal(t, 1, paramsWithLogModuli.PCount())
		require.Equal(t, ring.Standard, paramsWithLogModuli.RingType()) // Omitting the RingType field should result in a standard instance
		require.Equal(t, rlwe.DefaultXe, paramsWithLogModuli.Xe())      // Omitting Xe should result in Default being used
		require.Equal(t, float64(1<<30), paramsWithLogModuli.DefaultScale().Float64())

		// checks that ckks.Parameters can be unmarshalled with log-moduli definition with empty P without error
		dataWithLogModuliNoP := []byte(fmt.Sprintf(`{"LogN":%d,"LogQ":[50,50],"LogP":[], "RingType": "ConjugateInvariant"}`, tc.Params.LogN()))
		var paramsWithLogModuliNoP Parameters
		err = json.Unmarshal(dataWithLogModuliNoP, &paramsWithLogModuliNoP)
		require.Nil(t, err)
		require.Equal(t, 2, paramsWithLogModuliNoP.QCount())
		require.Equal(t, 0, paramsWithLogModuliNoP.PCount())
		require.Equal(t, ring.ConjugateInvariant, paramsWithLogModuliNoP.RingType())

		// checks that one can provide custom parameters for the secret-key and error distributions
		dataWithCustomSecrets := []byte(fmt.Sprintf(`{"LogN":%d,"LogQ":[50,50],"LogP":[60], "Xs": {"Type": "Ternary", "H": 192}, "Xe": {"Type": "DiscreteGaussian", "Sigma": 6.6, "Bound": 39.6}}`, tc.Params.LogN()))
		var paramsWithCustomSecrets Parameters
		err = json.Unmarshal(dataWithCustomSecrets, &paramsWithCustomSecrets)
		require.Nil(t, err)
		require.Equal(t, ring.DiscreteGaussian{Sigma: 6.6, Bound: 39.6}, paramsWithCustomSecrets.Xe())
		require.Equal(t, ring.Ternary{H: 192}, paramsWithCustomSecrets.Xs())
	})
}

func testEncoder(tc *TestContext, t *testing.T) {

	t.Run(name("Encoder/IsBatched=true", tc), func(t *testing.T) {

		values, plaintext, _ := tc.NewTestVector(-1-1i, 1+1i)

		VerifyTestVectors(tc.Params, tc.Ecd, nil, values, plaintext, tc.Params.LogDefaultScale(), 0, *printPrecisionStats, t)
	})

	logprec := float64(tc.Params.LogDefaultScale()) / 2

	t.Run(name("Encoder/IsBatched=true/DecodePublic/[]float64", tc), func(t *testing.T) {

		values, plaintext, _ := tc.NewTestVector(-1-1i, 1+1i)

		have := make([]float64, len(values))

		require.NoError(t, tc.Ecd.DecodePublic(plaintext, have, logprec))

		want := make([]float64, len(values))
		for i := range want {
			want[i], _ = values[i][0].Float64()
			want[i] -= have[i]
		}

		// Allows for a 10% error over the expected standard deviation of the error
		require.GreaterOrEqual(t, StandardDeviation(want, rlwe.NewScale(1)), math.Exp2(-logprec)/math.Sqrt(12)*0.9)
	})

	t.Run(name("Encoder/IsBatched=true/DecodePublic/[]complex128", tc), func(t *testing.T) {

		if tc.Params.RingType() == ring.ConjugateInvariant {
			t.Skip()
		}
		values, plaintext, _ := tc.NewTestVector(-1-1i, 1+1i)

		have := make([]complex128, len(values))
		require.NoError(t, tc.Ecd.DecodePublic(plaintext, have, logprec))

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

	t.Run(name("Encoder/IsBatched=true/DecodePublic/[]big.Float", tc), func(t *testing.T) {
		values, plaintext, _ := tc.NewTestVector(-1-1i, 1+1i)
		have := make([]*big.Float, len(values))
		require.NoError(t, tc.Ecd.DecodePublic(plaintext, have, logprec))

		want := make([]*big.Float, len(values))
		for i := range want {
			want[i] = values[i][0].Sub(values[i][0], have[i])
		}

		// Allows for a 10% error over the expected standard deviation of the error
		require.GreaterOrEqual(t, StandardDeviation(want, rlwe.NewScale(1)), math.Exp2(-logprec)/math.Sqrt(12)*0.9)
	})

	t.Run(name("Encoder/IsBatched=true/DecodePublic/[]bignum.Complex", tc), func(t *testing.T) {
		if tc.Params.RingType() == ring.ConjugateInvariant {
			t.Skip()
		}
		values, plaintext, _ := tc.NewTestVector(-1-1i, 1+1i)
		have := make([]*bignum.Complex, len(values))
		require.NoError(t, tc.Ecd.DecodePublic(plaintext, have, logprec))

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

	t.Run(name("Encoder/IsBatched=false", tc), func(t *testing.T) {

		slots := tc.Params.N()

		valuesWant := make([]float64, slots)

		for i := 0; i < slots; i++ {
			valuesWant[i] = sampling.RandFloat64(-1, 1)
		}

		valuesWant[0] = 0.607538

		pt := NewPlaintext(tc.Params, tc.Params.MaxLevel())
		pt.IsBatched = false

		tc.Ecd.Encode(valuesWant, pt)

		valuesTest := make([]float64, len(valuesWant))

		tc.Ecd.Decode(pt, valuesTest)

		var meanprec float64

		for i := range valuesWant {
			meanprec += math.Abs(valuesTest[i] - valuesWant[i])
		}

		meanprec /= float64(slots)

		if *printPrecisionStats {
			t.Logf("\nMean    precision : %.2f \n", math.Log2(1/meanprec))
		}

		minPrec := math.Log2(tc.Params.DefaultScale().Float64()) - float64(tc.Params.LogN()+2)
		if minPrec < 0 {
			minPrec = 0
		}

		require.GreaterOrEqual(t, math.Log2(1/meanprec), minPrec)

		// Also tests at level 0
		pt = NewPlaintext(tc.Params, tc.Params.LevelsConsumedPerRescaling()-1)
		pt.IsBatched = false

		tc.Ecd.Encode(valuesWant, pt)

		tc.Ecd.Decode(pt, valuesTest)

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

func testEvaluatorAdd(tc *TestContext, t *testing.T) {

	t.Run(name("Evaluator/AddNew/Ct", tc), func(t *testing.T) {

		values1, _, ciphertext1 := tc.NewTestVector(-1-1i, 1+1i)
		values2, _, ciphertext2 := tc.NewTestVector(-1-1i, 1+1i)

		for i := range values1 {
			values1[i].Add(values1[i], values2[i])
		}

		ciphertext3, err := tc.Evl.AddNew(ciphertext1, ciphertext2)
		require.NoError(t, err)

		VerifyTestVectors(tc.Params, tc.Ecd, tc.Dec, values1, ciphertext3, tc.Params.LogDefaultScale(), 0, *printPrecisionStats, t)
	})

	t.Run(name("Evaluator/Add/Ct", tc), func(t *testing.T) {

		values1, _, ciphertext1 := tc.NewTestVector(-1-1i, 1+1i)
		values2, _, ciphertext2 := tc.NewTestVector(-1-1i, 1+1i)

		for i := range values1 {
			values1[i].Add(values1[i], values2[i])
		}

		require.NoError(t, tc.Evl.Add(ciphertext1, ciphertext2, ciphertext1))

		VerifyTestVectors(tc.Params, tc.Ecd, tc.Dec, values1, ciphertext1, tc.Params.LogDefaultScale(), 0, *printPrecisionStats, t)
	})

	t.Run(name("Evaluator/Add/Pt", tc), func(t *testing.T) {

		values1, _, ciphertext1 := tc.NewTestVector(-1-1i, 1+1i)
		values2, plaintext2, _ := tc.NewTestVector(-1-1i, 1+1i)

		for i := range values1 {
			values1[i].Add(values1[i], values2[i])
		}

		require.NoError(t, tc.Evl.Add(ciphertext1, plaintext2, ciphertext1))

		VerifyTestVectors(tc.Params, tc.Ecd, tc.Dec, values1, ciphertext1, tc.Params.LogDefaultScale(), 0, *printPrecisionStats, t)
	})

	t.Run(name("Evaluator/Add/Scalar", tc), func(t *testing.T) {

		values, _, ciphertext := tc.NewTestVector(-1-1i, 1+1i)

		constant := randomConst(tc.Params.RingType(), tc.Ecd.Prec(), -1+1i, -1+1i)

		for i := range values {
			values[i].Add(values[i], constant)
		}

		require.NoError(t, tc.Evl.Add(ciphertext, constant, ciphertext))

		VerifyTestVectors(tc.Params, tc.Ecd, tc.Dec, values, ciphertext, tc.Params.LogDefaultScale(), 0, *printPrecisionStats, t)
	})

	t.Run(name("Evaluator/Add/Vector", tc), func(t *testing.T) {

		values1, _, ciphertext := tc.NewTestVector(-1-1i, 1+1i)
		values2, _, _ := tc.NewTestVector(-1-1i, 1+1i)

		for i := range values1 {
			values1[i].Add(values1[i], values2[i])
		}

		require.NoError(t, tc.Evl.Add(ciphertext, values2, ciphertext))

		VerifyTestVectors(tc.Params, tc.Ecd, tc.Dec, values1, ciphertext, tc.Params.LogDefaultScale(), 0, *printPrecisionStats, t)
	})
}

func testEvaluatorSub(tc *TestContext, t *testing.T) {

	t.Run(name("Evaluator/SubNew/Ct", tc), func(t *testing.T) {

		values1, _, ciphertext1 := tc.NewTestVector(-1-1i, 1+1i)
		values2, _, ciphertext2 := tc.NewTestVector(-1-1i, 1+1i)

		for i := range values1 {
			values1[i].Sub(values1[i], values2[i])
		}

		ciphertext3, err := tc.Evl.SubNew(ciphertext1, ciphertext2)
		require.NoError(t, err)

		VerifyTestVectors(tc.Params, tc.Ecd, tc.Dec, values1, ciphertext3, tc.Params.LogDefaultScale(), 0, *printPrecisionStats, t)
	})

	t.Run(name("Evaluator/Sub/Ct", tc), func(t *testing.T) {

		values1, _, ciphertext1 := tc.NewTestVector(-1-1i, 1+1i)
		values2, _, ciphertext2 := tc.NewTestVector(-1-1i, 1+1i)

		for i := range values1 {
			values1[i].Sub(values1[i], values2[i])
		}

		require.NoError(t, tc.Evl.Sub(ciphertext1, ciphertext2, ciphertext1))

		VerifyTestVectors(tc.Params, tc.Ecd, tc.Dec, values1, ciphertext1, tc.Params.LogDefaultScale(), 0, *printPrecisionStats, t)
	})

	t.Run(name("Evaluator/Sub/Pt", tc), func(t *testing.T) {

		values1, _, ciphertext1 := tc.NewTestVector(-1-1i, 1+1i)
		values2, plaintext2, ciphertext2 := tc.NewTestVector(-1-1i, 1+1i)

		valuesTest := make([]*bignum.Complex, len(values1))
		for i := range values1 {
			valuesTest[i] = bignum.NewComplex()
			valuesTest[i].Sub(values1[i], values2[i])
		}

		require.NoError(t, tc.Evl.Sub(ciphertext1, plaintext2, ciphertext2))

		VerifyTestVectors(tc.Params, tc.Ecd, tc.Dec, valuesTest, ciphertext2, tc.Params.LogDefaultScale(), 0, *printPrecisionStats, t)
	})

	t.Run(name("Evaluator/Sub/Scalar", tc), func(t *testing.T) {

		values, _, ciphertext := tc.NewTestVector(-1-1i, 1+1i)

		constant := randomConst(tc.Params.RingType(), tc.Ecd.Prec(), -1+1i, -1+1i)

		for i := range values {
			values[i].Sub(values[i], constant)
		}

		require.NoError(t, tc.Evl.Sub(ciphertext, constant, ciphertext))

		VerifyTestVectors(tc.Params, tc.Ecd, tc.Dec, values, ciphertext, tc.Params.LogDefaultScale(), 0, *printPrecisionStats, t)
	})

	t.Run(name("Evaluator/Sub/Vector", tc), func(t *testing.T) {

		values1, _, ciphertext := tc.NewTestVector(-1-1i, 1+1i)
		values2, _, _ := tc.NewTestVector(-1-1i, 1+1i)

		for i := range values1 {
			values1[i].Sub(values1[i], values2[i])
		}

		require.NoError(t, tc.Evl.Sub(ciphertext, values2, ciphertext))

		VerifyTestVectors(tc.Params, tc.Ecd, tc.Dec, values1, ciphertext, tc.Params.LogDefaultScale(), 0, *printPrecisionStats, t)
	})
}

func testEvaluatorRescale(tc *TestContext, t *testing.T) {

	t.Run(name("Evaluator/RescaleTo/Single", tc), func(t *testing.T) {

		if tc.Params.MaxLevel() < 2 {
			t.Skip("skipping test for params max level < 2")
		}

		values, _, ciphertext := tc.NewTestVector(-1-1i, 1+1i)

		constant := tc.Params.RingQ().SubRings[ciphertext.Level()].Modulus

		require.NoError(t, tc.Evl.Mul(ciphertext, constant, ciphertext))

		ciphertext.Scale = ciphertext.Scale.Mul(rlwe.NewScale(constant))

		if err := tc.Evl.RescaleTo(ciphertext, tc.Params.DefaultScale(), ciphertext); err != nil {
			t.Fatal(err)
		}

		VerifyTestVectors(tc.Params, tc.Ecd, tc.Dec, values, ciphertext, tc.Params.LogDefaultScale(), 0, *printPrecisionStats, t)
	})

	t.Run(name("Evaluator/RescaleTo/Many", tc), func(t *testing.T) {

		if tc.Params.MaxLevel() < 2 {
			t.Skip("skipping test for params max level < 2")
		}

		values, _, ciphertext := tc.NewTestVector(-1-1i, 1+1i)

		nbRescales := tc.Params.MaxLevel()
		if nbRescales > 5 {
			nbRescales = 5
		}

		for i := 0; i < nbRescales; i++ {
			constant := tc.Params.RingQ().SubRings[ciphertext.Level()-i].Modulus
			require.NoError(t, tc.Evl.Mul(ciphertext, constant, ciphertext))
			ciphertext.Scale = ciphertext.Scale.Mul(rlwe.NewScale(constant))
		}

		if err := tc.Evl.RescaleTo(ciphertext, tc.Params.DefaultScale(), ciphertext); err != nil {
			t.Fatal(err)
		}

		VerifyTestVectors(tc.Params, tc.Ecd, tc.Dec, values, ciphertext, tc.Params.LogDefaultScale(), 0, *printPrecisionStats, t)
	})
}

func testEvaluatorMul(tc *TestContext, t *testing.T) {

	t.Run(name("Evaluator/MulNew/Ct/Pt", tc), func(t *testing.T) {

		values1, plaintext1, ciphertext1 := tc.NewTestVector(-1-1i, 1+1i)

		mul := bignum.NewComplexMultiplier()

		for i := range values1 {
			mul.Mul(values1[i], values1[i], values1[i])
		}

		ciphertext2, err := tc.Evl.MulNew(ciphertext1, plaintext1)
		require.NoError(t, err)

		VerifyTestVectors(tc.Params, tc.Ecd, tc.Dec, values1, ciphertext2, tc.Params.LogDefaultScale(), 0, *printPrecisionStats, t)
	})

	t.Run(name("Evaluator/Mul/Ct/Scalar", tc), func(t *testing.T) {

		values, _, ciphertext := tc.NewTestVector(-1-1i, 1+1i)

		constant := randomConst(tc.Params.RingType(), tc.Ecd.Prec(), -1+1i, -1+1i)

		mul := bignum.NewComplexMultiplier()

		for i := range values {
			mul.Mul(values[i], constant, values[i])
		}

		require.NoError(t, tc.Evl.Mul(ciphertext, constant, ciphertext))

		VerifyTestVectors(tc.Params, tc.Ecd, tc.Dec, values, ciphertext, tc.Params.LogDefaultScale(), 0, *printPrecisionStats, t)
	})

	t.Run(name("Evaluator/Mul/Ct/Vector", tc), func(t *testing.T) {

		values1, _, ciphertext := tc.NewTestVector(-1-1i, 1+1i)
		values2, _, _ := tc.NewTestVector(-1-1i, 1+1i)

		mul := bignum.NewComplexMultiplier()

		for i := range values1 {
			mul.Mul(values1[i], values2[i], values1[i])
		}

		tc.Evl.Mul(ciphertext, values2, ciphertext)

		VerifyTestVectors(tc.Params, tc.Ecd, tc.Dec, values1, ciphertext, tc.Params.LogDefaultScale(), 0, *printPrecisionStats, t)
	})

	t.Run(name("Evaluator/Mul/Ct/Pt", tc), func(t *testing.T) {

		values1, plaintext1, ciphertext1 := tc.NewTestVector(-1-1i, 1+1i)

		mul := bignum.NewComplexMultiplier()

		for i := range values1 {
			mul.Mul(values1[i], values1[i], values1[i])
		}

		require.NoError(t, tc.Evl.MulRelin(ciphertext1, plaintext1, ciphertext1))

		VerifyTestVectors(tc.Params, tc.Ecd, tc.Dec, values1, ciphertext1, tc.Params.LogDefaultScale(), 0, *printPrecisionStats, t)
	})

	t.Run(name("Evaluator/Mul/Ct/Ct/Degree0", tc), func(t *testing.T) {

		values1, plaintext1, _ := tc.NewTestVector(-1-1i, 1+1i)
		values2, _, ciphertext2 := tc.NewTestVector(-1-1i, 1+1i)

		mul := bignum.NewComplexMultiplier()

		for i := range values1 {
			mul.Mul(values2[i], values1[i], values2[i])
		}

		ciphertext1 := &rlwe.Ciphertext{}
		ciphertext1.Value = []ring.Poly{plaintext1.Value}
		ciphertext1.MetaData = plaintext1.MetaData.CopyNew()

		require.NoError(t, tc.Evl.MulRelin(ciphertext1, ciphertext2, ciphertext1))

		VerifyTestVectors(tc.Params, tc.Ecd, tc.Dec, values2, ciphertext1, tc.Params.LogDefaultScale(), 0, *printPrecisionStats, t)
	})

	t.Run(name("Evaluator/MulRelin/Ct/Ct", tc), func(t *testing.T) {

		// op0 <- op0 * op1
		values1, _, ciphertext1 := tc.NewTestVector(-1-1i, 1+1i)
		values2, _, ciphertext2 := tc.NewTestVector(-1-1i, 1+1i)

		mul := bignum.NewComplexMultiplier()

		for i := range values1 {
			mul.Mul(values1[i], values2[i], values1[i])
		}

		require.NoError(t, tc.Evl.MulRelin(ciphertext1, ciphertext2, ciphertext1))
		require.Equal(t, ciphertext1.Degree(), 1)

		VerifyTestVectors(tc.Params, tc.Ecd, tc.Dec, values1, ciphertext1, tc.Params.LogDefaultScale(), 0, *printPrecisionStats, t)

		// op1 <- op0 * op1
		values1, _, ciphertext1 = tc.NewTestVector(-1-1i, 1+1i)
		values2, _, ciphertext2 = tc.NewTestVector(-1-1i, 1+1i)

		for i := range values1 {
			mul.Mul(values2[i], values1[i], values2[i])
		}

		require.NoError(t, tc.Evl.MulRelin(ciphertext1, ciphertext2, ciphertext2))
		require.Equal(t, ciphertext2.Degree(), 1)

		VerifyTestVectors(tc.Params, tc.Ecd, tc.Dec, values2, ciphertext2, tc.Params.LogDefaultScale(), 0, *printPrecisionStats, t)

		// op0 <- op0 * op0
		for i := range values1 {
			mul.Mul(values1[i], values1[i], values1[i])
		}

		require.NoError(t, tc.Evl.MulRelin(ciphertext1, ciphertext1, ciphertext1))
		require.Equal(t, ciphertext1.Degree(), 1)

		VerifyTestVectors(tc.Params, tc.Ecd, tc.Dec, values1, ciphertext1, tc.Params.LogDefaultScale(), 0, *printPrecisionStats, t)
	})
}

func testEvaluatorMulThenAdd(tc *TestContext, t *testing.T) {

	t.Run(name("Evaluator/MulThenAdd/Scalar", tc), func(t *testing.T) {

		values1, _, ciphertext1 := tc.NewTestVector(-1-1i, 1+1i)
		values2, _, ciphertext2 := tc.NewTestVector(-1-1i, 1+1i)

		constant := randomConst(tc.Params.RingType(), tc.Ecd.Prec(), -1+1i, -1+1i)

		mul := bignum.NewComplexMultiplier()

		tmp := new(bignum.Complex)
		tmp[0] = new(big.Float)
		tmp[1] = new(big.Float)

		for i := range values1 {
			mul.Mul(values1[i], constant, tmp)
			values2[i].Add(values2[i], tmp)
		}

		require.NoError(t, tc.Evl.MulThenAdd(ciphertext1, constant, ciphertext2))

		VerifyTestVectors(tc.Params, tc.Ecd, tc.Dec, values2, ciphertext2, tc.Params.LogDefaultScale(), 0, *printPrecisionStats, t)
	})

	t.Run(name("Evaluator/MulThenAdd/Vector", tc), func(t *testing.T) {

		values1, _, ciphertext1 := tc.NewTestVector(-1-1i, 1+1i)
		values2, _, ciphertext2 := tc.NewTestVector(-1-1i, 1+1i)

		require.NoError(t, tc.Evl.MulThenAdd(ciphertext2, values1, ciphertext1))

		mul := bignum.NewComplexMultiplier()

		tmp := new(bignum.Complex)
		tmp[0] = new(big.Float)
		tmp[1] = new(big.Float)

		for i := range values1 {
			mul.Mul(values2[i], values1[i], tmp)
			values1[i].Add(values1[i], tmp)
		}

		require.Equal(t, ciphertext1.Degree(), 1)

		VerifyTestVectors(tc.Params, tc.Ecd, tc.Dec, values1, ciphertext1, tc.Params.LogDefaultScale(), 0, *printPrecisionStats, t)
	})

	t.Run(name("Evaluator/MulThenAdd/Pt", tc), func(t *testing.T) {

		values1, plaintext1, ciphertext1 := tc.NewTestVector(-1, 1)
		values2, _, ciphertext2 := tc.NewTestVector(-1, 1)

		mul := bignum.NewComplexMultiplier()

		tmp := new(bignum.Complex)
		tmp[0] = new(big.Float)
		tmp[1] = new(big.Float)

		for i := range values1 {
			mul.Mul(values2[i], values1[i], tmp)
			values1[i].Add(values1[i], tmp)
		}

		require.NoError(t, tc.Evl.MulThenAdd(ciphertext2, plaintext1, ciphertext1))

		require.Equal(t, ciphertext1.Degree(), 1)

		VerifyTestVectors(tc.Params, tc.Ecd, tc.Dec, values1, ciphertext1, tc.Params.LogDefaultScale(), 0, *printPrecisionStats, t)
	})

	t.Run(name("Evaluator/MulRelinThenAdd/Ct", tc), func(t *testing.T) {

		// opOut = opOut + op1 * op0
		values1, _, ciphertext1 := tc.NewTestVector(-1-1i, 1+1i)
		values2, _, ciphertext2 := tc.NewTestVector(-1-1i, 1+1i)

		mul := bignum.NewComplexMultiplier()

		for i := range values1 {
			mul.Mul(values1[i], values2[i], values2[i])
		}

		ciphertext3 := NewCiphertext(tc.Params, 2, ciphertext1.Level())

		ciphertext3.Scale = ciphertext1.Scale.Mul(ciphertext2.Scale)

		require.NoError(t, tc.Evl.MulThenAdd(ciphertext1, ciphertext2, ciphertext3))

		require.Equal(t, ciphertext3.Degree(), 2)

		require.NoError(t, tc.Evl.Relinearize(ciphertext3, ciphertext3))

		require.Equal(t, ciphertext3.Degree(), 1)

		VerifyTestVectors(tc.Params, tc.Ecd, tc.Dec, values2, ciphertext3, tc.Params.LogDefaultScale(), 0, *printPrecisionStats, t)

		// op1 = op1 + op0*op0
		values1, _, ciphertext1 = tc.NewTestVector(-1-1i, 1+1i)
		values2, _, ciphertext2 = tc.NewTestVector(-1-1i, 1+1i)

		tmp := bignum.NewComplex()
		for i := range values1 {
			mul.Mul(values2[i], values2[i], tmp)
			values1[i].Add(values1[i], tmp)
		}

		require.NoError(t, tc.Evl.MulRelinThenAdd(ciphertext2, ciphertext2, ciphertext1))

		require.Equal(t, ciphertext1.Degree(), 1)

		VerifyTestVectors(tc.Params, tc.Ecd, tc.Dec, values1, ciphertext1, tc.Params.LogDefaultScale(), 0, *printPrecisionStats, t)
	})
}

func testBridge(tc *TestContext, t *testing.T) {

	t.Run(name("Bridge", tc), func(t *testing.T) {

		if tc.Params.RingType() != ring.ConjugateInvariant {
			t.Skip("only tested for params.RingType() == ring.ConjugateInvariant")
		}

		ciParams := tc.Params
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

		evkCtR, evkRtC := stdKeyGen.GenEvaluationKeysForRingSwapNew(stdSK, tc.Sk)

		switcher, err := NewDomainSwitcher(stdParams, evkCtR, evkRtC)
		if err != nil {
			t.Fatal(err)
		}

		evalStandar := NewEvaluator(stdParams, nil)

		values, _, ctCI := tc.NewTestVector(-1-1i, 1+1i)

		stdCTHave := NewCiphertext(stdParams, ctCI.Degree(), ctCI.Level())

		switcher.RealToComplex(evalStandar, ctCI, stdCTHave)

		VerifyTestVectors(stdParams, stdEncoder, stdDecryptor, values, stdCTHave, tc.Params.LogDefaultScale(), 0, *printPrecisionStats, t)

		stdCTImag, err := stdEvaluator.MulNew(stdCTHave, 1i)
		require.NoError(t, err)
		require.NoError(t, stdEvaluator.Add(stdCTHave, stdCTImag, stdCTHave))

		ciCTHave := NewCiphertext(ciParams, 1, stdCTHave.Level())
		switcher.ComplexToReal(evalStandar, stdCTHave, ciCTHave)

		VerifyTestVectors(tc.Params, tc.Ecd, tc.Dec, values, ciCTHave, tc.Params.LogDefaultScale(), 0, *printPrecisionStats, t)
	})
}

func name(opname string, tc *TestContext) string {

	var precMode string
	switch tc.Params.PrecisionMode() {
	case PREC64:
		precMode = "PREC64"
	case PREC128:
		precMode = "PREC128"
	}

	return fmt.Sprintf("%s/RingType=%s/logN=%d/logQP=%d/Qi=%d/Pi=%d/LogScale=%d/PrecMode=%s",
		opname,
		tc.Params.RingType(),
		tc.Params.LogN(),
		int(math.Round(tc.Params.LogQP())),
		tc.Params.QCount(),
		tc.Params.PCount(),
		int(math.Log2(tc.Params.DefaultScale().Float64())),
		precMode)
}
