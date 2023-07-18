package ckks

import (
	"encoding/json"
	"testing"

	"github.com/stretchr/testify/require"

	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/rlwe"
	"github.com/tuneinsight/lattigo/v4/utils/sampling"
)

func BenchmarkCKKSScheme(b *testing.B) {

	var err error

	var testParams []ParametersLiteral
	switch {
	case *flagParamString != "": // the custom test suite reads the parameters from the -params flag
		testParams = append(testParams, ParametersLiteral{})
		if err = json.Unmarshal([]byte(*flagParamString), &testParams[0]); err != nil {
			b.Fatal(err)
		}
	default:
		testParams = testParamsLiteral
	}

	for _, ringType := range []ring.Type{ring.Standard, ring.ConjugateInvariant} {

		for _, paramsLiteral := range testParams {

			paramsLiteral.RingType = ringType

			var params Parameters
			if params, err = NewParametersFromLiteral(paramsLiteral); err != nil {
				b.Fatal(err)
			}

			var tc *testContext
			if tc, err = genTestParams(params); err != nil {
				b.Fatal(err)
			}

			benchEncoder(tc, b)
			benchEvaluator(tc, b)
		}
	}
}

func benchEncoder(tc *testContext, b *testing.B) {

	encoder := tc.encoder

	b.Run(GetTestName(tc.params, "Encoder/Encode"), func(b *testing.B) {

		pt := NewPlaintext(tc.params, tc.params.MaxLevel())

		values := make([]complex128, 1<<pt.PlaintextLogDimensions.Cols)
		for i := range values {
			values[i] = sampling.RandComplex128(-1, 1)
		}

		b.ResetTimer()

		for i := 0; i < b.N; i++ {
			encoder.Encode(values, pt)
		}
	})

	b.Run(GetTestName(tc.params, "Encoder/Decode"), func(b *testing.B) {

		pt := NewPlaintext(tc.params, tc.params.MaxLevel())

		values := make([]complex128, 1<<pt.PlaintextLogDimensions.Cols)
		for i := range values {
			values[i] = sampling.RandComplex128(-1, 1)
		}

		encoder.Encode(values, pt)

		b.ResetTimer()

		for i := 0; i < b.N; i++ {
			encoder.Decode(pt, values)
		}
	})
}

func benchEvaluator(tc *testContext, b *testing.B) {

	plaintext := NewPlaintext(tc.params, tc.params.MaxLevel())
	ciphertext1 := rlwe.NewCiphertextRandom(tc.prng, tc.params.Parameters, 1, tc.params.MaxLevel())
	ciphertext2 := rlwe.NewCiphertextRandom(tc.prng, tc.params.Parameters, 1, tc.params.MaxLevel())
	receiver := rlwe.NewCiphertextRandom(tc.prng, tc.params.Parameters, 2, tc.params.MaxLevel())

	rlk, err := tc.kgen.GenRelinearizationKeyNew(tc.sk)
	require.NoError(b, err)

	eval := tc.evaluator.WithKey(rlwe.NewMemEvaluationKeySet(rlk))

	b.Run(GetTestName(tc.params, "Evaluator/Add/Scalar"), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			eval.Add(ciphertext1, 3.1415-1.4142i, ciphertext1)
		}
	})

	b.Run(GetTestName(tc.params, "Evaluator/Add/Pt"), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			eval.Add(ciphertext1, plaintext, ciphertext1)
		}
	})

	b.Run(GetTestName(tc.params, "Evaluator/Add/Ct"), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			eval.Add(ciphertext1, ciphertext2, ciphertext1)
		}
	})

	b.Run(GetTestName(tc.params, "Evaluator/Mul/Scalar"), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			eval.Mul(ciphertext1, 3.1415-1.4142i, ciphertext1)
		}
	})

	b.Run(GetTestName(tc.params, "Evaluator/Mul/Pt"), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			eval.Mul(ciphertext1, plaintext, ciphertext1)
		}
	})

	b.Run(GetTestName(tc.params, "Evaluator/Mul/Ct"), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			eval.Mul(ciphertext1, ciphertext2, receiver)
		}
	})

	b.Run(GetTestName(tc.params, "Evaluator/Square"), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			eval.Mul(ciphertext1, ciphertext1, receiver)
		}
	})

	b.Run(GetTestName(tc.params, "Evaluator/Rescale"), func(b *testing.B) {
		ciphertext1.PlaintextScale = tc.params.PlaintextScale().Mul(tc.params.PlaintextScale())

		for i := 0; i < b.N; i++ {
			eval.Rescale(ciphertext1, tc.params.PlaintextScale(), ciphertext2)
		}
	})
}
