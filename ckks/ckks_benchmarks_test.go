package ckks

import (
	"encoding/json"
	"testing"

	"github.com/tuneinsight/lattigo/v4/rlwe"
	"github.com/tuneinsight/lattigo/v4/utils/sampling"
)

func BenchmarkCKKSScheme(b *testing.B) {

	var err error

	defaultParams := append(DefaultParams, DefaultConjugateInvariantParams...)
	if testing.Short() {
		defaultParams = DefaultParams[:2]
	}

	if *flagParamString != "" {
		var jsonParams ParametersLiteral
		if err = json.Unmarshal([]byte(*flagParamString), &jsonParams); err != nil {
			b.Fatal(err)
		}
		defaultParams = []ParametersLiteral{jsonParams} // the custom test suite reads the parameters from the -params flag
	}

	for _, defaultParams := range defaultParams {
		var params Parameters
		if params, err = NewParametersFromLiteral(defaultParams); err != nil {
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

func benchEncoder(tc *testContext, b *testing.B) {

	encoder := tc.encoder
	logSlots := tc.params.LogSlots()

	b.Run(GetTestName(tc.params, "Encoder/Encode"), func(b *testing.B) {

		values := make([]complex128, 1<<logSlots)
		for i := 0; i < 1<<logSlots; i++ {
			values[i] = sampling.RandComplex128(-1, 1)
		}

		plaintext := NewPlaintext(tc.params, tc.params.MaxLevel())

		for i := 0; i < b.N; i++ {
			encoder.Encode(values, plaintext, logSlots)
		}
	})

	b.Run(GetTestName(tc.params, "Encoder/Decode"), func(b *testing.B) {

		values := make([]complex128, 1<<logSlots)
		for i := 0; i < 1<<logSlots; i++ {
			values[i] = sampling.RandComplex128(-1, 1)
		}

		plaintext := NewPlaintext(tc.params, tc.params.MaxLevel())
		encoder.Encode(values, plaintext, logSlots)

		for i := 0; i < b.N; i++ {
			encoder.Decode(plaintext, logSlots)
		}
	})
}

func benchEvaluator(tc *testContext, b *testing.B) {

	plaintext := NewPlaintext(tc.params, tc.params.MaxLevel())
	ciphertext1 := rlwe.NewCiphertextRandom(tc.prng, tc.params.Parameters, 1, tc.params.MaxLevel())
	ciphertext2 := rlwe.NewCiphertextRandom(tc.prng, tc.params.Parameters, 1, tc.params.MaxLevel())
	receiver := rlwe.NewCiphertextRandom(tc.prng, tc.params.Parameters, 2, tc.params.MaxLevel())

	eval := tc.evaluator.WithKey(&rlwe.EvaluationKeySet{RelinearizationKey: tc.kgen.GenRelinearizationKeyNew(tc.sk)})

	b.Run(GetTestName(tc.params, "Evaluator/Add"), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			eval.Add(ciphertext1, ciphertext2, ciphertext1)
		}
	})

	b.Run(GetTestName(tc.params, "Evaluator/AddScalar"), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			eval.AddConst(ciphertext1, 3.1415-1.4142i, ciphertext1)
		}
	})

	b.Run(GetTestName(tc.params, "Evaluator/MulScalar"), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			eval.MultByConst(ciphertext1, 3.1415-1.4142i, ciphertext1)
		}
	})

	b.Run(GetTestName(tc.params, "Evaluator/MulPlain"), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			eval.Mul(ciphertext1, plaintext, ciphertext1)
		}
	})

	b.Run(GetTestName(tc.params, "Evaluator/Mul"), func(b *testing.B) {
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
		ciphertext1.Scale = tc.params.DefaultScale().Mul(tc.params.DefaultScale())

		for i := 0; i < b.N; i++ {
			if err := eval.Rescale(ciphertext1, tc.params.DefaultScale(), ciphertext2); err != nil {
				panic(err)
			}
		}
	})
}
