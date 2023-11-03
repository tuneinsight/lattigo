package bgv

import (
	"encoding/json"
	"runtime"
	"testing"

	"github.com/tuneinsight/lattigo/v4/rlwe"
)

func BenchmarkBGV(b *testing.B) {

	var err error

	paramsLiterals := testParams

	if *flagParamString != "" {
		var jsonParams ParametersLiteral
		if err = json.Unmarshal([]byte(*flagParamString), &jsonParams); err != nil {
			b.Fatal(err)
		}
		paramsLiterals = []ParametersLiteral{jsonParams} // the custom test suite reads the parameters from the -params flag
	}

	for _, p := range paramsLiterals[:] {

		p.PlaintextModulus = testPlaintextModulus[1]

		var params Parameters
		if params, err = NewParametersFromLiteral(p); err != nil {
			b.Error(err)
			b.Fail()
		}

		var tc *testContext
		if tc, err = genTestParams(params); err != nil {
			b.Error(err)
			b.Fail()
		}

		for _, testSet := range []func(tc *testContext, b *testing.B){
			benchEncoder,
			benchEvaluator,
		} {
			testSet(tc, b)
			runtime.GC()
		}
	}
}

func benchEncoder(tc *testContext, b *testing.B) {

	params := tc.params

	poly := tc.uSampler.ReadNew()
	params.RingT().Reduce(poly, poly)
	coeffsUint64 := poly.Coeffs[0]
	coeffsInt64 := make([]int64, len(coeffsUint64))
	for i := range coeffsUint64 {
		coeffsInt64[i] = int64(coeffsUint64[i])
	}

	encoder := tc.encoder

	level := params.MaxLevel()
	plaintext := NewPlaintext(params, level)

	b.Run(GetTestName("Encoder/Encode/Uint", params, level), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			encoder.Encode(coeffsUint64, plaintext)
		}
	})

	b.Run(GetTestName("Encoder/Encode/Int", params, level), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			encoder.Encode(coeffsInt64, plaintext)
		}
	})

	b.Run(GetTestName("Encoder/Decode/Uint", params, level), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			encoder.Decode(plaintext, coeffsUint64)
		}
	})

	b.Run(GetTestName("Encoder/Decode/Int", params, level), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			encoder.Decode(plaintext, coeffsInt64)
		}
	})
}

func benchEvaluator(tc *testContext, b *testing.B) {

	params := tc.params
	eval := tc.evaluator
	scale := rlwe.NewScale(1)
	level := params.MaxLevel()

	ciphertext0 := rlwe.NewCiphertextRandom(tc.prng, params.Parameters, 1, level)
	ciphertext1 := rlwe.NewCiphertextRandom(tc.prng, params.Parameters, 1, level)
	ct := rlwe.NewCiphertextRandom(tc.prng, params.Parameters, 0, level)
	plaintext1 := &rlwe.Plaintext{Value: ct.Value[0]}
	plaintext1.Element.Value = ct.Value[:1]
	plaintext1.Scale = scale
	plaintext1.IsNTT = ciphertext0.IsNTT
	scalar := params.PlaintextModulus() >> 1

	b.Run(GetTestName("Evaluator/Add/Ct/Ct", params, level), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			eval.Add(ciphertext0, ciphertext1, ciphertext0)
		}
	})

	b.Run(GetTestName("Evaluator/Add/Ct/Pt", params, level), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			eval.Add(ciphertext0, plaintext1, ciphertext0)
		}
	})

	b.Run(GetTestName("Evaluator/Add/Ct/Scalar", params, level), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			eval.Add(ciphertext0, scalar, ciphertext0)
		}
	})

	b.Run(GetTestName("Evaluator/Mul/Ct/Ct", params, level), func(b *testing.B) {
		receiver := NewCiphertext(params, 2, level)
		for i := 0; i < b.N; i++ {
			eval.Mul(ciphertext0, ciphertext1, receiver)
		}
	})

	b.Run(GetTestName("Evaluator/MulInvariant/Ct/Ct", params, level), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			eval.MulScaleInvariant(ciphertext0, plaintext1.Value.Coeffs[0], ciphertext0)
		}
	})

	b.Run(GetTestName("Evaluator/Mul/Ct/Pt", params, level), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			eval.Mul(ciphertext0, plaintext1, ciphertext0)
		}
	})

	b.Run(GetTestName("Evaluator/Mul/Ct/Scalar", params, level), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			eval.Mul(ciphertext0, scalar, ciphertext0)
		}
	})

	b.Run(GetTestName("Evaluator/Mul/Ct/Vector", params, level), func(b *testing.B) {
		coeffs := plaintext1.Value.Coeffs[0][:params.MaxSlots()]
		for i := 0; i < b.N; i++ {
			eval.Mul(ciphertext0, coeffs, ciphertext0)
		}
	})

	b.Run(GetTestName("Evaluator/MulRelin/Ct/Ct", params, level), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			eval.MulRelin(ciphertext0, ciphertext1, ciphertext0)
		}
	})

	b.Run(GetTestName("Evaluator/MulRelinInvariant/Ct/Ct", params, level), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			eval.MulRelinScaleInvariant(ciphertext0, ciphertext1, ciphertext0)
		}
	})

	b.Run(GetTestName("Evaluator/MulRelinThenAdd/Ct/Ct", params, level), func(b *testing.B) {
		ciphertext2 := rlwe.NewCiphertextRandom(tc.prng, params.Parameters, 1, level)
		for i := 0; i < b.N; i++ {
			eval.MulRelinThenAdd(ciphertext0, ciphertext1, ciphertext2)
		}
	})

	b.Run(GetTestName("Evaluator/MulThenAdd/Ct/Pt", params, level), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			eval.MulThenAdd(ciphertext0, plaintext1, ciphertext1)
		}
	})

	b.Run(GetTestName("Evaluator/MulThenAdd/Ct/Scalar", params, level), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			eval.MulThenAdd(ciphertext0, scalar, ciphertext1)
		}
	})

	b.Run(GetTestName("Evaluator/MulThenAdd/Ct/Vector", params, level), func(b *testing.B) {
		coeffs := plaintext1.Value.Coeffs[0][:params.MaxSlots()]
		for i := 0; i < b.N; i++ {
			eval.MulThenAdd(ciphertext0, coeffs, ciphertext1)
		}
	})

	b.Run(GetTestName("Evaluator/Rescale", params, level), func(b *testing.B) {
		receiver := NewCiphertext(params, 1, level-1)
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			if err := eval.Rescale(ciphertext0, receiver); err != nil {
				b.Log(err)
				b.Fail()
			}
		}
	})
}
