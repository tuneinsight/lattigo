package bfv

import (
	"encoding/json"
	"runtime"
	"testing"

	"github.com/tuneinsight/lattigo/v5/core/rlwe"
)

func BenchmarkBFV(b *testing.B) {
	var err error

	var testParams []ParametersLiteral
	switch {
	case *flagParamString != "": // the custom test suite reads the parameters from the -params flag
		testParams = append(testParams, ParametersLiteral{})
		if err = json.Unmarshal([]byte(*flagParamString), &testParams[0]); err != nil {
			b.Fatal(err)
		}
	default:
		testParams = []ParametersLiteral{
			{
				LogN:             14,
				LogQ:             []int{50, 40, 40, 40, 40, 40, 40, 40},
				LogP:             []int{60},
				PlaintextModulus: 0x10001,
			},
		}
	}

	for _, paramsLiteral := range testParams {
		tc := NewTestContext(paramsLiteral)

		for _, testSet := range []func(tc *TestContext, b *testing.B){
			benchEncoder,
			benchEvaluator,
		} {
			testSet(tc, b)
			runtime.GC()
		}
	}
}

func benchEncoder(tc *TestContext, b *testing.B) {

	params := tc.Params

	poly := tc.Sampler.ReadNew()
	params.RingT().Reduce(poly, poly)
	coeffsUint64 := poly.Coeffs[0]
	coeffsInt64 := make([]int64, len(coeffsUint64))
	for i := range coeffsUint64 {
		coeffsInt64[i] = int64(coeffsUint64[i])
	}

	encoder := tc.Ecd

	level := params.MaxLevel()
	plaintext := NewPlaintext(params, level)

	b.Run(name("Encoder/Encode/Uint", tc, level), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			if err := encoder.Encode(coeffsUint64, plaintext); err != nil {
				b.Log(err)
				b.Fail()
			}
		}
	})

	b.Run(name("Encoder/Encode/Int", tc, level), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			if err := encoder.Encode(coeffsInt64, plaintext); err != nil {
				b.Log(err)
				b.Fail()
			}
		}
	})

	b.Run(name("Encoder/Decode/Uint", tc, level), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			if err := encoder.Decode(plaintext, coeffsUint64); err != nil {
				b.Log(err)
				b.Fail()
			}
		}
	})

	b.Run(name("Encoder/Decode/Int", tc, level), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			if err := encoder.Decode(plaintext, coeffsInt64); err != nil {
				b.Log(err)
				b.Fail()
			}
		}
	})
}

func benchEvaluator(tc *TestContext, b *testing.B) {

	params := tc.Params
	eval := tc.Evl

	level := params.MaxLevel()

	plaintext := NewPlaintext(params, level)
	plaintext.Value = rlwe.NewCiphertextRandom(tc.Prng, params.Parameters, 0, plaintext.Level()).Value[0]

	ciphertext1 := rlwe.NewCiphertextRandom(tc.Prng, params.Parameters, 1, level)
	ciphertext2 := rlwe.NewCiphertextRandom(tc.Prng, params.Parameters, 1, level)
	scalar := params.PlaintextModulus() >> 1

	*ciphertext1.MetaData = *plaintext.MetaData
	*ciphertext2.MetaData = *plaintext.MetaData

	vector := plaintext.Value.Coeffs[0][:params.MaxSlots()]

	b.Run(name("Evaluator/Add/Scalar", tc, level), func(b *testing.B) {
		receiver := NewCiphertext(params, 1, ciphertext1.Level())
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			if err := eval.Add(ciphertext1, scalar, receiver); err != nil {
				b.Log(err)
				b.Fail()
			}
		}
	})

	b.Run(name("Evaluator/Add/Vector", tc, level), func(b *testing.B) {
		receiver := NewCiphertext(params, 1, ciphertext1.Level())
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			if err := eval.Add(ciphertext1, vector, receiver); err != nil {
				b.Log(err)
				b.Fail()
			}
		}
	})

	b.Run(name("Evaluator/Add/Plaintext", tc, level), func(b *testing.B) {
		receiver := NewCiphertext(params, 1, ciphertext1.Level())
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			if err := eval.Add(ciphertext1, plaintext, receiver); err != nil {
				b.Log(err)
				b.Fail()
			}
		}
	})

	b.Run(name("Evaluator/Add/Ciphertext", tc, level), func(b *testing.B) {
		receiver := NewCiphertext(params, 1, ciphertext1.Level())
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			if err := eval.Add(ciphertext1, ciphertext2, receiver); err != nil {
				b.Log(err)
				b.Fail()
			}
		}
	})

	b.Run(name("Evaluator/Mul/Scalar", tc, level), func(b *testing.B) {
		receiver := NewCiphertext(params, 1, ciphertext1.Level())
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			if err := eval.Mul(ciphertext1, scalar, receiver); err != nil {
				b.Log(err)
				b.Fail()
			}
		}
	})

	b.Run(name("Evaluator/Mul/Plaintext", tc, level), func(b *testing.B) {
		receiver := NewCiphertext(params, 1, ciphertext1.Level())
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			if err := eval.Mul(ciphertext1, plaintext, receiver); err != nil {
				b.Log(err)
				b.Fail()
			}
		}
	})

	b.Run(name("Evaluator/Mul/Vector", tc, level), func(b *testing.B) {
		receiver := NewCiphertext(params, 1, ciphertext1.Level())
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			if err := eval.Mul(ciphertext1, vector, receiver); err != nil {
				b.Log(err)
				b.Fail()
			}
		}
	})

	b.Run(name("Evaluator/Mul/Ciphertext", tc, level), func(b *testing.B) {
		receiver := NewCiphertext(params, 2, ciphertext1.Level())
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			if err := eval.Mul(ciphertext1, ciphertext2, receiver); err != nil {
				b.Log(err)
				b.Fail()
			}
		}
	})

	b.Run(name("Evaluator/MulRelin/Ciphertext", tc, level), func(b *testing.B) {
		receiver := NewCiphertext(params, 1, ciphertext1.Level())
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			if err := eval.MulRelin(ciphertext1, ciphertext2, receiver); err != nil {
				b.Log(err)
				b.Fail()
			}
		}
	})

	b.Run(name("Evaluator/MulThenAdd/Scalar", tc, level), func(b *testing.B) {
		receiver := NewCiphertext(params, 1, ciphertext1.Level())
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			if err := eval.MulThenAdd(ciphertext1, scalar, receiver); err != nil {
				b.Log(err)
				b.Fail()
			}
		}
	})

	b.Run(name("Evaluator/MulThenAdd/Vector", tc, level), func(b *testing.B) {
		receiver := NewCiphertext(params, 1, ciphertext1.Level())
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			if err := eval.MulThenAdd(ciphertext1, vector, receiver); err != nil {
				b.Log(err)
				b.Fail()
			}
		}
	})

	b.Run(name("Evaluator/MulThenAdd/Plaintext", tc, level), func(b *testing.B) {
		receiver := NewCiphertext(params, 1, ciphertext1.Level())
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			if err := eval.MulThenAdd(ciphertext1, plaintext, receiver); err != nil {
				b.Log(err)
				b.Fail()
			}
		}
	})

	b.Run(name("Evaluator/MulThenAdd/Ciphertext", tc, level), func(b *testing.B) {
		receiver := NewCiphertext(params, 1, ciphertext1.Level())
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			if err := eval.MulThenAdd(ciphertext1, plaintext, receiver); err != nil {
				b.Log(err)
				b.Fail()
			}
		}
	})

	b.Run(name("Evaluator/MulRelinThenAdd/Ciphertext", tc, level), func(b *testing.B) {
		receiver := NewCiphertext(params, 2, ciphertext1.Level())
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			if err := eval.MulRelinThenAdd(ciphertext1, ciphertext2, receiver); err != nil {
				b.Log(err)
				b.Fail()
			}
		}
	})

	b.Run(name("Evaluator/Rotate", tc, level), func(b *testing.B) {
		gk := tc.Kgen.GenGaloisKeyNew(5, tc.Sk)
		evk := rlwe.NewMemEvaluationKeySet(nil, gk)
		eval := eval.WithKey(evk)
		receiver := NewCiphertext(params, 1, ciphertext2.Level())
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			if err := eval.RotateColumns(ciphertext2, 1, receiver); err != nil {
				b.Log(err)
				b.Fail()
			}
		}
	})
}
