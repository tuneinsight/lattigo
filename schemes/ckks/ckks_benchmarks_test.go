package ckks

import (
	"encoding/json"
	"runtime"
	"testing"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/utils/sampling"
)

func BenchmarkCKKS(b *testing.B) {
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
				LogN:            14,
				LogQ:            []int{50, 40, 40, 40, 40, 40, 40, 40},
				LogP:            []int{60},
				LogDefaultScale: 40,
				RingType:        ring.Standard,
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

	encoder := tc.Ecd

	b.Run(name("Encoder/Encode", tc), func(b *testing.B) {

		pt := NewPlaintext(tc.Params, tc.Params.MaxLevel())

		values := make([]complex128, 1<<pt.LogDimensions.Cols)
		for i := range values {
			values[i] = sampling.RandComplex128(-1, 1)
		}

		b.ResetTimer()

		for i := 0; i < b.N; i++ {
			if err := encoder.Encode(values, pt); err != nil {
				b.Log(err)
				b.Fail()
			}
		}
	})

	b.Run(name("Encoder/Decode", tc), func(b *testing.B) {

		pt := NewPlaintext(tc.Params, tc.Params.MaxLevel())

		values := make([]complex128, 1<<pt.LogDimensions.Cols)
		for i := range values {
			values[i] = sampling.RandComplex128(-1, 1)
		}

		encoder.Encode(values, pt)

		b.ResetTimer()

		for i := 0; i < b.N; i++ {
			if err := encoder.Decode(pt, values); err != nil {
				b.Log(err)
				b.Fail()
			}
		}
	})
}

func benchEvaluator(tc *TestContext, b *testing.B) {

	params := tc.Params
	plaintext := NewPlaintext(params, params.MaxLevel())
	plaintext.Value = rlwe.NewCiphertextRandom(tc.Prng, params.Parameters, 0, plaintext.Level()).Value[0]

	vector := make([]float64, params.MaxSlots())
	for i := range vector {
		vector[i] = 1
	}

	ciphertext1 := rlwe.NewCiphertextRandom(tc.Prng, params.Parameters, 1, params.MaxLevel())
	ciphertext2 := rlwe.NewCiphertextRandom(tc.Prng, params.Parameters, 1, params.MaxLevel())

	*ciphertext1.MetaData = *plaintext.MetaData
	*ciphertext2.MetaData = *plaintext.MetaData

	eval := tc.Evl.WithKey(rlwe.NewMemEvaluationKeySet(tc.Kgen.GenRelinearizationKeyNew(tc.Sk)))

	b.Run(name("Evaluator/Add/Scalar", tc), func(b *testing.B) {
		receiver := NewCiphertext(params, 1, ciphertext1.Level())
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			if err := eval.Add(ciphertext1, 3.1415-1.4142i, receiver); err != nil {
				b.Log(err)
				b.Fail()
			}
		}
	})

	b.Run(name("Evaluator/Add/Vector", tc), func(b *testing.B) {
		receiver := NewCiphertext(params, 1, ciphertext1.Level())
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			if err := eval.Add(ciphertext1, vector, receiver); err != nil {
				b.Log(err)
				b.Fail()
			}
		}
	})

	b.Run(name("Evaluator/Add/Plaintext", tc), func(b *testing.B) {
		receiver := NewCiphertext(params, 1, ciphertext1.Level())
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			if err := eval.Add(ciphertext1, plaintext, receiver); err != nil {
				b.Log(err)
				b.Fail()
			}
		}
	})

	b.Run(name("Evaluator/Add/Ciphertext", tc), func(b *testing.B) {
		receiver := NewCiphertext(params, 1, ciphertext1.Level())
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			if err := eval.Add(ciphertext1, ciphertext2, receiver); err != nil {
				b.Log(err)
				b.Fail()
			}
		}
	})

	b.Run(name("Evaluator/Mul/Scalar", tc), func(b *testing.B) {
		receiver := NewCiphertext(params, 1, ciphertext1.Level())
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			if err := eval.Mul(ciphertext1, 3.1415-1.4142i, receiver); err != nil {
				b.Log(err)
				b.Fail()
			}
		}
	})

	b.Run(name("Evaluator/Mul/Vector", tc), func(b *testing.B) {
		receiver := NewCiphertext(params, 1, ciphertext1.Level())
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			if err := eval.Mul(ciphertext1, vector, receiver); err != nil {
				b.Log(err)
				b.Fail()
			}
		}
	})

	b.Run(name("Evaluator/Mul/Plaintext", tc), func(b *testing.B) {
		receiver := NewCiphertext(params, 1, ciphertext1.Level())
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			if err := eval.Mul(ciphertext1, plaintext, receiver); err != nil {
				b.Log(err)
				b.Fail()
			}
		}
	})

	b.Run(name("Evaluator/Mul/Ciphertext", tc), func(b *testing.B) {
		receiver := NewCiphertext(params, 2, ciphertext1.Level())
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			if err := eval.Mul(ciphertext1, ciphertext2, receiver); err != nil {
				b.Log(err)
				b.Fail()
			}
		}
	})

	b.Run(name("Evaluator/MulRelin/Ciphertext", tc), func(b *testing.B) {
		receiver := NewCiphertext(params, 1, ciphertext1.Level())
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			if err := eval.MulRelin(ciphertext1, ciphertext2, receiver); err != nil {
				b.Log(err)
				b.Fail()
			}
		}
	})

	b.Run(name("Evaluator/MulThenAdd/Scalar", tc), func(b *testing.B) {
		receiver := NewCiphertext(params, 1, ciphertext1.Level())
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			if err := eval.MulThenAdd(ciphertext1, 3.1415-1.4142i, receiver); err != nil {
				b.Log(err)
				b.Fail()
			}
		}
	})

	b.Run(name("Evaluator/MulThenAdd/Vector", tc), func(b *testing.B) {
		receiver := NewCiphertext(params, 1, ciphertext1.Level())
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			if err := eval.MulThenAdd(ciphertext1, vector, receiver); err != nil {
				b.Log(err)
				b.Fail()
			}
		}
	})

	b.Run(name("Evaluator/MulThenAdd/Plaintext", tc), func(b *testing.B) {
		receiver := NewCiphertext(params, 1, ciphertext1.Level())
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			if err := eval.MulThenAdd(ciphertext1, plaintext, receiver); err != nil {
				b.Log(err)
				b.Fail()
			}
		}
	})

	b.Run(name("Evaluator/MulThenAdd/Ciphertext", tc), func(b *testing.B) {
		receiver := NewCiphertext(params, 2, ciphertext1.Level())
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			if err := eval.MulThenAdd(ciphertext1, ciphertext2, receiver); err != nil {
				b.Log(err)
				b.Fail()
			}
		}
	})

	b.Run(name("Evaluator/MulRelinThenAdd/Ciphertext", tc), func(b *testing.B) {
		receiver := NewCiphertext(params, 1, ciphertext1.Level())
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			if err := eval.MulRelinThenAdd(ciphertext1, ciphertext2, receiver); err != nil {
				b.Log(err)
				b.Fail()
			}
		}
	})

	b.Run(name("Evaluator/Rescale", tc), func(b *testing.B) {
		receiver := NewCiphertext(params, 1, ciphertext1.Level()-1)
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			if err := eval.Rescale(ciphertext1, receiver); err != nil {
				b.Log(err)
				b.Fail()
			}
		}
	})

	b.Run(name("Evaluator/Rotate", tc), func(b *testing.B) {
		gk := tc.Kgen.GenGaloisKeyNew(5, tc.Sk)
		evk := rlwe.NewMemEvaluationKeySet(nil, gk)
		eval := eval.WithKey(evk)
		receiver := NewCiphertext(params, 1, ciphertext1.Level())
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			if err := eval.Rotate(ciphertext1, 1, receiver); err != nil {
				b.Log(err)
				b.Fail()
			}
		}
	})
}
