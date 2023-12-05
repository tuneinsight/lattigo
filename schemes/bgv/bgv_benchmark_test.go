package bgv

import (
	"encoding/json"
	"fmt"
	"runtime"
	"testing"

	"github.com/tuneinsight/lattigo/v5/core/rlwe"
)

func GetBenchName(params Parameters, opname string) string {
	return fmt.Sprintf("%s/logN=%d/Qi=%d/Pi=%d/LogSlots=%d",
		opname,
		params.LogN(),
		params.QCount(),
		params.PCount(),
		params.LogMaxSlots())
}

func BenchmarkBGV(b *testing.B) {

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

		var params Parameters
		if params, err = NewParametersFromLiteral(paramsLiteral); err != nil {
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

	b.Run(GetBenchName(params, "Encoder/Encode/Uint"), func(b *testing.B) {
		plaintext := NewPlaintext(params, params.MaxLevel())
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			if err := encoder.Encode(coeffsUint64, plaintext); err != nil {
				b.Log(err)
				b.Fail()
			}
		}
	})

	b.Run(GetBenchName(params, "Encoder/Encode/Int"), func(b *testing.B) {
		plaintext := NewPlaintext(params, params.MaxLevel())
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			if err := encoder.Encode(coeffsInt64, plaintext); err != nil {
				b.Log(err)
				b.Fail()
			}
		}
	})

	b.Run(GetBenchName(params, "Encoder/Decode/Uint"), func(b *testing.B) {
		plaintext := NewPlaintext(params, params.MaxLevel())
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			if err := encoder.Decode(plaintext, coeffsUint64); err != nil {
				b.Log(err)
				b.Fail()
			}
		}
	})

	b.Run(GetBenchName(params, "Encoder/Decode/Int"), func(b *testing.B) {
		plaintext := NewPlaintext(params, params.MaxLevel())
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			if err := encoder.Decode(plaintext, coeffsInt64); err != nil {
				b.Log(err)
				b.Fail()
			}
		}
	})
}

func benchEvaluator(tc *testContext, b *testing.B) {

	params := tc.params
	eval := tc.evaluator

	plaintext := NewPlaintext(params, params.MaxLevel())
	plaintext.Value = rlwe.NewCiphertextRandom(tc.prng, params.Parameters, 0, plaintext.Level()).Value[0]

	ciphertext1 := rlwe.NewCiphertextRandom(tc.prng, params.Parameters, 1, params.MaxLevel())
	ciphertext2 := rlwe.NewCiphertextRandom(tc.prng, params.Parameters, 1, params.MaxLevel())
	scalar := params.PlaintextModulus() >> 1

	*ciphertext1.MetaData = *plaintext.MetaData
	*ciphertext2.MetaData = *plaintext.MetaData

	vector := plaintext.Value.Coeffs[0][:params.MaxSlots()]

	b.Run(GetBenchName(params, "Evaluator/Add/Scalar"), func(b *testing.B) {
		receiver := NewCiphertext(params, 1, ciphertext1.Level())
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			if err := eval.Add(ciphertext1, scalar, receiver); err != nil {
				b.Log(err)
				b.Fail()
			}
		}
	})

	b.Run(GetBenchName(params, "Evaluator/Add/Vector"), func(b *testing.B) {
		receiver := NewCiphertext(params, 1, ciphertext1.Level())
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			if err := eval.Add(ciphertext1, vector, receiver); err != nil {
				b.Log(err)
				b.Fail()
			}
		}
	})

	b.Run(GetBenchName(params, "Evaluator/Add/Plaintext"), func(b *testing.B) {
		receiver := NewCiphertext(params, 1, ciphertext1.Level())
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			if err := eval.Add(ciphertext1, plaintext, receiver); err != nil {
				b.Log(err)
				b.Fail()
			}
		}
	})

	b.Run(GetBenchName(params, "Evaluator/Add/Ciphertext"), func(b *testing.B) {
		receiver := NewCiphertext(params, 1, ciphertext1.Level())
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			if err := eval.Add(ciphertext1, ciphertext2, receiver); err != nil {
				b.Log(err)
				b.Fail()
			}
		}
	})

	b.Run(GetBenchName(params, "Evaluator/Mul/Scalar"), func(b *testing.B) {
		receiver := NewCiphertext(params, 1, ciphertext1.Level())
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			if err := eval.Mul(ciphertext1, scalar, receiver); err != nil {
				b.Log(err)
				b.Fail()
			}
		}
	})

	b.Run(GetBenchName(params, "Evaluator/Mul/Plaintext"), func(b *testing.B) {
		receiver := NewCiphertext(params, 1, ciphertext1.Level())
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			if err := eval.Mul(ciphertext1, plaintext, receiver); err != nil {
				b.Log(err)
				b.Fail()
			}
		}
	})

	b.Run(GetBenchName(params, "Evaluator/Mul/Vector"), func(b *testing.B) {
		receiver := NewCiphertext(params, 1, ciphertext1.Level())
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			if err := eval.Mul(ciphertext1, vector, receiver); err != nil {
				b.Log(err)
				b.Fail()
			}
		}
	})

	b.Run(GetBenchName(params, "Evaluator/Mul/Ciphertext"), func(b *testing.B) {
		receiver := NewCiphertext(params, 2, ciphertext1.Level())
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			if err := eval.Mul(ciphertext1, ciphertext2, receiver); err != nil {
				b.Log(err)
				b.Fail()
			}
		}
	})

	b.Run(GetBenchName(params, "Evaluator/MulRelin/Ciphertext"), func(b *testing.B) {
		receiver := NewCiphertext(params, 1, ciphertext1.Level())
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			if err := eval.MulRelin(ciphertext1, ciphertext2, receiver); err != nil {
				b.Log(err)
				b.Fail()
			}
		}
	})

	b.Run(GetBenchName(params, "Evaluator/MulInvariant/Ciphertext"), func(b *testing.B) {
		receiver := NewCiphertext(params, 2, ciphertext1.Level())
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			if err := eval.MulScaleInvariant(ciphertext1, ciphertext2, receiver); err != nil {
				b.Log(err)
				b.Fail()
			}
		}
	})

	b.Run(GetBenchName(params, "Evaluator/MulRelinInvariant/Ciphertext"), func(b *testing.B) {
		receiver := NewCiphertext(params, 1, ciphertext1.Level())
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			if err := eval.MulRelinScaleInvariant(ciphertext1, ciphertext2, receiver); err != nil {
				b.Log(err)
				b.Fail()
			}
		}
	})

	b.Run(GetBenchName(params, "Evaluator/MulThenAdd/Scalar"), func(b *testing.B) {
		receiver := NewCiphertext(params, 1, ciphertext1.Level())
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			if err := eval.MulThenAdd(ciphertext1, scalar, receiver); err != nil {
				b.Log(err)
				b.Fail()
			}
		}
	})

	b.Run(GetBenchName(params, "Evaluator/MulThenAdd/Vector"), func(b *testing.B) {
		receiver := NewCiphertext(params, 1, ciphertext1.Level())
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			if err := eval.MulThenAdd(ciphertext1, vector, receiver); err != nil {
				b.Log(err)
				b.Fail()
			}
		}
	})

	b.Run(GetBenchName(params, "Evaluator/MulThenAdd/Plaintext"), func(b *testing.B) {
		receiver := NewCiphertext(params, 1, ciphertext1.Level())
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			if err := eval.MulThenAdd(ciphertext1, plaintext, receiver); err != nil {
				b.Log(err)
				b.Fail()
			}
		}
	})

	b.Run(GetBenchName(params, "Evaluator/MulThenAdd/Ciphertext"), func(b *testing.B) {
		receiver := NewCiphertext(params, 1, ciphertext1.Level())
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			if err := eval.MulThenAdd(ciphertext1, plaintext, receiver); err != nil {
				b.Log(err)
				b.Fail()
			}
		}
	})

	b.Run(GetBenchName(params, "Evaluator/MulRelinThenAdd/Ciphertext"), func(b *testing.B) {
		receiver := NewCiphertext(params, 2, ciphertext1.Level())
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			if err := eval.MulRelinThenAdd(ciphertext1, ciphertext2, receiver); err != nil {
				b.Log(err)
				b.Fail()
			}
		}
	})

	b.Run(GetBenchName(params, "Evaluator/Rescale"), func(b *testing.B) {
		receiver := NewCiphertext(params, 1, ciphertext1.Level()-1)
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			if err := eval.Rescale(ciphertext1, receiver); err != nil {
				b.Log(err)
				b.Fail()
			}
		}
	})

	b.Run(GetBenchName(params, "Evaluator/Rotate"), func(b *testing.B) {
		gk := tc.kgen.GenGaloisKeyNew(5, tc.sk)
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
