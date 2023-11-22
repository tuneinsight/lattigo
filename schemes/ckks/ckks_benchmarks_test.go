package ckks

import (
	"encoding/json"
	"fmt"
	"runtime"
	"testing"

	"github.com/tuneinsight/lattigo/v5/core/rlwe"
	"github.com/tuneinsight/lattigo/v5/ring"
	"github.com/tuneinsight/lattigo/v5/utils/sampling"
)

func GetBenchName(params Parameters, opname string) string {

	var PrecisionMod string
	switch params.precisionMode {
	case PREC64:
		PrecisionMod = "PREC64"
	case PREC128:
		PrecisionMod = "PREC128"
	}

	return fmt.Sprintf("%s/RingType=%s/logN=%d/Qi=%d/Pi=%d/LogSlots=%d/%s",
		opname,
		params.RingType(),
		params.LogN(),
		params.QCount(),
		params.PCount(),
		params.LogMaxSlots(),
		PrecisionMod)
}

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

		var params Parameters
		if params, err = NewParametersFromLiteral(paramsLiteral); err != nil {
			b.Error(err)
			b.Fail()
		}

		var tc *testContext
		if tc, err = genTestParams(params); err != nil {
			b.Fatal(err)
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

	encoder := tc.encoder

	b.Run(GetBenchName(tc.params, "Encoder/Encode"), func(b *testing.B) {

		pt := NewPlaintext(tc.params, tc.params.MaxLevel())

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

	b.Run(GetBenchName(tc.params, "Encoder/Decode"), func(b *testing.B) {

		pt := NewPlaintext(tc.params, tc.params.MaxLevel())

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

func benchEvaluator(tc *testContext, b *testing.B) {

	params := tc.params
	plaintext := NewPlaintext(params, params.MaxLevel())
	plaintext.Value = rlwe.NewCiphertextRandom(tc.prng, params.Parameters, 0, plaintext.Level()).Value[0]

	vector := make([]float64, params.MaxSlots())
	for i := range vector {
		vector[i] = 1
	}

	ciphertext1 := rlwe.NewCiphertextRandom(tc.prng, params.Parameters, 1, params.MaxLevel())
	ciphertext2 := rlwe.NewCiphertextRandom(tc.prng, params.Parameters, 1, params.MaxLevel())

	*ciphertext1.MetaData = *plaintext.MetaData
	*ciphertext2.MetaData = *plaintext.MetaData

	eval := tc.evaluator.WithKey(rlwe.NewMemEvaluationKeySet(tc.kgen.GenRelinearizationKeyNew(tc.sk)))

	b.Run(GetBenchName(params, "Evaluator/Add/Scalar"), func(b *testing.B) {
		receiver := NewCiphertext(params, 1, ciphertext1.Level())
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			if err := eval.Add(ciphertext1, 3.1415-1.4142i, receiver); err != nil {
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
			if err := eval.Mul(ciphertext1, 3.1415-1.4142i, receiver); err != nil {
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

	b.Run(GetBenchName(params, "Evaluator/MulThenAdd/Scalar"), func(b *testing.B) {
		receiver := NewCiphertext(params, 1, ciphertext1.Level())
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			if err := eval.MulThenAdd(ciphertext1, 3.1415-1.4142i, receiver); err != nil {
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
		receiver := NewCiphertext(params, 2, ciphertext1.Level())
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			if err := eval.MulThenAdd(ciphertext1, ciphertext2, receiver); err != nil {
				b.Log(err)
				b.Fail()
			}
		}
	})

	b.Run(GetBenchName(params, "Evaluator/MulRelinThenAdd/Ciphertext"), func(b *testing.B) {
		receiver := NewCiphertext(params, 1, ciphertext1.Level())
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
