package heint_test

import (
	"encoding/json"
	"fmt"
	"runtime"
	"testing"

	"github.com/tuneinsight/lattigo/v5/core/rlwe"
	"github.com/tuneinsight/lattigo/v5/he/heint"
)

func GetBenchName(params heint.Parameters, opname string) string {
	return fmt.Sprintf("%s/logN=%d/Qi=%d/Pi=%d/LogSlots=%d",
		opname,
		params.LogN(),
		params.QCount(),
		params.PCount(),
		params.LogMaxSlots())
}

func BenchmarkHEInt(b *testing.B) {

	var err error

	var testParams []heint.ParametersLiteral
	switch {
	case *flagParamString != "": // the custom test suite reads the parameters from the -params flag
		testParams = append(testParams, heint.ParametersLiteral{})
		if err = json.Unmarshal([]byte(*flagParamString), &testParams[0]); err != nil {
			b.Fatal(err)
		}
	default:
		testParams = []heint.ParametersLiteral{
			{
				LogN:             14,
				LogQ:             []int{50, 40, 40, 40, 40, 40, 40, 40},
				LogP:             []int{60},
				PlaintextModulus: 0x10001,
			},
		}
	}

	for _, paramsLiteral := range testParams {

		var params heint.Parameters
		if params, err = heint.NewParametersFromLiteral(paramsLiteral); err != nil {
			b.Error(err)
			b.Fail()
		}

		var tc *testContext
		if tc, err = genTestParams(params); err != nil {
			b.Fatal(err)
		}

		for _, testSet := range []func(tc *testContext, b *testing.B){
			benchKeyGenerator,
			benchEncoder,
			benchEncryptor,
			benchEvaluator,
		} {
			testSet(tc, b)
			runtime.GC()
		}
	}
}

func benchKeyGenerator(tc *testContext, b *testing.B) {

	params := tc.params

	b.Run(GetBenchName(params, "KeyGenerator/GenSecretKey"), func(b *testing.B) {
		sk := rlwe.NewSecretKey(params)
		kgen := tc.kgen
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			kgen.GenSecretKey(sk)
		}
	})

	b.Run(GetBenchName(params, "KeyGenerator/GenPublicKey"), func(b *testing.B) {
		sk := tc.sk
		pk := rlwe.NewPublicKey(params)
		kgen := tc.kgen
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			kgen.GenPublicKey(sk, pk)
		}
	})

	b.Run(GetBenchName(params, "KeyGenerator/GenEvaluationKey"), func(b *testing.B) {
		sk := tc.sk
		kgen := tc.kgen
		evk := rlwe.NewEvaluationKey(params)
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			kgen.GenEvaluationKey(sk, sk, evk)
		}
	})
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
		plaintext := heint.NewPlaintext(params, params.MaxLevel())
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			if err := encoder.Encode(coeffsUint64, plaintext); err != nil {
				b.Log(err)
				b.Fail()
			}
		}
	})

	b.Run(GetBenchName(params, "Encoder/Encode/Int"), func(b *testing.B) {
		plaintext := heint.NewPlaintext(params, params.MaxLevel())
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			if err := encoder.Encode(coeffsInt64, plaintext); err != nil {
				b.Log(err)
				b.Fail()
			}
		}
	})

	b.Run(GetBenchName(params, "Encoder/Decode/Uint"), func(b *testing.B) {
		plaintext := heint.NewPlaintext(params, params.MaxLevel())
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			if err := encoder.Decode(plaintext, coeffsUint64); err != nil {
				b.Log(err)
				b.Fail()
			}
		}
	})

	b.Run(GetBenchName(params, "Encoder/Decode/Int"), func(b *testing.B) {
		plaintext := heint.NewPlaintext(params, params.MaxLevel())
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			if err := encoder.Decode(plaintext, coeffsInt64); err != nil {
				b.Log(err)
				b.Fail()
			}
		}
	})
}

func benchEncryptor(tc *testContext, b *testing.B) {

	params := tc.params

	b.Run(GetBenchName(params, "Encryptor/Encrypt/Sk"), func(b *testing.B) {

		pt := heint.NewPlaintext(params, params.MaxLevel())

		poly := tc.uSampler.ReadNew()
		params.RingT().Reduce(poly, poly)

		if err := tc.encoder.Encode(poly.Coeffs[0], pt); err != nil {
			b.Log(err)
			b.Fail()
		}

		ct := heint.NewCiphertext(params, 1, pt.Level())

		enc := tc.encryptorSk

		b.ResetTimer()

		for i := 0; i < b.N; i++ {
			if err := enc.Encrypt(pt, ct); err != nil {
				b.Log(err)
				b.Fail()
			}
		}
	})

	b.Run(GetBenchName(params, "Encryptor/Encrypt/Pk"), func(b *testing.B) {

		pt := heint.NewPlaintext(params, params.MaxLevel())

		poly := tc.uSampler.ReadNew()
		params.RingT().Reduce(poly, poly)

		if err := tc.encoder.Encode(poly.Coeffs[0], pt); err != nil {
			b.Log(err)
			b.Fail()
		}

		ct := heint.NewCiphertext(params, 1, pt.Level())

		enc := tc.encryptorPk

		b.ResetTimer()

		for i := 0; i < b.N; i++ {
			if err := enc.Encrypt(pt, ct); err != nil {
				b.Log(err)
				b.Fail()
			}
		}
	})

	b.Run(GetBenchName(params, "Decryptor/Decrypt"), func(b *testing.B) {

		pt := heint.NewPlaintext(params, params.MaxLevel())

		ct := rlwe.NewCiphertextRandom(tc.prng, params.Parameters, 1, params.MaxLevel())

		*ct.MetaData = *pt.MetaData

		dec := tc.decryptor

		b.ResetTimer()

		for i := 0; i < b.N; i++ {
			dec.Decrypt(ct, pt)
		}
	})
}

func benchEvaluator(tc *testContext, b *testing.B) {

	params := tc.params
	eval := tc.evaluator

	plaintext := heint.NewPlaintext(params, params.MaxLevel())
	plaintext.Value = rlwe.NewCiphertextRandom(tc.prng, params.Parameters, 0, plaintext.Level()).Value[0]

	ciphertext1 := rlwe.NewCiphertextRandom(tc.prng, params.Parameters, 1, params.MaxLevel())
	ciphertext2 := rlwe.NewCiphertextRandom(tc.prng, params.Parameters, 1, params.MaxLevel())
	scalar := params.PlaintextModulus() >> 1

	*ciphertext1.MetaData = *plaintext.MetaData
	*ciphertext2.MetaData = *plaintext.MetaData

	vector := plaintext.Value.Coeffs[0][:params.MaxSlots()]

	b.Run(GetBenchName(params, "Evaluator/Add/Scalar"), func(b *testing.B) {
		receiver := heint.NewCiphertext(params, 1, ciphertext1.Level())
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			if err := eval.Add(ciphertext1, scalar, receiver); err != nil {
				b.Log(err)
				b.Fail()
			}
		}
	})

	b.Run(GetBenchName(params, "Evaluator/Add/Vector"), func(b *testing.B) {
		receiver := heint.NewCiphertext(params, 1, ciphertext1.Level())
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			if err := eval.Add(ciphertext1, vector, receiver); err != nil {
				b.Log(err)
				b.Fail()
			}
		}
	})

	b.Run(GetBenchName(params, "Evaluator/Add/Plaintext"), func(b *testing.B) {
		receiver := heint.NewCiphertext(params, 1, ciphertext1.Level())
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			if err := eval.Add(ciphertext1, plaintext, receiver); err != nil {
				b.Log(err)
				b.Fail()
			}
		}
	})

	b.Run(GetBenchName(params, "Evaluator/Add/Ciphertext"), func(b *testing.B) {
		receiver := heint.NewCiphertext(params, 1, ciphertext1.Level())
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			if err := eval.Add(ciphertext1, ciphertext2, receiver); err != nil {
				b.Log(err)
				b.Fail()
			}
		}
	})

	b.Run(GetBenchName(params, "Evaluator/Mul/Scalar"), func(b *testing.B) {
		receiver := heint.NewCiphertext(params, 1, ciphertext1.Level())
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			if err := eval.Mul(ciphertext1, scalar, receiver); err != nil {
				b.Log(err)
				b.Fail()
			}
		}
	})

	b.Run(GetBenchName(params, "Evaluator/Mul/Plaintext"), func(b *testing.B) {
		receiver := heint.NewCiphertext(params, 1, ciphertext1.Level())
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			if err := eval.Mul(ciphertext1, plaintext, receiver); err != nil {
				b.Log(err)
				b.Fail()
			}
		}
	})

	b.Run(GetBenchName(params, "Evaluator/Mul/Vector"), func(b *testing.B) {
		receiver := heint.NewCiphertext(params, 1, ciphertext1.Level())
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			if err := eval.Mul(ciphertext1, vector, receiver); err != nil {
				b.Log(err)
				b.Fail()
			}
		}
	})

	b.Run(GetBenchName(params, "Evaluator/Mul/Ciphertext"), func(b *testing.B) {
		receiver := heint.NewCiphertext(params, 2, ciphertext1.Level())
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			if err := eval.Mul(ciphertext1, ciphertext2, receiver); err != nil {
				b.Log(err)
				b.Fail()
			}
		}
	})

	b.Run(GetBenchName(params, "Evaluator/MulRelin/Ciphertext"), func(b *testing.B) {
		receiver := heint.NewCiphertext(params, 1, ciphertext1.Level())
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			if err := eval.MulRelin(ciphertext1, ciphertext2, receiver); err != nil {
				b.Log(err)
				b.Fail()
			}
		}
	})

	b.Run(GetBenchName(params, "Evaluator/MulInvariant/Ciphertext"), func(b *testing.B) {
		receiver := heint.NewCiphertext(params, 2, ciphertext1.Level())
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			if err := eval.MulScaleInvariant(ciphertext1, ciphertext2, receiver); err != nil {
				b.Log(err)
				b.Fail()
			}
		}
	})

	b.Run(GetBenchName(params, "Evaluator/MulRelinInvariant/Ciphertext"), func(b *testing.B) {
		receiver := heint.NewCiphertext(params, 1, ciphertext1.Level())
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			if err := eval.MulRelinScaleInvariant(ciphertext1, ciphertext2, receiver); err != nil {
				b.Log(err)
				b.Fail()
			}
		}
	})

	b.Run(GetBenchName(params, "Evaluator/MulThenAdd/Scalar"), func(b *testing.B) {
		receiver := heint.NewCiphertext(params, 1, ciphertext1.Level())
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			if err := eval.MulThenAdd(ciphertext1, scalar, receiver); err != nil {
				b.Log(err)
				b.Fail()
			}
		}
	})

	b.Run(GetBenchName(params, "Evaluator/MulThenAdd/Vector"), func(b *testing.B) {
		receiver := heint.NewCiphertext(params, 1, ciphertext1.Level())
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			if err := eval.MulThenAdd(ciphertext1, vector, receiver); err != nil {
				b.Log(err)
				b.Fail()
			}
		}
	})

	b.Run(GetBenchName(params, "Evaluator/MulThenAdd/Plaintext"), func(b *testing.B) {
		receiver := heint.NewCiphertext(params, 1, ciphertext1.Level())
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			if err := eval.MulThenAdd(ciphertext1, plaintext, receiver); err != nil {
				b.Log(err)
				b.Fail()
			}
		}
	})

	b.Run(GetBenchName(params, "Evaluator/MulThenAdd/Ciphertext"), func(b *testing.B) {
		receiver := heint.NewCiphertext(params, 1, ciphertext1.Level())
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			if err := eval.MulThenAdd(ciphertext1, plaintext, receiver); err != nil {
				b.Log(err)
				b.Fail()
			}
		}
	})

	b.Run(GetBenchName(params, "Evaluator/MulRelinThenAdd/Ciphertext"), func(b *testing.B) {
		receiver := heint.NewCiphertext(params, 2, ciphertext1.Level())
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			if err := eval.MulRelinThenAdd(ciphertext1, ciphertext2, receiver); err != nil {
				b.Log(err)
				b.Fail()
			}
		}
	})

	b.Run(GetBenchName(params, "Evaluator/Rescale"), func(b *testing.B) {
		receiver := heint.NewCiphertext(params, 1, ciphertext1.Level()-1)
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
		receiver := heint.NewCiphertext(params, 1, ciphertext2.Level())
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			if err := eval.RotateColumns(ciphertext2, 1, receiver); err != nil {
				b.Log(err)
				b.Fail()
			}
		}
	})
}
