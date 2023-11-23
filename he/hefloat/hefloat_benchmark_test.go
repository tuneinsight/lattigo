package hefloat_test

import (
	"encoding/json"
	"fmt"
	"runtime"
	"testing"

	"github.com/tuneinsight/lattigo/v5/core/rlwe"
	"github.com/tuneinsight/lattigo/v5/he/hefloat"
	"github.com/tuneinsight/lattigo/v5/ring"
	"github.com/tuneinsight/lattigo/v5/schemes/ckks"
	"github.com/tuneinsight/lattigo/v5/utils/sampling"
)

func GetBenchName(params hefloat.Parameters, opname string) string {

	var PrecisionMod string
	switch params.PrecisionMode() {
	case ckks.PREC64:
		PrecisionMod = "PREC64"
	case ckks.PREC128:
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

func BenchmarkHEFloat(b *testing.B) {

	var err error

	var testParams []hefloat.ParametersLiteral
	switch {
	case *flagParamString != "": // the custom test suite reads the parameters from the -params flag
		testParams = append(testParams, hefloat.ParametersLiteral{})
		if err = json.Unmarshal([]byte(*flagParamString), &testParams[0]); err != nil {
			b.Fatal(err)
		}
	default:
		testParams = []hefloat.ParametersLiteral{
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

		var params hefloat.Parameters
		if params, err = hefloat.NewParametersFromLiteral(paramsLiteral); err != nil {
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

	encoder := tc.encoder

	b.Run(GetBenchName(tc.params, "Encoder/Encode"), func(b *testing.B) {

		pt := hefloat.NewPlaintext(tc.params, tc.params.MaxLevel())

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

		pt := hefloat.NewPlaintext(tc.params, tc.params.MaxLevel())

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

func benchEncryptor(tc *testContext, b *testing.B) {

	params := tc.params

	b.Run(GetBenchName(params, "Encryptor/Encrypt/Sk"), func(b *testing.B) {

		pt := hefloat.NewPlaintext(params, params.MaxLevel())

		values := make([]complex128, 1<<pt.LogDimensions.Cols)
		for i := range values {
			values[i] = sampling.RandComplex128(-1, 1)
		}

		if err := tc.encoder.Encode(values, pt); err != nil {
			b.Log(err)
			b.Fail()
		}

		ct := hefloat.NewCiphertext(params, 1, pt.Level())

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

		pt := hefloat.NewPlaintext(params, params.MaxLevel())

		values := make([]complex128, 1<<pt.LogDimensions.Cols)
		for i := range values {
			values[i] = sampling.RandComplex128(-1, 1)
		}

		if err := tc.encoder.Encode(values, pt); err != nil {
			b.Log(err)
			b.Fail()
		}

		ct := hefloat.NewCiphertext(params, 1, pt.Level())

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

		pt := hefloat.NewPlaintext(params, params.MaxLevel())

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
	plaintext := hefloat.NewPlaintext(params, params.MaxLevel())
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
		receiver := hefloat.NewCiphertext(params, 1, ciphertext1.Level())
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			if err := eval.Add(ciphertext1, 3.1415-1.4142i, receiver); err != nil {
				b.Log(err)
				b.Fail()
			}
		}
	})

	b.Run(GetBenchName(params, "Evaluator/Add/Vector"), func(b *testing.B) {
		receiver := hefloat.NewCiphertext(params, 1, ciphertext1.Level())
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			if err := eval.Add(ciphertext1, vector, receiver); err != nil {
				b.Log(err)
				b.Fail()
			}
		}
	})

	b.Run(GetBenchName(params, "Evaluator/Add/Plaintext"), func(b *testing.B) {
		receiver := hefloat.NewCiphertext(params, 1, ciphertext1.Level())
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			if err := eval.Add(ciphertext1, plaintext, receiver); err != nil {
				b.Log(err)
				b.Fail()
			}
		}
	})

	b.Run(GetBenchName(params, "Evaluator/Add/Ciphertext"), func(b *testing.B) {
		receiver := hefloat.NewCiphertext(params, 1, ciphertext1.Level())
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			if err := eval.Add(ciphertext1, ciphertext2, receiver); err != nil {
				b.Log(err)
				b.Fail()
			}
		}
	})

	b.Run(GetBenchName(params, "Evaluator/Mul/Scalar"), func(b *testing.B) {
		receiver := hefloat.NewCiphertext(params, 1, ciphertext1.Level())
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			if err := eval.Mul(ciphertext1, 3.1415-1.4142i, receiver); err != nil {
				b.Log(err)
				b.Fail()
			}
		}
	})

	b.Run(GetBenchName(params, "Evaluator/Mul/Vector"), func(b *testing.B) {
		receiver := hefloat.NewCiphertext(params, 1, ciphertext1.Level())
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			if err := eval.Mul(ciphertext1, vector, receiver); err != nil {
				b.Log(err)
				b.Fail()
			}
		}
	})

	b.Run(GetBenchName(params, "Evaluator/Mul/Plaintext"), func(b *testing.B) {
		receiver := hefloat.NewCiphertext(params, 1, ciphertext1.Level())
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			if err := eval.Mul(ciphertext1, plaintext, receiver); err != nil {
				b.Log(err)
				b.Fail()
			}
		}
	})

	b.Run(GetBenchName(params, "Evaluator/Mul/Ciphertext"), func(b *testing.B) {
		receiver := hefloat.NewCiphertext(params, 2, ciphertext1.Level())
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			if err := eval.Mul(ciphertext1, ciphertext2, receiver); err != nil {
				b.Log(err)
				b.Fail()
			}
		}
	})

	b.Run(GetBenchName(params, "Evaluator/MulRelin/Ciphertext"), func(b *testing.B) {
		receiver := hefloat.NewCiphertext(params, 1, ciphertext1.Level())
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			if err := eval.MulRelin(ciphertext1, ciphertext2, receiver); err != nil {
				b.Log(err)
				b.Fail()
			}
		}
	})

	b.Run(GetBenchName(params, "Evaluator/MulThenAdd/Scalar"), func(b *testing.B) {
		receiver := hefloat.NewCiphertext(params, 1, ciphertext1.Level())
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			if err := eval.MulThenAdd(ciphertext1, 3.1415-1.4142i, receiver); err != nil {
				b.Log(err)
				b.Fail()
			}
		}
	})

	b.Run(GetBenchName(params, "Evaluator/MulThenAdd/Vector"), func(b *testing.B) {
		receiver := hefloat.NewCiphertext(params, 1, ciphertext1.Level())
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			if err := eval.MulThenAdd(ciphertext1, vector, receiver); err != nil {
				b.Log(err)
				b.Fail()
			}
		}
	})

	b.Run(GetBenchName(params, "Evaluator/MulThenAdd/Plaintext"), func(b *testing.B) {
		receiver := hefloat.NewCiphertext(params, 1, ciphertext1.Level())
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			if err := eval.MulThenAdd(ciphertext1, plaintext, receiver); err != nil {
				b.Log(err)
				b.Fail()
			}
		}
	})

	b.Run(GetBenchName(params, "Evaluator/MulThenAdd/Ciphertext"), func(b *testing.B) {
		receiver := hefloat.NewCiphertext(params, 2, ciphertext1.Level())
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			if err := eval.MulThenAdd(ciphertext1, ciphertext2, receiver); err != nil {
				b.Log(err)
				b.Fail()
			}
		}
	})

	b.Run(GetBenchName(params, "Evaluator/MulRelinThenAdd/Ciphertext"), func(b *testing.B) {
		receiver := hefloat.NewCiphertext(params, 1, ciphertext1.Level())
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			if err := eval.MulRelinThenAdd(ciphertext1, ciphertext2, receiver); err != nil {
				b.Log(err)
				b.Fail()
			}
		}
	})

	b.Run(GetBenchName(params, "Evaluator/Rescale"), func(b *testing.B) {
		receiver := hefloat.NewCiphertext(params, 1, ciphertext1.Level()-1)
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
		receiver := hefloat.NewCiphertext(params, 1, ciphertext1.Level())
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			if err := eval.Rotate(ciphertext1, 1, receiver); err != nil {
				b.Log(err)
				b.Fail()
			}
		}
	})
}
