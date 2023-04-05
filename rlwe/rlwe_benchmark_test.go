package rlwe

import (
	"encoding/json"
	"runtime"
	"testing"
)

func BenchmarkRLWE(b *testing.B) {

	var err error

	defaultParamsLiteral := TestParamsLiteral[:]

	if *flagParamString != "" {
		var jsonParams ParametersLiteral
		if err = json.Unmarshal([]byte(*flagParamString), &jsonParams); err != nil {
			b.Fatal(err)
		}
		defaultParamsLiteral = []ParametersLiteral{jsonParams} // the custom test suite reads the parameters from the -params flag
	}

	for _, paramsLit := range defaultParamsLiteral {

		var params Parameters
		if params, err = NewParametersFromLiteral(paramsLit); err != nil {
			b.Fatal(err)
		}

		tc := NewTestContext(params)

		for _, testSet := range []func(tc *TestContext, b *testing.B){
			benchKeyGenerator,
			benchEncryptor,
			benchDecryptor,
			benchEvaluator,
		} {
			testSet(tc, b)
			runtime.GC()
		}
	}
}

func benchKeyGenerator(tc *TestContext, b *testing.B) {

	params := tc.params
	kgen := tc.kgen

	b.Run(testString(params, params.MaxLevel(), "KeyGenerator/GenSecretKey"), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			kgen.GenSecretKey(tc.sk)
		}
	})

	b.Run(testString(params, params.MaxLevel(), "KeyGenerator/GenPublicKey"), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			kgen.GenPublicKey(tc.sk, tc.pk)
		}

	})

	b.Run(testString(params, params.MaxLevel(), "KeyGenerator/GenEvaluationKey"), func(b *testing.B) {
		sk0, sk1 := tc.sk, kgen.GenSecretKeyNew()
		evk := NewEvaluationKey(params, params.MaxLevelQ(), params.MaxLevelP())
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			kgen.GenEvaluationKey(sk0, sk1, evk)
		}
	})
}

func benchEncryptor(tc *TestContext, b *testing.B) {

	params := tc.params

	b.Run(testString(params, params.MaxLevel(), "Encryptor/EncryptZero/SecretKey"), func(b *testing.B) {
		ct := NewCiphertext(params, 1, params.MaxLevel())
		enc := tc.enc.WithKey(tc.sk)
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			enc.EncryptZero(ct)
		}

	})

	b.Run(testString(params, params.MaxLevel(), "Encryptor/EncryptZero/PublicKey"), func(b *testing.B) {
		ct := NewCiphertext(params, 1, params.MaxLevel())
		enc := tc.enc.WithKey(tc.pk)
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			enc.EncryptZero(ct)
		}
	})
}

func benchDecryptor(tc *TestContext, b *testing.B) {

	params := tc.params

	b.Run(testString(params, params.MaxLevel(), "Decryptor/Decrypt"), func(b *testing.B) {
		dec := tc.dec
		ct := tc.enc.EncryptZeroNew(params.MaxLevel())
		pt := NewPlaintext(params, ct.Level())
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			dec.Decrypt(ct, pt)
		}
	})
}

func benchEvaluator(tc *TestContext, b *testing.B) {

	params := tc.params
	kgen := tc.kgen
	sk := tc.sk
	eval := tc.eval

	b.Run(testString(params, params.MaxLevel(), "Evaluator/GadgetProduct"), func(b *testing.B) {

		ct := NewEncryptor(params, sk).EncryptZeroNew(params.MaxLevel())
		evk := kgen.GenEvaluationKeyNew(sk, kgen.GenSecretKeyNew())

		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			eval.GadgetProduct(ct.Level(), ct.Value[1], &evk.GadgetCiphertext, ct)
		}
	})
}
