package rlwe

import (
	"encoding/json"
	"runtime"
	"testing"

	"github.com/stretchr/testify/require"
	"github.com/tuneinsight/lattigo/v6/utils"
)

func BenchmarkRLWE(b *testing.B) {

	var err error

	defaultParamsLiteral := testInsecure

	if *flagParamString != "" {
		var jsonParams TestParametersLiteral
		if err = json.Unmarshal([]byte(*flagParamString), &jsonParams); err != nil {
			b.Fatal(err)
		}
		defaultParamsLiteral = []TestParametersLiteral{jsonParams} // the custom test suite reads the parameters from the -params flag
	}

	for _, paramsLit := range defaultParamsLiteral[:] {

		var params Parameters
		if params, err = NewParametersFromLiteral(paramsLit.ParametersLiteral); err != nil {
			b.Fatal(err)
		}

		tc, err := NewTestContext(params)
		require.NoError(b, err)

		for _, testSet := range []func(tc *TestContext, BaseTwoDecomposition int, b *testing.B){
			benchKeyGenerator,
			benchEncryptor,
			benchDecryptor,
			benchEvaluator,
		} {
			testSet(tc, paramsLit.BaseTwoDecomposition, b)
			runtime.GC()
		}
	}
}

func benchKeyGenerator(tc *TestContext, bpw2 int, b *testing.B) {

	params := tc.params
	kgen := tc.kgen

	b.Run(testString(params, params.MaxLevelQ(), params.MaxLevelP(), bpw2, "KeyGenerator/GenSecretKey"), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			kgen.GenSecretKey(tc.sk)
		}
	})

	b.Run(testString(params, params.MaxLevelQ(), params.MaxLevelP(), bpw2, "KeyGenerator/GenPublicKey"), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			kgen.GenPublicKey(tc.sk, tc.pk)
		}

	})

	b.Run(testString(params, params.MaxLevelQ(), params.MaxLevelP(), bpw2, "KeyGenerator/GenEvaluationKey"), func(b *testing.B) {
		sk0, sk1 := tc.sk, kgen.GenSecretKeyNew()
		evk := NewEvaluationKey(params)
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			kgen.GenEvaluationKey(sk0, sk1, evk)
		}
	})
}

func benchEncryptor(tc *TestContext, bpw2 int, b *testing.B) {

	params := tc.params

	b.Run(testString(params, params.MaxLevelQ(), params.MaxLevelP(), bpw2, "Encryptor/EncryptZero/SecretKey"), func(b *testing.B) {
		ct := NewCiphertext(params, 1, params.MaxLevel())
		enc := tc.enc.WithKey(tc.sk)
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			enc.EncryptZero(ct)
		}

	})

	b.Run(testString(params, params.MaxLevelQ(), params.MaxLevelP(), bpw2, "Encryptor/EncryptZero/PublicKey"), func(b *testing.B) {
		ct := NewCiphertext(params, 1, params.MaxLevel())
		enc := tc.enc.WithKey(tc.pk)
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			enc.EncryptZero(ct)
		}
	})
}

func benchDecryptor(tc *TestContext, bpw2 int, b *testing.B) {

	params := tc.params

	b.Run(testString(params, params.MaxLevelQ(), params.MaxLevelP(), bpw2, "Decryptor/Decrypt"), func(b *testing.B) {
		dec := tc.dec
		ct := tc.enc.EncryptZeroNew(params.MaxLevel())
		pt := NewPlaintext(params, ct.Level())
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			dec.Decrypt(ct, pt)
		}
	})
}

func benchEvaluator(tc *TestContext, bpw2 int, b *testing.B) {

	params := tc.params
	kgen := tc.kgen
	sk := tc.sk
	eval := tc.eval

	levelsP := []int{0}

	if params.MaxLevelP() > 0 {
		levelsP = append(levelsP, params.MaxLevelP())
	}

	for _, levelP := range levelsP {

		b.Run(testString(params, params.MaxLevelQ(), levelP, bpw2, "Evaluator/GadgetProduct"), func(b *testing.B) {

			enc := NewEncryptor(params, sk)

			ct := enc.EncryptZeroNew(params.MaxLevel())

			evkParams := EvaluationKeyParameters{LevelQ: utils.Pointy(params.MaxLevelQ()), LevelP: utils.Pointy(levelP), BaseTwoDecomposition: utils.Pointy(bpw2)}

			evk := kgen.GenEvaluationKeyNew(sk, kgen.GenSecretKeyNew(), evkParams)

			b.ResetTimer()
			for i := 0; i < b.N; i++ {
				eval.GadgetProduct(ct.Level(), ct.Value[1], &evk.GadgetCiphertext, ct)
			}
		})
	}
}
