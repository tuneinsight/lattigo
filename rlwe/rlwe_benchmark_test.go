package rlwe

import (
	"encoding/json"
	"runtime"
	"testing"
)

func BenchmarkRLWE(b *testing.B) {
	defaultParams := TestParams
	if testing.Short() {
		defaultParams = TestParams[:2]
	}
	if *flagParamString != "" {
		var jsonParams ParametersLiteral
		if err := json.Unmarshal([]byte(*flagParamString), &jsonParams); err != nil {
			b.Fatal(err)
		}
		defaultParams = []ParametersLiteral{jsonParams} // the custom test suite reads the parameters from the -params flag
	}

	for _, defaultParam := range defaultParams {
		params, err := NewParametersFromLiteral(defaultParam)
		if err != nil {
			b.Fatal(err)
		}

		kgen := NewKeyGenerator(params)
		eval := NewEvaluator(params, nil)

		for _, testSet := range []func(kgen KeyGenerator, eval *Evaluator, b *testing.B){
			benchHoistedKeySwitch,
		} {
			testSet(kgen, eval, b)
			runtime.GC()
		}
	}
}

func benchHoistedKeySwitch(kgen KeyGenerator, eval *Evaluator, b *testing.B) {

	params := kgen.(*keyGenerator).params
	skIn := kgen.GenSecretKey()
	skOut := kgen.GenSecretKey()
	plaintext := NewPlaintext(params, params.MaxLevel())
	plaintext.Value.IsNTT = true
	encryptor := NewEncryptor(params, skIn)
	ciphertext := NewCiphertextNTT(params, 1, plaintext.Level())
	encryptor.Encrypt(plaintext, ciphertext)

	swk := kgen.GenSwitchingKey(skIn, skOut)

	b.Run(testString(params, "DecomposeNTT/"), func(b *testing.B) {
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			eval.DecomposeNTT(ciphertext.Level(), params.PCount()-1, params.PCount(), ciphertext.Value[1], eval.BuffDecompQP)
		}
	})

	b.Run(testString(params, "KeySwitchHoisted/"), func(b *testing.B) {
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			eval.KeyswitchHoisted(ciphertext.Level(), eval.BuffDecompQP, swk, ciphertext.Value[0], ciphertext.Value[1], eval.BuffQP[1].P, eval.BuffQP[2].P)
		}
	})
}
