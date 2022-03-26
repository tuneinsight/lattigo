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
		json.Unmarshal([]byte(*flagParamString), &jsonParams)
		defaultParams = []ParametersLiteral{jsonParams} // the custom test suite reads the parameters from the -params flag
	}

	for _, defaultParam := range defaultParams {
		params, err := NewParametersFromLiteral(defaultParam)
		if err != nil {
			panic(err)
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
<<<<<<< dev_bfv_poly
			keySwitcher.DecomposeNTT(ciphertext.Level(), params.PCount()-1, params.PCount(), ciphertext.Value[1], keySwitcher.BuffDecompQP)
=======
			eval.DecomposeNTT(ciphertext.Level(), params.PCount()-1, params.PCount(), ciphertext.Value[1], eval.PoolDecompQP)
>>>>>>> [rlwe]: complete refactoring
		}
	})

	b.Run(testString(params, "KeySwitchHoisted/"), func(b *testing.B) {
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
<<<<<<< dev_bfv_poly
			keySwitcher.KeyswitchHoisted(ciphertext.Level(), keySwitcher.BuffDecompQP, swk, ciphertext.Value[0], ciphertext.Value[1], keySwitcher.BuffQP[1].P, keySwitcher.BuffQP[2].P)
=======
			eval.KeyswitchHoisted(ciphertext.Level(), eval.PoolDecompQP, swk, ciphertext.Value[0], ciphertext.Value[1], eval.Pool[1].P, eval.Pool[2].P)
>>>>>>> [rlwe]: complete refactoring
		}
	})
}
