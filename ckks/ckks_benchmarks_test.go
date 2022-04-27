package ckks

import (
	"encoding/json"
	"testing"

	"github.com/tuneinsight/lattigo/v3/ring"
	"github.com/tuneinsight/lattigo/v3/rlwe"
	"github.com/tuneinsight/lattigo/v3/utils"
)

func BenchmarkCKKSScheme(b *testing.B) {

	var err error

	defaultParams := append(DefaultParams, DefaultConjugateInvariantParams...)
	if testing.Short() {
		defaultParams = DefaultParams[:2]
	}
	if *flagParamString != "" {
		var jsonParams ParametersLiteral
		if err = json.Unmarshal([]byte(*flagParamString), &jsonParams); err != nil {
			b.Error(err)
			b.Fail()
		}
		defaultParams = []ParametersLiteral{jsonParams} // the custom test suite reads the parameters from the -params flag
	}

	for _, defaultParams := range defaultParams {
		var params Parameters
		if params, err = NewParametersFromLiteral(defaultParams); err != nil {
			b.Error(err)
			b.Fail()
		}

		var tc *testContext
		if tc, err = genTestParams(params); err != nil {
			b.Error(err)
			b.Fail()
		}

		benchEncoder(tc, b)
		benchKeyGen(tc, b)
		benchEncrypt(tc, b)
		benchDecrypt(tc, b)
		benchEvaluator(tc, b)
		benchInnerSum(tc, b)
	}
}

func benchEncoder(tc *testContext, b *testing.B) {

	encoder := tc.encoder
	logSlots := tc.params.LogSlots()

	b.Run(GetTestName(tc.params, "Encoder/Encode"), func(b *testing.B) {

		values := make([]complex128, 1<<logSlots)
		for i := 0; i < 1<<logSlots; i++ {
			values[i] = utils.RandComplex128(-1, 1)
		}

		plaintext := NewPlaintext(tc.params, tc.params.MaxLevel(), tc.params.DefaultScale())

		for i := 0; i < b.N; i++ {
			encoder.Encode(values, plaintext, logSlots)
		}
	})

	b.Run(GetTestName(tc.params, "Encoder/Decode"), func(b *testing.B) {

		values := make([]complex128, 1<<logSlots)
		for i := 0; i < 1<<logSlots; i++ {
			values[i] = utils.RandComplex128(-1, 1)
		}

		plaintext := NewPlaintext(tc.params, tc.params.MaxLevel(), tc.params.DefaultScale())
		encoder.Encode(values, plaintext, logSlots)

		for i := 0; i < b.N; i++ {
			encoder.Decode(plaintext, logSlots)
		}
	})
}

func benchKeyGen(tc *testContext, b *testing.B) {

	kgen := tc.kgen
	sk := tc.sk

	b.Run(GetTestName(tc.params, "KeyGen/KeyPairGen"), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			kgen.GenKeyPair()
		}
	})

	b.Run(GetTestName(tc.params, "KeyGen/SwitchKeyGen"), func(b *testing.B) {

		if tc.params.PCount() == 0 {
			b.Skip("#Pi is empty")
		}

		for i := 0; i < b.N; i++ {
			kgen.GenRelinearizationKey(sk, 1)
		}
	})
}

func benchEncrypt(tc *testContext, b *testing.B) {

	encryptorPk := tc.encryptorPk
	encryptorSk := tc.encryptorSk

	plaintext := NewPlaintext(tc.params, tc.params.MaxLevel(), tc.params.DefaultScale())
	ciphertext := NewCiphertext(tc.params, 1, tc.params.MaxLevel(), tc.params.DefaultScale())

	b.Run(GetTestName(tc.params, "Encrypt/key=Pk"), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			encryptorPk.Encrypt(plaintext, ciphertext)
		}
	})

	b.Run(GetTestName(tc.params, "Encrypt/key=Sk"), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			encryptorSk.Encrypt(plaintext, ciphertext)
		}
	})
}

func benchDecrypt(tc *testContext, b *testing.B) {

	decryptor := tc.decryptor

	plaintext := NewPlaintext(tc.params, tc.params.MaxLevel(), tc.params.DefaultScale())
	ciphertext := NewCiphertextRandom(tc.prng, tc.params, 1, tc.params.MaxLevel(), tc.params.DefaultScale())

	b.Run(GetTestName(tc.params, "Decrypt"), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			decryptor.Decrypt(ciphertext, plaintext)
		}
	})
}

func benchEvaluator(tc *testContext, b *testing.B) {

	plaintext := NewPlaintext(tc.params, tc.params.MaxLevel(), tc.params.DefaultScale())
	ciphertext1 := NewCiphertextRandom(tc.prng, tc.params, 1, tc.params.MaxLevel(), tc.params.DefaultScale())
	ciphertext2 := NewCiphertextRandom(tc.prng, tc.params, 1, tc.params.MaxLevel(), tc.params.DefaultScale())
	receiver := NewCiphertextRandom(tc.prng, tc.params, 2, tc.params.MaxLevel(), tc.params.DefaultScale())

	var rlk *rlwe.RelinearizationKey
	var rotkey *rlwe.RotationKeySet
	if tc.params.PCount() != 0 {
		rlk = tc.kgen.GenRelinearizationKey(tc.sk, 1)
		rotkey = tc.kgen.GenRotationKeysForRotations([]int{1}, tc.params.RingType() == ring.Standard, tc.sk)
	}

	eval := tc.evaluator.WithKey(rlwe.EvaluationKey{Rlk: rlk, Rtks: rotkey})

	b.Run(GetTestName(tc.params, "Evaluator/Add"), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			eval.Add(ciphertext1, ciphertext2, ciphertext1)
		}
	})

	b.Run(GetTestName(tc.params, "Evaluator/AddScalar"), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			eval.AddConst(ciphertext1, ciphertext2, ciphertext1)
		}
	})

	b.Run(GetTestName(tc.params, "Evaluator/MulScalar"), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			eval.MultByConst(ciphertext1, complex(3.1415, -1.4142), ciphertext1)
		}
	})

	b.Run(GetTestName(tc.params, "Evaluator/MulPlain"), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			eval.Mul(ciphertext1, plaintext, ciphertext1)
		}
	})

	b.Run(GetTestName(tc.params, "Evaluator/Mul"), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			eval.Mul(ciphertext1, ciphertext2, receiver)
		}
	})

	b.Run(GetTestName(tc.params, "Evaluator/Square"), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			eval.Mul(ciphertext1, ciphertext1, receiver)
		}
	})

	b.Run(GetTestName(tc.params, "Evaluator/Rescale"), func(b *testing.B) {

		if tc.params.PCount() == 0 {
			b.Skip("#Pi is empty")
		}

		ciphertext1.Scale = tc.params.DefaultScale() * tc.params.DefaultScale()

		for i := 0; i < b.N; i++ {
			if err := eval.Rescale(ciphertext1, tc.params.DefaultScale(), ciphertext2); err != nil {
				panic(err)
			}
		}
	})

	b.Run(GetTestName(tc.params, "Evaluator/PermuteNTTWithIndexLvl"), func(b *testing.B) {

		if tc.params.PCount() == 0 {
			b.Skip("#Pi is empty")
		}

		galEL := tc.params.GaloisElementForColumnRotationBy(1)
		for i := 0; i < b.N; i++ {
			tc.params.RingQ().PermuteNTTWithIndexLvl(ciphertext1.Level(), ciphertext1.Value[0], eval.(*evaluator).PermuteNTTIndex[galEL], ciphertext1.Value[0])
			tc.params.RingQ().PermuteNTTWithIndexLvl(ciphertext1.Level(), ciphertext1.Value[1], eval.(*evaluator).PermuteNTTIndex[galEL], ciphertext1.Value[1])
		}
	})

	b.Run(GetTestName(tc.params, "Evaluator/Conjugate"), func(b *testing.B) {

		if tc.params.RingType() != ring.Standard {
			b.Skip("#Pi is empty")
		}

		if tc.params.PCount() == 0 {
			b.Skip("#Pi is empty")
		}

		for i := 0; i < b.N; i++ {
			eval.Conjugate(ciphertext1, ciphertext1)
		}
	})

	b.Run(GetTestName(tc.params, "Evaluator/Relin"), func(b *testing.B) {

		if tc.params.PCount() == 0 {
			b.Skip("#Pi is empty")
		}

		for i := 0; i < b.N; i++ {
			eval.Relinearize(receiver, ciphertext1)
		}
	})

	b.Run(GetTestName(tc.params, "Evaluator/Rotate"), func(b *testing.B) {

		if tc.params.PCount() == 0 {
			b.Skip("#Pi is empty")
		}

		for i := 0; i < b.N; i++ {
			eval.Rotate(ciphertext1, 1, ciphertext1)
		}
	})
}

func benchInnerSum(tc *testContext, b *testing.B) {

	ciphertext1 := NewCiphertextRandom(tc.prng, tc.params, 1, tc.params.MaxLevel(), tc.params.DefaultScale())

	batch := 1
	n := 4

	b.Run(GetTestName(tc.params, "InnerSum"), func(b *testing.B) {

		if tc.params.PCount() == 0 {
			b.Skip("#Pi is empty")
		}

		rotKey := tc.kgen.GenRotationKeysForRotations(tc.params.RotationsForInnerSum(batch, n), false, tc.sk)
		eval := tc.evaluator.WithKey(rlwe.EvaluationKey{Rlk: tc.rlk, Rtks: rotKey})

		b.ResetTimer()

		for i := 0; i < b.N; i++ {
			eval.InnerSum(ciphertext1, batch, n, ciphertext1)
		}
	})

	b.Run(GetTestName(tc.params, "InnerSumLog"), func(b *testing.B) {

		if tc.params.PCount() == 0 {
			b.Skip("#Pi is empty")
		}

		rotKey := tc.kgen.GenRotationKeysForRotations(tc.params.RotationsForInnerSumLog(batch, n), false, tc.sk)
		eval := tc.evaluator.WithKey(rlwe.EvaluationKey{Rlk: tc.rlk, Rtks: rotKey})

		b.ResetTimer()

		for i := 0; i < b.N; i++ {
			eval.InnerSumLog(ciphertext1, batch, n, ciphertext1)
		}
	})

}
