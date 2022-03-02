package bfv

import (
	"encoding/json"
	"testing"

	"github.com/tuneinsight/lattigo/v3/rlwe"
)

func BenchmarkBFV(b *testing.B) {

	defaultParams := DefaultParams
	if testing.Short() {
		defaultParams = DefaultParams[:2]
	} else {

	}

	if *flagParamString != "" {
		var jsonParams ParametersLiteral
		json.Unmarshal([]byte(*flagParamString), &jsonParams)
		defaultParams = []ParametersLiteral{jsonParams} // the custom test suite reads the parameters from the -params flag
	}

	for _, p := range defaultParams {

		params, err := NewParametersFromLiteral(p)
		var tc *testContext
		if tc, err = genTestParams(params); err != nil {
			panic(err)
		}

		benchEncoder(tc, b)
		benchKeyGen(tc, b)
		benchEncrypt(tc, b)
		benchDecrypt(tc, b)
		benchEvaluator(tc, b)
	}
}

func benchEncoder(tc *testContext, b *testing.B) {

	encoder := tc.encoder
	coeffs := tc.uSampler.ReadNew()
	coeffsOut := make([]uint64, tc.params.N())

	plaintext := NewPlaintext(tc.params)
	plaintextRingT := NewPlaintextRingT(tc.params)
	plaintextMul := NewPlaintextMul(tc.params)

	b.Run(testString("Encoder/EncodeUint", tc.params, tc.params.MaxLevel()), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			encoder.Encode(coeffs.Coeffs[0], plaintext)
		}
	})

	b.Run(testString("Encoder/DecodeUint/pt=Plaintext", tc.params, tc.params.MaxLevel()), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			encoder.DecodeUint(plaintext, coeffsOut)
		}
	})

	b.Run(testString("Encoder/EncodeUintRingT", tc.params, tc.params.MaxLevel()), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			encoder.EncodeRingT(coeffs.Coeffs[0], plaintextRingT)
		}
	})

	b.Run(testString("Encoder/DecodeUint/pt=PlaintextRingT", tc.params, tc.params.MaxLevel()), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			encoder.DecodeUint(plaintextRingT, coeffsOut)
		}
	})

	b.Run(testString("Encoder/EncodeUintMul", tc.params, tc.params.MaxLevel()), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			encoder.EncodeMul(coeffs.Coeffs[0], plaintextMul)
		}
	})
}

func benchKeyGen(tc *testContext, b *testing.B) {

	kgen := tc.kgen
	sk := tc.sk

	b.Run(testString("KeyGen/KeyPairGen", tc.params, tc.params.MaxLevel()), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			kgen.GenKeyPair()
		}
	})

	b.Run(testString("KeyGen/SwitchKeyGen", tc.params, tc.params.MaxLevel()), func(b *testing.B) {

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

	plaintext := NewPlaintext(tc.params)
	ciphertext := NewCiphertextRandom(tc.prng, tc.params, 1)

	b.Run(testString("Encrypt/key=Pk", tc.params, tc.params.MaxLevel()), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			encryptorPk.Encrypt(plaintext, ciphertext)
		}
	})

	b.Run(testString("Encrypt/key=Sk", tc.params, tc.params.MaxLevel()), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			encryptorSk.Encrypt(plaintext, ciphertext)
		}
	})
}

func benchDecrypt(tc *testContext, b *testing.B) {

	decryptor := tc.decryptor
	ciphertext := NewCiphertextRandom(tc.prng, tc.params, 1)

	b.Run(testString("Decrypt/", tc.params, tc.params.MaxLevel()), func(b *testing.B) {
		plaintext := NewPlaintext(tc.params)
		for i := 0; i < b.N; i++ {
			decryptor.Decrypt(ciphertext, plaintext)
		}
	})
}

func benchEvaluator(tc *testContext, b *testing.B) {

	encoder := tc.encoder

	plaintext := NewPlaintext(tc.params)
	plaintextRingT := NewPlaintextRingT(tc.params)
	plaintextMul := NewPlaintextMul(tc.params)

	coeffs := tc.uSampler.ReadNew()
	encoder.EncodeRingT(coeffs.Coeffs[0], plaintextRingT)
	encoder.Encode(coeffs.Coeffs[0], plaintext)
	encoder.EncodeMul(coeffs.Coeffs[0], plaintextMul)

	ciphertext1 := NewCiphertextRandom(tc.prng, tc.params, 1)
	ciphertext2 := NewCiphertextRandom(tc.prng, tc.params, 1)
	receiver := NewCiphertextRandom(tc.prng, tc.params, 2)

	var rotkey *rlwe.RotationKeySet
<<<<<<< dev_bfv_poly
	if tc.params.PCount() != 0 {
		rotkey = tc.kgen.GenRotationKeysForRotations([]int{1}, true, tc.sk)
=======
	if testctx.params.PCount() != 0 {
<<<<<<< dev_bfv_poly
		rotkey = testctx.kgen.GenRotationKeysForRotations([]int{1}, true, testctx.sk, 0)
>>>>>>> First step for adding bit-decomp
=======
		rotkey = testctx.kgen.GenRotationKeysForRotations([]int{1}, true, testctx.sk)
>>>>>>> fixed bfv & ckks
	}
	evaluator := tc.evaluator.WithKey(rlwe.EvaluationKey{Rlk: tc.rlk, Rtks: rotkey})

	b.Run(testString("Evaluator/Add/Ct", tc.params, tc.params.MaxLevel()), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			evaluator.Add(ciphertext1, ciphertext2, ciphertext1)
		}
	})

	b.Run(testString("Evaluator/Add/op1=Ciphertext/op2=PlaintextRingT", tc.params, tc.params.MaxLevel()), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			evaluator.Add(ciphertext1, plaintextRingT, ciphertext1)
		}
	})

	b.Run(testString("Evaluator/Add/op1=Ciphertext/op2=Plaintext", tc.params, tc.params.MaxLevel()), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			evaluator.Add(ciphertext1, plaintext, ciphertext1)
		}
	})

	b.Run(testString("Evaluator/MulScalar", tc.params, tc.params.MaxLevel()), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			evaluator.MulScalar(ciphertext1, 5, ciphertext1)
		}
	})

	b.Run(testString("Evaluator/Mul/op1=Ciphertext/op2=Ciphertext", tc.params, tc.params.MaxLevel()), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			evaluator.Mul(ciphertext1, ciphertext2, receiver)
		}
	})

	b.Run(testString("Evaluator/Mul/op1=Ciphertext/op2=Plaintext/", tc.params, tc.params.MaxLevel()), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			evaluator.Mul(ciphertext1, plaintext, ciphertext1)
		}
	})

	b.Run(testString("Evaluator/Mul/op1=Ciphertext/op2=PlaintextRingT", tc.params, tc.params.MaxLevel()), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			evaluator.Mul(ciphertext1, plaintextRingT, ciphertext1)
		}
	})

	b.Run(testString("Evaluator/Mul/op1=Ciphertext/op2=PlaintextMul", tc.params, tc.params.MaxLevel()), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			evaluator.Mul(ciphertext1, plaintextMul, ciphertext1)
		}
	})

	b.Run(testString("Evaluator/Square", tc.params, tc.params.MaxLevel()), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			evaluator.Mul(ciphertext1, ciphertext1, receiver)
		}
	})

	b.Run(testString("Evaluator/Relin", tc.params, tc.params.MaxLevel()), func(b *testing.B) {

		if tc.params.PCount() == 0 {
			b.Skip("#Pi is empty")
		}

		for i := 0; i < b.N; i++ {
			evaluator.Relinearize(receiver, ciphertext1)
		}
	})

	b.Run(testString("Evaluator/RotateRows", tc.params, tc.params.MaxLevel()), func(b *testing.B) {

		if tc.params.PCount() == 0 {
			b.Skip("#Pi is empty")
		}

		for i := 0; i < b.N; i++ {
			evaluator.RotateRows(ciphertext1, ciphertext1)
		}
	})

	b.Run(testString("Evaluator/RotateCols", tc.params, tc.params.MaxLevel()), func(b *testing.B) {

		if tc.params.PCount() == 0 {
			b.Skip("#Pi is empty")
		}

		for i := 0; i < b.N; i++ {
			evaluator.RotateColumns(ciphertext1, 1, ciphertext1)
		}
	})
}
