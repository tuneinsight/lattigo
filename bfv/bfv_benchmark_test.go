package bfv

import (
	"testing"
)

func BenchmarkBFV(b *testing.B) {
	for _, p := range defaultParams {

		if err = genTestParams(p); err != nil {
			panic(err)
		}

		b.Run("Encoder", benchEncoder)
		b.Run("KeyGen", benchKeyGen)
		b.Run("Encrypt", benchEncrypt)
		b.Run("Decrypt", benchDecrypt)
		b.Run("Evaluator", benchEvaluator)
	}
}

func benchEncoder(b *testing.B) {

	encoder := params.encoder
	coeffs := params.uSampler.ReadNew()
	plaintext := NewPlaintext(params.params)

	b.Run(testString("Encode/", params.params), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			encoder.EncodeUint(coeffs.Coeffs[0], plaintext)
		}
	})

	b.Run(testString("Decode/", params.params), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			params.encoder.DecodeUint(plaintext)
		}
	})
}

func benchKeyGen(b *testing.B) {

	kgen := params.kgen
	sk := params.sk

	b.Run(testString("KeyPairGen/", params.params), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			kgen.GenKeyPair()
		}
	})

	b.Run(testString("SwitchKeyGen/", params.params), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			kgen.GenRelinKey(sk, 1)
		}
	})
}

func benchEncrypt(b *testing.B) {

	encryptorPk := params.encryptorPk
	encryptorSk := params.encryptorSk

	plaintext := NewPlaintext(params.params)
	ciphertext := NewCiphertextRandom(params.prng, params.params, 1)

	b.Run(testString("Sk/", params.params), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			encryptorPk.Encrypt(plaintext, ciphertext)
		}
	})

	b.Run(testString("Pk/", params.params), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			encryptorSk.Encrypt(plaintext, ciphertext)
		}
	})
}

func benchDecrypt(b *testing.B) {

	decryptor := params.decryptor
	plaintext := NewPlaintext(params.params)
	ciphertext := NewCiphertextRandom(params.prng, params.params, 1)

	b.Run(testString("", params.params), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			decryptor.Decrypt(ciphertext, plaintext)
		}
	})
}

func benchEvaluator(b *testing.B) {

	evaluator := params.evaluator

	ciphertext1 := NewCiphertextRandom(params.prng, params.params, 1)
	ciphertext2 := NewCiphertextRandom(params.prng, params.params, 1)
	receiver := NewCiphertextRandom(params.prng, params.params, 2)

	rotkey := NewRotationKeys()
	params.kgen.GenRot(RotationLeft, params.sk, 1, rotkey)
	params.kgen.GenRot(RotationRow, params.sk, 0, rotkey)

	b.Run(testString("Add/", params.params), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			evaluator.Add(ciphertext1, ciphertext2, ciphertext1)
		}
	})

	b.Run(testString("MulScalar/", params.params), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			evaluator.MulScalar(ciphertext1, 5, ciphertext1)
		}
	})

	b.Run(testString("Mul/", params.params), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			evaluator.Mul(ciphertext1, ciphertext2, receiver)
		}
	})

	b.Run(testString("Square/", params.params), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			evaluator.Mul(ciphertext1, ciphertext1, receiver)
		}
	})

	b.Run(testString("Relin/", params.params), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			evaluator.Relinearize(receiver, params.rlk, ciphertext1)
		}
	})

	b.Run(testString("RotateRows/", params.params), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			evaluator.RotateRows(ciphertext1, rotkey, ciphertext1)
		}
	})

	b.Run(testString("RotateCols/", params.params), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			evaluator.RotateColumns(ciphertext1, 1, rotkey, ciphertext1)
		}
	})
}
