package bfv

import (
	"testing"
)

func BenchmarkBFV(b *testing.B) {
	var err error
	var testctx = new(testContext)
	var defaultParams []*Parameters

	if testing.Short() {
		defaultParams = DefaultParams[PN12QP109 : PN12QP109+3]
	} else {
		defaultParams = DefaultParams
	}

	for _, p := range defaultParams {

		if testctx, err = genTestParams(p); err != nil {
			panic(err)
		}

		benchEncoder(testctx, b)
		benchKeyGen(testctx, b)
		benchEncrypt(testctx, b)
		benchDecrypt(testctx, b)
		benchEvaluator(testctx, b)
	}
}

func benchEncoder(testctx *testContext, b *testing.B) {

	encoder := testctx.encoder
	coeffs := testctx.uSampler.ReadNew()
	plaintext := NewPlaintext(testctx.params)

	b.Run(testString("Encode/", testctx.params), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			encoder.EncodeUint(coeffs.Coeffs[0], plaintext)
		}
	})

	b.Run(testString("Decode/", testctx.params), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			testctx.encoder.DecodeUint(plaintext)
		}
	})
}

func benchKeyGen(testctx *testContext, b *testing.B) {

	kgen := testctx.kgen
	sk := testctx.sk

	b.Run(testString("KeyPairGen/", testctx.params), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			kgen.GenKeyPair()
		}
	})

	b.Run(testString("SwitchKeyGen/", testctx.params), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			kgen.GenRelinKey(sk, 1)
		}
	})
}

func benchEncrypt(testctx *testContext, b *testing.B) {

	encryptorPk := testctx.encryptorPk
	encryptorSk := testctx.encryptorSk

	plaintext := NewPlaintext(testctx.params)
	ciphertext := NewCiphertextRandom(testctx.prng, testctx.params, 1)

	b.Run(testString("Sk/", testctx.params), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			encryptorPk.Encrypt(plaintext, ciphertext)
		}
	})

	b.Run(testString("Pk/", testctx.params), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			encryptorSk.Encrypt(plaintext, ciphertext)
		}
	})
}

func benchDecrypt(testctx *testContext, b *testing.B) {

	decryptor := testctx.decryptor
	plaintext := NewPlaintext(testctx.params)
	ciphertext := NewCiphertextRandom(testctx.prng, testctx.params, 1)

	b.Run(testString("", testctx.params), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			decryptor.Decrypt(ciphertext, plaintext)
		}
	})
}

func benchEvaluator(testctx *testContext, b *testing.B) {

	evaluator := testctx.evaluator

	ciphertext1 := NewCiphertextRandom(testctx.prng, testctx.params, 1)
	ciphertext2 := NewCiphertextRandom(testctx.prng, testctx.params, 1)
	receiver := NewCiphertextRandom(testctx.prng, testctx.params, 2)

	rotkey := NewRotationKeys()
	testctx.kgen.GenRot(RotationLeft, testctx.sk, 1, rotkey)
	testctx.kgen.GenRot(RotationRow, testctx.sk, 0, rotkey)

	b.Run(testString("Add/", testctx.params), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			evaluator.Add(ciphertext1, ciphertext2, ciphertext1)
		}
	})

	b.Run(testString("MulScalar/", testctx.params), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			evaluator.MulScalar(ciphertext1, 5, ciphertext1)
		}
	})

	b.Run(testString("Mul/", testctx.params), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			evaluator.Mul(ciphertext1, ciphertext2, receiver)
		}
	})

	b.Run(testString("Square/", testctx.params), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			evaluator.Mul(ciphertext1, ciphertext1, receiver)
		}
	})

	b.Run(testString("Relin/", testctx.params), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			evaluator.Relinearize(receiver, testctx.rlk, ciphertext1)
		}
	})

	b.Run(testString("RotateRows/", testctx.params), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			evaluator.RotateRows(ciphertext1, rotkey, ciphertext1)
		}
	})

	b.Run(testString("RotateCols/", testctx.params), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			evaluator.RotateColumns(ciphertext1, 1, rotkey, ciphertext1)
		}
	})
}
