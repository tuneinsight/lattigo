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

	plaintextZT := NewPlaintextZT(testctx.params)
	plaintextZQ := NewPlaintextZQ(testctx.params)

	b.Run(testString("Encoder/Encode/ZT/", testctx.params), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			encoder.EncodeUint(coeffs.Coeffs[0], plaintextZT)
		}
	})

	b.Run(testString("Encoder/Decode/ZT/", testctx.params), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			testctx.encoder.DecodeUint(plaintextZT)
		}
	})

	b.Run(testString("Encoder/Encode/ZQ/", testctx.params), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			encoder.EncodeUint(coeffs.Coeffs[0], plaintextZQ)
		}
	})

	b.Run(testString("Encoder/Decode/ZQ/", testctx.params), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			testctx.encoder.DecodeUint(plaintextZQ)
		}
	})
}

func benchKeyGen(testctx *testContext, b *testing.B) {

	kgen := testctx.kgen
	sk := testctx.sk

	b.Run(testString("KeyGen/KeyPairGen/", testctx.params), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			kgen.GenKeyPair()
		}
	})

	b.Run(testString("KeyGen/SwitchKeyGen/", testctx.params), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			kgen.GenRelinKey(sk, 1)
		}
	})
}

func benchEncrypt(testctx *testContext, b *testing.B) {

	encryptorPk := testctx.encryptorPk
	encryptorSk := testctx.encryptorSk

	plaintextZQ := NewPlaintextZQ(testctx.params)
	plaintextZT := NewPlaintextZT(testctx.params)
	ciphertext := NewCiphertextRandom(testctx.prng, testctx.params, 1)

	b.Run(testString("Encrypt/key=Pk/ZT/", testctx.params), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			encryptorPk.Encrypt(plaintextZT, ciphertext)
		}
	})

	b.Run(testString("Encrypt/key=Pk/Fast/ZT/", testctx.params), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			encryptorPk.EncryptFast(plaintextZT, ciphertext)
		}
	})

	b.Run(testString("Encrypt/key=Sk/ZT/", testctx.params), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			encryptorSk.Encrypt(plaintextZT, ciphertext)
		}
	})

	b.Run(testString("Encrypt/key=Pk/ZQ/", testctx.params), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			encryptorPk.Encrypt(plaintextZQ, ciphertext)
		}
	})

	b.Run(testString("Encrypt/key=Pk/Fast/ZQ/", testctx.params), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			encryptorPk.EncryptFast(plaintextZQ, ciphertext)
		}
	})

	b.Run(testString("Encrypt/key=Sk/ZQ/", testctx.params), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			encryptorSk.Encrypt(plaintextZQ, ciphertext)
		}
	})
}

func benchDecrypt(testctx *testContext, b *testing.B) {

	decryptor := testctx.decryptor
	plaintextZQ := NewPlaintextZQ(testctx.params)
	plaintextZT := NewPlaintextZQ(testctx.params)
	ciphertext := NewCiphertextRandom(testctx.prng, testctx.params, 1)

	b.Run(testString("Decrypt/ZT", testctx.params), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			decryptor.Decrypt(ciphertext, plaintextZT)
		}
	})

	b.Run(testString("Decrypt/ZQ", testctx.params), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			decryptor.Decrypt(ciphertext, plaintextZQ)
		}
	})
}

func benchEvaluator(testctx *testContext, b *testing.B) {

	evaluator := testctx.evaluator
	encoder := testctx.encoder

	plaintextZQ := NewPlaintextZQ(testctx.params)
	plaintextZT := NewPlaintextZT(testctx.params)

	coeffs := testctx.uSampler.ReadNew()
	encoder.EncodeUint(coeffs.Coeffs[0], plaintextZT)
	encoder.EncodeUint(coeffs.Coeffs[0], plaintextZQ)

	ciphertext1 := NewCiphertextRandom(testctx.prng, testctx.params, 1)
	ciphertext2 := NewCiphertextRandom(testctx.prng, testctx.params, 1)
	receiver := NewCiphertextRandom(testctx.prng, testctx.params, 2)

	rotkey := NewRotationKeys()
	testctx.kgen.GenRot(RotationLeft, testctx.sk, 1, rotkey)
	testctx.kgen.GenRot(RotationRow, testctx.sk, 0, rotkey)

	b.Run(testString("Evaluator/Add/Ct/", testctx.params), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			evaluator.Add(ciphertext1, ciphertext2, ciphertext1)
		}
	})

	b.Run(testString("Evaluator/Add/Pt/ZT/", testctx.params), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			evaluator.Add(ciphertext1, plaintextZT, ciphertext1)
		}
	})

	b.Run(testString("Evaluator/Add/Pt/ZQ/", testctx.params), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			evaluator.Add(ciphertext1, plaintextZQ, ciphertext1)
		}
	})

	b.Run(testString("Evaluator/MulScalar/", testctx.params), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			evaluator.MulScalar(ciphertext1, 5, ciphertext1)
		}
	})

	b.Run(testString("Evaluator/Mul/Ct/", testctx.params), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			evaluator.Mul(ciphertext1, ciphertext2, receiver)
		}
	})

	b.Run(testString("Evaluator/Mul/Pt/ZT/", testctx.params), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			evaluator.Mul(ciphertext1, plaintextZT, ciphertext1)
		}
	})

	b.Run(testString("Evaluator/Mul/Pt/ZQ/", testctx.params), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			evaluator.Mul(ciphertext1, plaintextZQ, ciphertext1)
		}
	})

	b.Run(testString("Evaluator/Square/", testctx.params), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			evaluator.Mul(ciphertext1, ciphertext1, receiver)
		}
	})

	b.Run(testString("Evaluator/Relin/", testctx.params), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			evaluator.Relinearize(receiver, testctx.rlk, ciphertext1)
		}
	})

	b.Run(testString("Evaluator/RotateRows/", testctx.params), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			evaluator.RotateRows(ciphertext1, rotkey, ciphertext1)
		}
	})

	b.Run(testString("Evaluator/RotateCols/", testctx.params), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			evaluator.RotateColumns(ciphertext1, 1, rotkey, ciphertext1)
		}
	})
}
