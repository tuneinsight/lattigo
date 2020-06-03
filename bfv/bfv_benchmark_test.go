package bfv

import (
	"github.com/ldsec/lattigo/ring"
	"github.com/ldsec/lattigo/utils"
	"testing"
)

type benchParams struct {
	params Parameters
}

func BenchmarkBFV(b *testing.B) {
	b.Run("Encoder", benchEncoder)
	b.Run("KeyGen", benchKeyGen)
	b.Run("Encrypt", benchEncrypt)
	b.Run("Decrypt", benchDecrypt)
	b.Run("Evaluator", benchEvaluator)
}

func benchEncoder(b *testing.B) {

	for _, parameters := range testParams.bfvParameters {

		params := genBfvParams(parameters)

		encoder := params.encoder
		prng, err := utils.NewPRNG()
		if err != nil {
			panic(err)
		}
		uniformSampler := ring.NewUniformSampler(prng, params.bfvContext.contextT)

		coeffs := uniformSampler.NewUniformPoly()
		plaintext := NewPlaintext(params.params)

		b.Run(testString("Encode/", parameters), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				encoder.EncodeUint(coeffs.Coeffs[0], plaintext)
			}
		})

		b.Run(testString("Decode/", parameters), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				params.encoder.DecodeUint(plaintext)
			}
		})
	}
}

func benchKeyGen(b *testing.B) {

	for _, parameters := range testParams.bfvParameters {

		params := genBfvParams(parameters)
		kgen := params.kgen
		sk := params.sk

		b.Run(testString("KeyPairGen/", parameters), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				kgen.GenKeyPair()
			}
		})

		b.Run(testString("SwitchKeyGen/", parameters), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				kgen.GenRelinKey(sk, 1)
			}
		})
	}
}

func benchEncrypt(b *testing.B) {
	for _, parameters := range testParams.bfvParameters {

		params := genBfvParams(parameters)
		encryptorPk := params.encryptorPk
		encryptorSk := params.encryptorSk

		plaintext := NewPlaintext(parameters)
		ciphertext := NewCiphertextRandom(parameters, 1)

		b.Run(testString("Sk/", parameters), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				encryptorPk.Encrypt(plaintext, ciphertext)
			}
		})

		b.Run(testString("Pk/", parameters), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				encryptorSk.Encrypt(plaintext, ciphertext)
			}
		})
	}
}

func benchDecrypt(b *testing.B) {
	for _, parameters := range testParams.bfvParameters {

		params := genBfvParams(parameters)
		decryptor := params.decryptor

		plaintext := NewPlaintext(parameters)
		ciphertext := NewCiphertextRandom(parameters, 1)

		b.Run(testString("", parameters), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				decryptor.Decrypt(ciphertext, plaintext)
			}
		})
	}
}

func benchEvaluator(b *testing.B) {
	for _, parameters := range testParams.bfvParameters {
		params := genBfvParams(parameters)
		evaluator := params.evaluator

		ciphertext1 := NewCiphertextRandom(parameters, 1)
		ciphertext2 := NewCiphertextRandom(parameters, 1)
		receiver := NewCiphertextRandom(parameters, 2)

		rlk := params.kgen.GenRelinKey(params.sk, 1)
		rotkey := NewRotationKeys()
		params.kgen.GenRot(RotationLeft, params.sk, 1, rotkey)
		params.kgen.GenRot(RotationRow, params.sk, 0, rotkey)

		b.Run(testString("Add/", parameters), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				evaluator.Add(ciphertext1, ciphertext2, ciphertext1)
			}
		})

		b.Run(testString("MulScalar/", parameters), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				evaluator.MulScalar(ciphertext1, 5, ciphertext1)
			}
		})

		b.Run(testString("Mul/", parameters), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				evaluator.Mul(ciphertext1, ciphertext2, receiver)
			}
		})

		b.Run(testString("Square/", parameters), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				evaluator.Mul(ciphertext1, ciphertext1, receiver)
			}
		})

		b.Run(testString("Relin/", parameters), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				evaluator.Relinearize(receiver, rlk, ciphertext1)
			}
		})

		b.Run(testString("RotateRows/", parameters), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				evaluator.RotateRows(ciphertext1, rotkey, ciphertext1)
			}
		})

		b.Run(testString("RotateCols/", parameters), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				evaluator.RotateColumns(ciphertext1, 1, rotkey, ciphertext1)
			}
		})

	}
}
