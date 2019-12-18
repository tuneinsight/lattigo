package ckks

import (
	"github.com/ldsec/lattigo/ring"
	"testing"
)

func BenchmarkCKKSScheme(b *testing.B) {
	b.Run("Encoder", benchEncoder)
	b.Run("KeyGen", benchKeyGen)
	b.Run("Encrypt", benchEncrypt)
	b.Run("Decrypt", benchDecrypt)
	b.Run("Evaluator", benchEvaluator)
	b.Run("HoistedRotations", benchHoistedRotations)
}

func benchEncoder(b *testing.B) {

	for _, parameters := range testParams.ckksParameters {

		params := genCkksParams(parameters)

		encoder := params.encoder
		slots := uint64(1 << parameters.LogSlots)

		b.Run(testString("Encode/", parameters), func(b *testing.B) {

			values := make([]complex128, slots)
			for i := uint64(0); i < slots; i++ {
				values[i] = complex(randomFloat(0.1, 1), 0)
			}

			plaintext := NewPlaintext(parameters, parameters.MaxLevel(), parameters.Scale)

			for i := 0; i < b.N; i++ {
				encoder.Encode(plaintext, values, slots)
			}
		})

		b.Run(testString("Decode/", parameters), func(b *testing.B) {

			values := make([]complex128, slots)
			for i := uint64(0); i < slots; i++ {
				values[i] = complex(randomFloat(0.1, 1), 0)
			}

			plaintext := NewPlaintext(parameters, parameters.MaxLevel(), parameters.Scale)
			encoder.Encode(plaintext, values, slots)

			for i := 0; i < b.N; i++ {
				encoder.Decode(plaintext, slots)
			}
		})
	}
}

func benchKeyGen(b *testing.B) {

	for _, parameters := range testParams.ckksParameters {

		params := genCkksParams(parameters)
		kgen := params.kgen
		sk := params.sk

		b.Run(testString("KeyPairGen/", parameters), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				kgen.GenKeyPair()
			}
		})

		b.Run(testString("SwitchKeyGen/", parameters), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				kgen.GenRelinKey(sk)
			}
		})
	}
}

func benchEncrypt(b *testing.B) {
	for _, parameters := range testParams.ckksParameters {

		params := genCkksParams(parameters)
		encryptorPk := params.encryptorPk
		encryptorSk := params.encryptorSk

		plaintext := NewPlaintext(parameters, parameters.MaxLevel(), parameters.Scale)
		ciphertext := NewCiphertext(parameters, 1, parameters.MaxLevel(), parameters.Scale)

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
	for _, parameters := range testParams.ckksParameters {

		params := genCkksParams(parameters)
		decryptor := params.decryptor

		plaintext := NewPlaintext(parameters, parameters.MaxLevel(), parameters.Scale)
		ciphertext := NewCiphertextRandom(parameters, 1, parameters.MaxLevel(), parameters.Scale)

		b.Run(testString("", parameters), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				decryptor.Decrypt(ciphertext, plaintext)
			}
		})
	}
}

func benchEvaluator(b *testing.B) {
	for _, parameters := range testParams.ckksParameters {
		params := genCkksParams(parameters)
		evaluator := params.evaluator

		ciphertext1 := NewCiphertextRandom(parameters, 1, parameters.MaxLevel(), parameters.Scale)
		ciphertext2 := NewCiphertextRandom(parameters, 1, parameters.MaxLevel(), parameters.Scale)
		receiver := NewCiphertextRandom(parameters, 2, parameters.MaxLevel(), parameters.Scale)

		rlk := params.kgen.GenRelinKey(params.sk)
		rotkey := NewRotationKeys()
		params.kgen.GenRot(RotationLeft, params.sk, 1, rotkey)
		params.kgen.GenRot(Conjugate, params.sk, 0, rotkey)

		b.Run(testString("Add/", parameters), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				evaluator.Add(ciphertext1, ciphertext2, ciphertext1)
			}
		})

		b.Run(testString("AddScalar/", parameters), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				evaluator.AddConst(ciphertext1, ciphertext2, ciphertext1)
			}
		})

		b.Run(testString("MulScalar/", parameters), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				evaluator.MultByConst(ciphertext1, complex(3.1415, -1.4142), ciphertext1)
			}
		})

		b.Run(testString("Rescale/", parameters), func(b *testing.B) {

			contextQ := params.ckkscontext.contextQ

			for i := 0; i < b.N; i++ {
				contextQ.DivRoundByLastModulusNTT(ciphertext1.Value()[0])
				contextQ.DivRoundByLastModulusNTT(ciphertext1.Value()[1])

				b.StopTimer()
				ciphertext1 = NewCiphertextRandom(parameters, 1, parameters.MaxLevel(), parameters.Scale)
				b.StartTimer()
			}
		})

		b.Run(testString("Mul/", parameters), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				evaluator.MulRelin(ciphertext1, ciphertext2, nil, receiver)
			}
		})

		b.Run(testString("Square/", parameters), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				evaluator.MulRelin(ciphertext1, ciphertext1, nil, receiver)
			}
		})

		b.Run(testString("Relin/", parameters), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				evaluator.Relinearize(receiver, rlk, ciphertext1)
			}
		})

		b.Run(testString("Conjugate/", parameters), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				evaluator.Conjugate(ciphertext1, rotkey, ciphertext1)
			}
		})

		b.Run(testString("Rotate/", parameters), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				evaluator.RotateColumns(ciphertext1, 1, rotkey, ciphertext1)
			}
		})

	}
}

func benchHoistedRotations(b *testing.B) {

	for _, parameters := range testParams.ckksParameters {

		params := genCkksParams(parameters)
		evaluator := params.evaluator

		rotkey := NewRotationKeys()
		params.kgen.GenRot(RotationLeft, params.sk, 5, rotkey)

		ciphertext := NewCiphertextRandom(parameters, 1, parameters.MaxLevel(), parameters.Scale)

		contextQ := params.ckkscontext.contextQ
		contextP := params.ckkscontext.contextP

		c2NTT := ciphertext.value[1]
		c2InvNTT := contextQ.NewPoly()
		contextQ.InvNTTLvl(ciphertext.Level(), c2NTT, c2InvNTT)

		c2QiQDecomp := make([]*ring.Poly, parameters.Beta())
		c2QiPDecomp := make([]*ring.Poly, parameters.Beta())

		for i := uint64(0); i < parameters.Beta(); i++ {
			c2QiQDecomp[i] = contextQ.NewPoly()
			c2QiPDecomp[i] = contextP.NewPoly()
		}

		b.Run(testString("DecomposeNTT/", parameters), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				for j := uint64(0); j < parameters.Beta(); j++ {
					evaluator.decomposeAndSplitNTT(ciphertext.Level(), j, c2NTT, c2InvNTT, c2QiQDecomp[j], c2QiPDecomp[j])
				}
			}
		})

		b.Run(testString("RotateHoisted/", parameters), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				evaluator.switchKeyHoisted(ciphertext, c2QiQDecomp, c2QiPDecomp, 5, rotkey, ciphertext)
			}
		})
	}
}
