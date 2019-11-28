package ckks

import (
	"github.com/ldsec/lattigo/ring"
	"testing"
)

type benchParams struct {
	params *Parameters
}

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

		ckkscontext := params.ckkscontext
		encoder := params.encoder
		slots := params.ckkscontext.Slots()

		b.Run(testString("Encode/", params), func(b *testing.B) {

			values := make([]complex128, slots)
			for i := uint64(0); i < slots; i++ {
				values[i] = complex(randomFloat(0.1, 1), 0)
			}

			plaintext := ckkscontext.NewPlaintext(ckkscontext.levels-1, ckkscontext.scale)

			for i := 0; i < b.N; i++ {
				encoder.Encode(plaintext, values, slots)
			}
		})

		b.Run(testString("Decode/", params), func(b *testing.B) {

			values := make([]complex128, slots)
			for i := uint64(0); i < slots; i++ {
				values[i] = complex(randomFloat(0.1, 1), 0)
			}

			plaintext := ckkscontext.NewPlaintext(ckkscontext.levels-1, ckkscontext.scale)
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

		b.Run(testString("KeyPairGen/", params), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				kgen.NewKeyPair()
			}
		})

		b.Run(testString("SwitchKeyGen/", params), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				kgen.NewRelinKey(sk)
			}
		})
	}
}

func benchEncrypt(b *testing.B) {
	for _, parameters := range testParams.ckksParameters {

		params := genCkksParams(parameters)
		ckkscontext := params.ckkscontext
		encryptorPk := params.encryptorPk
		encryptorSk := params.encryptorSk

		plaintext := ckkscontext.NewPlaintext(ckkscontext.levels-1, ckkscontext.scale)
		ciphertext := ckkscontext.NewCiphertext(1, ckkscontext.levels-1, ckkscontext.scale)

		b.Run(testString("Sk/", params), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				encryptorPk.Encrypt(plaintext, ciphertext)
			}
		})

		b.Run(testString("Pk/", params), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				encryptorSk.Encrypt(plaintext, ciphertext)
			}
		})
	}
}

func benchDecrypt(b *testing.B) {
	for _, parameters := range testParams.ckksParameters {

		params := genCkksParams(parameters)
		ckkscontext := params.ckkscontext
		decryptor := params.decryptor

		plaintext := ckkscontext.NewPlaintext(ckkscontext.levels-1, ckkscontext.scale)
		ciphertext := ckkscontext.NewRandomCiphertext(1, ckkscontext.levels-1, ckkscontext.scale)

		b.Run(testString("", params), func(b *testing.B) {
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
		ckkscontext := params.ckkscontext

		ciphertext1 := ckkscontext.NewRandomCiphertext(1, ckkscontext.levels-1, ckkscontext.scale)
		ciphertext2 := ckkscontext.NewRandomCiphertext(1, ckkscontext.levels-1, ckkscontext.scale)
		receiver := ckkscontext.NewRandomCiphertext(2, ckkscontext.levels-1, ckkscontext.scale)

		rlk := params.kgen.NewRelinKey(params.sk)
		rotkey := ckkscontext.NewRotationKeys()
		params.kgen.GenRot(RotationLeft, params.sk, 1, rotkey)
		params.kgen.GenRot(Conjugate, params.sk, 0, rotkey)

		b.Run(testString("Add/", params), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				evaluator.Add(ciphertext1, ciphertext2, ciphertext1)
			}
		})

		b.Run(testString("AddScalar/", params), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				evaluator.AddConst(ciphertext1, ciphertext2, ciphertext1)
			}
		})

		b.Run(testString("MulScalar/", params), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				evaluator.MultByConst(ciphertext1, complex(3.1415, -1.4142), ciphertext1)
			}
		})

		b.Run(testString("Rescale/", params), func(b *testing.B) {

			contextQ := ckkscontext.contextQ

			for i := 0; i < b.N; i++ {
				contextQ.DivRoundByLastModulusNTT(ciphertext1.Value()[0])
				contextQ.DivRoundByLastModulusNTT(ciphertext1.Value()[1])

				b.StopTimer()
				ciphertext1 = ckkscontext.NewRandomCiphertext(1, ckkscontext.levels-1, ckkscontext.scale)
				b.StartTimer()
			}
		})

		b.Run(testString("Mul/", params), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				evaluator.MulRelin(ciphertext1, ciphertext2, nil, receiver)
			}
		})

		b.Run(testString("Square/", params), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				evaluator.MulRelin(ciphertext1, ciphertext1, nil, receiver)
			}
		})

		b.Run(testString("Relin/", params), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				evaluator.Relinearize(receiver, rlk, ciphertext1)
			}
		})

		b.Run(testString("Conjugate/", params), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				evaluator.Conjugate(ciphertext1, rotkey, ciphertext1)
			}
		})

		b.Run(testString("Rotate/", params), func(b *testing.B) {
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
		ckkscontext := params.ckkscontext

		rotkey := ckkscontext.NewRotationKeys()
		params.kgen.GenRot(RotationLeft, params.sk, 5, rotkey)

		ciphertext := ckkscontext.NewRandomCiphertext(1, ckkscontext.Levels()-1, ckkscontext.Scale())

		contextQ := ckkscontext.contextQ
		contextP := ckkscontext.contextP

		c2NTT := ciphertext.value[1]
		c2InvNTT := contextQ.NewPoly()
		contextQ.InvNTTLvl(ciphertext.Level(), c2NTT, c2InvNTT)

		c2QiQDecomp := make([]*ring.Poly, ckkscontext.beta)
		c2QiPDecomp := make([]*ring.Poly, ckkscontext.beta)

		for i := uint64(0); i < ckkscontext.beta; i++ {
			c2QiQDecomp[i] = contextQ.NewPoly()
			c2QiPDecomp[i] = contextP.NewPoly()
		}

		b.Run(testString("DecomposeNTT/", params), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				for j := uint64(0); j < ckkscontext.beta; j++ {
					evaluator.decomposeAndSplitNTT(ciphertext.Level(), j, c2NTT, c2InvNTT, c2QiQDecomp[j], c2QiPDecomp[j])
				}
			}
		})

		b.Run(testString("RotateHoisted/", params), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				evaluator.switchKeyHoisted(ciphertext, c2QiQDecomp, c2QiPDecomp, 5, rotkey, ciphertext)
			}
		})
	}
}
