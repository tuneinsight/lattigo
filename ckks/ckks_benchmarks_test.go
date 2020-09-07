package ckks

import (
	"testing"

	"github.com/ldsec/lattigo/ring"
)

func BenchmarkCKKSScheme(b *testing.B) {

	for _, defaultParams := range defaultParams {

		if err = genTestParams(defaultParams); err != nil {
			panic(err)
		}

		b.Run("Encoder", benchEncoder)
		b.Run("KeyGen", benchKeyGen)
		b.Run("Encrypt", benchEncrypt)
		b.Run("Decrypt", benchDecrypt)
		b.Run("Evaluator", benchEvaluator)
		b.Run("HoistedRotations", benchHoistedRotations)
	}
}

func benchEncoder(b *testing.B) {

	encoder := params.encoder
	slots := params.params.Slots()

	b.Run(testString("Encode/"), func(b *testing.B) {

		values := make([]complex128, slots)
		for i := uint64(0); i < slots; i++ {
			values[i] = complex(randomFloat(-1, 1), randomFloat(-1, 1))
		}

		plaintext := NewPlaintext(params.params, params.params.MaxLevel(), params.params.Scale())

		for i := 0; i < b.N; i++ {
			encoder.Encode(plaintext, values, slots)
		}
	})

	b.Run(testString("Decode/"), func(b *testing.B) {

		values := make([]complex128, slots)
		for i := uint64(0); i < slots; i++ {
			values[i] = complex(randomFloat(-1, 1), randomFloat(-1, 1))
		}

		plaintext := NewPlaintext(params.params, params.params.MaxLevel(), params.params.Scale())
		encoder.Encode(plaintext, values, slots)

		for i := 0; i < b.N; i++ {
			encoder.Decode(plaintext, slots)
		}
	})

}

func benchKeyGen(b *testing.B) {

	kgen := params.kgen
	sk := params.sk

	b.Run(testString("KeyPairGen/"), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			kgen.GenKeyPair()
		}
	})

	b.Run(testString("SwitchKeyGen/"), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			kgen.GenRelinKey(sk)
		}
	})
}

func benchEncrypt(b *testing.B) {

	encryptorPk := params.encryptorPk
	encryptorSk := params.encryptorSk

	plaintext := NewPlaintext(params.params, params.params.MaxLevel(), params.params.Scale())
	ciphertext := NewCiphertext(params.params, 1, params.params.MaxLevel(), params.params.Scale())

	b.Run(testString("Pk/Slow"), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			encryptorPk.Encrypt(plaintext, ciphertext)
		}
	})

	b.Run(testString("Pk/Fast"), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			encryptorPk.EncryptFast(plaintext, ciphertext)
		}
	})

	b.Run(testString("Sk/"), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			encryptorSk.Encrypt(plaintext, ciphertext)
		}
	})
}

func benchDecrypt(b *testing.B) {

	decryptor := params.decryptor

	plaintext := NewPlaintext(params.params, params.params.MaxLevel(), params.params.Scale())
	ciphertext := NewCiphertextRandom(params.prng, params.params, 1, params.params.MaxLevel(), params.params.Scale())

	b.Run(testString(""), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			decryptor.Decrypt(ciphertext, plaintext)
		}
	})
}

func benchEvaluator(b *testing.B) {

	evaluator := params.evaluator

	plaintext := NewPlaintext(params.params, params.params.MaxLevel(), params.params.scale)
	ciphertext1 := NewCiphertextRandom(params.prng, params.params, 1, params.params.MaxLevel(), params.params.Scale())
	ciphertext2 := NewCiphertextRandom(params.prng, params.params, 1, params.params.MaxLevel(), params.params.Scale())
	receiver := NewCiphertextRandom(params.prng, params.params, 2, params.params.MaxLevel(), params.params.Scale())

	rlk := params.kgen.GenRelinKey(params.sk)
	rotkey := NewRotationKeys()
	params.kgen.GenRot(RotationLeft, params.sk, 1, rotkey)
	params.kgen.GenRot(Conjugate, params.sk, 0, rotkey)

	b.Run(testString("Add/"), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			evaluator.Add(ciphertext1, ciphertext2, ciphertext1)
		}
	})

	b.Run(testString("AddScalar/"), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			evaluator.AddConst(ciphertext1, ciphertext2, ciphertext1)
		}
	})

	b.Run(testString("MulScalar/"), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			evaluator.MultByConst(ciphertext1, complex(3.1415, -1.4142), ciphertext1)
		}
	})

	b.Run(testString("MulPlain/"), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			evaluator.MulRelin(ciphertext1, plaintext, nil, receiver)
		}
	})

	b.Run(testString("Mul/"), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			evaluator.MulRelin(ciphertext1, ciphertext2, nil, receiver)
		}
	})

	b.Run(testString("Square/"), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			evaluator.MulRelin(ciphertext1, ciphertext1, nil, receiver)
		}
	})

	b.Run(testString("Rescale/"), func(b *testing.B) {

		ringQ := params.ringQ

		for i := 0; i < b.N; i++ {
			ringQ.DivRoundByLastModulusNTT(ciphertext1.Value()[0])
			ringQ.DivRoundByLastModulusNTT(ciphertext1.Value()[1])

			b.StopTimer()
			ciphertext1 = NewCiphertextRandom(params.prng, params.params, 1, params.params.MaxLevel(), params.params.Scale())
			b.StartTimer()
		}
	})

	b.Run(testString("PermuteNTT/"), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			ring.PermuteNTTWithIndexLvl(ciphertext1.Level(), ciphertext1.value[0], rotkey.permuteNTTLeftIndex[1], ciphertext1.value[0])
			ring.PermuteNTTWithIndexLvl(ciphertext1.Level(), ciphertext1.value[1], rotkey.permuteNTTLeftIndex[1], ciphertext1.value[1])
		}
	})

	b.Run(testString("Conjugate/"), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			evaluator.Conjugate(ciphertext1, rotkey, ciphertext1)
		}
	})

	b.Run(testString("Relin/"), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			evaluator.Relinearize(receiver, rlk, ciphertext1)
		}
	})

	b.Run(testString("Conjugate/"), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			evaluator.Conjugate(ciphertext1, rotkey, ciphertext1)
		}
	})

	b.Run(testString("Rotate/"), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			evaluator.RotateColumns(ciphertext1, 1, rotkey, ciphertext1)
		}
	})
}

func benchHoistedRotations(b *testing.B) {

	evaluator := params.evaluator.(*evaluator)

	rotkey := NewRotationKeys()
	params.kgen.GenRot(RotationLeft, params.sk, 5, rotkey)

	ciphertext := NewCiphertextRandom(params.prng, params.params, 1, params.params.MaxLevel(), params.params.Scale())

	ringQ := params.ringQ
	ringP := params.ringP

	c2NTT := ciphertext.value[1]
	c2InvNTT := ringQ.NewPoly()
	ringQ.InvNTTLvl(ciphertext.Level(), c2NTT, c2InvNTT)

	c2QiQDecomp := make([]*ring.Poly, params.params.Beta())
	c2QiPDecomp := make([]*ring.Poly, params.params.Beta())

	for i := uint64(0); i < params.params.Beta(); i++ {
		c2QiQDecomp[i] = ringQ.NewPoly()
		c2QiPDecomp[i] = ringP.NewPoly()
	}

	b.Run(testString("DecomposeNTT/"), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			for j := uint64(0); j < params.params.Beta(); j++ {
				evaluator.decomposeAndSplitNTT(ciphertext.Level(), j, c2NTT, c2InvNTT, c2QiQDecomp[j], c2QiPDecomp[j])
			}
		}
	})

	b.Run(testString("RotateHoisted/"), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			evaluator.permuteNTTHoisted(ciphertext, c2QiQDecomp, c2QiPDecomp, 5, rotkey, ciphertext)
		}
	})
}
