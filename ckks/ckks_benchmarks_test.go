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

	b.Run(testString("Encode/", params.params), func(b *testing.B) {

		values := make([]complex128, slots)
		for i := uint64(0); i < slots; i++ {
			values[i] = complex(randomFloat(-1, 1), randomFloat(-1, 1))
		}

		plaintext := NewPlaintext(params.params, params.params.MaxLevel(), params.params.Scale())

		for i := 0; i < b.N; i++ {
			encoder.Encode(plaintext, values, slots)
		}
	})

	b.Run(testString("Decode/", params.params), func(b *testing.B) {

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

	b.Run(testString("KeyPairGen/", params.params), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			kgen.GenKeyPair()
		}
	})

	b.Run(testString("SwitchKeyGen/", params.params), func(b *testing.B) {
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

	plaintext := NewPlaintext(params.params, params.params.MaxLevel(), params.params.Scale())
	ciphertext := NewCiphertextRandom(params.prng, params.params, 1, params.params.MaxLevel(), params.params.Scale())

	b.Run(testString("", params.params), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			decryptor.Decrypt(ciphertext, plaintext)
		}
	})
}

func benchEvaluator(b *testing.B) {

	evaluator := params.evaluator

	ciphertext1 := NewCiphertextRandom(params.prng, params.params, 1, params.params.MaxLevel(), params.params.Scale())
	ciphertext2 := NewCiphertextRandom(params.prng, params.params, 1, params.params.MaxLevel(), params.params.Scale())
	receiver := NewCiphertextRandom(params.prng, params.params, 2, params.params.MaxLevel(), params.params.Scale())

	rlk := params.kgen.GenRelinKey(params.sk)
	rotkey := NewRotationKeys()
	params.kgen.GenRot(RotationLeft, params.sk, 1, rotkey)
	params.kgen.GenRot(Conjugate, params.sk, 0, rotkey)

	b.Run(testString("Add/", params.params), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			evaluator.Add(ciphertext1, ciphertext2, ciphertext1)
		}
	})

	b.Run(testString("AddScalar/", params.params), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			evaluator.AddConst(ciphertext1, ciphertext2, ciphertext1)
		}
	})

	b.Run(testString("MulScalar/", params.params), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			evaluator.MultByConst(ciphertext1, complex(3.1415, -1.4142), ciphertext1)
		}
	})

	b.Run(testString("Rescale/", params.params), func(b *testing.B) {

		contextQ := params.ckkscontext.contextQ

		for i := 0; i < b.N; i++ {
			contextQ.DivRoundByLastModulusNTT(ciphertext1.Value()[0])
			contextQ.DivRoundByLastModulusNTT(ciphertext1.Value()[1])

			b.StopTimer()
			ciphertext1 = NewCiphertextRandom(params.prng, params.params, 1, params.params.MaxLevel(), params.params.Scale())
			b.StartTimer()
		}
	})

	b.Run(testString("Mul/", params.params), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			evaluator.MulRelin(ciphertext1, ciphertext2, nil, receiver)
		}
	})

	b.Run(testString("Square/", params.params), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			evaluator.MulRelin(ciphertext1, ciphertext1, nil, receiver)
		}
	})

	b.Run(testString("Relin/", params.params), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			evaluator.Relinearize(receiver, rlk, ciphertext1)
		}
	})

	b.Run(testString("Conjugate/", params.params), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			evaluator.Conjugate(ciphertext1, rotkey, ciphertext1)
		}
	})

	b.Run(testString("Rotate/", params.params), func(b *testing.B) {
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

	contextQ := params.ckkscontext.contextQ
	contextP := params.ckkscontext.contextP

	c2NTT := ciphertext.value[1]
	c2InvNTT := contextQ.NewPoly()
	contextQ.InvNTTLvl(ciphertext.Level(), c2NTT, c2InvNTT)

	c2QiQDecomp := make([]*ring.Poly, params.params.Beta())
	c2QiPDecomp := make([]*ring.Poly, params.params.Beta())

	for i := uint64(0); i < params.params.Beta(); i++ {
		c2QiQDecomp[i] = contextQ.NewPoly()
		c2QiPDecomp[i] = contextP.NewPoly()
	}

	b.Run(testString("DecomposeNTT/", params.params), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			for j := uint64(0); j < params.params.Beta(); j++ {
				evaluator.decomposeAndSplitNTT(ciphertext.Level(), j, c2NTT, c2InvNTT, c2QiQDecomp[j], c2QiPDecomp[j])
			}
		}
	})

	b.Run(testString("RotateHoisted/", params.params), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			evaluator.permuteNTTHoisted(ciphertext, c2QiQDecomp, c2QiPDecomp, 5, rotkey, ciphertext)
		}
	})

}
