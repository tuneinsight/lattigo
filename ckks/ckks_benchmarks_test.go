package ckks

import (
	"testing"

	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/utils"
)

func BenchmarkCKKSScheme(b *testing.B) {

	var err error

	var defaultParams []*Parameters

	if testing.Short() {
		defaultParams = DefaultParams[PN12QP109 : PN12QP109+3]
	} else {
		defaultParams = DefaultParams
	}

	for _, defaultParams := range defaultParams {
		var testContext *testParams
		if testContext, err = genTestParams(defaultParams, 0); err != nil {
			panic(err)
		}

		benchEncoder(testContext, b)
		benchKeyGen(testContext, b)
		benchEncrypt(testContext, b)
		benchDecrypt(testContext, b)
		benchEvaluator(testContext, b)
		benchInnerSum(testContext, b)
		benchHoistedRotations(testContext, b)
	}
}

func benchEncoder(testContext *testParams, b *testing.B) {

	encoder := testContext.encoder
	logSlots := testContext.params.LogSlots()

	b.Run(testString(testContext, "Encoder/Encode/"), func(b *testing.B) {

		values := make([]complex128, 1<<logSlots)
		for i := 0; i < 1<<logSlots; i++ {
			values[i] = utils.RandComplex128(-1, 1)
		}

		plaintext := NewPlaintext(testContext.params, testContext.params.MaxLevel(), testContext.params.Scale())

		for i := 0; i < b.N; i++ {
			encoder.Encode(plaintext, values, logSlots)
		}
	})

	b.Run(testString(testContext, "Encoder/Decode/"), func(b *testing.B) {

		values := make([]complex128, 1<<logSlots)
		for i := 0; i < 1<<logSlots; i++ {
			values[i] = utils.RandComplex128(-1, 1)
		}

		plaintext := NewPlaintext(testContext.params, testContext.params.MaxLevel(), testContext.params.Scale())
		encoder.Encode(plaintext, values, logSlots)

		for i := 0; i < b.N; i++ {
			encoder.Decode(plaintext, logSlots)
		}
	})
}

func benchKeyGen(testContext *testParams, b *testing.B) {

	kgen := testContext.kgen
	sk := testContext.sk

	b.Run(testString(testContext, "KeyGen/KeyPairGen/"), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			kgen.GenKeyPair()
		}
	})

	b.Run(testString(testContext, "KeyGen/SwitchKeyGen/"), func(b *testing.B) {

		if testContext.params.PiCount() == 0 {
			b.Skip("#Pi is empty")
		}

		for i := 0; i < b.N; i++ {
			kgen.GenRelinearizationKey(sk)
		}
	})
}

func benchEncrypt(testContext *testParams, b *testing.B) {

	encryptorPk := testContext.encryptorPk
	encryptorSk := testContext.encryptorSk

	plaintext := NewPlaintext(testContext.params, testContext.params.MaxLevel(), testContext.params.Scale())
	ciphertext := NewCiphertext(testContext.params, 1, testContext.params.MaxLevel(), testContext.params.Scale())

	b.Run(testString(testContext, "Encrypt/key=Pk/"), func(b *testing.B) {

		if testContext.params.PiCount() == 0 {
			b.Skip("#Pi is empty")
		}

		for i := 0; i < b.N; i++ {
			encryptorPk.Encrypt(plaintext, ciphertext)
		}
	})

	b.Run(testString(testContext, "EncryptFast/key=Pk/"), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			encryptorPk.EncryptFast(plaintext, ciphertext)
		}
	})

	b.Run(testString(testContext, "Encrypt/key=Sk/"), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			encryptorSk.Encrypt(plaintext, ciphertext)
		}
	})
}

func benchDecrypt(testContext *testParams, b *testing.B) {

	decryptor := testContext.decryptor

	plaintext := NewPlaintext(testContext.params, testContext.params.MaxLevel(), testContext.params.Scale())
	ciphertext := NewCiphertextRandom(testContext.prng, testContext.params, 1, testContext.params.MaxLevel(), testContext.params.Scale())

	b.Run(testString(testContext, "Decrypt/"), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			decryptor.Decrypt(ciphertext, plaintext)
		}
	})
}

func benchEvaluator(testContext *testParams, b *testing.B) {

	plaintext := NewPlaintext(testContext.params, testContext.params.MaxLevel(), testContext.params.scale)
	ciphertext1 := NewCiphertextRandom(testContext.prng, testContext.params, 1, testContext.params.MaxLevel(), testContext.params.Scale())
	ciphertext2 := NewCiphertextRandom(testContext.prng, testContext.params, 1, testContext.params.MaxLevel(), testContext.params.Scale())
	receiver := NewCiphertextRandom(testContext.prng, testContext.params, 2, testContext.params.MaxLevel(), testContext.params.Scale())

	var rlk *RelinearizationKey
	var rotkey *RotationKeySet
	if testContext.params.PiCount() != 0 {
		rlk = testContext.kgen.GenRelinearizationKey(testContext.sk)
		rotkey = testContext.kgen.GenRotationKeysForRotations([]int{1}, true, testContext.sk)
	}

	eval := testContext.evaluator.WithKey(EvaluationKey{rlk, rotkey})

	b.Run(testString(testContext, "Evaluator/Add/"), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			eval.Add(ciphertext1, ciphertext2, ciphertext1)
		}
	})

	b.Run(testString(testContext, "Evaluator/AddScalar/"), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			eval.AddConst(ciphertext1, ciphertext2, ciphertext1)
		}
	})

	b.Run(testString(testContext, "Evaluator/MulScalar/"), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			eval.MultByConst(ciphertext1, complex(3.1415, -1.4142), ciphertext1)
		}
	})

	b.Run(testString(testContext, "Evaluator/MulPlain/"), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			eval.Mul(ciphertext1, plaintext, ciphertext1)
		}
	})

	b.Run(testString(testContext, "Evaluator/Mul/"), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			eval.Mul(ciphertext1, ciphertext2, receiver)
		}
	})

	b.Run(testString(testContext, "Evaluator/Square/"), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			eval.Mul(ciphertext1, ciphertext1, receiver)
		}
	})

	b.Run(testString(testContext, "Evaluator/Rescale/"), func(b *testing.B) {

		if testContext.params.PiCount() == 0 {
			b.Skip("#Pi is empty")
		}

		ciphertext1.SetScale(testContext.params.Scale() * testContext.params.Scale())

		for i := 0; i < b.N; i++ {
			if err := eval.Rescale(ciphertext1, testContext.params.Scale(), ciphertext2); err != nil {
				panic(err)
			}
		}
	})

	b.Run(testString(testContext, "Evaluator/PermuteNTT/"), func(b *testing.B) {

		if testContext.params.PiCount() == 0 {
			b.Skip("#Pi is empty")
		}

		galEL := testContext.params.GaloisElementForColumnRotationBy(1)
		for i := 0; i < b.N; i++ {
			ring.PermuteNTTWithIndexLvl(ciphertext1.Level(), ciphertext1.value[0], eval.(*evaluator).permuteNTTIndex[galEL], ciphertext1.value[0])
			ring.PermuteNTTWithIndexLvl(ciphertext1.Level(), ciphertext1.value[1], eval.(*evaluator).permuteNTTIndex[galEL], ciphertext1.value[1])
		}
	})

	b.Run(testString(testContext, "Evaluator/Conjugate/"), func(b *testing.B) {

		if testContext.params.PiCount() == 0 {
			b.Skip("#Pi is empty")
		}

		for i := 0; i < b.N; i++ {
			eval.Conjugate(ciphertext1, ciphertext1)
		}
	})

	b.Run(testString(testContext, "Evaluator/Relin/"), func(b *testing.B) {

		if testContext.params.PiCount() == 0 {
			b.Skip("#Pi is empty")
		}

		for i := 0; i < b.N; i++ {
			eval.Relinearize(receiver, ciphertext1)
		}
	})

	b.Run(testString(testContext, "Evaluator/Rotate/"), func(b *testing.B) {

		if testContext.params.PiCount() == 0 {
			b.Skip("#Pi is empty")
		}

		for i := 0; i < b.N; i++ {
			eval.Rotate(ciphertext1, 1, ciphertext1)
		}
	})
}

func benchInnerSum(testContext *testParams, b *testing.B) {

	ciphertext1 := NewCiphertextRandom(testContext.prng, testContext.params, 1, testContext.params.MaxLevel(), testContext.params.Scale())

	batch := 1
	n := 7

	b.Run(testString(testContext, "InnerSum/Naive"), func(b *testing.B) {

		if testContext.params.PiCount() == 0 {
			b.Skip("#Pi is empty")
		}

		rots := []int{}
		for i := 1; i < n; i++ {
			rots = append(rots, i*batch)
		}

		rotKey := testContext.kgen.GenRotationKeysForRotations(rots, false, testContext.sk)
		eval := testContext.evaluator.WithKey(EvaluationKey{testContext.rlk, rotKey})

		for i := 0; i < b.N; i++ {
			eval.InnerSumNaive(ciphertext1, batch, n, ciphertext1)
		}
	})

	b.Run(testString(testContext, "InnerSum/Log"), func(b *testing.B) {

		if testContext.params.PiCount() == 0 {
			b.Skip("#Pi is empty")
		}

		rots := []int{}
		var rot int
		for i := 1; i < n; i <<= 1 {

			rot = i
			rot *= batch

			if !utils.IsInSliceInt(rot, rots) && rot != 0 {
				rots = append(rots, rot)
			}

			rot = n - (n & ((i << 1) - 1))

			rot *= batch

			if !utils.IsInSliceInt(rot, rots) && rot != 0 {
				rots = append(rots, rot)
			}
		}

		rotKey := testContext.kgen.GenRotationKeysForRotations(rots, false, testContext.sk)
		eval := testContext.evaluator.WithKey(EvaluationKey{testContext.rlk, rotKey})

		for i := 0; i < b.N; i++ {
			eval.InnerSum(ciphertext1, batch, n, ciphertext1)
		}
	})

}

func benchHoistedRotations(testContext *testParams, b *testing.B) {

	b.Run(testString(testContext, "HoistedRotations/"), func(b *testing.B) {

		if testContext.params.PiCount() == 0 {
			b.Skip("#Pi is empty")
		}

		rotkey := testContext.kgen.GenRotationKeysForRotations([]int{5}, false, testContext.sk)
		evaluator := testContext.evaluator.WithKey(EvaluationKey{testContext.rlk, rotkey}).(*evaluator)

		ciphertext := NewCiphertextRandom(testContext.prng, testContext.params, 1, testContext.params.MaxLevel(), testContext.params.Scale())

		ringQ := testContext.ringQ
		ringP := testContext.ringP

		c2NTT := ciphertext.value[1]
		c2InvNTT := ringQ.NewPoly()
		ringQ.InvNTTLvl(ciphertext.Level(), c2NTT, c2InvNTT)

		c2QiQDecomp := make([]*ring.Poly, testContext.params.Beta())
		c2QiPDecomp := make([]*ring.Poly, testContext.params.Beta())

		for i := 0; i < testContext.params.Beta(); i++ {
			c2QiQDecomp[i] = ringQ.NewPoly()
			c2QiPDecomp[i] = ringP.NewPoly()
		}

		b.Run(testString(testContext, "/DecomposeNTT/"), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				for j := 0; j < testContext.params.Beta(); j++ {
					evaluator.decomposeAndSplitNTT(ciphertext.Level(), j, c2NTT, c2InvNTT, c2QiQDecomp[j], c2QiPDecomp[j])
				}
			}
		})

		b.Run(testString(testContext, "RotateHoisted/"), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				evaluator.permuteNTTHoisted(ciphertext.Level(), ciphertext.value[0], ciphertext.value[1], c2QiQDecomp, c2QiPDecomp, 5, ciphertext.value[0], ciphertext.value[1])
			}
		})
	})
}
