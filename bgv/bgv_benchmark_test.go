package bgv

import (
	"runtime"
	"testing"

	"github.com/tuneinsight/lattigo/v3/rlwe"
)

func BenchmarkBGV(b *testing.B) {

	var err error

	for _, p := range TestParams[:] {

		var params Parameters
		if params, err = NewParametersFromLiteral(p); err != nil {
			b.Error(err)
			b.Fail()
		}

		var tc *testContext
		if tc, err = genTestParams(params); err != nil {
			b.Error(err)
			b.Fail()
		}

		for _, testSet := range []func(tc *testContext, b *testing.B){
			benchEncoder,
			benchKeyGenerator,
			benchEncryptor,
			benchEvaluator,
		} {
			testSet(tc, b)
			runtime.GC()
		}
	}
}

func benchEncoder(tc *testContext, b *testing.B) {

	poly := tc.uSampler.ReadNew()

	tc.params.RingT().Reduce(poly, poly)

	scale := NewScale(tc.params, 1)

	coeffsUint64 := poly.Coeffs[0]

	coeffsInt64 := make([]int64, len(coeffsUint64))
	for i := range coeffsUint64 {
		coeffsInt64[i] = int64(coeffsUint64[i])
	}

	encoder := tc.encoder

	for _, lvl := range tc.testLevel {
		plaintext := NewPlaintext(tc.params, lvl, scale)
		b.Run(GetTestName("Encoder/Encode/Uint", tc.params, lvl), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				encoder.Encode(coeffsUint64, plaintext)
			}
		})
	}

	for _, lvl := range tc.testLevel {
		plaintext := NewPlaintext(tc.params, lvl, scale)
		b.Run(GetTestName("Encoder/Encode/Int", tc.params, lvl), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				encoder.Encode(coeffsInt64, plaintext)
			}
		})
	}

	for _, lvl := range tc.testLevel {
		plaintext := NewPlaintext(tc.params, lvl, scale)
		b.Run(GetTestName("Encoder/Decode/Uint", tc.params, lvl), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				encoder.DecodeUint(plaintext, coeffsUint64)
			}
		})
	}

	for _, lvl := range tc.testLevel {
		plaintext := NewPlaintext(tc.params, lvl, scale)
		b.Run(GetTestName("Encoder/Decode/Int", tc.params, lvl), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				encoder.DecodeInt(plaintext, coeffsInt64)
			}
		})
	}
}

func benchKeyGenerator(tc *testContext, b *testing.B) {

	kgen := tc.kgen

	b.Run(GetTestName("KeyGen/KeyPairGen", tc.params, tc.params.MaxLevel()), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			kgen.GenKeyPair()
		}
	})

	b.Run(GetTestName("KeyGen/SwitchKeyGen", tc.params, tc.params.MaxLevel()), func(b *testing.B) {
		sk := tc.sk
		for i := 0; i < b.N; i++ {
			kgen.GenRelinearizationKey(sk, 1)
		}
	})
}

func benchEncryptor(tc *testContext, b *testing.B) {

	scale := NewScale(tc.params, 1)

	for _, lvl := range tc.testLevel {
		b.Run(GetTestName("Encrypt/key=Pk", tc.params, lvl), func(b *testing.B) {
			plaintext := NewPlaintext(tc.params, lvl, scale)
			ciphertext := NewCiphertext(tc.params, 1, lvl)
			encryptorPk := tc.encryptorPk
			for i := 0; i < b.N; i++ {
				encryptorPk.Encrypt(plaintext, ciphertext)
			}
		})
	}

	for _, lvl := range tc.testLevel {
		b.Run(GetTestName("Encrypt/key=Sk", tc.params, lvl), func(b *testing.B) {
			plaintext := NewPlaintext(tc.params, lvl, scale)
			ciphertext := NewCiphertext(tc.params, 1, lvl)
			encryptorSk := tc.encryptorSk
			for i := 0; i < b.N; i++ {
				encryptorSk.Encrypt(plaintext, ciphertext)
			}
		})
	}
}

func benchEvaluator(tc *testContext, b *testing.B) {

	scale := NewScale(tc.params, 1)

	eval := tc.evaluator

	for _, lvl := range tc.testLevel {
		ciphertext0 := NewCiphertextRandom(tc.prng, tc.params, 1, lvl)
		ciphertext1 := NewCiphertextRandom(tc.prng, tc.params, 1, lvl)
		b.Run(GetTestName("Evaluator/Add/op0=ct/op1=ct", tc.params, lvl), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				eval.Add(ciphertext0, ciphertext1, ciphertext0)
			}
		})
	}

	for _, lvl := range tc.testLevel {
		ciphertext0 := NewCiphertextRandom(tc.prng, tc.params, 1, lvl)
		plaintext1 := &Plaintext{Plaintext: &rlwe.Plaintext{Value: NewCiphertextRandom(tc.prng, tc.params, 0, lvl).Value[0], Scale: scale.CopyNew()}}
		b.Run(GetTestName("Evaluator/Add/op0=ct/op1=pt", tc.params, lvl), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				eval.Add(ciphertext0, plaintext1, ciphertext0)
			}
		})
	}

	for _, lvl := range tc.testLevel {
		ciphertext0 := NewCiphertextRandom(tc.prng, tc.params, 1, lvl)
		scalar := tc.params.T() >> 1
		b.Run(GetTestName("Evaluator/AddScalar/op0=ct", tc.params, lvl), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				eval.AddScalar(ciphertext0, scalar, ciphertext0)
			}
		})
	}

	for _, lvl := range tc.testLevel {
		ciphertext0 := NewCiphertextRandom(tc.prng, tc.params, 1, lvl)
		scalar := tc.params.T() >> 1
		b.Run(GetTestName("Evaluator/MulScalar/op0=ct", tc.params, lvl), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				eval.MulScalar(ciphertext0, scalar, ciphertext0)
			}
		})
	}

	for _, lvl := range tc.testLevel {
		ciphertext0 := NewCiphertextRandom(tc.prng, tc.params, 1, lvl)
		ciphertext1 := NewCiphertextRandom(tc.prng, tc.params, 1, lvl)
		scalar := tc.params.T() >> 1
		b.Run(GetTestName("Evaluator/MulScalarAndAdd/op0=ct", tc.params, lvl), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				eval.MulScalarAndAdd(ciphertext0, scalar, ciphertext1)
			}
		})
	}

	for _, lvl := range tc.testLevel {
		ciphertext0 := NewCiphertextRandom(tc.prng, tc.params, 1, lvl)
		ciphertext1 := NewCiphertextRandom(tc.prng, tc.params, 1, lvl)
		receiver := NewCiphertext(tc.params, 2, lvl)
		b.Run(GetTestName("Evaluator/Mul/op0=ct/op1=ct", tc.params, lvl), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				eval.Mul(ciphertext0, ciphertext1, receiver)
			}
		})
	}

	for _, lvl := range tc.testLevel {
		ciphertext0 := NewCiphertextRandom(tc.prng, tc.params, 1, lvl)
		plaintext1 := &Plaintext{Plaintext: &rlwe.Plaintext{Value: NewCiphertextRandom(tc.prng, tc.params, 0, lvl).Value[0], Scale: scale.CopyNew()}}
		b.Run(GetTestName("Evaluator/Mul/op0=ct/op1=pt", tc.params, lvl), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				eval.Mul(ciphertext0, plaintext1, ciphertext0)
			}
		})
	}

	for _, lvl := range tc.testLevel {
		ciphertext0 := NewCiphertextRandom(tc.prng, tc.params, 1, lvl)
		ciphertext1 := NewCiphertextRandom(tc.prng, tc.params, 1, lvl)
		b.Run(GetTestName("Evaluator/MulRelin/op0=ct/op1=ct", tc.params, lvl), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				eval.MulRelin(ciphertext0, ciphertext1, ciphertext0)
			}
		})
	}

	for _, lvl := range tc.testLevel {
		ciphertext0 := NewCiphertextRandom(tc.prng, tc.params, 1, lvl)
		ciphertext1 := NewCiphertextRandom(tc.prng, tc.params, 1, lvl)
		ciphertext2 := NewCiphertextRandom(tc.prng, tc.params, 1, lvl)
		b.Run(GetTestName("Evaluator/MulRelinAndAdd/op0=ct/op1=ct", tc.params, lvl), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				eval.MulRelinAndAdd(ciphertext0, ciphertext1, ciphertext2)
			}
		})
	}

	for _, lvl := range tc.testLevel[1:] {
		ciphertext0 := NewCiphertextRandom(tc.prng, tc.params, 1, lvl)
		receiver := NewCiphertext(tc.params, 1, lvl-1)
		b.Run(GetTestName("Evaluator/Rescale/op0=ct", tc.params, lvl), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				if err := eval.Rescale(ciphertext0, receiver); err != nil {
					b.Log(err)
					b.Fail()
				}
			}
		})
	}

	for _, lvl := range tc.testLevel[1:] {
		ciphertext0 := NewCiphertextRandom(tc.prng, tc.params, 1, lvl)
		receiver := NewCiphertext(tc.params, 1, lvl-1)
		b.Run(GetTestName("Evaluator/Rescale(OLD)/op0=ct", tc.params, lvl), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				if err := eval.(*evaluator).rescaleOriginal(ciphertext0, receiver); err != nil {
					b.Log(err)
					b.Fail()
				}
			}
		})
	}

	for _, lvl := range tc.testLevel {
		ciphertext0 := NewCiphertextRandom(tc.prng, tc.params, 2, lvl)
		receiver := NewCiphertext(tc.params, 1, lvl)
		b.Run(GetTestName("Evaluator/Relin/op0=ct", tc.params, lvl), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				eval.Relinearize(ciphertext0, receiver)
			}
		})
	}

	for _, lvl := range tc.testLevel {
		rotkey := tc.kgen.GenRotationKeysForRotations([]int{}, true, tc.sk)
		eval := eval.WithKey(rlwe.EvaluationKey{Rtks: rotkey})
		ciphertext0 := NewCiphertextRandom(tc.prng, tc.params, 1, lvl)
		b.Run(GetTestName("Evaluator/RotateRwos/op0=ct", tc.params, lvl), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				eval.RotateRows(ciphertext0, ciphertext0)
			}
		})
	}

	for _, lvl := range tc.testLevel {
		rotkey := tc.kgen.GenRotationKeysForRotations([]int{1}, false, tc.sk)
		eval := eval.WithKey(rlwe.EvaluationKey{Rtks: rotkey})
		ciphertext0 := NewCiphertextRandom(tc.prng, tc.params, 1, lvl)
		b.Run(GetTestName("Evaluator/RotateColumns/op0=ct", tc.params, lvl), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				eval.RotateColumns(ciphertext0, 1, ciphertext0)
			}
		})
	}
}
