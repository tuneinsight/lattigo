package ckks

import (
	"fmt"
	"testing"
)

type benchParams struct {
	params *Parameters
}

func Benchmark_CKKSScheme(b *testing.B) {

	var err error
	var ckkscontext *CkksContext
	var encoder *Encoder
	var kgen *KeyGenerator
	var sk *SecretKey
	var pk *PublicKey
	var rlk *EvaluationKey
	var rotkey *RotationKey
	var encryptorPk *Encryptor
	var encryptorSk *Encryptor
	var decryptor *Decryptor
	var evaluator *Evaluator
	var plaintext *Plaintext
	var ciphertext1 *Ciphertext
	var ciphertext2 *Ciphertext

	params := []benchParams{
		//{params: DefaultParams[12]},
		//{params: DefaultParams[13]},
		{params: DefaultParams[14]},
		//{params: DefaultParams[15]},
		//{params: DefaultParams[16]},
	}

	var logN, slots, levels uint64
	var scale, sigma float64

	for _, param := range params {

		logN = uint64(param.params.LogN)
		scale = param.params.Scale
		sigma = param.params.Sigma
		levels = uint64(len(param.params.Modulichain))
		slots = uint64(1 << (logN - 1))

		if ckkscontext, err = NewCkksContext(param.params); err != nil {
			b.Error(err)
		}

		encoder = ckkscontext.NewEncoder()

		kgen = ckkscontext.NewKeyGenerator()

		sk, pk = kgen.NewKeyPair()

		rlk = kgen.NewRelinKey(sk)

		rotkey, err = kgen.NewRotationKeys(sk, []uint64{1}, nil, true)
		if err != nil {
			b.Error(err)
		}

		if encryptorPk, err = ckkscontext.NewEncryptorFromPk(pk); err != nil {
			b.Error(err)
		}

		if encryptorSk, err = ckkscontext.NewEncryptorFromSk(sk); err != nil {
			b.Error(err)
		}

		if decryptor, err = ckkscontext.NewDecryptor(sk); err != nil {
			b.Error(err)
		}

		evaluator = ckkscontext.NewEvaluator()

		_ = rotkey

		ciphertext1 = ckkscontext.NewCiphertext(1, levels-1, scale)
		ciphertext2 = ckkscontext.NewCiphertext(1, levels-1, scale)

		var values []complex128
		b.Run(fmt.Sprintf("logN=%d/logQ=%d/Levels=%d/sigma=%.2f/Encode", logN, ckkscontext.LogQ(), levels, sigma), func(b *testing.B) {

			values = make([]complex128, slots)
			for i := uint64(0); i < slots; i++ {
				values[i] = complex(randomFloat(0.1, 1), 0)
			}

			plaintext = ckkscontext.NewPlaintext(levels-1, scale)

			for i := 0; i < b.N; i++ {
				if err = encoder.Encode(plaintext, values, slots); err != nil {
					b.Error(err)
				}
			}
		})

		b.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/sigma=%.2f/Decode", logN, ckkscontext.LogQ(), levels, sigma), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				if values = encoder.Decode(plaintext, slots); err != nil {
					b.Error(err)
				}
			}
		})

		// Key Pair Generation
		b.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/sigma=%.2f/KeyPairGen", logN, ckkscontext.LogQ(), levels, sigma), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				sk, pk = kgen.NewKeyPair()
			}
		})

		// SwitchKeyGen
		b.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/sigma=%.2f/SwitchKeyGen", logN, ckkscontext.LogQ(), levels, sigma), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				rlk = kgen.NewRelinKey(sk)
			}
		})

		// Encrypt
		b.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/sigma=%.2f/EncryptPk", logN, ckkscontext.LogQ(), levels, sigma), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				if err = encryptorPk.Encrypt(plaintext, ciphertext1); err != nil {
					b.Error(err)
				}
			}
		})

		// Encrypt
		b.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/sigma=%.2f/EncryptSk", logN, ckkscontext.LogQ(), levels, sigma), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				if err = encryptorSk.Encrypt(plaintext, ciphertext1); err != nil {
					b.Error(err)
				}
			}
		})

		// Decrypt
		b.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/sigma=%.2f/Decrypt", logN, ckkscontext.LogQ(), levels, sigma), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				decryptor.Decrypt(ciphertext1, plaintext)
			}
		})

		// Add
		b.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/sigma=%.2f/Add", logN, ckkscontext.LogQ(), levels, sigma), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				if err = evaluator.Add(ciphertext1, ciphertext2, ciphertext1); err != nil {
					b.Error(err)
				}
			}
		})

		// Add Scalar
		b.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/sigma=%.2f/AddScalar", logN, ckkscontext.LogQ(), levels, sigma), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				if err = evaluator.AddConst(ciphertext1, complex(3.1415, -1.4142), ciphertext1); err != nil {
					b.Error(err)
				}
			}
		})

		// Mult Scalar
		b.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/sigma=%.2f/MultScalar", logN, ckkscontext.LogQ(), levels, sigma), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				if err = evaluator.MultConst(ciphertext1, complex(3.1415, -1.4142), ciphertext1); err != nil {
					b.Error(err)
				}
			}
		})

		// Rescale
		receiver := ckkscontext.NewRandomCiphertext(2, levels-1, scale)
		b.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/sigma=%.2f/Rescale", logN, ckkscontext.LogQ(), levels, sigma), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				rescale(evaluator, ciphertext1.Value()[0], receiver.Value()[0])
				rescale(evaluator, ciphertext1.Value()[1], receiver.Value()[1])
			}
		})

		b.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/sigma=%.2f/RescaleAll", logN, ckkscontext.LogQ(), levels, sigma), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				rescale(evaluator, ciphertext1.Value()[0], receiver.Value()[0])
				rescale(evaluator, ciphertext1.Value()[1], receiver.Value()[1])

				for receiver.Level() != 0 {
					rescale(evaluator, receiver.Value()[0], receiver.Value()[0])
					rescale(evaluator, receiver.Value()[1], receiver.Value()[1])
				}

				b.StopTimer()
				receiver = ckkscontext.NewRandomCiphertext(2, levels-1, scale)
				b.StartTimer()
			}
		})

		b.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/sigma=%.2f/RescaleMultipleAll", logN, ckkscontext.LogQ(), levels, sigma), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				rescaleMany(evaluator, ciphertext1.Level(), ciphertext1.Value()[0], receiver.Value()[0])
				rescaleMany(evaluator, ciphertext1.Level(), ciphertext1.Value()[1], receiver.Value()[1])

				b.StopTimer()
				receiver = ckkscontext.NewRandomCiphertext(2, levels-1, scale)
				b.StartTimer()
			}
		})

		// Mul
		receiver = ckkscontext.NewRandomCiphertext(2, levels-1, scale)
		b.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/sigma=%.2f/Multiply", logN, ckkscontext.LogQ(), levels, sigma), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				if err = evaluator.MulRelin(ciphertext1, ciphertext2, nil, receiver); err != nil {
					b.Error(err)
				}
			}
		})

		// Square
		b.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/sigma=%.2f/Square", logN, ckkscontext.LogQ(), levels, sigma), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				if err = evaluator.MulRelin(ciphertext1, ciphertext1, nil, receiver); err != nil {
					b.Error(err)
				}
			}
		})

		// Mul Relin
		b.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/sigma=%.2f/MulRelin", logN, ckkscontext.LogQ(), levels, sigma), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				if err = evaluator.MulRelin(ciphertext1, ciphertext2, rlk, ciphertext1); err != nil {
					b.Error(err)
				}
			}
		})

		// Mul Relin
		b.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/sigma=%.2f/SquareRelin", logN, ckkscontext.LogQ(), levels, sigma), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				if err = evaluator.MulRelin(ciphertext1, ciphertext1, rlk, ciphertext1); err != nil {
					b.Error(err)
				}
			}
		})

		// Relin
		b.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/sigma=%.2f/Relin", logN, ckkscontext.LogQ(), levels, sigma), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				if err = evaluator.Relinearize(receiver, rlk, ciphertext1); err != nil {
					b.Error(err)
				}
			}
		})

		// Conjugate / Rotate
		b.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/sigma=%.2f/Conjugate", logN, ckkscontext.LogQ(), levels, sigma), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				if err = evaluator.Conjugate(ciphertext1, rotkey, ciphertext1); err != nil {
					b.Error(err)
				}
			}
		})

		// Rotate Cols
		b.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/sigma=%.2f/RotateCols", logN, ckkscontext.LogQ(), levels, sigma), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				if err = evaluator.RotateColumns(ciphertext1, 1, rotkey, ciphertext1); err != nil {
					b.Error(err)
				}
			}
		})
	}
}
