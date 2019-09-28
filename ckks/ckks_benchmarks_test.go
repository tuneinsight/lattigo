package ckks

import (
	"fmt"
	"testing"
)

type benchParams struct {
	params *Parameters
	bdc    uint64
}

func BenchmarkCKKSScheme(b *testing.B) {

	var err error
	var ckkscontext *CkksContext
	var encoder *Encoder
	var kgen *keygenerator
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
		{params: DefaultParams[12], bdc: 60},
		{params: DefaultParams[13], bdc: 60},
		{params: DefaultParams[14], bdc: 60},
		{params: DefaultParams[15], bdc: 60},
	}

	var logN, logScale, levels, bdc uint64
	var sigma float64

	for _, param := range params {

		logN = uint64(param.params.LogN)
		logScale = uint64(param.params.Logscale)
		sigma = param.params.Sigma
		bdc = param.bdc
		levels = uint64(len(param.params.Modulichain))

		if ckkscontext, err = NewCkksContext(param.params); err != nil {
			b.Error(err)
		}

		encoder = ckkscontext.NewEncoder()

		kgen = ckkscontext.NewKeyGenerator()

		if sk, pk, err = kgen.NewKeyPair(1.0 / 3); err != nil {
			b.Error(err)
		}

		if rlk, err = kgen.NewRelinKey(sk, bdc); err != nil {
			b.Error(err)
		}

		if rotkey, err = kgen.NewRotationKeysPow2(sk, bdc, true); err != nil {
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

		ciphertext1 = ckkscontext.NewCiphertext(1, levels-1, logScale)
		ciphertext2 = ckkscontext.NewCiphertext(1, levels-1, logScale)

		var values []complex128
		b.Run(fmt.Sprintf("logN=%d/logQ=%d/Levels=%d/decomp=%d/sigma=%.2f/Encode", logN, ckkscontext.LogQ(), levels, bdc, sigma), func(b *testing.B) {
			slots := 1 << (logN - 2)
			values = make([]complex128, slots)
			for i := 0; i < slots; i++ {
				values[i] = complex(randomFloat(0.1, 1), 0)
			}

			plaintext = ckkscontext.NewPlaintext(levels-1, logScale)

			for i := 0; i < b.N; i++ {
				if err = encoder.EncodeComplex(plaintext, values); err != nil {
					b.Error(err)
				}
			}
		})

		b.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/decomp=%d/sigma=%.2f/Decode", logN, ckkscontext.LogQ(), levels, bdc, sigma), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				if values = encoder.DecodeComplex(plaintext); err != nil {
					b.Error(err)
				}
			}
		})

		// Key Pair Generation
		b.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/decomp=%d/sigma=%.2f/KeyPairGen", logN, ckkscontext.LogQ(), levels, bdc, sigma), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				sk, pk, err = kgen.NewKeyPair(1.0 / 3)
				if err != nil {
					b.Error(err)
				}
			}
		})

		// SwitchKeyGen
		b.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/decomp=%d/sigma=%.2f/SwitchKeyGen", logN, ckkscontext.LogQ(), levels, bdc, sigma), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				if rlk, err = kgen.NewRelinKey(sk, bdc); err != nil {
					b.Error(err)
				}
			}
		})

		// Encrypt
		b.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/decomp=%d/sigma=%.2f/EncryptPk", logN, ckkscontext.LogQ(), levels, bdc, sigma), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				if err = encryptorPk.Encrypt(plaintext, ciphertext1); err != nil {
					b.Error(err)
				}
			}
		})

		// Encrypt
		b.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/decomp=%d/sigma=%.2f/EncryptSk", logN, ckkscontext.LogQ(), levels, bdc, sigma), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				if err = encryptorSk.Encrypt(plaintext, ciphertext1); err != nil {
					b.Error(err)
				}
			}
		})

		// Decrypt
		b.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/decomp=%d/sigma=%.2f/Decrypt", logN, ckkscontext.LogQ(), levels, bdc, sigma), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				if err = decryptor.Decrypt(ciphertext1, plaintext); err != nil {
					b.Error(err)
				}
			}
		})

		// Add
		b.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/decomp=%d/sigma=%.2f/Add", logN, ckkscontext.LogQ(), levels, bdc, sigma), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				if err = evaluator.Add(ciphertext1, ciphertext2, ciphertext1); err != nil {
					b.Error(err)
				}
			}
		})

		// Add Scalar
		b.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/decomp=%d/sigma=%.2f/AddScalar", logN, ckkscontext.LogQ(), levels, bdc, sigma), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				if err = evaluator.AddConst(ciphertext1, complex(3.1415, -1.4142), ciphertext1); err != nil {
					b.Error(err)
				}
			}
		})

		// Mult Scalar
		b.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/decomp=%d/sigma=%.2f/MultScalar", logN, ckkscontext.LogQ(), levels, bdc, sigma), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				if err = evaluator.MultConst(ciphertext1, complex(3.1415, -1.4142), ciphertext1); err != nil {
					b.Error(err)
				}
			}
		})

		// Mul
		receiver := ckkscontext.NewRandomCiphertext(2, levels-1, logScale)
		b.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/decomp=%d/sigma=%.2f/Multiply", logN, ckkscontext.LogQ(), levels, bdc, sigma), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				if err = evaluator.MulRelin(ciphertext1, ciphertext2, nil, receiver); err != nil {
					b.Error(err)
				}
			}
		})

		// Square
		b.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/decomp=%d/sigma=%.2f/Square", logN, ckkscontext.LogQ(), levels, bdc, sigma), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				if err = evaluator.MulRelin(ciphertext1, ciphertext1, nil, receiver); err != nil {
					b.Error(err)
				}
			}
		})

		// Mul Relin
		b.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/decomp=%d/sigma=%.2f/MulRelin", logN, ckkscontext.LogQ(), levels, bdc, sigma), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				if err = evaluator.MulRelin(ciphertext1, ciphertext2, rlk, ciphertext1); err != nil {
					b.Error(err)
				}
			}
		})

		// Mul Relin
		b.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/decomp=%d/sigma=%.2f/SquareRelin", logN, ckkscontext.LogQ(), levels, bdc, sigma), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				if err = evaluator.MulRelin(ciphertext1, ciphertext1, rlk, ciphertext1); err != nil {
					b.Error(err)
				}
			}
		})

		// Relin
		b.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/decomp=%d/sigma=%.2f/Relin", logN, ckkscontext.LogQ(), levels, bdc, sigma), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				if err = evaluator.Relinearize(receiver, rlk, ciphertext1); err != nil {
					b.Error(err)
				}
			}
		})

		// Rescale
		b.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/decomp=%d/sigma=%.2f/Rescale", logN, ckkscontext.LogQ(), levels, bdc, sigma), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				rescale(evaluator, ciphertext1.Value()[0], receiver.Value()[0])
				rescale(evaluator, ciphertext1.Value()[1], receiver.Value()[1])
			}
		})

		// Conjugate / Rotate
		b.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/decomp=%d/sigma=%.2f/Conjugate", logN, ckkscontext.LogQ(), levels, bdc, sigma), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				if err = evaluator.Conjugate(ciphertext1, rotkey, ciphertext1); err != nil {
					b.Error(err)
				}
			}
		})

		// Rotate Cols
		b.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/decomp=%d/sigma=%.2f/RotateCols", logN, ckkscontext.LogQ(), levels, bdc, sigma), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				if err = evaluator.RotateColumns(ciphertext1, 1, rotkey, ciphertext1); err != nil {
					b.Error(err)
				}
			}
		})
	}
}
