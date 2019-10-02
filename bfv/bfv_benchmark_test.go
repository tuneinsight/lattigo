package bfv

import (
	"fmt"
	"testing"
)

func Benchmark_BFVScheme(b *testing.B) {

	paramSets := DefaultParams

	bitDecomps := []uint64{60}

	for _, params := range paramSets {

		bfvContext := NewBfvContext()
		if err := bfvContext.SetParameters(&params); err != nil {
			b.Error(err)
		}

		var sk *SecretKey
		var pk *PublicKey
		var err error

		kgen := bfvContext.NewKeyGenerator()

		// Public Key Generation
		b.Run(fmt.Sprintf("params=%d/KeyGen", params.N), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				sk, pk = kgen.NewKeyPair()
				if err != nil {
					b.Error(err)
				}
			}
		})

		// Encryption
		encryptorPk, err := bfvContext.NewEncryptorFromPk(pk)
		encryptorSk, err := bfvContext.NewEncryptorFromSk(sk)

		if err != nil {
			b.Error(err)
		}

		ptcoeffs := bfvContext.NewRandomPlaintextCoeffs()
		pt := bfvContext.NewPlaintext()
		pt.setCoefficientsUint64(bfvContext, ptcoeffs)
		b.Run(fmt.Sprintf("params=%d/EncryptFromPkNew", params.N), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				_, _ = encryptorPk.EncryptNew(pt)
			}
		})

		ctd1 := bfvContext.NewCiphertext(1)
		b.Run(fmt.Sprintf("params=%d/EncryptFromPk", params.N), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				_ = encryptorPk.Encrypt(pt, ctd1)
			}
		})

		b.Run(fmt.Sprintf("params=%d/EncryptFromSkNew", params.N), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				_, _ = encryptorSk.EncryptNew(pt)
			}
		})

		b.Run(fmt.Sprintf("params=%d/EncryptFromSk", params.N), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				_ = encryptorSk.Encrypt(pt, ctd1)
			}
		})

		// Decryption
		decryptor, err := bfvContext.NewDecryptor(sk)
		if err != nil {
			b.Error(err)
		}
		ptp := bfvContext.NewPlaintext()
		b.Run(fmt.Sprintf("params=%d/Decrypt", params.N), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				decryptor.Decrypt(ctd1, ptp)
			}
			_ = ptp
		})

		evaluator := bfvContext.NewEvaluator()

		ct1, err := encryptorSk.EncryptNew(pt)
		if err != nil {
			b.Error(err)
		}

		ct2, err := encryptorSk.EncryptNew(pt)
		if err != nil {
			b.Error(err)
		}

		// Addition
		b.Run(fmt.Sprintf("params=%d/Add", params.N), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				if err := evaluator.Add(ct1, ct2, ctd1); err != nil {
					b.Error(err)
				}
			}
		})

		// Subtraction
		b.Run(fmt.Sprintf("params=%d/Sub", params.N), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				if err := evaluator.Sub(ct1, ct2, ctd1); err != nil {
					b.Error(err)
				}
			}
		})

		// Multiplication
		ctd2 := bfvContext.NewCiphertext(2)
		b.Run(fmt.Sprintf("params=%d/Multiply", params.N), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				_ = evaluator.Mul(ct1, ct2, ctd2)
			}
		})

		// Square is Mul(ct, ct) for now
		b.Run(fmt.Sprintf("params=%d/Square", params.N), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				_ = evaluator.Mul(ct1, ct1, ctd2)
			}
		})

		for _, bitDecomp := range bitDecomps {

			// Relinearization Key Generation not becnhmarked (no inplace gen)
			rlk := kgen.NewRelinKey(sk, 2, bitDecomp)

			// Relinearization
			b.Run(fmt.Sprintf("params=%d/decomp=%d/Relin", params.N, bitDecomp), func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					if err := evaluator.Relinearize(ctd2, rlk, ctd1); err != nil {
						b.Error(err)
					}
				}
			})

			// Rotation Key Generation not benchmarked (no inplace gen)
			rtk := kgen.NewRotationKeysPow2(sk, bitDecomp, true)

			// Rotation Rows
			b.Run(fmt.Sprintf("params=%d/decomp=%d/RotateRows", params.N, bitDecomp), func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					if err := evaluator.RotateRows(ct1, rtk, ctd1); err != nil {
						b.Error(err)
					}
				}
			})

			// Rotation Cols
			b.Run(fmt.Sprintf("params=%d/decomp=%d/RotateCols", params.N, bitDecomp), func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					if err := evaluator.RotateColumns(ct1, 1, rtk, ctd1); err != nil {
						b.Error(err)
					}
				}
			})
		}
	}
}
