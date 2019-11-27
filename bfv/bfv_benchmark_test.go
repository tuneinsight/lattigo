package bfv

import (
	"fmt"
	"testing"
)

func Benchmark_BFV(b *testing.B) {

	paramSets := DefaultParams

	for _, params := range paramSets {

		bfvContext := NewContext()
		bfvContext.SetParameters(&params)

		var sk *SecretKey
		var pk *PublicKey
		var err error

		kgen := NewKeyGenerator(&params)

		// Public Key Generation
		b.Run(testString("KeyGen", &params), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				sk, pk = kgen.NewKeyPair()
				if err != nil {
					b.Error(err)
				}
			}
		})

		// Encryption
		encryptorPk := NewEncryptorFromPk(pk, &params)
		encryptorSk := NewEncryptorFromSk(sk, &params)

		ringCtx := NewRingContext(&params)

		ptcoeffs := bfvContext.NewRandomPlaintextCoeffs()
		pt := NewPlaintext(ringCtx)
		pt.setCoefficientsUint64(bfvContext, ptcoeffs)

		ctd1 := NewCiphertext(1, ringCtx)
		b.Run(testString("EncryptFromPk", &params), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				encryptorPk.Encrypt(pt, ctd1)
			}
		})

		b.Run(testString("EncryptFromSk", &params), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				encryptorSk.Encrypt(pt, ctd1)
			}
		})

		// Decryption
		decryptor := NewDecryptor(sk, &params)
		ptp := NewPlaintext(ringCtx)
		b.Run(testString("Decrypt", &params), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				decryptor.Decrypt(ctd1, ptp)
			}
			_ = ptp
		})

		evaluator := NewEvaluator(&params)

		ct1 := encryptorSk.EncryptNew(pt)
		ct2 := encryptorSk.EncryptNew(pt)

		// Addition
		b.Run(testString("Add", &params), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				evaluator.Add(ct1, ct2, ctd1)
			}
		})

		// Subtraction
		b.Run(testString("Sub", &params), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				evaluator.Sub(ct1, ct2, ctd1)
			}
		})

		// Multiplication
		ctd2 := NewCiphertext(2, ringCtx)
		b.Run(testString("Multiply", &params), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				evaluator.Mul(ct1, ct2, ctd2)
			}
		})

		// Square is Mul(ct, ct) for now
		b.Run(testString("Square", &params), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				evaluator.Mul(ct1, ct1, ctd2)
			}
		})

		// Relinearization Key Generation not becnhmarked (no inplace gen)
		rlk := kgen.NewRelinKey(sk, 2)

		// Relinearization
		b.Run(testString("Relin", &params), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				evaluator.Relinearize(ctd2, rlk, ctd1)
			}
		})

		// Rotation Key Generation not benchmarked (no inplace gen)
		rtk := kgen.NewRotationKeysPow2(sk)

		// Rotation Rows
		b.Run(testString("RotateRows", &params), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				evaluator.RotateRows(ct1, rtk, ctd1)
			}
		})

		// Rotation Cols
		b.Run(testString("RotateCols", &params), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				evaluator.RotateColumns(ct1, 1, rtk, ctd1)
			}
		})

	}
}

func testString(opname string, params *Parameters) string {
	return fmt.Sprintf("%s/params=%d", opname, params.N)
}
