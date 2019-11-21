package ckks

import (
	"fmt"
	"github.com/ldsec/lattigo/ring"
	"testing"
)

type benchParams struct {
	params *Parameters
}

func Benchmark_CKKSScheme(b *testing.B) {

	var ckkscontext *CkksContext
	var encoder *Encoder
	var kgen *KeyGenerator
	var sk *SecretKey
	var pk *PublicKey
	var rlk *EvaluationKey
	var rotkey *RotationKeys
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

		ckkscontext = NewCkksContext(param.params)

		encoder = ckkscontext.NewEncoder()

		kgen = ckkscontext.NewKeyGenerator()

		sk, pk = kgen.NewKeyPair()

		rlk = kgen.NewRelinKey(sk)

		rotkey = ckkscontext.NewRotationKeys()
		kgen.GenRot(RotationLeft, sk, 1, rotkey)

		encryptorPk = ckkscontext.NewEncryptorFromPk(pk)
		encryptorSk = ckkscontext.NewEncryptorFromSk(sk)
		decryptor = ckkscontext.NewDecryptor(sk)

		evaluator = ckkscontext.NewEvaluator()

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
				encoder.Encode(plaintext, values, slots)
			}
		})

		b.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/sigma=%.2f/Decode", logN, ckkscontext.LogQ(), levels, sigma), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				encoder.Decode(plaintext, slots)
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
				encryptorPk.Encrypt(plaintext, ciphertext1)
			}
		})

		// Encrypt
		b.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/sigma=%.2f/EncryptSk", logN, ckkscontext.LogQ(), levels, sigma), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				encryptorSk.Encrypt(plaintext, ciphertext1)
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
				evaluator.Add(ciphertext1, ciphertext2, ciphertext1)
			}
		})

		// Add Scalar
		b.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/sigma=%.2f/AddScalar", logN, ckkscontext.LogQ(), levels, sigma), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				evaluator.AddConst(ciphertext1, complex(3.1415, -1.4142), ciphertext1)
			}
		})

		// Mult Scalar
		b.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/sigma=%.2f/MultScalar", logN, ckkscontext.LogQ(), levels, sigma), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				evaluator.MultConst(ciphertext1, complex(3.1415, -1.4142), ciphertext1)
			}
		})

		// Rescale
		contextQ := ckkscontext.contextQ
		receiver := ckkscontext.NewRandomCiphertext(2, levels-1, scale)
		b.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/sigma=%.2f/Rescale", logN, ckkscontext.LogQ(), levels, sigma), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				contextQ.DivRoundByLastModulusNTT(ciphertext1.Value()[0])
				contextQ.DivRoundByLastModulusNTT(ciphertext1.Value()[1])

				b.StopTimer()
				ciphertext1 = ckkscontext.NewRandomCiphertext(1, levels-1, scale)
				b.StartTimer()
			}
		})

		b.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/sigma=%.2f/RescaleAll", logN, ckkscontext.LogQ(), levels, sigma), func(b *testing.B) {
			for i := 0; i < b.N; i++ {

				for ciphertext1.Level() != 0 {
					contextQ.DivRoundByLastModulusNTT(ciphertext1.Value()[0])
					contextQ.DivRoundByLastModulusNTT(ciphertext1.Value()[1])
				}

				b.StopTimer()
				ciphertext1 = ckkscontext.NewRandomCiphertext(1, levels-1, scale)
				b.StartTimer()
			}
		})

		b.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/sigma=%.2f/RescaleMultipleAll", logN, ckkscontext.LogQ(), levels, sigma), func(b *testing.B) {
			for i := 0; i < b.N; i++ {

				ckkscontext.contextQ.DivRoundByLastModulusManyNTT(ciphertext1.Value()[0], ciphertext1.Level()+1)
				ckkscontext.contextQ.DivRoundByLastModulusManyNTT(ciphertext1.Value()[1], ciphertext1.Level()+1)

				b.StopTimer()
				ciphertext1 = ckkscontext.NewRandomCiphertext(1, levels-1, scale)
				b.StartTimer()
			}
		})

		// Mul
		receiver = ckkscontext.NewRandomCiphertext(2, levels-1, scale)
		b.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/sigma=%.2f/Multiply", logN, ckkscontext.LogQ(), levels, sigma), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				evaluator.MulRelin(ciphertext1, ciphertext2, nil, receiver)
			}
		})

		// Square
		b.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/sigma=%.2f/Square", logN, ckkscontext.LogQ(), levels, sigma), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				evaluator.MulRelin(ciphertext1, ciphertext1, nil, receiver)
			}
		})

		// Mul Relin
		b.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/sigma=%.2f/MulRelin", logN, ckkscontext.LogQ(), levels, sigma), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				evaluator.MulRelin(ciphertext1, ciphertext2, rlk, ciphertext1)
			}
		})

		// Mul Relin
		b.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/sigma=%.2f/SquareRelin", logN, ckkscontext.LogQ(), levels, sigma), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				evaluator.MulRelin(ciphertext1, ciphertext1, rlk, ciphertext1)
			}
		})

		// Relin
		b.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/sigma=%.2f/Relin", logN, ckkscontext.LogQ(), levels, sigma), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				evaluator.Relinearize(receiver, rlk, ciphertext1)
			}
		})

		// Conjugate / Rotate
		b.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/sigma=%.2f/Conjugate", logN, ckkscontext.LogQ(), levels, sigma), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				evaluator.Conjugate(ciphertext1, rotkey, ciphertext1)
			}
		})

		// Rotate Cols
		b.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/sigma=%.2f/RotateCols", logN, ckkscontext.LogQ(), levels, sigma), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				evaluator.RotateColumns(ciphertext1, 1, rotkey, ciphertext1)
			}
		})
	}
}

func BenchmarkRotationHoisting(b *testing.B) {

	var ckkscontext *CkksContext
	var kgen *KeyGenerator
	var sk *SecretKey
	var evaluator *Evaluator
	var ciphertext *Ciphertext

	params := []benchParams{
		{params: &Parameters{16, []uint8{55, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 45, 45, 55, 55, 55, 55, 55, 55, 55, 55, 55, 45, 45, 45}, []uint8{55, 55, 55, 55}, 1 << 40, 3.2}},
	}

	var logN, logSlots uint64

	for _, param := range params {

		logN = uint64(param.params.LogN)

		ckkscontext = NewCkksContext(param.params)
		kgen = ckkscontext.NewKeyGenerator()

		sk = kgen.NewSecretKey()

		rotkeys := ckkscontext.NewRotationKeys()
		kgen.GenRot(RotationLeft, sk, 5, rotkeys)

		evaluator = ckkscontext.NewEvaluator()

		ciphertext = ckkscontext.NewRandomCiphertext(1, ckkscontext.Levels()-1, ckkscontext.Scale())

		contextQ := ckkscontext.contextQ
		contextP := ckkscontext.contextP

		c2NTT := ciphertext.value[1]
		c2InvNTT := contextQ.NewPoly()
		contextQ.InvNTTLvl(ciphertext.Level(), c2NTT, c2InvNTT)

		c2_qiQDecomp := make([]*ring.Poly, ckkscontext.beta)
		c2_qiPDecomp := make([]*ring.Poly, ckkscontext.beta)

		for i := uint64(0); i < ckkscontext.beta; i++ {
			c2_qiQDecomp[i] = contextQ.NewPoly()
			c2_qiPDecomp[i] = contextP.NewPoly()
		}

		b.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/logSlots=%d/DecomposeNTT", logN, ckkscontext.LogQ(), ckkscontext.Levels(), logSlots), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				for j := uint64(0); j < ckkscontext.beta; j++ {
					evaluator.decomposeAndSplitNTT(ciphertext.Level(), j, c2NTT, c2InvNTT, c2_qiQDecomp[j], c2_qiPDecomp[j])
				}
			}
		})

		b.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/logSlots=%d/RotateHoisted", logN, ckkscontext.LogQ(), ckkscontext.Levels(), logSlots), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				evaluator.switchKeyHoisted(ciphertext, c2_qiQDecomp, c2_qiPDecomp, 5, rotkeys, ciphertext)
			}
		})

		b.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/logSlots=%d/RotateNormal", logN, ckkscontext.LogQ(), ckkscontext.Levels(), logSlots), func(b *testing.B) {

			ciphertext = ckkscontext.NewRandomCiphertext(1, ckkscontext.Levels()-1, ckkscontext.Scale())

			for i := 0; i < b.N; i++ {
				evaluator.RotateColumns(ciphertext, 5, rotkeys, ciphertext)
			}
		})
	}
}

func BenchmarkBootstrapp(b *testing.B) {

	var err error
	var bootcontext *BootContext
	var ckkscontext *CkksContext
	var kgen *KeyGenerator
	var sk *SecretKey
	var evaluator *Evaluator
	var ciphertext *Ciphertext

	params := []benchParams{
		{params: &Parameters{16, []uint8{55, 40, 40, 40, 40, 40, 40, 40, 40, 45, 45, 45, 55, 55, 55, 55, 55, 55, 55, 55, 55, 45, 45, 45, 45}, []uint8{55, 55, 55, 55}, 1 << 40, 3.2}},
	}

	var logN, logSlots, levels, ctsDepth, stcDepth uint64

	for _, param := range params {

		logN = uint64(param.params.LogN)

		logSlots = 14
		ctsDepth = 4
		stcDepth = 3

		levels = uint64(len(param.params.Modulichain))

		ckkscontext = NewCkksContext(param.params)

		kgen = ckkscontext.NewKeyGenerator()

		sk = kgen.NewSecretKey()

		evaluator = ckkscontext.NewEvaluator()

		if bootcontext, err = ckkscontext.NewBootContext(1<<logSlots, sk, ctsDepth, stcDepth); err != nil {
			b.Error()
		}

		// Coeffs To Slots
		var ct0, ct1 *Ciphertext
		b.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/logSlots=%d/CoeffsToSlots", logN, ckkscontext.LogQ(), ckkscontext.Levels(), logSlots), func(b *testing.B) {

			for i := 0; i < b.N; i++ {
				b.StopTimer()
				ciphertext = ckkscontext.NewRandomCiphertext(1, ckkscontext.Levels()-1, ckkscontext.Scale())
				b.StartTimer()

				ct0, ct1 = bootcontext.coeffsToSlots(evaluator, ciphertext)
			}
		})

		// Sine evaluation
		var ct2, ct3 *Ciphertext
		b.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/logSlots=%d/EvaluateSine", logN, ckkscontext.LogQ(), levels, logSlots), func(b *testing.B) {

			for i := 0; i < b.N; i++ {

				b.StopTimer()
				ct0 = ckkscontext.NewRandomCiphertext(1, ct0.Level(), ckkscontext.Scale())
				if logSlots == logN-1 {
					ct1 = ckkscontext.NewRandomCiphertext(1, ct1.Level(), ckkscontext.Scale())
				} else {
					ct1 = nil
				}
				b.StartTimer()

				ct2, ct3 = bootcontext.evaluateSine(ct0, ct1, evaluator)
			}
		})

		// Slots To Coeffs
		b.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/logSlots=%d/SlotsToCoeffs", logN, ckkscontext.LogQ(), levels, logSlots), func(b *testing.B) {

			for i := 0; i < b.N; i++ {

				b.StopTimer()
				ct0 = ckkscontext.NewRandomCiphertext(1, ct2.Level(), ckkscontext.Scale())
				if logSlots == logN-1 {
					ct1 = ckkscontext.NewRandomCiphertext(1, ct3.Level(), ckkscontext.Scale())
				} else {
					ct1 = nil
				}
				b.StartTimer()

				bootcontext.slotsToCoeffs(evaluator, ct2, ct3)
			}
		})
	}
}
