package ckks

import (
	"testing"
)

func BenchmarkBootstrappSine(b *testing.B) {

	var bootcontext *BootContext
	var kgen KeyGenerator
	var sk *SecretKey
	var ciphertext *Ciphertext

	var LTScale float64

	LTScale = 1 << 45
	//SineScale = 1 << 55

	bootparams := BootstrappParams[4]

	parameters := &bootparams.Parameters

	bootparams.Gen()

	ctsDepth := bootparams.CtSDepth
	sinDepth := bootparams.SinDepth

	testString("Params/", parameters)

	kgen = NewKeyGenerator(parameters)

	sk = kgen.GenSecretKey()

	bootcontext = NewBootContext(bootparams, sk)

	b.Run(testString("ModUp/", parameters), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			b.StopTimer()
			ciphertext = NewCiphertextRandom(parameters, 1, 0, LTScale)
			b.StartTimer()

			ciphertext = bootcontext.modUp(ciphertext)
		}
	})

	b.Run(testString("SubSum/", parameters), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			ciphertext = bootcontext.subSum(ciphertext)
		}
	})

	// Coeffs To Slots
	var ct0, ct1 *Ciphertext
	b.Run(testString("CoeffsToSlots/", parameters), func(b *testing.B) {

		for i := 0; i < b.N; i++ {

			b.StopTimer()
			ciphertext = NewCiphertextRandom(parameters, 1, parameters.MaxLevel, LTScale)
			b.StartTimer()

			ct0, ct1 = bootcontext.coeffsToSlots(ciphertext)
		}
	})

	// Sine evaluation
	var ct2, ct3 *Ciphertext
	b.Run(testString("EvalSine/", parameters), func(b *testing.B) {

		for i := 0; i < b.N; i++ {

			b.StopTimer()
			ct0 = NewCiphertextRandom(parameters, 1, parameters.MaxLevel-ctsDepth, LTScale)
			if parameters.LogSlots == parameters.LogN-1 {
				ct1 = NewCiphertextRandom(parameters, 1, parameters.MaxLevel-ctsDepth, LTScale)
			} else {
				ct1 = nil
			}
			b.StartTimer()

			ct2, ct3 = bootcontext.evaluateSine(ct0, ct1)

			if ct2.Level() != parameters.MaxLevel-ctsDepth-sinDepth {
				panic("scaling error during eval sinebetter bench")
			}

			if ct3 != nil {
				if ct3.Level() != parameters.MaxLevel-ctsDepth-sinDepth {
					panic("scaling error during eval sinebetter bench")
				}
			}
		}
	})

	// Slots To Coeffs
	b.Run(testString("SlotsToCoeffs/", parameters), func(b *testing.B) {

		for i := 0; i < b.N; i++ {

			b.StopTimer()
			ct2 = NewCiphertextRandom(parameters, 1, parameters.MaxLevel-ctsDepth-sinDepth, LTScale)
			if parameters.LogSlots == parameters.LogN-1 {
				ct3 = NewCiphertextRandom(parameters, 1, parameters.MaxLevel-ctsDepth-sinDepth, LTScale)
			} else {
				ct3 = nil
			}
			b.StartTimer()

			bootcontext.slotsToCoeffs(ct2, ct3)

		}
	})

	b.Run(testString("Bootstrapp/", parameters), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			b.StopTimer()
			ct := NewCiphertext(parameters, 1, 0, parameters.Scale)
			b.StartTimer()

			bootcontext.Bootstrapp(ct)
		}
	})
}

/*
func BenchmarkBootstrappMultiplications(b *testing.B) {

	var kgen KeyGenerator
	var sk *SecretKey
	var rlk *EvaluationKey
	var eval Evaluator

	var DefaultScale, LTScale float64

	DefaultScale = 1 << 40
	LTScale = 1 << 45
	//SineScale = 1 << 55

	bootParams := new(Parameters)
	bootParams.LogN = 16
	bootParams.LogSlots = 10
	bootParams.Scale = DefaultScale
	bootParams.LogQi = []uint64{55, 40, 40, 40, 40, 40, 40, 40, 40, 40, 45, 45, 45, 55, 55, 55, 55, 55, 55, 55, 55, 55, 45, 45, 45}
	bootParams.LogPi = []uint64{55, 55, 55, 55}
	bootParams.Sigma = 3.2

	bootParams.GenFromLogModuli()

	var ctsDepth uint64

	ctsDepth = 3

	kgen = NewKeyGenerator(bootParams)

	sk = kgen.GenSecretKey()

	rlk = kgen.GenRelinKey(sk)

	eval = NewEvaluator(bootParams)


	ct0 := NewCiphertextRandom(bootParams, 1, bootParams.MaxLevel - ctsDepth, LTScale)
	ct1 := NewCiphertextRandom(bootParams, 1, bootParams.MaxLevel - ctsDepth, LTScale)




	b.Run(testString("Mul20/", bootParams), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			eval.MulRelinNew(ct0, ct1, rlk)
		}
	})

	eval.DropLevel(ct0, 1)
	eval.DropLevel(ct1, 1)

	b.Run(testString("Mul19/", bootParams), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			eval.MulRelinNew(ct0, ct1, rlk)
		}
	})

	eval.DropLevel(ct0, 1)
	eval.DropLevel(ct1, 1)

	b.Run(testString("Mul18/", bootParams), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			eval.MulRelinNew(ct0, ct1, rlk)
		}
	})

	eval.DropLevel(ct0, 1)
	eval.DropLevel(ct1, 1)

	b.Run(testString("Mul17/", bootParams), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			eval.MulRelinNew(ct0, ct1, rlk)
		}
	})

	eval.DropLevel(ct0, 1)
	eval.DropLevel(ct1, 1)

	b.Run(testString("Mul16/", bootParams), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			eval.MulRelinNew(ct0, ct1, rlk)
		}
	})

	eval.DropLevel(ct0, 1)
	eval.DropLevel(ct1, 1)

	b.Run(testString("Mul15/", bootParams), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			eval.MulRelinNew(ct0, ct1, rlk)
		}
	})

	eval.DropLevel(ct0, 1)
	eval.DropLevel(ct1, 1)

	b.Run(testString("Mul14/", bootParams), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			eval.MulRelinNew(ct0, ct1, rlk)
		}
	})

	eval.DropLevel(ct0, 1)
	eval.DropLevel(ct1, 1)

	b.Run(testString("Mul13/", bootParams), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			eval.MulRelinNew(ct0, ct1, rlk)
		}
	})

	eval.DropLevel(ct0, 1)
	eval.DropLevel(ct1, 1)

	b.Run(testString("Mul12/", bootParams), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			eval.MulRelinNew(ct0, ct1, rlk)
		}
	})
}
*/
