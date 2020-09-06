package ckks

import (
	"github.com/ldsec/lattigo/utils"
	"testing"
)

func BenchmarkBootstrapp(b *testing.B) {

	var bootcontext *BootContext
	var kgen KeyGenerator
	var sk *SecretKey
	var ciphertext *Ciphertext

	var LTScale float64

	LTScale = 1 << 45
	//SineScale = 1 << 55

	bootparams := BootstrappParams[0]

	parameters := &bootparams.Parameters

	bootparams.Gen()

	prng, err := utils.NewPRNG()
	if err != nil {
		panic(err)
	}

	ctsDepth := uint64(len(bootparams.CtSLevel))
	sinDepth := bootparams.SinDepth

	testString("Params/", parameters)

	kgen = NewKeyGenerator(parameters)

	sk = kgen.GenSecretKey()

	bootcontext = NewBootContext(bootparams)
	bootcontext.GenBootKeys(sk)

	b.Run(testString("ModUp/", parameters), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			b.StopTimer()
			ciphertext = NewCiphertextRandom(prng, parameters, 1, 0, LTScale)
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
			ciphertext = NewCiphertextRandom(prng, parameters, 1, parameters.MaxLevel(), LTScale)
			b.StartTimer()

			ct0, ct1 = bootcontext.coeffsToSlots(ciphertext)
		}
	})

	// Sine evaluation
	var ct2, ct3 *Ciphertext
	b.Run(testString("EvalSine/", parameters), func(b *testing.B) {

		for i := 0; i < b.N; i++ {

			b.StopTimer()
			ct0 = NewCiphertextRandom(prng, parameters, 1, parameters.MaxLevel()-ctsDepth, LTScale)
			if parameters.logSlots == parameters.logN-1 {
				ct1 = NewCiphertextRandom(prng, parameters, 1, parameters.MaxLevel()-ctsDepth, LTScale)
			} else {
				ct1 = nil
			}
			b.StartTimer()

			ct2, ct3 = bootcontext.evaluateSine(ct0, ct1)

			if ct2.Level() != parameters.MaxLevel()-ctsDepth-sinDepth {
				panic("scaling error during eval sinebetter bench")
			}

			if ct3 != nil {
				if ct3.Level() != parameters.MaxLevel()-ctsDepth-sinDepth {
					panic("scaling error during eval sinebetter bench")
				}
			}
		}
	})

	// Slots To Coeffs
	b.Run(testString("SlotsToCoeffs/", parameters), func(b *testing.B) {

		for i := 0; i < b.N; i++ {

			b.StopTimer()
			ct2 = NewCiphertextRandom(prng, parameters, 1, parameters.MaxLevel()-ctsDepth-sinDepth, LTScale)
			if parameters.logSlots == parameters.logN-1 {
				ct3 = NewCiphertextRandom(prng, parameters, 1, parameters.MaxLevel()-ctsDepth-sinDepth, LTScale)
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
			ct := NewCiphertext(parameters, 1, 0, parameters.scale)
			b.StartTimer()

			bootcontext.Bootstrapp(ct)
		}
	})
}
