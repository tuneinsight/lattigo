package ckks

import (
	"math"
	"testing"

	"github.com/ldsec/lattigo/utils"
)

func BenchmarkBootstrapp(b *testing.B) {
	var err error
	var testContext = new(testParams)
	var btp *Bootstrapper
	var ciphertext *Ciphertext

	var LTScale float64

	LTScale = 1 << 45
	//SineScale = 1 << 55

	paramSet := uint64(4)

	btpParams := DefaultBootstrappParams[paramSet]
	if testContext, err = genTestParams(DefaultBootstrappSchemeParams[paramSet], btpParams.H); err != nil {
		panic(err)
	}

	prng, err := utils.NewPRNG()
	if err != nil {
		panic(err)
	}

	ctsDepth := uint64(len(btpParams.CtSLevel))
	sinDepth := uint64(math.Ceil(math.Log2(float64(btpParams.SinDeg))) + float64(btpParams.SinRescal))

	btpKey := testContext.kgen.GenBootstrappingKey(testContext.params.logSlots, btpParams, testContext.sk)
	btp, err = NewBootstrapper(testContext.params, btpParams, btpKey)
	if err != nil {
		panic(err)
	}

	b.Run(testString(testContext, "ModUp/"), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			b.StopTimer()
			ciphertext = NewCiphertextRandom(prng, testContext.params, 1, 0, LTScale)
			b.StartTimer()

			ciphertext = btp.modUp(ciphertext)
		}
	})

	b.Run(testString(testContext, "SubSum/"), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			ciphertext = btp.subSum(ciphertext)
		}
	})

	// Coeffs To Slots
	var ct0, ct1 *Ciphertext
	b.Run(testString(testContext, "CoeffsToSlots/"), func(b *testing.B) {

		for i := 0; i < b.N; i++ {

			b.StopTimer()
			ciphertext = NewCiphertextRandom(prng, testContext.params, 1, testContext.params.MaxLevel(), LTScale)
			b.StartTimer()

			ct0, ct1 = btp.coeffsToSlots(ciphertext)
		}
	})

	// Sine evaluation
	var ct2, ct3 *Ciphertext
	b.Run(testString(testContext, "EvalSine/"), func(b *testing.B) {

		for i := 0; i < b.N; i++ {

			b.StopTimer()
			ct0 = NewCiphertextRandom(prng, testContext.params, 1, testContext.params.MaxLevel()-ctsDepth, LTScale)
			if testContext.params.logSlots == testContext.params.MaxLogSlots() {
				ct1 = NewCiphertextRandom(prng, testContext.params, 1, testContext.params.MaxLevel()-ctsDepth, LTScale)
			} else {
				ct1 = nil
			}
			b.StartTimer()

			ct2, ct3 = btp.evaluateSine(ct0, ct1)

			if ct2.Level() != testContext.params.MaxLevel()-ctsDepth-sinDepth {
				panic("scaling error during eval sinebetter bench")
			}

			if ct3 != nil {
				if ct3.Level() != testContext.params.MaxLevel()-ctsDepth-sinDepth {
					panic("scaling error during eval sinebetter bench")
				}
			}
		}
	})

	// Slots To Coeffs
	b.Run(testString(testContext, "SlotsToCoeffs/"), func(b *testing.B) {

		for i := 0; i < b.N; i++ {

			b.StopTimer()
			ct2 = NewCiphertextRandom(prng, testContext.params, 1, testContext.params.MaxLevel()-ctsDepth-sinDepth, LTScale)
			if testContext.params.logSlots == testContext.params.MaxLogSlots() {
				ct3 = NewCiphertextRandom(prng, testContext.params, 1, testContext.params.MaxLevel()-ctsDepth-sinDepth, LTScale)
			} else {
				ct3 = nil
			}
			b.StartTimer()

			btp.slotsToCoeffs(ct2, ct3)

		}
	})

	b.Run(testString(testContext, "Bootstrapp/"), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			b.StopTimer()
			ct := NewCiphertext(testContext.params, 1, 0, testContext.params.scale)
			b.StartTimer()

			btp.Bootstrapp(ct)
		}
	})
}
