package ckks

import (
	"github.com/ldsec/lattigo/utils"
	"math"
	"testing"
)

func BenchmarkBootstrapp(b *testing.B) {

	var btp *Bootstrapper
	var ciphertext *Ciphertext

	var LTScale float64

	LTScale = 1 << 45
	//SineScale = 1 << 55

	paramSet := uint64(4)

	btpParams := DefaultBootstrappParams[paramSet]
	if err := genTestParams(DefaultBootstrappSchemeParams[paramSet], btpParams.H); err != nil {
		panic(err)
	}

	prng, err := utils.NewPRNG()
	if err != nil {
		panic(err)
	}

	ctsDepth := uint64(len(btpParams.CtSLevel))
	sinDepth := uint64(math.Ceil(math.Log2(float64(btpParams.SinDeg))) + float64(btpParams.SinRescal))

	testString("Params/")

	btp, err = NewBootstrapper(params.params, btpParams)
	if err != nil {
		panic(err)
	}
	btp.GenKeys(params.sk)

	b.Run(testString("ModUp/"), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			b.StopTimer()
			ciphertext = NewCiphertextRandom(prng, params.params, 1, 0, LTScale)
			b.StartTimer()

			ciphertext = btp.modUp(ciphertext)
		}
	})

	b.Run(testString("SubSum/"), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			ciphertext = btp.subSum(ciphertext)
		}
	})

	// Coeffs To Slots
	var ct0, ct1 *Ciphertext
	b.Run(testString("CoeffsToSlots/"), func(b *testing.B) {

		for i := 0; i < b.N; i++ {

			b.StopTimer()
			ciphertext = NewCiphertextRandom(prng, params.params, 1, params.params.MaxLevel(), LTScale)
			b.StartTimer()

			ct0, ct1 = btp.coeffsToSlots(ciphertext)
		}
	})

	// Sine evaluation
	var ct2, ct3 *Ciphertext
	b.Run(testString("EvalSine/"), func(b *testing.B) {

		for i := 0; i < b.N; i++ {

			b.StopTimer()
			ct0 = NewCiphertextRandom(prng, params.params, 1, params.params.MaxLevel()-ctsDepth, LTScale)
			if params.params.logSlots == params.params.LogMaxSlots() {
				ct1 = NewCiphertextRandom(prng, params.params, 1, params.params.MaxLevel()-ctsDepth, LTScale)
			} else {
				ct1 = nil
			}
			b.StartTimer()

			ct2, ct3 = btp.evaluateSine(ct0, ct1)

			if ct2.Level() != params.params.MaxLevel()-ctsDepth-sinDepth {
				panic("scaling error during eval sinebetter bench")
			}

			if ct3 != nil {
				if ct3.Level() != params.params.MaxLevel()-ctsDepth-sinDepth {
					panic("scaling error during eval sinebetter bench")
				}
			}
		}
	})

	// Slots To Coeffs
	b.Run(testString("SlotsToCoeffs/"), func(b *testing.B) {

		for i := 0; i < b.N; i++ {

			b.StopTimer()
			ct2 = NewCiphertextRandom(prng, params.params, 1, params.params.MaxLevel()-ctsDepth-sinDepth, LTScale)
			if params.params.logSlots == params.params.LogMaxSlots() {
				ct3 = NewCiphertextRandom(prng, params.params, 1, params.params.MaxLevel()-ctsDepth-sinDepth, LTScale)
			} else {
				ct3 = nil
			}
			b.StartTimer()

			btp.slotsToCoeffs(ct2, ct3)

		}
	})

	b.Run(testString("Bootstrapp/"), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			b.StopTimer()
			ct := NewCiphertext(params.params, 1, 0, params.params.scale)
			b.StartTimer()

			btp.Bootstrapp(ct)
		}
	})
}
