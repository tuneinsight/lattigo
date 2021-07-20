package ckks

import (
	"math"
	"testing"
	"time"
)

func BenchmarkBootstrapp(b *testing.B) {

	if !*testBootstrapping {
		b.Skip("skipping bootstrapping test")
	}

	var err error
	var testContext = new(testParams)
	var btp *Bootstrapper

	paramSet := 2

	btpParams := DefaultBootstrapParams[paramSet]

	params, err := btpParams.Params()
	if err != nil {
		panic(err)
	}
	if testContext, err = genTestParams(params, btpParams.H); err != nil {
		panic(err)
	}

	rotations := btpParams.RotationsForBootstrapping(testContext.params.LogSlots())

	rotkeys := testContext.kgen.GenRotationKeysForRotations(rotations, true, testContext.sk)

	btpKey := BootstrappingKey{testContext.rlk, rotkeys}

	if btp, err = NewBootstrapper(testContext.params, *btpParams, btpKey); err != nil {
		panic(err)
	}

	b.Run(testString(testContext, "Bootstrapp/"), func(b *testing.B) {
		for i := 0; i < b.N; i++ {

			bootstrappingScale := math.Exp2(math.Round(math.Log2(btp.params.QiFloat64(0) / btp.MessageRatio)))

			b.StopTimer()
			ct := NewCiphertextRandom(testContext.prng, testContext.params, 1, 0, bootstrappingScale)
			b.StartTimer()

			var t time.Time
			var ct0, ct1 *Ciphertext

			// ModUp ct_{Q_0} -> ct_{Q_L}
			t = time.Now()
			ct = btp.modUpFromQ0(ct)
			b.Log("After ModUp  :", time.Since(t), ct.Level(), ct.Scale)

			//SubSum X -> (N/dslots) * Y^dslots
			t = time.Now()
			ct = btp.evaluator.Trace(ct, btp.params.LogSlots())
			b.Log("After SubSum :", time.Since(t), ct.Level(), ct.Scale)

			// Part 1 : Coeffs to slots
			t = time.Now()
			ct0, ct1 = btp.evaluator.CoeffsToSlots(ct, btp.ctsMatrices)
			b.Log("After CtS    :", time.Since(t), ct0.Level(), ct0.Scale)

			// Part 2 : SineEval
			t = time.Now()
			ct0 = btp.evaluator.EvalMod(ct0, btp.evalModPoly)
			ct0.Scale = btp.params.Scale()

			if ct1 != nil {
				ct1 = btp.evaluator.EvalMod(ct1, btp.evalModPoly)
				ct1.Scale = btp.params.Scale()
			}
			b.Log("After Sine   :", time.Since(t), ct0.Level(), ct0.Scale)

			// Part 3 : Slots to coeffs
			t = time.Now()
			ct0 = btp.evaluator.SlotsToCoeffs(ct0, ct1, btp.stcMatrices)
			ct0.Scale = math.Exp2(math.Round(math.Log2(ct0.Scale)))
			b.Log("After StC    :", time.Since(t), ct0.Level(), ct0.Scale)
		}
	})
}
