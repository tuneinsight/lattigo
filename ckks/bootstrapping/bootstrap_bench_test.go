package bootstrapping

import (
	"math"
	"testing"
	"time"

	"github.com/stretchr/testify/require"
	"github.com/tuneinsight/lattigo/v3/ckks"
)

func BenchmarkBootstrapp(b *testing.B) {

	var err error
	var btp *Bootstrapper

	paramSet := DefaultParametersSparse[0]

	ckksParamsLit, btpParams, err := NewParametersFromLiteral(paramSet.SchemeParams, paramSet.BootstrappingParams)
	require.Nil(b, err)

	params, err := ckks.NewParametersFromLiteral(ckksParamsLit)
	if err != nil {
		panic(err)
	}

	kgen := ckks.NewKeyGenerator(params)
	sk := kgen.GenSecretKey()

	evk := GenEvaluationKeys(btpParams, params, sk)

	if btp, err = NewBootstrapper(params, btpParams, evk); err != nil {
		panic(err)
	}

	b.Run(ParamsToString(params, "Bootstrapp/"), func(b *testing.B) {
		for i := 0; i < b.N; i++ {

			bootstrappingScale := math.Exp2(math.Round(math.Log2(btp.params.QiFloat64(0) / btp.evalModPoly.MessageRatio())))

			b.StopTimer()
			ct := ckks.NewCiphertext(params, 1, 0, bootstrappingScale)
			b.StartTimer()

			var t time.Time
			var ct0, ct1 *ckks.Ciphertext

			// ModUp ct_{Q_0} -> ct_{Q_L}
			t = time.Now()
			ct = btp.modUpFromQ0(ct)
			b.Log("After ModUp  :", time.Since(t), ct.Level(), ct.Scale)

			//SubSum X -> (N/dslots) * Y^dslots
			t = time.Now()
			btp.Trace(ct, btp.params.LogSlots(), ct)
			b.Log("After SubSum :", time.Since(t), ct.Level(), ct.Scale)

			// Part 1 : Coeffs to slots
			t = time.Now()
			ct0, ct1 = btp.CoeffsToSlotsNew(ct, btp.ctsMatrices)
			b.Log("After CtS    :", time.Since(t), ct0.Level(), ct0.Scale)

			// Part 2 : SineEval
			t = time.Now()
			ct0 = btp.EvalModNew(ct0, btp.evalModPoly)
			ct0.Scale = btp.params.DefaultScale()

			if ct1 != nil {
				ct1 = btp.EvalModNew(ct1, btp.evalModPoly)
				ct1.Scale = btp.params.DefaultScale()
			}
			b.Log("After Sine   :", time.Since(t), ct0.Level(), ct0.Scale)

			// Part 3 : Slots to coeffs
			t = time.Now()
			ct0 = btp.SlotsToCoeffsNew(ct0, ct1, btp.stcMatrices)
			ct0.Scale = math.Exp2(math.Round(math.Log2(ct0.Scale)))
			b.Log("After StC    :", time.Since(t), ct0.Level(), ct0.Scale)
		}
	})
}
