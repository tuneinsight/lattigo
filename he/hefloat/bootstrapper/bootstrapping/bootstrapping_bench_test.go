package bootstrapping

import (
	"math"
	"testing"
	"time"

	"github.com/stretchr/testify/require"
	"github.com/tuneinsight/lattigo/v4/core/rlwe"
	"github.com/tuneinsight/lattigo/v4/he/hefloat"
)

func BenchmarkBootstrap(b *testing.B) {

	paramSet := DefaultParametersDense[0]

	params, err := hefloat.NewParametersFromLiteral(paramSet.SchemeParams)
	require.NoError(b, err)

	btpParams, err := NewParametersFromLiteral(params, paramSet.BootstrappingParams)
	require.Nil(b, err)

	kgen := rlwe.NewKeyGenerator(params)
	sk := kgen.GenSecretKeyNew()

	btp, err := NewBootstrapper(btpParams, btpParams.GenEvaluationKeySetNew(sk))
	require.NoError(b, err)

	b.Run(ParamsToString(params, btpParams.LogMaxDimensions().Cols, "Bootstrap/"), func(b *testing.B) {

		var err error

		for i := 0; i < b.N; i++ {

			bootstrappingScale := rlwe.NewScale(math.Exp2(math.Round(math.Log2(float64(btp.params.Q()[0]) / btp.mod1Parameters.MessageRatio()))))

			b.StopTimer()
			ct := hefloat.NewCiphertext(params, 1, 0)
			ct.Scale = bootstrappingScale
			b.StartTimer()

			var t time.Time
			var ct0, ct1 *rlwe.Ciphertext

			// ModUp ct_{Q_0} -> ct_{Q_L}
			t = time.Now()
			ct, err = btp.modUpFromQ0(ct)
			require.NoError(b, err)
			b.Log("After ModUp  :", time.Since(t), ct.Level(), ct.Scale.Float64())

			//SubSum X -> (N/dslots) * Y^dslots
			t = time.Now()
			require.NoError(b, btp.Trace(ct, ct.LogDimensions.Cols, ct))
			b.Log("After SubSum :", time.Since(t), ct.Level(), ct.Scale.Float64())

			// Part 1 : Coeffs to slots
			t = time.Now()
			ct0, ct1, err = btp.CoeffsToSlotsNew(ct, btp.ctsMatrices)
			require.NoError(b, err)
			b.Log("After CtS    :", time.Since(t), ct0.Level(), ct0.Scale.Float64())

			// Part 2 : SineEval
			t = time.Now()
			ct0, err = btp.Mod1Evaluator.EvaluateNew(ct0)
			require.NoError(b, err)
			ct0.Scale = btp.params.DefaultScale()

			if ct1 != nil {
				ct1, err = btp.Mod1Evaluator.EvaluateNew(ct1)
				require.NoError(b, err)
				ct1.Scale = btp.params.DefaultScale()
			}
			b.Log("After Sine   :", time.Since(t), ct0.Level(), ct0.Scale.Float64())

			// Part 3 : Slots to coeffs
			t = time.Now()
			ct0, err = btp.SlotsToCoeffsNew(ct0, ct1, btp.stcMatrices)
			require.NoError(b, err)
			ct0.Scale = rlwe.NewScale(math.Exp2(math.Round(math.Log2(ct0.Scale.Float64()))))
			b.Log("After StC    :", time.Since(t), ct0.Level(), ct0.Scale.Float64())
		}
	})
}
