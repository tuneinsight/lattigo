package bootstrapping

import (
	"testing"
	"time"

	"github.com/stretchr/testify/require"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
)

func BenchmarkBootstrap(b *testing.B) {

	paramSet := DefaultParametersDense[0]

	params, err := ckks.NewParametersFromLiteral(paramSet.SchemeParams)
	require.NoError(b, err)

	btpParams, err := NewParametersFromLiteral(params, paramSet.BootstrappingParams)
	require.Nil(b, err)

	kgen := rlwe.NewKeyGenerator(params)
	sk := kgen.GenSecretKeyNew()

	evk, _, err := btpParams.GenEvaluationKeys(sk)
	require.NoError(b, err)

	eval, err := NewEvaluator(btpParams, evk)
	require.NoError(b, err)

	b.Run(ParamsToString(params, btpParams.LogMaxDimensions().Cols, "Bootstrap/"), func(b *testing.B) {

		var err error

		for i := 0; i < b.N; i++ {

			b.StopTimer()
			ct := ckks.NewCiphertext(params, 1, 0)
			b.StartTimer()

			var t time.Time
			var ct0, ct1 *rlwe.Ciphertext

			// ScaleDown
			t = time.Now()
			ct, _, err = eval.ScaleDown(ct)
			require.NoError(b, err)
			b.Log("ScaleDown:", time.Since(t), ct.Level(), ct.Scale.Float64())

			// ModUp ct_{Q_0} -> ct_{Q_L}
			t = time.Now()
			ct, err = eval.ModUp(ct)
			require.NoError(b, err)
			b.Log("ModUp    :", time.Since(t), ct.Level(), ct.Scale.Float64())

			// Part 1 : Coeffs to slots
			t = time.Now()
			ct0, ct1, err = eval.CoeffsToSlots(ct)
			require.NoError(b, err)
			b.Log("CtS      :", time.Since(t), ct0.Level(), ct0.Scale.Float64())

			// Part 2 : SineEval
			t = time.Now()
			ct0, err = eval.EvalMod(ct0)
			require.NoError(b, err)
			if ct1 != nil {
				ct1, err = eval.EvalMod(ct1)
				require.NoError(b, err)
			}
			b.Log("EvalMod  :", time.Since(t), ct0.Level(), ct0.Scale.Float64())

			// Part 3 : Slots to coeffs
			t = time.Now()
			ct0, err = eval.SlotsToCoeffs(ct0, ct1)
			require.NoError(b, err)
			b.Log("StC      :", time.Since(t), ct0.Level(), ct0.Scale.Float64())
		}
	})
}
