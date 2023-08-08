package float

import (
	"math"
	"testing"

	"github.com/tuneinsight/lattigo/v4/ckks"
	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/rlwe"

	"github.com/stretchr/testify/require"
)

func TestPieceWiseFunction(t *testing.T) {

	paramsLiteral := testPrec45

	for _, ringType := range []ring.Type{ring.Standard, ring.ConjugateInvariant} {

		paramsLiteral.RingType = ringType

		if testing.Short() {
			paramsLiteral.LogN = 10
		}

		params, err := ckks.NewParametersFromLiteral(paramsLiteral)
		require.NoError(t, err)

		var tc *ckksTestContext
		if tc, err = genCKKSTestParams(params); err != nil {
			t.Fatal(err)
		}

		enc := tc.encryptorSk
		sk := tc.sk
		ecd := tc.encoder
		dec := tc.decryptor
		kgen := tc.kgen

		btp := ckks.NewSecretKeyBootstrapper(params, sk)

		var galKeys []*rlwe.GaloisKey
		if params.RingType() == ring.Standard {
			galKeys = append(galKeys, kgen.GenGaloisKeyNew(params.GaloisElementForComplexConjugation(), sk))
		}

		PWFEval := NewPieceWiseFunctionEvaluator(params, tc.evaluator.WithKey(rlwe.NewMemEvaluationKeySet(kgen.GenRelinearizationKeyNew(sk), galKeys...)))

		prec := 30
		threshold := math.Exp2(-float64(prec))
		t.Run(GetTestName(params, "Sign"), func(t *testing.T) {

			values, _, ct := newCKKSTestVectors(tc, enc, complex(-1, 0), complex(1, 0), t)

			ct, err = PWFEval.EvaluateSign(ct, prec, btp)
			require.NoError(t, err)

			have := make([]complex128, params.MaxSlots())

			require.NoError(t, ecd.Decode(dec.DecryptNew(ct), have))

			want := make([]complex128, params.MaxSlots())
			for i := range have {

				vc128 := real(values[i].Complex128())

				if math.Abs(vc128) < threshold {
					want[i] = have[i] // Ignores values outside of the interval
					t.Log(vc128, have[i])
				} else {

					if vc128 < 0 {
						want[i] = -1
					} else {
						want[i] = 1
					}
				}
			}

			stats := ckks.GetPrecisionStats(params, ecd, nil, want, have, nil, false)

			if *printPrecisionStats {
				t.Log(stats.String())
			}

			rf64, _ := stats.MeanPrecision.Real.Float64()
			if64, _ := stats.MeanPrecision.Imag.Float64()

			require.Greater(t, rf64, 25.0)
			require.Greater(t, if64, 25.0)
		})

		t.Run(GetTestName(params, "Step"), func(t *testing.T) {
			values, _, ct := newCKKSTestVectors(tc, enc, complex(0.5, 0), complex(1, 0), t)

			ct, err = PWFEval.EvaluateStep(ct, 30, btp)
			require.NoError(t, err)

			have := make([]complex128, params.MaxSlots())

			require.NoError(t, ecd.Decode(dec.DecryptNew(ct), have))

			want := make([]complex128, params.MaxSlots())

			for i := range have {

				vc128 := real(values[i].Complex128())

				if math.Abs(vc128) < threshold {
					want[i] = have[i] // Ignores values outside of the interval
				} else {

					if vc128 < 0.5 {
						want[i] = 0
					} else {
						want[i] = 1
					}
				}
			}

			stats := ckks.GetPrecisionStats(params, ecd, nil, want, have, nil, false)

			if *printPrecisionStats {
				t.Log(stats.String())
			}

			rf64, _ := stats.MeanPrecision.Real.Float64()
			if64, _ := stats.MeanPrecision.Imag.Float64()

			require.Greater(t, rf64, 25.0)
			require.Greater(t, if64, 25.0)
		})
	}
}
