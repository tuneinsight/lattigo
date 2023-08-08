package float

import (
	"math"
	"math/big"
	"testing"

	"github.com/tuneinsight/lattigo/v4/ckks"
	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/rlwe"

	"github.com/stretchr/testify/require"
)

func testGoldschmidtDivisionNew(tc *ckksTestContext, t *testing.T) {

	params := tc.params

	t.Run(GetTestName(params, "GoldschmidtDivisionNew"), func(t *testing.T) {

		min := 0.1

		values, _, ciphertext := newCKKSTestVectors(tc, tc.encryptorSk, complex(min, 0), complex(2-min, 0), t)

		one := new(big.Float).SetInt64(1)
		for i := range values {
			values[i][0].Quo(one, values[i][0])
		}

		btp := ckks.NewSecretKeyBootstrapper(params, tc.sk)

		var err error
		if ciphertext, err = NewInverseEvaluator(params, tc.evaluator, nil).GoldschmidtDivisionNew(ciphertext, min, btp); err != nil {
			t.Fatal(err)
		}

		ckks.VerifyTestVectors(params, tc.encoder, tc.decryptor, values, ciphertext, nil, *printPrecisionStats, t)
	})
}

func TestInverse(t *testing.T) {

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

		t.Run(GetTestName(params, "FullDomain"), func(t *testing.T) {

			r := 10

			// 2^{-r}
			min := math.Exp2(-float64(r))

			// 2^{r}
			max := math.Exp2(float64(r))

			require.NoError(t, err)

			var galKeys []*rlwe.GaloisKey
			if params.RingType() == ring.Standard {
				galKeys = append(galKeys, kgen.GenGaloisKeyNew(params.GaloisElementForComplexConjugation(), sk))
			}

			evk := rlwe.NewMemEvaluationKeySet(kgen.GenRelinearizationKeyNew(sk), galKeys...)

			eval := tc.evaluator.WithKey(evk)

			values, _, ct := newCKKSTestVectors(tc, enc, complex(-max, 0), complex(max, 0), t)

			btp := ckks.NewSecretKeyBootstrapper(params, sk)

			invEval := NewInverseEvaluator(params, eval, eval)

			canBeNegative := true

			cInv, err := invEval.EvaluateNew(ct, min, max, canBeNegative, btp)
			require.NoError(t, err)

			have := make([]complex128, params.MaxSlots())

			require.NoError(t, ecd.Decode(dec.DecryptNew(cInv), have))

			want := make([]complex128, params.MaxSlots())

			for i := range have {

				vc128 := values[i].Complex128()

				have[i] *= vc128

				if math.Abs(real(vc128)) < min {
					want[i] = have[i] // Ignores values outside of the interval
				} else {
					want[i] = 1.0
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
