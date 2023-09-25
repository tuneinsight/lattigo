package float_test

import (
	"math/big"
	"testing"

	"github.com/tuneinsight/lattigo/v4/circuits/float"
	"github.com/tuneinsight/lattigo/v4/circuits/float/bootstrapper"
	"github.com/tuneinsight/lattigo/v4/ckks"
	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/rlwe"

	"github.com/stretchr/testify/require"
)

func TestComparisons(t *testing.T) {

	paramsLiteral := float.TestPrec90

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

		btp := bootstrapper.NewSecretKeyBootstrapper(params, sk)

		var galKeys []*rlwe.GaloisKey
		if params.RingType() == ring.Standard {
			galKeys = append(galKeys, kgen.GenGaloisKeyNew(params.GaloisElementForComplexConjugation(), sk))
		}

		eval := tc.evaluator.WithKey(rlwe.NewMemEvaluationKeySet(kgen.GenRelinearizationKeyNew(sk), galKeys...))

		MCPEval := float.NewMinimaxCompositePolynomialEvaluator(params, eval, btp)

		polys := float.NewMinimaxCompositePolynomial(float.DefaultMinimaxCompositePolynomialForSign)

		CmpEval := float.NewComparisonEvaluator(MCPEval, polys)

		t.Run(GetTestName(params, "Sign"), func(t *testing.T) {

			values, _, ct := newCKKSTestVectors(tc, enc, complex(-1, 0), complex(1, 0), t)

			var sign *rlwe.Ciphertext
			sign, err = CmpEval.Sign(ct)
			require.NoError(t, err)

			have := make([]*big.Float, params.MaxSlots())

			require.NoError(t, ecd.Decode(dec.DecryptNew(sign), have))

			want := make([]*big.Float, params.MaxSlots())

			for i := range have {
				want[i] = polys.Evaluate(values[i])[0]
			}

			ckks.VerifyTestVectors(params, ecd, nil, want, have, params.LogDefaultScale(), nil, *printPrecisionStats, t)
		})

		t.Run(GetTestName(params, "Step"), func(t *testing.T) {

			values, _, ct := newCKKSTestVectors(tc, enc, complex(-1, 0), complex(1, 0), t)

			var step *rlwe.Ciphertext
			step, err = CmpEval.Step(ct)
			require.NoError(t, err)

			have := make([]*big.Float, params.MaxSlots())

			require.NoError(t, ecd.Decode(dec.DecryptNew(step), have))

			want := make([]*big.Float, params.MaxSlots())

			half := new(big.Float).SetFloat64(0.5)

			for i := range have {
				want[i] = polys.Evaluate(values[i])[0]
				want[i].Mul(want[i], half)
				want[i].Add(want[i], half)
			}

			ckks.VerifyTestVectors(params, ecd, nil, want, have, params.LogDefaultScale(), nil, *printPrecisionStats, t)
		})

		t.Run(GetTestName(params, "Max"), func(t *testing.T) {

			values0, _, ct0 := newCKKSTestVectors(tc, enc, complex(-0.5, 0), complex(0.5, 0), t)
			values1, _, ct1 := newCKKSTestVectors(tc, enc, complex(-0.5, 0), complex(0.5, 0), t)

			var max *rlwe.Ciphertext
			max, err = CmpEval.Max(ct0, ct1)
			require.NoError(t, err)

			have := make([]*big.Float, params.MaxSlots())

			require.NoError(t, ecd.Decode(dec.DecryptNew(max), have))

			want := make([]*big.Float, params.MaxSlots())

			for i := range have {

				if values0[i][0].Cmp(values1[i][0]) == -1 {
					want[i] = values1[i][0]
				} else {
					want[i] = values0[i][0]
				}
			}

			ckks.VerifyTestVectors(params, ecd, nil, want, have, params.LogDefaultScale(), nil, *printPrecisionStats, t)
		})

		t.Run(GetTestName(params, "Min"), func(t *testing.T) {

			values0, _, ct0 := newCKKSTestVectors(tc, enc, complex(-0.5, 0), complex(0.5, 0), t)
			values1, _, ct1 := newCKKSTestVectors(tc, enc, complex(-0.5, 0), complex(0.5, 0), t)

			var max *rlwe.Ciphertext
			max, err = CmpEval.Min(ct0, ct1)
			require.NoError(t, err)

			have := make([]*big.Float, params.MaxSlots())

			require.NoError(t, ecd.Decode(dec.DecryptNew(max), have))

			want := make([]*big.Float, params.MaxSlots())

			for i := range have {

				if values0[i][0].Cmp(values1[i][0]) == 1 {
					want[i] = values1[i][0]
				} else {
					want[i] = values0[i][0]
				}
			}

			ckks.VerifyTestVectors(params, ecd, nil, want, have, params.LogDefaultScale(), nil, *printPrecisionStats, t)
		})
	}
}
