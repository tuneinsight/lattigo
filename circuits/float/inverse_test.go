package float_test

import (
	"math"
	"math/big"
	"testing"

	"github.com/stretchr/testify/require"
	"github.com/tuneinsight/lattigo/v4/circuits/float"
	"github.com/tuneinsight/lattigo/v4/circuits/float/bootstrapper"
	"github.com/tuneinsight/lattigo/v4/ckks"
	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/rlwe"
	"github.com/tuneinsight/lattigo/v4/utils/bignum"
)

func TestInverse(t *testing.T) {

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

		logmin := -30.0
		logmax := 10.0

		// 2^{-r}
		min := math.Exp2(float64(logmin))

		// 2^{r}
		max := math.Exp2(float64(logmax))

		var galKeys []*rlwe.GaloisKey
		if params.RingType() == ring.Standard {
			galKeys = append(galKeys, kgen.GenGaloisKeyNew(params.GaloisElementForComplexConjugation(), sk))
		}

		evk := rlwe.NewMemEvaluationKeySet(kgen.GenRelinearizationKeyNew(sk), galKeys...)

		eval := tc.evaluator.WithKey(evk)

		t.Run(GetTestName(params, "GoldschmidtDivisionNew"), func(t *testing.T) {

			values, _, ciphertext := newCKKSTestVectors(tc, tc.encryptorSk, complex(min, 0), complex(2-min, 0), t)

			one := new(big.Float).SetInt64(1)
			for i := range values {
				values[i][0].Quo(one, values[i][0])
			}

			invEval := float.NewInverseEvaluator(params, eval, btp)

			var err error
			if ciphertext, err = invEval.GoldschmidtDivisionNew(ciphertext, logmin); err != nil {
				t.Fatal(err)
			}

			ckks.VerifyTestVectors(params, tc.encoder, tc.decryptor, values, ciphertext, 70, nil, *printPrecisionStats, t)
		})

		t.Run(GetTestName(params, "PositiveDomain"), func(t *testing.T) {

			values, _, ct := newCKKSTestVectors(tc, enc, complex(0, 0), complex(max, 0), t)

			invEval := float.NewInverseEvaluator(params, eval, btp)

			cInv, err := invEval.EvaluatePositiveDomainNew(ct, logmin, logmax)
			require.NoError(t, err)

			have := make([]*big.Float, params.MaxSlots())

			require.NoError(t, ecd.Decode(dec.DecryptNew(cInv), have))

			want := make([]*big.Float, params.MaxSlots())

			threshold := bignum.NewFloat(min, params.EncodingPrecision())
			for i := range have {
				if new(big.Float).Abs(values[i][0]).Cmp(threshold) == -1 {
					want[i] = have[i] // Ignores values outside of the interval
				} else {
					want[i] = new(big.Float).Quo(bignum.NewFloat(1, params.EncodingPrecision()), values[i][0])
				}
			}

			ckks.VerifyTestVectors(params, tc.encoder, nil, want, have, 70, nil, *printPrecisionStats, t)
		})

		t.Run(GetTestName(params, "NegativeDomain"), func(t *testing.T) {

			values, _, ct := newCKKSTestVectors(tc, enc, complex(-max, 0), complex(0, 0), t)

			invEval := float.NewInverseEvaluator(params, eval, btp)

			cInv, err := invEval.EvaluateNegativeDomainNew(ct, logmin, logmax)
			require.NoError(t, err)

			have := make([]*big.Float, params.MaxSlots())

			require.NoError(t, ecd.Decode(dec.DecryptNew(cInv), have))

			want := make([]*big.Float, params.MaxSlots())

			threshold := bignum.NewFloat(min, params.EncodingPrecision())
			for i := range have {
				if new(big.Float).Abs(values[i][0]).Cmp(threshold) == -1 {
					want[i] = have[i] // Ignores values outside of the interval
				} else {
					want[i] = new(big.Float).Quo(bignum.NewFloat(1, params.EncodingPrecision()), values[i][0])
				}
			}

			ckks.VerifyTestVectors(params, tc.encoder, nil, want, have, 70, nil, *printPrecisionStats, t)
		})

		t.Run(GetTestName(params, "FullDomain"), func(t *testing.T) {

			values, _, ct := newCKKSTestVectors(tc, enc, complex(-max, 0), complex(max, 0), t)

			invEval := float.NewInverseEvaluator(params, eval, btp)

			cInv, err := invEval.EvaluateFullDomainNew(ct, logmin, logmax, float.NewMinimaxCompositePolynomial(float.DefaultMinimaxCompositePolynomialForSign))
			require.NoError(t, err)

			have := make([]*big.Float, params.MaxSlots())

			require.NoError(t, ecd.Decode(dec.DecryptNew(cInv), have))

			want := make([]*big.Float, params.MaxSlots())

			threshold := bignum.NewFloat(min, params.EncodingPrecision())
			for i := range have {
				if new(big.Float).Abs(values[i][0]).Cmp(threshold) == -1 {
					want[i] = have[i] // Ignores values outside of the interval
				} else {
					want[i] = new(big.Float).Quo(bignum.NewFloat(1, params.EncodingPrecision()), values[i][0])
				}
			}

			ckks.VerifyTestVectors(params, tc.encoder, nil, want, have, 70, nil, *printPrecisionStats, t)
		})
	}
}
