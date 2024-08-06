package inverse

import (
	"flag"
	"fmt"
	"math"
	"math/big"
	"testing"

	"github.com/stretchr/testify/require"

	"github.com/tuneinsight/lattigo/v6/circuits/ckks/bootstrapping"
	"github.com/tuneinsight/lattigo/v6/circuits/ckks/comparison"
	"github.com/tuneinsight/lattigo/v6/circuits/ckks/minimax"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
	"github.com/tuneinsight/lattigo/v6/utils/bignum"
)

var printPrecisionStats = flag.Bool("print-precision", false, "print precision stats")

func TestInverse(t *testing.T) {

	paramsLiteral := testInsecurePrec90

	for _, ringType := range []ring.Type{ring.Standard, ring.ConjugateInvariant} {

		paramsLiteral.RingType = ringType

		if testing.Short() {
			paramsLiteral.LogN = 10
		}

		tc := ckks.NewTestContext(paramsLiteral)

		btp := bootstrapping.NewSecretKeyBootstrapper(tc.Params, tc.Sk)

		logmin := -30.0
		logmax := 10.0

		// 2^{-r}
		min := math.Exp2(float64(logmin))

		// 2^{r}
		max := math.Exp2(float64(logmax))

		var galKeys []*rlwe.GaloisKey
		if tc.Params.RingType() == ring.Standard {
			galKeys = append(galKeys, tc.Kgen.GenGaloisKeyNew(tc.Params.GaloisElementForComplexConjugation(), tc.Sk))
		}

		evk := rlwe.NewMemEvaluationKeySet(tc.Kgen.GenRelinearizationKeyNew(tc.Sk), galKeys...)

		eval := tc.Evl.WithKey(evk)

		t.Run(name("GoldschmidtDivisionNew", tc), func(t *testing.T) {

			values, _, ciphertext := tc.NewTestVector(complex(min, 0), complex(2-min, 0))

			one := new(big.Float).SetInt64(1)
			for i := range values {
				values[i][0].Quo(one, values[i][0])
			}

			minEvl := minimax.NewEvaluator(tc.Params, eval, btp)

			invEval := NewEvaluator(tc.Params, minEvl)

			var err error
			if ciphertext, err = invEval.GoldschmidtDivisionNew(ciphertext, logmin); err != nil {
				t.Fatal(err)
			}

			ckks.VerifyTestVectors(tc.Params, tc.Ecd, tc.Dec, values, ciphertext, 70, 0, *printPrecisionStats, t)
		})

		t.Run(name("PositiveDomain", tc), func(t *testing.T) {

			values, _, ct := tc.NewTestVector(complex(0, 0), complex(max, 0))

			minEvl := minimax.NewEvaluator(tc.Params, eval, btp)

			invEval := NewEvaluator(tc.Params, minEvl)

			cInv, err := invEval.EvaluatePositiveDomainNew(ct, logmin, logmax)
			require.NoError(t, err)

			have := make([]*big.Float, tc.Params.MaxSlots())

			require.NoError(t, tc.Ecd.Decode(tc.Dec.DecryptNew(cInv), have))

			want := make([]*big.Float, tc.Params.MaxSlots())

			threshold := bignum.NewFloat(min, tc.Params.EncodingPrecision())
			for i := range have {
				if new(big.Float).Abs(values[i][0]).Cmp(threshold) == -1 {
					want[i] = have[i] // Ignores values outside of the interval
				} else {
					want[i] = new(big.Float).Quo(bignum.NewFloat(1, tc.Params.EncodingPrecision()), values[i][0])
				}
			}

			ckks.VerifyTestVectors(tc.Params, tc.Ecd, nil, want, have, 70, 0, *printPrecisionStats, t)
		})

		t.Run(name("NegativeDomain", tc), func(t *testing.T) {

			values, _, ct := tc.NewTestVector(complex(-max, 0), complex(0, 0))

			minEvl := minimax.NewEvaluator(tc.Params, eval, btp)

			invEval := NewEvaluator(tc.Params, minEvl)

			cInv, err := invEval.EvaluateNegativeDomainNew(ct, logmin, logmax)
			require.NoError(t, err)

			have := make([]*big.Float, tc.Params.MaxSlots())

			require.NoError(t, tc.Ecd.Decode(tc.Dec.DecryptNew(cInv), have))

			want := make([]*big.Float, tc.Params.MaxSlots())

			threshold := bignum.NewFloat(min, tc.Params.EncodingPrecision())
			for i := range have {
				if new(big.Float).Abs(values[i][0]).Cmp(threshold) == -1 {
					want[i] = have[i] // Ignores values outside of the interval
				} else {
					want[i] = new(big.Float).Quo(bignum.NewFloat(1, tc.Params.EncodingPrecision()), values[i][0])
				}
			}

			ckks.VerifyTestVectors(tc.Params, tc.Ecd, nil, want, have, 70, 0, *printPrecisionStats, t)
		})

		t.Run(name("FullDomain", tc), func(t *testing.T) {

			values, _, ct := tc.NewTestVector(complex(-max, 0), complex(max, 0))

			minEvl := minimax.NewEvaluator(tc.Params, eval, btp)

			invEval := NewEvaluator(tc.Params, minEvl)

			cInv, err := invEval.EvaluateFullDomainNew(ct, logmin, logmax, minimax.NewPolynomial(comparison.DefaultCompositePolynomialForSign))
			require.NoError(t, err)

			have := make([]*big.Float, tc.Params.MaxSlots())

			require.NoError(t, tc.Ecd.Decode(tc.Dec.DecryptNew(cInv), have))

			want := make([]*big.Float, tc.Params.MaxSlots())

			threshold := bignum.NewFloat(min, tc.Params.EncodingPrecision())
			for i := range have {
				if new(big.Float).Abs(values[i][0]).Cmp(threshold) == -1 {
					want[i] = have[i] // Ignores values outside of the interval
				} else {
					want[i] = new(big.Float).Quo(bignum.NewFloat(1, tc.Params.EncodingPrecision()), values[i][0])
				}
			}

			ckks.VerifyTestVectors(tc.Params, tc.Ecd, nil, want, have, 70, 0, *printPrecisionStats, t)
		})
	}
}

func name(opname string, tc *ckks.TestContext) string {
	return fmt.Sprintf("%s/RingType=%s/logN=%d/logQP=%d/Qi=%d/Pi=%d/LogScale=%d",
		opname,
		tc.Params.RingType(),
		tc.Params.LogN(),
		int(math.Round(tc.Params.LogQP())),
		tc.Params.QCount(),
		tc.Params.PCount(),
		int(math.Log2(tc.Params.DefaultScale().Float64())))
}

// testInsecurePrec90 are insecure parameters used for the sole purpose of fast testing.
var testInsecurePrec90 = ckks.ParametersLiteral{
	LogN:            10,
	LogQ:            []int{55, 55, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45},
	LogP:            []int{60, 60},
	LogDefaultScale: 90,
}
