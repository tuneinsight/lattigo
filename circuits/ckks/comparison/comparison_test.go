package comparison

import (
	"flag"
	"fmt"
	"math"
	"math/big"
	"testing"

	"github.com/stretchr/testify/require"

	"github.com/tuneinsight/lattigo/v6/circuits/ckks/bootstrapping"
	"github.com/tuneinsight/lattigo/v6/circuits/ckks/minimax"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
)

var printPrecisionStats = flag.Bool("print-precision", false, "print precision stats")

func TestComparisons(t *testing.T) {
	var err error

	paramsLiteral := testInsecurePrec90

	for _, ringType := range []ring.Type{ring.Standard, ring.ConjugateInvariant} {

		paramsLiteral.RingType = ringType

		if testing.Short() {
			paramsLiteral.LogN = 10
		}

		tc := ckks.NewTestContext(paramsLiteral)

		btp := bootstrapping.NewSecretKeyBootstrapper(tc.Params, tc.Sk)

		var galKeys []*rlwe.GaloisKey
		if tc.Params.RingType() == ring.Standard {
			galKeys = append(galKeys, tc.Kgen.GenGaloisKeyNew(tc.Params.GaloisElementForComplexConjugation(), tc.Sk))
		}

		eval := tc.Evl.WithKey(rlwe.NewMemEvaluationKeySet(tc.Kgen.GenRelinearizationKeyNew(tc.Sk), galKeys...))

		polys := minimax.NewPolynomial(DefaultCompositePolynomialForSign)

		minimaxEvl := minimax.NewEvaluator(tc.Params, eval, btp)

		CmpEval := NewEvaluator(tc.Params, minimaxEvl, polys)

		t.Run(name("Sign", tc), func(t *testing.T) {

			values, _, ct := tc.NewTestVector(complex(-1, 0), complex(1, 0))

			var sign *rlwe.Ciphertext
			sign, err = CmpEval.Sign(ct)
			require.NoError(t, err)

			have := make([]*big.Float, tc.Params.MaxSlots())

			require.NoError(t, tc.Ecd.Decode(tc.Dec.DecryptNew(sign), have))

			want := make([]*big.Float, tc.Params.MaxSlots())

			for i := range have {
				want[i] = polys.Evaluate(values[i])[0]
			}

			ckks.VerifyTestVectors(tc.Params, tc.Ecd, nil, want, have, tc.Params.LogDefaultScale(), 0, *printPrecisionStats, t)
		})

		t.Run(name("Step", tc), func(t *testing.T) {

			values, _, ct := tc.NewTestVector(complex(-1, 0), complex(1, 0))

			var step *rlwe.Ciphertext
			step, err = CmpEval.Step(ct)
			require.NoError(t, err)

			have := make([]*big.Float, tc.Params.MaxSlots())

			require.NoError(t, tc.Ecd.Decode(tc.Dec.DecryptNew(step), have))

			want := make([]*big.Float, tc.Params.MaxSlots())

			half := new(big.Float).SetFloat64(0.5)

			for i := range have {
				want[i] = polys.Evaluate(values[i])[0]
				want[i].Mul(want[i], half)
				want[i].Add(want[i], half)
			}

			ckks.VerifyTestVectors(tc.Params, tc.Ecd, nil, want, have, tc.Params.LogDefaultScale(), 0, *printPrecisionStats, t)
		})

		t.Run(name("Max", tc), func(t *testing.T) {

			values0, _, ct0 := tc.NewTestVector(complex(-0.5, 0), complex(0.5, 0))
			values1, _, ct1 := tc.NewTestVector(complex(-0.5, 0), complex(0.5, 0))

			var max *rlwe.Ciphertext
			max, err = CmpEval.Max(ct0, ct1)
			require.NoError(t, err)

			have := make([]*big.Float, tc.Params.MaxSlots())

			require.NoError(t, tc.Ecd.Decode(tc.Dec.DecryptNew(max), have))

			want := make([]*big.Float, tc.Params.MaxSlots())

			for i := range have {

				if values0[i][0].Cmp(values1[i][0]) == -1 {
					want[i] = values1[i][0]
				} else {
					want[i] = values0[i][0]
				}
			}

			ckks.VerifyTestVectors(tc.Params, tc.Ecd, nil, want, have, tc.Params.LogDefaultScale(), 0, *printPrecisionStats, t)
		})

		t.Run(name("Min", tc), func(t *testing.T) {

			values0, _, ct0 := tc.NewTestVector(complex(-0.5, 0), complex(0.5, 0))
			values1, _, ct1 := tc.NewTestVector(complex(-0.5, 0), complex(0.5, 0))

			var max *rlwe.Ciphertext
			max, err = CmpEval.Min(ct0, ct1)
			require.NoError(t, err)

			have := make([]*big.Float, tc.Params.MaxSlots())

			require.NoError(t, tc.Ecd.Decode(tc.Dec.DecryptNew(max), have))

			want := make([]*big.Float, tc.Params.MaxSlots())

			for i := range have {

				if values0[i][0].Cmp(values1[i][0]) == 1 {
					want[i] = values1[i][0]
				} else {
					want[i] = values0[i][0]
				}
			}

			ckks.VerifyTestVectors(tc.Params, tc.Ecd, nil, want, have, tc.Params.LogDefaultScale(), 0, *printPrecisionStats, t)
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
