package polynomial

import (
	"encoding/json"
	"flag"
	"runtime"
	"testing"

	"github.com/stretchr/testify/require"

	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/schemes/bgv"
	"github.com/tuneinsight/lattigo/v6/utils/bignum"
)

var flagParamString = flag.String("params", "", "specify the test cryptographic parameters as a JSON string. Overrides -short.")

func TestPolynomialEvaluator(t *testing.T) {
	var err error

	paramsLiterals := testParams

	if *flagParamString != "" {
		var jsonParams bgv.ParametersLiteral
		if err = json.Unmarshal([]byte(*flagParamString), &jsonParams); err != nil {
			t.Fatal(err)
		}
		paramsLiterals = []bgv.ParametersLiteral{jsonParams} // the custom test suite reads the parameters from the -params flag
	}

	for _, p := range paramsLiterals[:] {

		for _, plaintextModulus := range testPlaintextModulus[:] {
			p.PlaintextModulus = plaintextModulus

			tc := bgv.NewTestContext(p, false)

			for _, testSet := range []func(tc *bgv.TestContext, t *testing.T){
				run,
			} {
				testSet(tc, t)
				runtime.GC()
			}
		}
	}
}

func run(tc *bgv.TestContext, t *testing.T) {
	t.Run("Single", func(t *testing.T) {

		if tc.Params.MaxLevel() < 4 {
			t.Skip("MaxLevel() to low")
		}

		values, _, ciphertext := bgv.NewTestVector(tc.Params, tc.Ecd, tc.Enc, tc.Params.MaxLevel(), tc.Params.DefaultScale())

		coeffs := []uint64{0, 0, 1}

		T := tc.Params.PlaintextModulus()
		for i := range values {
			values[i] = ring.EvalPolyModP(values[i], coeffs, T)
		}

		poly := bignum.NewPolynomial(bignum.Monomial, coeffs, nil)

		t.Run("Standard"+tc.String(), func(t *testing.T) {

			polyEval := NewEvaluator(tc.Params, tc.Evl)

			res, err := polyEval.Evaluate(ciphertext, poly, tc.Params.DefaultScale())
			require.NoError(t, err)

			require.Equal(t, res.Scale.Cmp(tc.Params.DefaultScale()), 0)

			bgv.VerifyTestVectors(tc.Params, tc.Ecd, tc.Dec, res, values, t)
		})

		t.Run("Invariant"+tc.String(), func(t *testing.T) {
			tc.Evl.ScaleInvariant = true

			polyEval := NewEvaluator(tc.Params, tc.Evl)

			res, err := polyEval.Evaluate(ciphertext, poly, tc.Params.DefaultScale())
			require.NoError(t, err)

			require.Equal(t, res.Level(), ciphertext.Level())
			require.Equal(t, res.Scale.Cmp(tc.Params.DefaultScale()), 0)

			bgv.VerifyTestVectors(tc.Params, tc.Ecd, tc.Dec, res, values, t)
		})
	})
}

var (
	// testInsecure are insecure parameters used for the sole purpose of fast testing.
	testInsecure = bgv.ParametersLiteral{
		LogN: 10,
		Q:    []uint64{0x3fffffa8001, 0x1000090001, 0x10000c8001, 0x10000f0001, 0xffff00001},
		P:    []uint64{0x7fffffd8001},
	}

	testPlaintextModulus = []uint64{0x101, 0xffc001}

	testParams = []bgv.ParametersLiteral{testInsecure}
)
