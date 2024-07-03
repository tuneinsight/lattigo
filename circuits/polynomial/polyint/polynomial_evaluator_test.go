package polyint

import (
	"encoding/json"
	"flag"
	"runtime"
	"testing"

	"github.com/stretchr/testify/require"

	"github.com/tuneinsight/lattigo/v5/ring"
	"github.com/tuneinsight/lattigo/v5/schemes/bgv"
	"github.com/tuneinsight/lattigo/v5/utils/bignum"
)

var flagParamString = flag.String("params", "", "specify the test cryptographic parameters as a JSON string. Overrides -short.")

func TestPolynomialEvaluator(t *testing.T) {
	var err error

	paramsLiterals := bgv.TestParams

	if *flagParamString != "" {
		var jsonParams bgv.ParametersLiteral
		if err = json.Unmarshal([]byte(*flagParamString), &jsonParams); err != nil {
			t.Fatal(err)
		}
		paramsLiterals = []bgv.ParametersLiteral{jsonParams} // the custom test suite reads the parameters from the -params flag
	}

	for _, p := range paramsLiterals[:] {

		for _, plaintextModulus := range bgv.TestPlaintextModulus[:] {
			p.PlaintextModulus = plaintextModulus

			tc := bgv.NewTestContext(p)

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

			polyEval := NewPolynomialEvaluator(tc.Params, tc.Evl, false)

			res, err := polyEval.Evaluate(ciphertext, poly, tc.Params.DefaultScale())
			require.NoError(t, err)

			require.Equal(t, res.Scale.Cmp(tc.Params.DefaultScale()), 0)

			bgv.VerifyTestVectors(tc.Params, tc.Ecd, tc.Dec, res, values, t)
		})

		t.Run("Invariant"+tc.String(), func(t *testing.T) {

			polyEval := NewPolynomialEvaluator(tc.Params, tc.Evl, true)

			res, err := polyEval.Evaluate(ciphertext, poly, tc.Params.DefaultScale())
			require.NoError(t, err)

			require.Equal(t, res.Level(), ciphertext.Level())
			require.Equal(t, res.Scale.Cmp(tc.Params.DefaultScale()), 0)

			bgv.VerifyTestVectors(tc.Params, tc.Ecd, tc.Dec, res, values, t)
		})
	})
}
