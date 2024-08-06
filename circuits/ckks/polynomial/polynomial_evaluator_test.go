package polynomial

import (
	"encoding/json"
	"flag"
	"fmt"
	"math"
	"math/big"
	"runtime"
	"testing"

	"github.com/stretchr/testify/require"

	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
	"github.com/tuneinsight/lattigo/v6/utils/bignum"
)

var flagParamString = flag.String("params", "", "specify the test cryptographic parameters as a JSON string. Overrides -short.")
var printPrecisionStats = flag.Bool("print-precision", false, "print precision stats")

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

func TestPolynomialEvaluator(t *testing.T) {
	var err error

	var testParams []ckks.ParametersLiteral
	switch {
	case *flagParamString != "": // the custom test suite reads the parameters from the -params flag
		testParams = append(testParams, ckks.ParametersLiteral{})
		if err = json.Unmarshal([]byte(*flagParamString), &testParams[0]); err != nil {
			t.Fatal(err)
		}
	default:
		testParams = testParametersLiteral
	}

	for _, ringType := range []ring.Type{ring.Standard, ring.ConjugateInvariant} {

		for _, paramsLiteral := range testParams {

			paramsLiteral.RingType = ringType

			if testing.Short() {
				paramsLiteral.LogN = 10
			}

			tc := ckks.NewTestContext(paramsLiteral)

			for _, testSet := range []func(tc *ckks.TestContext, t *testing.T){
				run,
			} {
				testSet(tc, t)
				runtime.GC()
			}
		}
	}
}

func run(tc *ckks.TestContext, t *testing.T) {

	params := tc.Params

	var err error

	polyEval := NewEvaluator(params, tc.Evl)

	t.Run(name("EvaluatePoly/PolySingle/Exp", tc), func(t *testing.T) {

		if params.MaxLevel() < 3 {
			t.Skip("skipping test for params max level < 3")
		}

		values, _, ciphertext := tc.NewTestVector(-1, 1)

		prec := params.EncodingPrecision()

		coeffs := []*big.Float{
			bignum.NewFloat(1, prec),
			bignum.NewFloat(1, prec),
			new(big.Float).Quo(bignum.NewFloat(1, prec), bignum.NewFloat(2, prec)),
			new(big.Float).Quo(bignum.NewFloat(1, prec), bignum.NewFloat(6, prec)),
			new(big.Float).Quo(bignum.NewFloat(1, prec), bignum.NewFloat(24, prec)),
			new(big.Float).Quo(bignum.NewFloat(1, prec), bignum.NewFloat(120, prec)),
			new(big.Float).Quo(bignum.NewFloat(1, prec), bignum.NewFloat(720, prec)),
			new(big.Float).Quo(bignum.NewFloat(1, prec), bignum.NewFloat(5040, prec)),
		}

		poly := bignum.NewPolynomial(bignum.Monomial, coeffs, nil)

		for i := range values {
			values[i] = poly.Evaluate(values[i])
		}

		if ciphertext, err = polyEval.Evaluate(ciphertext, poly, ciphertext.Scale); err != nil {
			t.Fatal(err)
		}

		ckks.VerifyTestVectors(params, tc.Ecd, tc.Dec, values, ciphertext, params.LogDefaultScale(), 0, *printPrecisionStats, t)
	})

	t.Run(name("Polynomial/PolyVector/Exp", tc), func(t *testing.T) {

		if params.MaxLevel() < 3 {
			t.Skip("skipping test for params max level < 3")
		}

		values, _, ciphertext := tc.NewTestVector(-1, 1)

		prec := params.EncodingPrecision()

		coeffs := []*big.Float{
			bignum.NewFloat(1, prec),
			bignum.NewFloat(1, prec),
			new(big.Float).Quo(bignum.NewFloat(1, prec), bignum.NewFloat(2, prec)),
			new(big.Float).Quo(bignum.NewFloat(1, prec), bignum.NewFloat(6, prec)),
			new(big.Float).Quo(bignum.NewFloat(1, prec), bignum.NewFloat(24, prec)),
			new(big.Float).Quo(bignum.NewFloat(1, prec), bignum.NewFloat(120, prec)),
			new(big.Float).Quo(bignum.NewFloat(1, prec), bignum.NewFloat(720, prec)),
			new(big.Float).Quo(bignum.NewFloat(1, prec), bignum.NewFloat(5040, prec)),
		}

		poly := bignum.NewPolynomial(bignum.Monomial, coeffs, nil)

		slots := ciphertext.Slots()

		slotIndex := make(map[int][]int)
		idx := make([]int, slots>>1)
		for i := 0; i < slots>>1; i++ {
			idx[i] = 2 * i
		}

		slotIndex[0] = idx

		valuesWant := make([]*bignum.Complex, slots)
		for _, j := range idx {
			valuesWant[j] = poly.Evaluate(values[j])
		}

		polyVector, err := NewPolynomialVector([]bignum.Polynomial{poly}, slotIndex)
		require.NoError(t, err)

		if ciphertext, err = polyEval.Evaluate(ciphertext, polyVector, ciphertext.Scale); err != nil {
			t.Fatal(err)
		}

		ckks.VerifyTestVectors(params, tc.Ecd, tc.Dec, valuesWant, ciphertext, params.LogDefaultScale(), 0, *printPrecisionStats, t)
	})
}

var (
	// testInsecurePrec45 are insecure parameters used for the sole purpose of fast testing.
	testInsecurePrec45 = ckks.ParametersLiteral{
		LogN:            10,
		LogQ:            []int{55, 45, 45, 45, 45, 45, 45},
		LogP:            []int{60},
		LogDefaultScale: 45,
	}

	// testInsecurePrec90 are insecure parameters used for the sole purpose of fast testing.
	testInsecurePrec90 = ckks.ParametersLiteral{
		LogN:            10,
		LogQ:            []int{55, 55, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45},
		LogP:            []int{60, 60},
		LogDefaultScale: 90,
	}

	testParametersLiteral = []ckks.ParametersLiteral{testInsecurePrec45, testInsecurePrec90}
)
