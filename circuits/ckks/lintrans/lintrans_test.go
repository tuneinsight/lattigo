package lintrans

import (
	"encoding/json"
	"flag"
	"fmt"
	"math"
	"math/big"
	"math/rand"
	"runtime"
	"testing"

	"github.com/stretchr/testify/require"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
	"github.com/tuneinsight/lattigo/v6/utils"
	"github.com/tuneinsight/lattigo/v6/utils/bignum"
	"github.com/tuneinsight/lattigo/v6/utils/sampling"
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

	mulCmplx := bignum.NewComplexMultiplier().Mul

	add := func(a, b, c []*bignum.Complex) {
		for i := range c {
			if a[i] != nil && b[i] != nil {
				c[i].Add(a[i], b[i])
			}
		}
	}

	muladd := func(a, b, c []*bignum.Complex) {
		tmp := &bignum.Complex{new(big.Float), new(big.Float)}
		for i := range c {
			if a[i] != nil && b[i] != nil {
				mulCmplx(a[i], b[i], tmp)
				c[i].Add(c[i], tmp)
			}
		}
	}

	prec := tc.Ecd.Prec()

	newVec := func(size int) (vec []*bignum.Complex) {
		vec = make([]*bignum.Complex, size)
		for i := range vec {
			vec[i] = &bignum.Complex{new(big.Float).SetPrec(prec), new(big.Float).SetPrec(0)}
		}
		return
	}

	t.Run(name("Average", tc), func(t *testing.T) {

		values, _, ciphertext := tc.NewTestVector(-1-1i, 1+1i)

		slots := ciphertext.Slots()

		logBatch := 9
		batch := 1 << logBatch
		n := slots / batch

		eval := tc.Evl.WithKey(rlwe.NewMemEvaluationKeySet(nil, tc.Kgen.GenGaloisKeysNew(rlwe.GaloisElementsForInnerSum(params, batch, n), tc.Sk)...))

		require.NoError(t, eval.Average(ciphertext, logBatch, ciphertext))

		tmp0 := make([]*bignum.Complex, len(values))
		for i := range tmp0 {
			tmp0[i] = values[i].Clone()
		}

		for i := 1; i < n; i++ {

			tmp1 := utils.RotateSlice(tmp0, i*batch)

			for j := range values {
				values[j].Add(values[j], tmp1[j])
			}
		}

		nB := new(big.Float).SetFloat64(float64(n))

		for i := range values {
			values[i][0].Quo(values[i][0], nB)
			values[i][1].Quo(values[i][1], nB)
		}

		ckks.VerifyTestVectors(params, tc.Ecd, tc.Dec, values, ciphertext, params.LogDefaultScale(), 0, *printPrecisionStats, t)
	})

	t.Run(name("LinearTransform/BSGS=True", tc), func(t *testing.T) {

		values, _, ciphertext := tc.NewTestVector(-1-1i, 1+1i)

		slots := ciphertext.Slots()

		nonZeroDiags := []int{-15, -4, -1, 0, 1, 2, 3, 4, 15}

		one := new(big.Float).SetInt64(1)
		zero := new(big.Float)

		diagonals := make(Diagonals[*bignum.Complex])
		for _, i := range nonZeroDiags {
			diagonals[i] = make([]*bignum.Complex, slots)
			for j := 0; j < slots>>1; j++ {
				diagonals[i][j] = &bignum.Complex{one, zero}
			}
		}

		ltparams := Parameters{
			DiagonalsIndexList:        diagonals.DiagonalsIndexList(),
			LevelQ:                    ciphertext.Level(),
			LevelP:                    params.MaxLevelP(),
			Scale:                     params.GetOptimalScalingFactor(ciphertext.Scale, params.DefaultScale(), ciphertext.Level()),
			LogDimensions:             ciphertext.LogDimensions,
			LogBabyStepGiantStepRatio: 1,
		}

		// Allocate the linear transformation
		linTransf := NewTransformation(params, ltparams)

		// Encode on the linear transformation
		require.NoError(t, Encode(tc.Ecd, diagonals, linTransf))

		galEls := linTransf.GaloisElements(params)

		evk := rlwe.NewMemEvaluationKeySet(nil, tc.Kgen.GenGaloisKeysNew(galEls, tc.Sk)...)

		ltEval := NewEvaluator(tc.Evl.WithKey(evk))

		require.NoError(t, ltEval.Evaluate(ciphertext, linTransf, ciphertext))

		values = diagonals.Evaluate(values, newVec, add, muladd)

		ckks.VerifyTestVectors(params, tc.Ecd, tc.Dec, values, ciphertext, params.LogDefaultScale(), 0, *printPrecisionStats, t)
	})

	t.Run(name("LinearTransform/BSGS=False", tc), func(t *testing.T) {

		values, _, ciphertext := tc.NewTestVector(-1-1i, 1+1i)

		slots := ciphertext.Slots()

		nonZeroDiags := []int{-15, -4, -1, 0, 1, 2, 3, 4, 15}

		one := new(big.Float).SetInt64(1)
		zero := new(big.Float)

		diagonals := make(Diagonals[*bignum.Complex])
		for _, i := range nonZeroDiags {
			diagonals[i] = make([]*bignum.Complex, slots)

			for j := 0; j < slots; j++ {
				diagonals[i][j] = &bignum.Complex{one, zero}
			}
		}

		ltparams := Parameters{
			DiagonalsIndexList:        diagonals.DiagonalsIndexList(),
			LevelQ:                    ciphertext.Level(),
			LevelP:                    params.MaxLevelP(),
			Scale:                     rlwe.NewScale(params.Q()[ciphertext.Level()]),
			LogDimensions:             ciphertext.LogDimensions,
			LogBabyStepGiantStepRatio: -1,
		}

		// Allocate the linear transformation
		linTransf := NewTransformation(params, ltparams)

		// Encode on the linear transformation
		require.NoError(t, Encode(tc.Ecd, diagonals, linTransf))

		galEls := linTransf.GaloisElements(params)

		evk := rlwe.NewMemEvaluationKeySet(nil, tc.Kgen.GenGaloisKeysNew(galEls, tc.Sk)...)

		ltEval := NewEvaluator(tc.Evl.WithKey(evk))

		require.NoError(t, ltEval.Evaluate(ciphertext, linTransf, ciphertext))

		values = diagonals.Evaluate(values, newVec, add, muladd)

		ckks.VerifyTestVectors(params, tc.Ecd, tc.Dec, values, ciphertext, params.LogDefaultScale(), 0, *printPrecisionStats, t)
	})

	t.Run(name("LinearTransform/Permutation", tc), func(t *testing.T) {
		idx := make([]int, params.MaxSlots())
		for i := range idx {
			idx[i] = i
		}
		r := rand.New(rand.NewSource(0))
		r.Shuffle(len(idx), func(i, j int) {
			idx[i], idx[j] = idx[j], idx[i]
		})

		idx = idx[:len(idx)>>1] // truncates to test partial permutation

		permutation := make([]PermutationMapping[*bignum.Complex], len(idx))

		for i := range permutation {
			permutation[i] = PermutationMapping[*bignum.Complex]{
				From: i,
				To:   idx[i],
				Scaling: &bignum.Complex{
					bignum.NewFloat(sampling.RandFloat64(real(-1), real(1)), prec),
					bignum.NewFloat(sampling.RandFloat64(imag(-1), imag(1)), prec),
				},
			}
		}

		diagonals := Permutation[*bignum.Complex](permutation).GetDiagonals(params.LogMaxSlots())

		values, _, ciphertext := tc.NewTestVector(-1-1i, 1+1i)

		ltparams := Parameters{
			DiagonalsIndexList:        diagonals.DiagonalsIndexList(),
			LevelQ:                    ciphertext.Level(),
			LevelP:                    params.MaxLevelP(),
			Scale:                     params.GetOptimalScalingFactor(ciphertext.Scale, params.DefaultScale(), ciphertext.Level()),
			LogDimensions:             ciphertext.LogDimensions,
			LogBabyStepGiantStepRatio: 1,
		}

		// Allocate the linear transformation
		linTransf := NewTransformation(params, ltparams)

		// Encode on the linear transformation
		require.NoError(t, Encode(tc.Ecd, diagonals, linTransf))

		galEls := linTransf.GaloisElements(params)

		evk := rlwe.NewMemEvaluationKeySet(nil, tc.Kgen.GenGaloisKeysNew(galEls, tc.Sk)...)

		ltEval := NewEvaluator(tc.Evl.WithKey(evk))

		require.NoError(t, ltEval.Evaluate(ciphertext, linTransf, ciphertext))

		values = diagonals.Evaluate(values, newVec, add, muladd)

		ckks.VerifyTestVectors(params, tc.Ecd, tc.Dec, values, ciphertext, params.LogDefaultScale(), 0, *printPrecisionStats, t)
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
