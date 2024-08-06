package lintrans

import (
	"encoding/json"
	"flag"
	"math/rand"
	"runtime"
	"testing"

	"github.com/stretchr/testify/require"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/schemes/bgv"
	"github.com/tuneinsight/lattigo/v6/utils/sampling"
)

var flagParamString = flag.String("params", "", "specify the test cryptographic parameters as a JSON string. Overrides -short.")

func TestLinearTransformation(t *testing.T) {
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
	rT := tc.Params.RingT().SubRings[0]

	add := func(a, b, c []uint64) {
		rT.Add(a, b, c)
	}
	muladd := func(a, b, c []uint64) {
		rT.MulCoeffsBarrettThenAdd(a, b, c)
	}
	newVec := func(size int) (vec []uint64) {
		return make([]uint64, size)
	}

	params := tc.Params

	T := params.PlaintextModulus()

	t.Run("Evaluator/LinearTransformationBSGS=true/"+tc.String(), func(t *testing.T) {

		values, _, ciphertext := bgv.NewTestVector(params, tc.Ecd, tc.Enc, params.MaxLevel(), params.DefaultScale())

		slots := ciphertext.Slots()

		nonZeroDiags := []int{-15, -4, -1, 0, 1, 2, 3, 4, 15}

		diagonals := make(Diagonals[uint64])
		for _, i := range nonZeroDiags {
			diagonals[i] = make([]uint64, slots)
			for j := 0; j < slots>>1; j++ {
				diagonals[i][j] = sampling.RandUint64() % T
			}
		}

		ltparams := Parameters{
			DiagonalsIndexList:        diagonals.DiagonalsIndexList(),
			LevelQ:                    ciphertext.Level(),
			LevelP:                    params.MaxLevelP(),
			Scale:                     tc.Params.DefaultScale(),
			LogDimensions:             ciphertext.LogDimensions,
			LogBabyStepGiantStepRatio: 1,
		}

		// Allocate the linear transformation
		linTransf := NewLinearTransformation(params, ltparams)

		// Encode on the linear transformation
		require.NoError(t, Encode(tc.Ecd, diagonals, linTransf))

		galEls := linTransf.GaloisElements(params)

		eval := tc.Evl.WithKey(rlwe.NewMemEvaluationKeySet(nil, tc.Kgen.GenGaloisKeysNew(galEls, tc.Sk)...))
		ltEval := NewEvaluator(eval)

		require.NoError(t, ltEval.Evaluate(ciphertext, linTransf, ciphertext))

		values = diagonals.Evaluate(values, newVec, add, muladd)

		bgv.VerifyTestVectors(params, tc.Ecd, tc.Dec, ciphertext, values, t)
	})

	t.Run("Evaluator/LinearTransformationBSGS=false"+tc.String(), func(t *testing.T) {

		values, _, ciphertext := bgv.NewTestVector(params, tc.Ecd, tc.Enc, params.MaxLevel(), params.DefaultScale())

		slots := ciphertext.Slots()

		nonZeroDiags := []int{-15, -4, -1, 0, 1, 2, 3, 4, 15}

		diagonals := make(Diagonals[uint64])
		for _, i := range nonZeroDiags {
			diagonals[i] = make([]uint64, slots)
			for j := 0; j < slots>>1; j++ {
				diagonals[i][j] = sampling.RandUint64() % T
			}
		}

		ltparams := Parameters{
			DiagonalsIndexList:        diagonals.DiagonalsIndexList(),
			LevelQ:                    ciphertext.Level(),
			LevelP:                    params.MaxLevelP(),
			Scale:                     params.DefaultScale(),
			LogDimensions:             ciphertext.LogDimensions,
			LogBabyStepGiantStepRatio: -1,
		}

		// Allocate the linear transformation
		linTransf := NewLinearTransformation(params, ltparams)

		// Encode on the linear transformation
		require.NoError(t, Encode(tc.Ecd, diagonals, linTransf))

		galEls := linTransf.GaloisElements(params)

		eval := tc.Evl.WithKey(rlwe.NewMemEvaluationKeySet(nil, tc.Kgen.GenGaloisKeysNew(galEls, tc.Sk)...))
		ltEval := NewEvaluator(eval)

		require.NoError(t, ltEval.Evaluate(ciphertext, linTransf, ciphertext))

		values = diagonals.Evaluate(values, newVec, add, muladd)

		bgv.VerifyTestVectors(params, tc.Ecd, tc.Dec, ciphertext, values, t)
	})

	t.Run("Evaluator/LinearTransformation/Permutation"+tc.String(), func(t *testing.T) {

		idx := [2][]int{
			make([]int, params.MaxSlots()>>1),
			make([]int, params.MaxSlots()>>1),
		}

		for i := range idx[0] {
			idx[0][i] = i
			idx[1][i] = i
		}

		r := rand.New(rand.NewSource(0))
		r.Shuffle(len(idx[0]), func(i, j int) {
			idx[0][i], idx[0][j] = idx[0][j], idx[0][i]
		})
		r.Shuffle(len(idx[1]), func(i, j int) {
			idx[1][i], idx[1][j] = idx[1][j], idx[1][i]
		})

		idx[0] = idx[0][:len(idx[0])>>1]
		idx[1] = idx[1][:len(idx[1])>>1] // truncates to test partial permutation

		permutation := [2][]PermutationMapping[uint64]{
			make([]PermutationMapping[uint64], len(idx[0])),
			make([]PermutationMapping[uint64], len(idx[1])),
		}

		for i := range permutation {
			for j := range permutation[i] {
				permutation[i][j] = PermutationMapping[uint64]{
					From:    j,
					To:      idx[i][j],
					Scaling: sampling.RandUint64() % T,
				}
				permutation[i][j] = PermutationMapping[uint64]{
					From:    j,
					To:      idx[i][j],
					Scaling: sampling.RandUint64() % T,
				}
			}
		}

		diagonals := Permutation[uint64](permutation).GetDiagonals(params.LogMaxSlots())

		values, _, ciphertext := bgv.NewTestVector(params, tc.Ecd, tc.Enc, params.MaxLevel(), params.DefaultScale())

		ltparams := Parameters{
			DiagonalsIndexList:        diagonals.DiagonalsIndexList(),
			LevelQ:                    ciphertext.Level(),
			LevelP:                    params.MaxLevelP(),
			Scale:                     params.DefaultScale(),
			LogDimensions:             ciphertext.LogDimensions,
			LogBabyStepGiantStepRatio: 1,
		}

		// Allocate the linear transformation
		linTransf := NewLinearTransformation(params, ltparams)

		// Encode on the linear transformation
		require.NoError(t, Encode(tc.Ecd, diagonals, linTransf))

		galEls := linTransf.GaloisElements(params)

		evk := rlwe.NewMemEvaluationKeySet(nil, tc.Kgen.GenGaloisKeysNew(galEls, tc.Sk)...)

		ltEval := NewEvaluator(tc.Evl.WithKey(evk))

		require.NoError(t, ltEval.Evaluate(ciphertext, linTransf, ciphertext))

		values = diagonals.Evaluate(values, newVec, add, muladd)

		bgv.VerifyTestVectors(params, tc.Ecd, tc.Dec, ciphertext, values, t)
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
