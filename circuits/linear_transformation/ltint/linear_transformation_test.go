package ltint

import (
	"encoding/json"
	"flag"
	"fmt"
	"math"
	"math/rand"
	"runtime"
	"slices"
	"testing"

	"github.com/stretchr/testify/require"
	"github.com/tuneinsight/lattigo/v5/core/rlwe"
	"github.com/tuneinsight/lattigo/v5/ring"
	"github.com/tuneinsight/lattigo/v5/schemes"
	"github.com/tuneinsight/lattigo/v5/schemes/bgv"
	"github.com/tuneinsight/lattigo/v5/utils/sampling"
)

type testContext struct {
	params      bgv.Parameters
	ringQ       *ring.Ring
	ringT       *ring.Ring
	prng        sampling.PRNG
	uSampler    *ring.UniformSampler
	encoder     schemes.Encoder
	kgen        *rlwe.KeyGenerator
	sk          *rlwe.SecretKey
	pk          *rlwe.PublicKey
	encryptorPk *rlwe.Encryptor
	encryptorSk *rlwe.Encryptor
	decryptor   *rlwe.Decryptor
	evaluator   *bgv.Evaluator
	testLevel   []int
}

func genTestParams(params bgv.Parameters) (tc *testContext, err error) {

	tc = new(testContext)
	tc.params = params

	if tc.prng, err = sampling.NewPRNG(); err != nil {
		return nil, err
	}

	tc.ringQ = params.RingQ()
	tc.ringT = params.RingT()

	tc.uSampler = ring.NewUniformSampler(tc.prng, tc.ringT)
	tc.kgen = rlwe.NewKeyGenerator(tc.params)
	tc.sk, tc.pk = tc.kgen.GenKeyPairNew()
	tc.encoder = bgv.NewEncoder(tc.params)

	tc.encryptorPk = rlwe.NewEncryptor(tc.params, tc.pk)
	tc.encryptorSk = rlwe.NewEncryptor(tc.params, tc.sk)
	tc.decryptor = rlwe.NewDecryptor(tc.params, tc.sk)
	tc.evaluator = bgv.NewEvaluator(tc.params, rlwe.NewMemEvaluationKeySet(tc.kgen.GenRelinearizationKeyNew(tc.sk)))

	tc.testLevel = []int{0, params.MaxLevel()}

	return
}

var flagPrintNoise = flag.Bool("print-noise", false, "print the residual noise")
var flagParamString = flag.String("params", "", "specify the test cryptographic parameters as a JSON string. Overrides -short.")

func GetTestName(opname string, p bgv.Parameters, lvl int) string {
	return fmt.Sprintf("%s/LogN=%d/logQ=%d/logP=%d/LogSlots=%dx%d/logT=%d/Qi=%d/Pi=%d/lvl=%d",
		opname,
		p.LogN(),
		int(math.Round(p.LogQ())),
		int(math.Round(p.LogP())),
		p.LogMaxDimensions().Rows,
		p.LogMaxDimensions().Cols,
		int(math.Round(p.LogT())),
		p.QCount(),
		p.PCount(),
		lvl)
}

func TestPolynomialEvaluator(t *testing.T) {
	var err error

	paramsLiterals := schemes.BgvTestParams

	if *flagParamString != "" {
		var jsonParams bgv.ParametersLiteral
		if err = json.Unmarshal([]byte(*flagParamString), &jsonParams); err != nil {
			t.Fatal(err)
		}
		paramsLiterals = []bgv.ParametersLiteral{jsonParams} // the custom test suite reads the parameters from the -params flag
	}

	for _, p := range paramsLiterals[:] {

		for _, plaintextModulus := range schemes.BgvTestPlaintextModulus[:] {

			p.PlaintextModulus = plaintextModulus

			var params bgv.Parameters
			if params, err = bgv.NewParametersFromLiteral(p); err != nil {
				t.Error(err)
				t.Fail()
			}

			var tc *testContext
			if tc, err = genTestParams(params); err != nil {
				t.Error(err)
				t.Fail()
			}

			for _, testSet := range []func(tc *testContext, t *testing.T){
				run,
			} {
				testSet(tc, t)
				runtime.GC()
			}
		}
	}
}

func run(tc *testContext, t *testing.T) {
	rT := tc.params.RingT().SubRings[0]

	add := func(a, b, c []uint64) {
		rT.Add(a, b, c)
	}

	muladd := func(a, b, c []uint64) {
		rT.MulCoeffsBarrettThenAdd(a, b, c)
	}

	newVec := func(size int) (vec []uint64) {
		return make([]uint64, size)
	}

	params := tc.params

	T := params.PlaintextModulus()

	level := tc.params.MaxLevel()

	t.Run(GetTestName("Evaluator/LinearTransformationBSGS=true", params, level), func(t *testing.T) {

		values, _, ciphertext := newTestVectorsLvl(level, params.DefaultScale(), tc, tc.encryptorSk)

		slots := ciphertext.Slots()

		nonZeroDiags := []int{-15, -4, -1, 0, 1, 2, 3, 4, 15}

		diagonals := make(Diagonals[uint64])
		for _, i := range nonZeroDiags {
			diagonals[i] = make([]uint64, slots)
			for j := 0; j < slots>>1; j++ {
				diagonals[i][j] = sampling.RandUint64() % T
			}
		}

		ltparams := LinearTransformationParameters{
			DiagonalsIndexList:        diagonals.DiagonalsIndexList(),
			LevelQ:                    ciphertext.Level(),
			LevelP:                    params.MaxLevelP(),
			Scale:                     tc.params.DefaultScale(),
			LogDimensions:             ciphertext.LogDimensions,
			LogBabyStepGiantStepRatio: 1,
		}

		// Allocate the linear transformation
		linTransf := NewLinearTransformation(params, ltparams)

		// Encode on the linear transformation
		require.NoError(t, EncodeLinearTransformation(tc.encoder, diagonals, linTransf))

		galEls := linTransf.GaloisElements(params)

		eval := tc.evaluator.WithKey(rlwe.NewMemEvaluationKeySet(nil, tc.kgen.GenGaloisKeysNew(galEls, tc.sk)...))
		ltEval := NewLinearTransformationEvaluator(eval)

		require.NoError(t, ltEval.Evaluate(ciphertext, linTransf, ciphertext))

		values.Coeffs[0] = diagonals.Evaluate(values.Coeffs[0], newVec, add, muladd)

		verifyTestVectors(tc, tc.decryptor, values, ciphertext, t)
	})

	t.Run(GetTestName("Evaluator/LinearTransformationBSGS=false", params, level), func(t *testing.T) {

		values, _, ciphertext := newTestVectorsLvl(level, params.DefaultScale(), tc, tc.encryptorSk)

		slots := ciphertext.Slots()

		nonZeroDiags := []int{-15, -4, -1, 0, 1, 2, 3, 4, 15}

		diagonals := make(Diagonals[uint64])
		for _, i := range nonZeroDiags {
			diagonals[i] = make([]uint64, slots)
			for j := 0; j < slots>>1; j++ {
				diagonals[i][j] = sampling.RandUint64() % T
			}
		}

		ltparams := LinearTransformationParameters{
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
		require.NoError(t, EncodeLinearTransformation(tc.encoder, diagonals, linTransf))

		galEls := linTransf.GaloisElements(params)

		eval := tc.evaluator.WithKey(rlwe.NewMemEvaluationKeySet(nil, tc.kgen.GenGaloisKeysNew(galEls, tc.sk)...))
		ltEval := NewLinearTransformationEvaluator(eval)

		require.NoError(t, ltEval.Evaluate(ciphertext, linTransf, ciphertext))

		values.Coeffs[0] = diagonals.Evaluate(values.Coeffs[0], newVec, add, muladd)

		verifyTestVectors(tc, tc.decryptor, values, ciphertext, t)
	})

	t.Run(GetTestName("Evaluator/LinearTransformation/Permutation", params, level), func(t *testing.T) {

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

		values, _, ciphertext := newTestVectorsLvl(level, tc.params.NewScale(1), tc, tc.encryptorSk)

		ltparams := LinearTransformationParameters{
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
		require.NoError(t, EncodeLinearTransformation(tc.encoder, diagonals, linTransf))

		galEls := linTransf.GaloisElements(params)

		evk := rlwe.NewMemEvaluationKeySet(nil, tc.kgen.GenGaloisKeysNew(galEls, tc.sk)...)

		ltEval := NewLinearTransformationEvaluator(tc.evaluator.WithKey(evk))

		require.NoError(t, ltEval.Evaluate(ciphertext, linTransf, ciphertext))

		values.Coeffs[0] = diagonals.Evaluate(values.Coeffs[0], newVec, add, muladd)

		verifyTestVectors(tc, tc.decryptor, values, ciphertext, t)
	})
}

func newTestVectorsLvl(level int, scale rlwe.Scale, tc *testContext, encryptor *rlwe.Encryptor) (coeffs ring.Poly, plaintext *rlwe.Plaintext, ciphertext *rlwe.Ciphertext) {
	coeffs = tc.uSampler.ReadNew()
	for i := range coeffs.Coeffs[0] {
		coeffs.Coeffs[0][i] = uint64(i)
	}

	plaintext = bgv.NewPlaintext(tc.params, level)
	plaintext.Scale = scale
	tc.encoder.Encode(coeffs.Coeffs[0], plaintext)
	if encryptor != nil {
		var err error
		ciphertext, err = encryptor.EncryptNew(plaintext)
		if err != nil {
			panic(err)
		}
	}

	return coeffs, plaintext, ciphertext
}

func verifyTestVectors(tc *testContext, decryptor *rlwe.Decryptor, coeffs ring.Poly, element rlwe.ElementInterface[ring.Poly], t *testing.T) {

	coeffsTest := make([]uint64, tc.params.MaxSlots())

	switch el := element.(type) {
	case *rlwe.Plaintext:
		require.NoError(t, tc.encoder.Decode(el, coeffsTest))
	case *rlwe.Ciphertext:

		pt := decryptor.DecryptNew(el)

		require.NoError(t, tc.encoder.Decode(pt, coeffsTest))

		if *flagPrintNoise {
			require.NoError(t, tc.encoder.Encode(coeffsTest, pt))
			ct, err := tc.evaluator.SubNew(el, pt)
			require.NoError(t, err)
			vartmp, _, _ := rlwe.Norm(ct, decryptor)
			t.Logf("STD(noise): %f\n", vartmp)
		}

	default:
		t.Error("invalid test object to verify")
	}

	require.True(t, slices.Equal(coeffs.Coeffs[0], coeffsTest))
}
