package ltfloat

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

	"github.com/tuneinsight/lattigo/v5/core/rlwe"
	"github.com/tuneinsight/lattigo/v5/ring"
	"github.com/tuneinsight/lattigo/v5/schemes"
	"github.com/tuneinsight/lattigo/v5/schemes/ckks"
	"github.com/tuneinsight/lattigo/v5/utils"
	"github.com/tuneinsight/lattigo/v5/utils/bignum"
	"github.com/tuneinsight/lattigo/v5/utils/sampling"
)

var flagParamString = flag.String("params", "", "specify the test cryptographic parameters as a JSON string. Overrides -short.")
var printPrecisionStats = flag.Bool("print-precision", false, "print precision stats")

func GetTestName(params ckks.Parameters, opname string) string {
	return fmt.Sprintf("%s/RingType=%s/logN=%d/logQP=%d/Qi=%d/Pi=%d/LogScale=%d",
		opname,
		params.RingType(),
		params.LogN(),
		int(math.Round(params.LogQP())),
		params.QCount(),
		params.PCount(),
		int(math.Log2(params.DefaultScale().Float64())))
}

type testContext struct {
	params      ckks.Parameters
	ringQ       *ring.Ring
	ringP       *ring.Ring
	prng        sampling.PRNG
	encoder     ckks.Encoder
	kgen        *rlwe.KeyGenerator
	sk          *rlwe.SecretKey
	pk          *rlwe.PublicKey
	encryptorPk *rlwe.Encryptor
	encryptorSk *rlwe.Encryptor
	decryptor   *rlwe.Decryptor
	evaluator   ckks.Evaluator
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
		testParams = schemes.CkksTestParametersLiteral
	}

	for _, ringType := range []ring.Type{ring.Standard, ring.ConjugateInvariant} {

		for _, paramsLiteral := range testParams {

			paramsLiteral.RingType = ringType

			if testing.Short() {
				paramsLiteral.LogN = 10
			}

			var params ckks.Parameters
			if params, err = ckks.NewParametersFromLiteral(paramsLiteral); err != nil {
				t.Fatal(err)
			}

			var tc *testContext
			if tc, err = genTestParams(params); err != nil {
				t.Fatal(err)
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

	params := tc.params

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

	prec := tc.encoder.Prec()

	newVec := func(size int) (vec []*bignum.Complex) {
		vec = make([]*bignum.Complex, size)
		for i := range vec {
			vec[i] = &bignum.Complex{new(big.Float).SetPrec(prec), new(big.Float).SetPrec(0)}
		}
		return
	}

	t.Run(GetTestName(params, "Average"), func(t *testing.T) {

		values, _, ciphertext := newTestVectors(tc, tc.encryptorSk, -1-1i, 1+1i, t)

		slots := ciphertext.Slots()

		logBatch := 9
		batch := 1 << logBatch
		n := slots / batch

		eval := tc.evaluator.WithKey(rlwe.NewMemEvaluationKeySet(nil, tc.kgen.GenGaloisKeysNew(rlwe.GaloisElementsForInnerSum(params, batch, n), tc.sk)...))

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

		ckks.VerifyTestVectors(params, &tc.encoder, tc.decryptor, values, ciphertext, params.LogDefaultScale(), 0, *printPrecisionStats, t)
	})

	t.Run(GetTestName(params, "LinearTransform/BSGS=True"), func(t *testing.T) {

		values, _, ciphertext := newTestVectors(tc, tc.encryptorSk, -1-1i, 1+1i, t)

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

		ltparams := LinearTransformationParameters{
			DiagonalsIndexList:        diagonals.DiagonalsIndexList(),
			LevelQ:                    ciphertext.Level(),
			LevelP:                    params.MaxLevelP(),
			Scale:                     params.GetOptimalScalingFactor(ciphertext.Scale, params.DefaultScale(), ciphertext.Level()),
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

		values = diagonals.Evaluate(values, newVec, add, muladd)

		ckks.VerifyTestVectors(params, &tc.encoder, tc.decryptor, values, ciphertext, params.LogDefaultScale(), 0, *printPrecisionStats, t)
	})

	t.Run(GetTestName(params, "LinearTransform/BSGS=False"), func(t *testing.T) {

		values, _, ciphertext := newTestVectors(tc, tc.encryptorSk, -1-1i, 1+1i, t)

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

		ltparams := LinearTransformationParameters{
			DiagonalsIndexList:        diagonals.DiagonalsIndexList(),
			LevelQ:                    ciphertext.Level(),
			LevelP:                    params.MaxLevelP(),
			Scale:                     rlwe.NewScale(params.Q()[ciphertext.Level()]),
			LogDimensions:             ciphertext.LogDimensions,
			LogBabyStepGiantStepRatio: -1,
		}

		// Allocate the linear transformation
		linTransf := NewLinearTransformation(params, ltparams)

		// Encode on the linear transformation
		require.NoError(t, EncodeLinearTransformation(tc.encoder, diagonals, linTransf))

		galEls := linTransf.GaloisElements(params)

		evk := rlwe.NewMemEvaluationKeySet(nil, tc.kgen.GenGaloisKeysNew(galEls, tc.sk)...)

		ltEval := NewLinearTransformationEvaluator(tc.evaluator.WithKey(evk))

		require.NoError(t, ltEval.Evaluate(ciphertext, linTransf, ciphertext))

		values = diagonals.Evaluate(values, newVec, add, muladd)

		ckks.VerifyTestVectors(params, &tc.encoder, tc.decryptor, values, ciphertext, params.LogDefaultScale(), 0, *printPrecisionStats, t)
	})

	t.Run(GetTestName(params, "LinearTransform/Permutation"), func(t *testing.T) {
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

		values, _, ciphertext := newTestVectors(tc, tc.encryptorSk, -1-1i, 1+1i, t)

		ltparams := LinearTransformationParameters{
			DiagonalsIndexList:        diagonals.DiagonalsIndexList(),
			LevelQ:                    ciphertext.Level(),
			LevelP:                    params.MaxLevelP(),
			Scale:                     params.GetOptimalScalingFactor(ciphertext.Scale, params.DefaultScale(), ciphertext.Level()),
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

		values = diagonals.Evaluate(values, newVec, add, muladd)

		ckks.VerifyTestVectors(params, &tc.encoder, tc.decryptor, values, ciphertext, params.LogDefaultScale(), 0, *printPrecisionStats, t)
	})
}

func genTestParams(defaultParam ckks.Parameters) (tc *testContext, err error) {

	tc = new(testContext)

	tc.params = defaultParam

	tc.kgen = rlwe.NewKeyGenerator(tc.params)

	tc.sk, tc.pk = tc.kgen.GenKeyPairNew()

	tc.ringQ = defaultParam.RingQ()
	if tc.params.PCount() != 0 {
		tc.ringP = defaultParam.RingP()
	}

	if tc.prng, err = sampling.NewPRNG(); err != nil {
		return nil, err
	}

	tc.encoder = *ckks.NewEncoder(tc.params)

	tc.encryptorPk = rlwe.NewEncryptor(tc.params, tc.pk)
	tc.encryptorSk = rlwe.NewEncryptor(tc.params, tc.sk)
	tc.decryptor = rlwe.NewDecryptor(tc.params, tc.sk)
	tc.evaluator = *ckks.NewEvaluator(tc.params, rlwe.NewMemEvaluationKeySet(tc.kgen.GenRelinearizationKeyNew(tc.sk)))

	return tc, nil

}

func newTestVectors(tc *testContext, encryptor *rlwe.Encryptor, a, b complex128, t *testing.T) (values []*bignum.Complex, pt *rlwe.Plaintext, ct *rlwe.Ciphertext) {

	var err error

	prec := tc.params.EncodingPrecision()

	pt = ckks.NewPlaintext(tc.params, tc.params.MaxLevel())

	values = make([]*bignum.Complex, pt.Slots())

	switch tc.params.RingType() {
	case ring.Standard:
		for i := range values {
			values[i] = &bignum.Complex{
				bignum.NewFloat(sampling.RandFloat64(real(a), real(b)), prec),
				bignum.NewFloat(sampling.RandFloat64(imag(a), imag(b)), prec),
			}
		}
	case ring.ConjugateInvariant:
		for i := range values {
			values[i] = &bignum.Complex{
				bignum.NewFloat(sampling.RandFloat64(real(a), real(b)), prec),
				new(big.Float),
			}
		}
	default:
		t.Fatal("invalid ring type")
	}

	tc.encoder.Encode(values, pt)

	if encryptor != nil {
		ct, err = encryptor.EncryptNew(pt)
		require.NoError(t, err)
	}

	return values, pt, ct
}
