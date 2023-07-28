package circuits

import (
	"encoding/json"
	"flag"
	"fmt"
	"math"
	"math/big"
	"runtime"
	"testing"

	"github.com/stretchr/testify/require"
	"github.com/tuneinsight/lattigo/v4/ckks"
	"github.com/tuneinsight/lattigo/v4/ring"

	"github.com/tuneinsight/lattigo/v4/rlwe"
	"github.com/tuneinsight/lattigo/v4/utils"
	"github.com/tuneinsight/lattigo/v4/utils/bignum"
	"github.com/tuneinsight/lattigo/v4/utils/sampling"
)

var printPrecisionStats = flag.Bool("print-precision", false, "print precision stats")

func GetCKKSTestName(params ckks.Parameters, opname string) string {
	return fmt.Sprintf("%s/RingType=%s/logN=%d/logQP=%d/Qi=%d/Pi=%d/LogScale=%d",
		opname,
		params.RingType(),
		params.LogN(),
		int(math.Round(params.LogQP())),
		params.QCount(),
		params.PCount(),
		int(math.Log2(params.DefaultScale().Float64())))
}

type ckksTestContext struct {
	params      ckks.Parameters
	ringQ       *ring.Ring
	ringP       *ring.Ring
	prng        sampling.PRNG
	encoder     *ckks.Encoder
	kgen        *rlwe.KeyGenerator
	sk          *rlwe.SecretKey
	pk          *rlwe.PublicKey
	encryptorPk *rlwe.Encryptor
	encryptorSk *rlwe.Encryptor
	decryptor   *rlwe.Decryptor
	evaluator   *ckks.Evaluator
}

func TestCKKS(t *testing.T) {

	var err error

	var testParams []ckks.ParametersLiteral
	switch {
	case *flagParamString != "": // the custom test suite reads the parameters from the -params flag
		testParams = append(testParams, ckks.ParametersLiteral{})
		if err = json.Unmarshal([]byte(*flagParamString), &testParams[0]); err != nil {
			t.Fatal(err)
		}
	default:
		testParams = testParamsLiteral
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

			var tc *ckksTestContext
			if tc, err = genCKKSTestParams(params); err != nil {
				t.Fatal(err)
			}

			for _, testSet := range []func(tc *ckksTestContext, t *testing.T){
				testCKKSLinearTransformation,
			} {
				testSet(tc, t)
				runtime.GC()
			}
		}
	}

}

func genCKKSTestParams(defaultParam ckks.Parameters) (tc *ckksTestContext, err error) {

	tc = new(ckksTestContext)

	tc.params = defaultParam

	tc.kgen = ckks.NewKeyGenerator(tc.params)

	tc.sk, tc.pk = tc.kgen.GenKeyPairNew()

	tc.ringQ = defaultParam.RingQ()
	if tc.params.PCount() != 0 {
		tc.ringP = defaultParam.RingP()
	}

	if tc.prng, err = sampling.NewPRNG(); err != nil {
		return nil, err
	}

	tc.encoder = ckks.NewEncoder(tc.params)

	if tc.encryptorPk, err = ckks.NewEncryptor(tc.params, tc.pk); err != nil {
		return
	}

	if tc.encryptorSk, err = ckks.NewEncryptor(tc.params, tc.sk); err != nil {
		return
	}

	if tc.decryptor, err = ckks.NewDecryptor(tc.params, tc.sk); err != nil {
		return
	}

	rlk, err := tc.kgen.GenRelinearizationKeyNew(tc.sk)
	if err != nil {
		return nil, err
	}

	tc.evaluator = ckks.NewEvaluator(tc.params, rlwe.NewMemEvaluationKeySet(rlk))

	return tc, nil

}

func newTestVectors(tc *ckksTestContext, encryptor *rlwe.Encryptor, a, b complex128, t *testing.T) (values []*bignum.Complex, pt *rlwe.Plaintext, ct *rlwe.Ciphertext) {

	var err error

	prec := tc.encoder.Prec()

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

func verifyCKKSTestVectors(params ckks.Parameters, encoder *ckks.Encoder, decryptor *rlwe.Decryptor, valuesWant, valuesHave interface{}, noise ring.DistributionParameters, t *testing.T) {

	precStats := ckks.GetPrecisionStats(params, encoder, decryptor, valuesWant, valuesHave, noise, false)

	if *printPrecisionStats {
		t.Log(precStats.String())
	}

	rf64, _ := precStats.MeanPrecision.Real.Float64()
	if64, _ := precStats.MeanPrecision.Imag.Float64()

	minPrec := math.Log2(params.DefaultScale().Float64()) - float64(params.LogN()+2)
	if minPrec < 0 {
		minPrec = 0
	}

	require.GreaterOrEqual(t, rf64, minPrec)
	require.GreaterOrEqual(t, if64, minPrec)
}

func testCKKSLinearTransformation(tc *ckksTestContext, t *testing.T) {

	t.Run(GetCKKSTestName(tc.params, "Average"), func(t *testing.T) {

		values, _, ciphertext := newTestVectors(tc, tc.encryptorSk, -1-1i, 1+1i, t)

		slots := ciphertext.Slots()

		logBatch := 9
		batch := 1 << logBatch
		n := slots / batch

		gks, err := tc.kgen.GenGaloisKeysNew(rlwe.GaloisElementsForInnerSum(tc.params, batch, n), tc.sk)
		require.NoError(t, err)
		evk := rlwe.NewMemEvaluationKeySet(nil, gks...)

		eval := tc.evaluator.WithKey(evk)

		eval.Average(ciphertext, logBatch, ciphertext)

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

		verifyCKKSTestVectors(tc.params, tc.encoder, tc.decryptor, values, ciphertext, nil, t)
	})

	t.Run(GetCKKSTestName(tc.params, "LinearTransform/BSGS=True"), func(t *testing.T) {

		params := tc.params

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
			DiagonalsIndexList:       nonZeroDiags,
			Level:                    ciphertext.Level(),
			Scale:                    rlwe.NewScale(params.Q()[ciphertext.Level()]),
			LogDimensions:            ciphertext.LogDimensions,
			LogBabyStepGianStepRatio: 1,
		}

		// Allocate the linear transformation
		linTransf := NewLinearTransformation(params, ltparams)

		// Encode on the linear transformation
		require.NoError(t, EncodeFloatLinearTransformation[*bignum.Complex](ltparams, tc.encoder, diagonals, linTransf))

		galEls := GaloisElementsForLinearTransformation(params, ltparams)

		gks, err := tc.kgen.GenGaloisKeysNew(galEls, tc.sk)
		require.NoError(t, err)
		evk := rlwe.NewMemEvaluationKeySet(nil, gks...)

		ltEval := NewEvaluator(tc.evaluator.WithKey(evk))

		require.NoError(t, ltEval.LinearTransformation(ciphertext, []*rlwe.Ciphertext{ciphertext}, linTransf))

		tmp := make([]*bignum.Complex, len(values))
		for i := range tmp {
			tmp[i] = values[i].Clone()
		}

		for i := 0; i < slots; i++ {
			values[i].Add(values[i], tmp[(i-15+slots)%slots])
			values[i].Add(values[i], tmp[(i-4+slots)%slots])
			values[i].Add(values[i], tmp[(i-1+slots)%slots])
			values[i].Add(values[i], tmp[(i+1)%slots])
			values[i].Add(values[i], tmp[(i+2)%slots])
			values[i].Add(values[i], tmp[(i+3)%slots])
			values[i].Add(values[i], tmp[(i+4)%slots])
			values[i].Add(values[i], tmp[(i+15)%slots])
		}

		verifyCKKSTestVectors(tc.params, tc.encoder, tc.decryptor, values, ciphertext, nil, t)
	})

	t.Run(GetCKKSTestName(tc.params, "LinearTransform/BSGS=False"), func(t *testing.T) {

		params := tc.params

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
			DiagonalsIndexList:       nonZeroDiags,
			Level:                    ciphertext.Level(),
			Scale:                    rlwe.NewScale(params.Q()[ciphertext.Level()]),
			LogDimensions:            ciphertext.LogDimensions,
			LogBabyStepGianStepRatio: -1,
		}

		// Allocate the linear transformation
		linTransf := NewLinearTransformation(params, ltparams)

		// Encode on the linear transformation
		require.NoError(t, EncodeFloatLinearTransformation[*bignum.Complex](ltparams, tc.encoder, diagonals, linTransf))

		galEls := GaloisElementsForLinearTransformation(params, ltparams)

		gks, err := tc.kgen.GenGaloisKeysNew(galEls, tc.sk)
		require.NoError(t, err)
		evk := rlwe.NewMemEvaluationKeySet(nil, gks...)

		ltEval := NewEvaluator(tc.evaluator.WithKey(evk))

		require.NoError(t, ltEval.LinearTransformation(ciphertext, []*rlwe.Ciphertext{ciphertext}, linTransf))

		tmp := make([]*bignum.Complex, len(values))
		for i := range tmp {
			tmp[i] = values[i].Clone()
		}

		for i := 0; i < slots; i++ {
			values[i].Add(values[i], tmp[(i-15+slots)%slots])
			values[i].Add(values[i], tmp[(i-4+slots)%slots])
			values[i].Add(values[i], tmp[(i-1+slots)%slots])
			values[i].Add(values[i], tmp[(i+1)%slots])
			values[i].Add(values[i], tmp[(i+2)%slots])
			values[i].Add(values[i], tmp[(i+3)%slots])
			values[i].Add(values[i], tmp[(i+4)%slots])
			values[i].Add(values[i], tmp[(i+15)%slots])
		}

		verifyCKKSTestVectors(tc.params, tc.encoder, tc.decryptor, values, ciphertext, nil, t)
	})
}

var (
	testPrec45 = ckks.ParametersLiteral{
		LogN: 10,
		Q: []uint64{
			0x80000000080001,
			0x2000000a0001,
			0x2000000e0001,
			0x2000001d0001,
			0x1fffffcf0001,
			0x1fffffc20001,
			0x200000440001,
		},
		P: []uint64{
			0x80000000130001,
			0x7fffffffe90001,
		},
		LogDefaultScale: 45,
	}

	testPrec90 = ckks.ParametersLiteral{
		LogN: 10,
		Q: []uint64{
			0x80000000080001,
			0x80000000440001,
			0x2000000a0001,
			0x2000000e0001,
			0x1fffffc20001,
			0x200000440001,
			0x200000500001,
			0x200000620001,
			0x1fffff980001,
			0x2000006a0001,
			0x1fffff7e0001,
			0x200000860001,
		},
		P: []uint64{
			0xffffffffffc0001,
			0x10000000006e0001,
		},
		LogDefaultScale: 90,
	}

	testParamsLiteral = []ckks.ParametersLiteral{testPrec45, testPrec90}
)
