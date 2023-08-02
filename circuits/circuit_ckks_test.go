package circuits

import (
	"encoding/json"
	"flag"
	"fmt"
	"math"
	"math/big"
	"math/bits"
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
				testDecryptPublic,
				testEvaluatePoly,
				testChebyshevInterpolator,
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

func newCKKSTestVectors(tc *ckksTestContext, encryptor *rlwe.Encryptor, a, b complex128, t *testing.T) (values []*bignum.Complex, pt *rlwe.Plaintext, ct *rlwe.Ciphertext) {

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

func VerifyCKKSTestVectors(params ckks.Parameters, encoder *ckks.Encoder, decryptor *rlwe.Decryptor, valuesWant, valuesHave interface{}, noise ring.DistributionParameters, printPrecisionStats bool, t *testing.T) {

	precStats := ckks.GetPrecisionStats(params, encoder, decryptor, valuesWant, valuesHave, noise, false)

	if printPrecisionStats {
		t.Log(precStats.String())
	}

	rf64, _ := precStats.MeanPrecision.Real.Float64()
	if64, _ := precStats.MeanPrecision.Imag.Float64()

	minPrec := math.Log2(params.DefaultScale().Float64())

	switch params.RingType() {
	case ring.Standard:
		minPrec -= float64(params.LogN()) + 2 // Z[X]/(X^{N} + 1)
	case ring.ConjugateInvariant:
		minPrec -= float64(params.LogN()) + 2.5 // Z[X + X^1]/(X^{2N} + 1)
	}
	if minPrec < 0 {
		minPrec = 0
	}

	require.GreaterOrEqual(t, rf64, minPrec)
	require.GreaterOrEqual(t, if64, minPrec)
}

func testCKKSLinearTransformation(tc *ckksTestContext, t *testing.T) {

	t.Run(GetCKKSTestName(tc.params, "Average"), func(t *testing.T) {

		values, _, ciphertext := newCKKSTestVectors(tc, tc.encryptorSk, -1-1i, 1+1i, t)

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

		values, _, ciphertext := newCKKSTestVectors(tc, tc.encryptorSk, -1-1i, 1+1i, t)

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

		require.NoError(t, ltEval.LinearTransformation(ciphertext, linTransf, ciphertext))

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

		values, _, ciphertext := newCKKSTestVectors(tc, tc.encryptorSk, -1-1i, 1+1i, t)

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

		require.NoError(t, ltEval.LinearTransformation(ciphertext, linTransf, ciphertext))

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

func testEvaluatePoly(tc *ckksTestContext, t *testing.T) {

	var err error

	polyEval := NewCKKSPolynomialEvaluator(tc.params, tc.evaluator)

	t.Run(GetCKKSTestName(tc.params, "EvaluatePoly/PolySingle/Exp"), func(t *testing.T) {

		if tc.params.MaxLevel() < 3 {
			t.Skip("skipping test for params max level < 3")
		}

		values, _, ciphertext := newCKKSTestVectors(tc, tc.encryptorSk, -1, 1, t)

		prec := tc.encoder.Prec()

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

		if ciphertext, err = polyEval.Polynomial(ciphertext, poly, ciphertext.Scale); err != nil {
			t.Fatal(err)
		}

		VerifyCKKSTestVectors(tc.params, tc.encoder, tc.decryptor, values, ciphertext, nil, *printPrecisionStats, t)
	})

	t.Run(GetCKKSTestName(tc.params, "Polynomial/PolyVector/Exp"), func(t *testing.T) {

		if tc.params.MaxLevel() < 3 {
			t.Skip("skipping test for params max level < 3")
		}

		values, _, ciphertext := newCKKSTestVectors(tc, tc.encryptorSk, -1, 1, t)

		prec := tc.encoder.Prec()

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

		polyVector, err := NewPolynomialVector([]Polynomial{NewPolynomial(poly)}, slotIndex)
		require.NoError(t, err)

		if ciphertext, err = polyEval.Polynomial(ciphertext, polyVector, ciphertext.Scale); err != nil {
			t.Fatal(err)
		}

		VerifyCKKSTestVectors(tc.params, tc.encoder, tc.decryptor, valuesWant, ciphertext, nil, *printPrecisionStats, t)
	})
}

func testChebyshevInterpolator(tc *ckksTestContext, t *testing.T) {

	var err error

	polyEval := NewCKKSPolynomialEvaluator(tc.params, tc.evaluator)

	t.Run(GetCKKSTestName(tc.params, "ChebyshevInterpolator/Sin"), func(t *testing.T) {

		degree := 13

		if tc.params.MaxDepth() < bits.Len64(uint64(degree)) {
			t.Skip("skipping test: not enough levels")
		}

		eval := tc.evaluator

		values, _, ciphertext := newCKKSTestVectors(tc, tc.encryptorSk, -1, 1, t)

		prec := tc.params.EncodingPrecision()

		interval := bignum.Interval{
			Nodes: degree,
			A:     *new(big.Float).SetPrec(prec).SetFloat64(-8),
			B:     *new(big.Float).SetPrec(prec).SetFloat64(8),
		}

		poly := NewPolynomial(bignum.ChebyshevApproximation(math.Sin, interval))

		scalar, constant := poly.ChangeOfBasis()
		eval.Mul(ciphertext, scalar, ciphertext)
		eval.Add(ciphertext, constant, ciphertext)
		if err = eval.RescaleTo(ciphertext, tc.params.DefaultScale(), ciphertext); err != nil {
			t.Fatal(err)
		}

		if ciphertext, err = polyEval.Polynomial(ciphertext, poly, ciphertext.Scale); err != nil {
			t.Fatal(err)
		}

		for i := range values {
			values[i] = poly.Evaluate(values[i])
		}

		VerifyCKKSTestVectors(tc.params, tc.encoder, tc.decryptor, values, ciphertext, nil, *printPrecisionStats, t)
	})
}

func testDecryptPublic(tc *ckksTestContext, t *testing.T) {

	var err error

	t.Run(GetCKKSTestName(tc.params, "DecryptPublic/Sin"), func(t *testing.T) {

		degree := 7
		a, b := -1.5, 1.5

		if tc.params.MaxDepth() < bits.Len64(uint64(degree)) {
			t.Skip("skipping test: not enough levels")
		}

		eval := tc.evaluator

		values, _, ciphertext := newCKKSTestVectors(tc, tc.encryptorSk, complex(a, 0), complex(b, 0), t)

		prec := tc.params.EncodingPrecision()

		sin := func(x *bignum.Complex) (y *bignum.Complex) {
			xf64, _ := x[0].Float64()
			y = bignum.NewComplex()
			y.SetPrec(prec)
			y[0].SetFloat64(math.Sin(xf64))
			return
		}

		interval := bignum.Interval{
			Nodes: degree,
			A:     *new(big.Float).SetPrec(prec).SetFloat64(a),
			B:     *new(big.Float).SetPrec(prec).SetFloat64(b),
		}

		poly := bignum.ChebyshevApproximation(sin, interval)

		for i := range values {
			values[i] = poly.Evaluate(values[i])
		}

		scalar, constant := poly.ChangeOfBasis()

		require.NoError(t, eval.Mul(ciphertext, scalar, ciphertext))
		require.NoError(t, eval.Add(ciphertext, constant, ciphertext))
		if err := eval.RescaleTo(ciphertext, tc.params.DefaultScale(), ciphertext); err != nil {
			t.Fatal(err)
		}

		polyEval := NewCKKSPolynomialEvaluator(tc.params, tc.evaluator)

		if ciphertext, err = polyEval.Polynomial(ciphertext, poly, ciphertext.Scale); err != nil {
			t.Fatal(err)
		}

		plaintext := tc.decryptor.DecryptNew(ciphertext)

		valuesHave := make([]*big.Float, plaintext.Slots())

		require.NoError(t, tc.encoder.Decode(plaintext, valuesHave))

		VerifyCKKSTestVectors(tc.params, tc.encoder, nil, values, valuesHave, nil, *printPrecisionStats, t)

		for i := range valuesHave {
			valuesHave[i].Sub(valuesHave[i], values[i][0])
		}

		// This should make it lose at most ~0.5 bit or precision.
		sigma := ckks.StandardDeviation(valuesHave, rlwe.NewScale(plaintext.Scale.Float64()/math.Sqrt(float64(len(values)))))

		tc.encoder.DecodePublic(plaintext, valuesHave, ring.DiscreteGaussian{Sigma: sigma, Bound: 2.5066282746310002 * sigma})

		VerifyCKKSTestVectors(tc.params, tc.encoder, nil, values, valuesHave, nil, *printPrecisionStats, t)
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
