package polyint

import (
	"encoding/json"
	"flag"
	"fmt"
	"math"
	"runtime"
	"slices"
	"testing"

	"github.com/stretchr/testify/require"

	"github.com/tuneinsight/lattigo/v5/core/rlwe"
	"github.com/tuneinsight/lattigo/v5/ring"
	"github.com/tuneinsight/lattigo/v5/schemes"
	"github.com/tuneinsight/lattigo/v5/schemes/bgv"
	"github.com/tuneinsight/lattigo/v5/utils/bignum"
	"github.com/tuneinsight/lattigo/v5/utils/sampling"
)

type testContext struct {
	params      bgv.Parameters
	ringQ       *ring.Ring
	ringT       *ring.Ring
	prng        sampling.PRNG
	uSampler    *ring.UniformSampler
	encoder     *bgv.Encoder
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
	t.Run("Single", func(t *testing.T) {

		if tc.params.MaxLevel() < 4 {
			t.Skip("MaxLevel() to low")
		}

		values, _, ciphertext := newTestVectorsLvl(tc.params.MaxLevel(), tc.params.NewScale(1), tc, tc.encryptorSk)

		coeffs := []uint64{0, 0, 1}

		T := tc.params.PlaintextModulus()
		for i := range values.Coeffs[0] {
			values.Coeffs[0][i] = ring.EvalPolyModP(values.Coeffs[0][i], coeffs, T)
		}

		poly := bignum.NewPolynomial(bignum.Monomial, coeffs, nil)

		t.Run(GetTestName("Standard", tc.params, tc.params.MaxLevel()), func(t *testing.T) {

			polyEval := NewPolynomialEvaluator(tc.params, tc.evaluator, false)

			res, err := polyEval.Evaluate(ciphertext, poly, tc.params.DefaultScale())
			require.NoError(t, err)

			require.Equal(t, res.Scale.Cmp(tc.params.DefaultScale()), 0)

			verifyTestVectors(tc, tc.decryptor, values, res, t)
		})

		t.Run(GetTestName("Invariant", tc.params, tc.params.MaxLevel()), func(t *testing.T) {

			polyEval := NewPolynomialEvaluator(tc.params, tc.evaluator, true)

			res, err := polyEval.Evaluate(ciphertext, poly, tc.params.DefaultScale())
			require.NoError(t, err)

			require.Equal(t, res.Level(), ciphertext.Level())
			require.Equal(t, res.Scale.Cmp(tc.params.DefaultScale()), 0)

			verifyTestVectors(tc, tc.decryptor, values, res, t)
		})
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
