package bgv

import (
	"encoding/json"
	"flag"
	"fmt"
	"math"
	"math/big"
	"runtime"
	"testing"

	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/rlwe"
	"github.com/tuneinsight/lattigo/v4/utils"

	"github.com/stretchr/testify/require"
	"github.com/tuneinsight/lattigo/v4/utils/bignum"
	"github.com/tuneinsight/lattigo/v4/utils/sampling"
)

var flagPrintNoise = flag.Bool("print-noise", false, "print the residual noise")
var flagParamString = flag.String("params", "", "specify the test cryptographic parameters as a JSON string. Overrides -short.")

func GetTestName(opname string, p Parameters, lvl int) string {
	return fmt.Sprintf("%s/LogN=%d/logQ=%d/logP=%d/LogSlots=%dx%d/logT=%d/Qi=%d/Pi=%d/lvl=%d",
		opname,
		p.LogN(),
		int(math.Round(p.LogQ())),
		int(math.Round(p.LogP())),
		p.PlaintextLogDimensions()[0],
		p.PlaintextLogDimensions()[1],
		int(math.Round(p.LogT())),
		p.QCount(),
		p.PCount(),
		lvl)
}

func TestBGV(t *testing.T) {

	var err error

	paramsLiterals := testParams

	if *flagParamString != "" {
		var jsonParams ParametersLiteral
		if err = json.Unmarshal([]byte(*flagParamString), &jsonParams); err != nil {
			t.Fatal(err)
		}
		paramsLiterals = []ParametersLiteral{jsonParams} // the custom test suite reads the parameters from the -params flag
	}

	for _, p := range paramsLiterals[:] {

		for _, plaintextModulus := range testPlaintextModulus[:] {

			p.T = plaintextModulus

			var params Parameters
			if params, err = NewParametersFromLiteral(p); err != nil {
				t.Error(err)
				t.Fail()
			}

			var tc *testContext
			if tc, err = genTestParams(params); err != nil {
				t.Error(err)
				t.Fail()
			}

			for _, testSet := range []func(tc *testContext, t *testing.T){
				testParameters,
				testEncoder,
				testEvaluator,
				testLinearTransformation,
			} {
				testSet(tc, t)
				runtime.GC()
			}
		}
	}
}

type testContext struct {
	params      Parameters
	ringQ       *ring.Ring
	ringT       *ring.Ring
	prng        sampling.PRNG
	uSampler    *ring.UniformSampler
	encoder     *Encoder
	kgen        *rlwe.KeyGenerator
	sk          *rlwe.SecretKey
	pk          *rlwe.PublicKey
	encryptorPk *rlwe.Encryptor
	encryptorSk *rlwe.Encryptor
	decryptor   *rlwe.Decryptor
	evaluator   *Evaluator
	testLevel   []int
}

func genTestParams(params Parameters) (tc *testContext, err error) {

	tc = new(testContext)
	tc.params = params

	if tc.prng, err = sampling.NewPRNG(); err != nil {
		return nil, err
	}

	tc.ringQ = params.RingQ()
	tc.ringT = params.RingT()

	tc.uSampler = ring.NewUniformSampler(tc.prng, tc.ringT)
	tc.kgen = NewKeyGenerator(tc.params)
	tc.sk, tc.pk = tc.kgen.GenKeyPairNew()
	tc.encoder = NewEncoder(tc.params)

	if tc.encryptorPk, err = NewEncryptor(tc.params, tc.pk); err != nil {
		return
	}

	if tc.encryptorSk, err = NewEncryptor(tc.params, tc.sk); err != nil {
		return
	}

	if tc.decryptor, err = NewDecryptor(tc.params, tc.sk); err != nil {
		return
	}

	var rlk *rlwe.RelinearizationKey
	if rlk, err = tc.kgen.GenRelinearizationKeyNew(tc.sk); err != nil {
		return
	}

	evk := rlwe.NewMemEvaluationKeySet(rlk)
	tc.evaluator = NewEvaluator(tc.params, evk)

	tc.testLevel = []int{0, params.MaxLevel()}

	return
}

func newTestVectorsLvl(level int, scale rlwe.Scale, tc *testContext, encryptor *rlwe.Encryptor) (coeffs ring.Poly, plaintext *rlwe.Plaintext, ciphertext *rlwe.Ciphertext) {
	coeffs = tc.uSampler.ReadNew()
	for i := range coeffs.Coeffs[0] {
		coeffs.Coeffs[0][i] = uint64(i)
	}

	plaintext = NewPlaintext(tc.params, level)
	plaintext.PlaintextScale = scale
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

func verifyTestVectors(tc *testContext, decryptor *rlwe.Decryptor, coeffs ring.Poly, element rlwe.OperandInterface[ring.Poly], t *testing.T) {

	coeffsTest := make([]uint64, tc.params.PlaintextSlots())

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

	require.True(t, utils.EqualSlice(coeffs.Coeffs[0], coeffsTest))
}

func testParameters(tc *testContext, t *testing.T) {
	t.Run(GetTestName("Parameters/Binary", tc.params, 0), func(t *testing.T) {

		bytes, err := tc.params.MarshalBinary()
		require.Nil(t, err)
		var p Parameters
		require.Nil(t, p.UnmarshalBinary(bytes))
		require.True(t, tc.params.Equal(p))

	})

	t.Run(GetTestName("Parameters/JSON", tc.params, 0), func(t *testing.T) {
		// checks that parameters can be marshalled without error
		data, err := json.Marshal(tc.params)
		require.Nil(t, err)
		require.NotNil(t, data)

		// checks that ckks.Parameters can be unmarshalled without error
		var paramsRec Parameters
		err = json.Unmarshal(data, &paramsRec)
		require.Nil(t, err)
		require.True(t, tc.params.Equal(paramsRec))

		// checks that ckks.Parameters can be unmarshalled with log-moduli definition without error
		dataWithLogModuli := []byte(fmt.Sprintf(`{"LogN":%d,"LogQ":[50,50],"LogP":[60], "T":65537}`, tc.params.LogN()))
		var paramsWithLogModuli Parameters
		err = json.Unmarshal(dataWithLogModuli, &paramsWithLogModuli)
		require.Nil(t, err)
		require.Equal(t, 2, paramsWithLogModuli.QCount())
		require.Equal(t, 1, paramsWithLogModuli.PCount())
		require.Equal(t, rlwe.DefaultXe, paramsWithLogModuli.Xe()) // Omitting Xe should result in Default being used
		require.Equal(t, rlwe.DefaultXs, paramsWithLogModuli.Xs()) // Omitting Xe should result in Default being used

		// checks that one can provide custom parameters for the secret-key and error distributions
		dataWithCustomSecrets := []byte(fmt.Sprintf(`{"LogN":%d,"LogQ":[50,50],"LogP":[60], "T":65537, "Xs": {"Type": "Ternary", "H": 192}, "Xe": {"Type": "DiscreteGaussian", "Sigma": 6.6, "Bound": 39.6}}`, tc.params.LogN()))
		var paramsWithCustomSecrets Parameters
		err = json.Unmarshal(dataWithCustomSecrets, &paramsWithCustomSecrets)
		require.Nil(t, err)
		require.Equal(t, ring.DiscreteGaussian{Sigma: 6.6, Bound: 39.6}, paramsWithCustomSecrets.Xe())
		require.Equal(t, ring.Ternary{H: 192}, paramsWithCustomSecrets.Xs())
	})
}

func testEncoder(tc *testContext, t *testing.T) {

	for _, lvl := range tc.testLevel {
		t.Run(GetTestName("Encoder/Uint", tc.params, lvl), func(t *testing.T) {
			values, plaintext, _ := newTestVectorsLvl(lvl, tc.params.PlaintextScale(), tc, nil)
			verifyTestVectors(tc, nil, values, plaintext, t)
		})
	}

	for _, lvl := range tc.testLevel {
		t.Run(GetTestName("Encoder/Int", tc.params, lvl), func(t *testing.T) {

			T := tc.params.T()
			THalf := T >> 1
			coeffs := tc.uSampler.ReadNew()
			coeffsInt := make([]int64, coeffs.N())
			for i, c := range coeffs.Coeffs[0] {
				c %= T
				if c >= THalf {
					coeffsInt[i] = -int64(T - c)
				} else {
					coeffsInt[i] = int64(c)
				}
			}

			plaintext := NewPlaintext(tc.params, lvl)
			tc.encoder.Encode(coeffsInt, plaintext)
			have := make([]int64, tc.params.PlaintextSlots())
			tc.encoder.Decode(plaintext, have)
			require.True(t, utils.EqualSlice(coeffsInt, have))
		})
	}
}

func testEvaluator(tc *testContext, t *testing.T) {

	for _, lvl := range tc.testLevel {
		t.Run(GetTestName("Evaluator/Add/Ct/Ct/New", tc.params, lvl), func(t *testing.T) {

			values0, _, ciphertext0 := newTestVectorsLvl(lvl, tc.params.NewScale(3), tc, tc.encryptorSk)
			values1, _, ciphertext1 := newTestVectorsLvl(lvl, tc.params.NewScale(7), tc, tc.encryptorSk)

			require.True(t, ciphertext0.PlaintextScale.Cmp(ciphertext1.PlaintextScale) != 0)

			ciphertext2, err := tc.evaluator.AddNew(ciphertext0, ciphertext1)
			require.NoError(t, err)
			tc.ringT.Add(values0, values1, values0)

			verifyTestVectors(tc, tc.decryptor, values0, ciphertext2, t)

		})
	}

	for _, lvl := range tc.testLevel {
		t.Run(GetTestName("Evaluator/Add/Ct/Ct/Inplace", tc.params, lvl), func(t *testing.T) {

			values0, _, ciphertext0 := newTestVectorsLvl(lvl, tc.params.NewScale(3), tc, tc.encryptorSk)
			values1, _, ciphertext1 := newTestVectorsLvl(lvl, tc.params.NewScale(7), tc, tc.encryptorSk)

			require.True(t, ciphertext0.PlaintextScale.Cmp(ciphertext1.PlaintextScale) != 0)

			require.NoError(t, tc.evaluator.Add(ciphertext0, ciphertext1, ciphertext0))
			tc.ringT.Add(values0, values1, values0)

			verifyTestVectors(tc, tc.decryptor, values0, ciphertext0, t)

		})
	}

	for _, lvl := range tc.testLevel {
		t.Run(GetTestName("Evaluator/Add/Ct/Pt/Inplace", tc.params, lvl), func(t *testing.T) {

			values0, _, ciphertext0 := newTestVectorsLvl(lvl, tc.params.NewScale(3), tc, tc.encryptorSk)
			values1, plaintext, _ := newTestVectorsLvl(lvl, tc.params.NewScale(7), tc, tc.encryptorSk)

			require.True(t, ciphertext0.PlaintextScale.Cmp(plaintext.PlaintextScale) != 0)

			require.NoError(t, tc.evaluator.Add(ciphertext0, plaintext, ciphertext0))
			tc.ringT.Add(values0, values1, values0)

			verifyTestVectors(tc, tc.decryptor, values0, ciphertext0, t)

		})
	}

	for _, lvl := range tc.testLevel {
		t.Run(GetTestName("Evaluator/Add/Ct/Scalar/Inplace", tc.params, lvl), func(t *testing.T) {

			values, _, ciphertext := newTestVectorsLvl(lvl, tc.params.PlaintextScale(), tc, tc.encryptorSk)

			scalar := tc.params.T() >> 1

			require.NoError(t, tc.evaluator.Add(ciphertext, scalar, ciphertext))
			tc.ringT.AddScalar(values, scalar, values)

			verifyTestVectors(tc, tc.decryptor, values, ciphertext, t)

		})
	}

	for _, lvl := range tc.testLevel {
		t.Run(GetTestName("Evaluator/Add/Ct/Vector/Inplace", tc.params, lvl), func(t *testing.T) {

			values, _, ciphertext := newTestVectorsLvl(lvl, tc.params.PlaintextScale(), tc, tc.encryptorSk)

			require.NoError(t, tc.evaluator.Add(ciphertext, values.Coeffs[0], ciphertext))
			tc.ringT.Add(values, values, values)

			verifyTestVectors(tc, tc.decryptor, values, ciphertext, t)

		})
	}

	for _, lvl := range tc.testLevel {
		t.Run(GetTestName("Evaluator/Sub/Ct/Ct/New", tc.params, lvl), func(t *testing.T) {

			values0, _, ciphertext0 := newTestVectorsLvl(lvl, tc.params.NewScale(3), tc, tc.encryptorSk)
			values1, _, ciphertext1 := newTestVectorsLvl(lvl, tc.params.NewScale(7), tc, tc.encryptorSk)

			require.True(t, ciphertext0.PlaintextScale.Cmp(ciphertext1.PlaintextScale) != 0)

			ciphertext0, err := tc.evaluator.SubNew(ciphertext0, ciphertext1)
			require.NoError(t, err)
			tc.ringT.Sub(values0, values1, values0)

			verifyTestVectors(tc, tc.decryptor, values0, ciphertext0, t)

		})
	}

	for _, lvl := range tc.testLevel {
		t.Run(GetTestName("Evaluator/Sub/Ct/Ct/Inplace", tc.params, lvl), func(t *testing.T) {

			values0, _, ciphertext0 := newTestVectorsLvl(lvl, tc.params.NewScale(3), tc, tc.encryptorSk)
			values1, _, ciphertext1 := newTestVectorsLvl(lvl, tc.params.NewScale(7), tc, tc.encryptorSk)

			require.True(t, ciphertext0.PlaintextScale.Cmp(ciphertext1.PlaintextScale) != 0)

			require.NoError(t, tc.evaluator.Sub(ciphertext0, ciphertext1, ciphertext0))
			tc.ringT.Sub(values0, values1, values0)

			verifyTestVectors(tc, tc.decryptor, values0, ciphertext0, t)

		})
	}

	for _, lvl := range tc.testLevel {
		t.Run(GetTestName("Evaluator/Sub/Ct/Pt/Inplace", tc.params, lvl), func(t *testing.T) {

			values0, _, ciphertext0 := newTestVectorsLvl(lvl, tc.params.NewScale(3), tc, tc.encryptorSk)
			values1, plaintext, _ := newTestVectorsLvl(lvl, tc.params.NewScale(7), tc, tc.encryptorSk)

			require.True(t, ciphertext0.PlaintextScale.Cmp(plaintext.PlaintextScale) != 0)

			require.NoError(t, tc.evaluator.Sub(ciphertext0, plaintext, ciphertext0))
			tc.ringT.Sub(values0, values1, values0)

			verifyTestVectors(tc, tc.decryptor, values0, ciphertext0, t)

		})
	}

	for _, lvl := range tc.testLevel {
		t.Run(GetTestName("Evaluator/Sub/Ct/Scalar/Inplace", tc.params, lvl), func(t *testing.T) {

			values, _, ciphertext := newTestVectorsLvl(lvl, tc.params.PlaintextScale(), tc, tc.encryptorSk)

			scalar := tc.params.T() >> 1

			require.NoError(t, tc.evaluator.Sub(ciphertext, scalar, ciphertext))
			tc.ringT.SubScalar(values, scalar, values)

			verifyTestVectors(tc, tc.decryptor, values, ciphertext, t)

		})
	}

	for _, lvl := range tc.testLevel {
		t.Run(GetTestName("Evaluator/Sub/Ct/Vector/Inplace", tc.params, lvl), func(t *testing.T) {

			values, _, ciphertext := newTestVectorsLvl(lvl, tc.params.PlaintextScale(), tc, tc.encryptorSk)

			require.NoError(t, tc.evaluator.Sub(ciphertext, values.Coeffs[0], ciphertext))
			tc.ringT.Sub(values, values, values)

			verifyTestVectors(tc, tc.decryptor, values, ciphertext, t)

		})
	}

	for _, lvl := range tc.testLevel {
		t.Run(GetTestName("Evaluator/Mul/Ct/Ct/Inplace", tc.params, lvl), func(t *testing.T) {

			if lvl == 0 {
				t.Skip("Level = 0")
			}

			values0, _, ciphertext0 := newTestVectorsLvl(lvl, tc.params.NewScale(3), tc, tc.encryptorSk)
			values1, _, ciphertext1 := newTestVectorsLvl(lvl, tc.params.NewScale(7), tc, tc.encryptorSk)

			require.True(t, ciphertext0.PlaintextScale.Cmp(ciphertext1.PlaintextScale) != 0)

			require.NoError(t, tc.evaluator.Mul(ciphertext0, ciphertext1, ciphertext0))
			tc.ringT.MulCoeffsBarrett(values0, values1, values0)

			verifyTestVectors(tc, tc.decryptor, values0, ciphertext0, t)

		})
	}

	for _, lvl := range tc.testLevel {
		t.Run(GetTestName("Evaluator/Mul/Ct/Pt/Inplace", tc.params, lvl), func(t *testing.T) {

			if lvl == 0 {
				t.Skip("Level = 0")
			}

			values0, _, ciphertext0 := newTestVectorsLvl(lvl, tc.params.NewScale(3), tc, tc.encryptorSk)
			values1, plaintext, _ := newTestVectorsLvl(lvl, tc.params.NewScale(7), tc, tc.encryptorSk)

			require.True(t, ciphertext0.PlaintextScale.Cmp(plaintext.PlaintextScale) != 0)

			require.NoError(t, tc.evaluator.Mul(ciphertext0, plaintext, ciphertext0))
			tc.ringT.MulCoeffsBarrett(values0, values1, values0)

			verifyTestVectors(tc, tc.decryptor, values0, ciphertext0, t)

		})
	}

	for _, lvl := range tc.testLevel {
		t.Run(GetTestName("Evaluator/Mul/Ct/Scalar/Inplace", tc.params, lvl), func(t *testing.T) {

			if lvl == 0 {
				t.Skip("Level = 0")
			}

			values, _, ciphertext := newTestVectorsLvl(lvl, tc.params.PlaintextScale(), tc, tc.encryptorSk)

			scalar := tc.params.T() >> 1

			require.NoError(t, tc.evaluator.Mul(ciphertext, scalar, ciphertext))
			tc.ringT.MulScalar(values, scalar, values)

			verifyTestVectors(tc, tc.decryptor, values, ciphertext, t)
		})
	}

	for _, lvl := range tc.testLevel {
		t.Run(GetTestName("Evaluator/Mul/Ct/Vector/Inplace", tc.params, lvl), func(t *testing.T) {

			if lvl == 0 {
				t.Skip("Level = 0")
			}

			values, _, ciphertext := newTestVectorsLvl(lvl, tc.params.PlaintextScale(), tc, tc.encryptorSk)

			require.NoError(t, tc.evaluator.Mul(ciphertext, values.Coeffs[0], ciphertext))
			tc.ringT.MulCoeffsBarrett(values, values, values)

			verifyTestVectors(tc, tc.decryptor, values, ciphertext, t)
		})
	}

	for _, lvl := range tc.testLevel {
		t.Run(GetTestName("Evaluator/Square/Ct/Ct/Inplace", tc.params, lvl), func(t *testing.T) {

			if lvl == 0 {
				t.Skip("Level = 0")
			}

			values0, _, ciphertext0 := newTestVectorsLvl(lvl, tc.params.NewScale(3), tc, tc.encryptorSk)

			require.NoError(t, tc.evaluator.Mul(ciphertext0, ciphertext0, ciphertext0))
			tc.ringT.MulCoeffsBarrett(values0, values0, values0)

			verifyTestVectors(tc, tc.decryptor, values0, ciphertext0, t)
		})
	}

	for _, lvl := range tc.testLevel {
		t.Run(GetTestName("Evaluator/MulRelin/Ct/Ct/Inplace", tc.params, lvl), func(t *testing.T) {

			if lvl == 0 {
				t.Skip("Level = 0")
			}

			values0, _, ciphertext0 := newTestVectorsLvl(lvl, tc.params.NewScale(3), tc, tc.encryptorSk)
			values1, _, ciphertext1 := newTestVectorsLvl(lvl, tc.params.NewScale(7), tc, tc.encryptorSk)

			tc.ringT.MulCoeffsBarrett(values0, values1, values0)

			require.True(t, ciphertext0.PlaintextScale.Cmp(ciphertext1.PlaintextScale) != 0)

			receiver := NewCiphertext(tc.params, 1, lvl)

			require.NoError(t, tc.evaluator.MulRelin(ciphertext0, ciphertext1, receiver))

			require.NoError(t, tc.evaluator.Rescale(receiver, receiver))

			verifyTestVectors(tc, tc.decryptor, values0, receiver, t)
		})
	}

	for _, lvl := range tc.testLevel {
		t.Run(GetTestName("Evaluator/MulThenAdd/Ct/Ct/Inplace", tc.params, lvl), func(t *testing.T) {

			if lvl == 0 {
				t.Skip("Level = 0")
			}

			values0, _, ciphertext0 := newTestVectorsLvl(lvl, tc.params.PlaintextScale(), tc, tc.encryptorSk)
			values1, _, ciphertext1 := newTestVectorsLvl(lvl, rlwe.NewScale(2), tc, tc.encryptorSk)
			values2, _, ciphertext2 := newTestVectorsLvl(lvl, tc.params.NewScale(7), tc, tc.encryptorSk)

			require.True(t, ciphertext0.PlaintextScale.Cmp(ciphertext1.PlaintextScale) != 0)
			require.True(t, ciphertext0.PlaintextScale.Cmp(ciphertext2.PlaintextScale) != 0)

			require.NoError(t, tc.evaluator.MulThenAdd(ciphertext0, ciphertext1, ciphertext2))
			tc.ringT.MulCoeffsBarrettThenAdd(values0, values1, values2)

			verifyTestVectors(tc, tc.decryptor, values2, ciphertext2, t)
		})
	}

	for _, lvl := range tc.testLevel {
		t.Run(GetTestName("Evaluator/MulThenAdd/Ct/Pt/Inplace", tc.params, lvl), func(t *testing.T) {

			if lvl == 0 {
				t.Skip("Level = 0")
			}

			values0, _, ciphertext0 := newTestVectorsLvl(lvl, tc.params.PlaintextScale(), tc, tc.encryptorSk)
			values1, plaintext1, _ := newTestVectorsLvl(lvl, rlwe.NewScale(2), tc, tc.encryptorSk)
			values2, _, ciphertext2 := newTestVectorsLvl(lvl, tc.params.NewScale(7), tc, tc.encryptorSk)

			require.True(t, ciphertext0.PlaintextScale.Cmp(plaintext1.PlaintextScale) != 0)
			require.True(t, ciphertext0.PlaintextScale.Cmp(ciphertext2.PlaintextScale) != 0)

			require.NoError(t, tc.evaluator.MulThenAdd(ciphertext0, plaintext1, ciphertext2))
			tc.ringT.MulCoeffsBarrettThenAdd(values0, values1, values2)

			verifyTestVectors(tc, tc.decryptor, values2, ciphertext2, t)
		})
	}

	for _, lvl := range tc.testLevel {
		t.Run(GetTestName("Evaluator/MulThenAdd/Ct/Scalar/Inplace", tc.params, lvl), func(t *testing.T) {

			if lvl == 0 {
				t.Skip("Level = 0")
			}

			values0, _, ciphertext0 := newTestVectorsLvl(lvl, tc.params.NewScale(7), tc, tc.encryptorSk)
			values1, _, ciphertext1 := newTestVectorsLvl(lvl, tc.params.NewScale(3), tc, tc.encryptorSk)

			require.True(t, ciphertext0.PlaintextScale.Cmp(ciphertext1.PlaintextScale) != 0)

			scalar := tc.params.T() >> 1

			require.NoError(t, tc.evaluator.MulThenAdd(ciphertext0, scalar, ciphertext1))
			tc.ringT.MulScalarThenAdd(values0, scalar, values1)

			verifyTestVectors(tc, tc.decryptor, values1, ciphertext1, t)
		})
	}

	for _, lvl := range tc.testLevel {
		t.Run(GetTestName("Evaluator/MulThenAdd/Ct/Vector/Inplace", tc.params, lvl), func(t *testing.T) {

			if lvl == 0 {
				t.Skip("Level = 0")
			}

			values0, _, ciphertext0 := newTestVectorsLvl(lvl, tc.params.NewScale(7), tc, tc.encryptorSk)
			values1, _, ciphertext1 := newTestVectorsLvl(lvl, tc.params.NewScale(3), tc, tc.encryptorSk)

			require.True(t, ciphertext0.PlaintextScale.Cmp(ciphertext1.PlaintextScale) != 0)

			scale := ciphertext1.PlaintextScale

			require.NoError(t, tc.evaluator.MulThenAdd(ciphertext0, values1.Coeffs[0], ciphertext1))
			tc.ringT.MulCoeffsBarrettThenAdd(values0, values1, values1)

			// Checks that output scale isn't changed
			require.True(t, scale.Equal(ciphertext1.PlaintextScale))

			verifyTestVectors(tc, tc.decryptor, values1, ciphertext1, t)
		})
	}

	for _, lvl := range tc.testLevel {
		t.Run(GetTestName("Evaluator/MulRelinThenAdd/Ct/Ct/Inplace", tc.params, lvl), func(t *testing.T) {

			if lvl == 0 {
				t.Skip("Level = 0")
			}

			values0, _, ciphertext0 := newTestVectorsLvl(lvl, tc.params.PlaintextScale(), tc, tc.encryptorSk)
			values1, _, ciphertext1 := newTestVectorsLvl(lvl, rlwe.NewScale(2), tc, tc.encryptorSk)
			values2, _, ciphertext2 := newTestVectorsLvl(lvl, tc.params.NewScale(7), tc, tc.encryptorSk)

			require.True(t, ciphertext0.PlaintextScale.Cmp(ciphertext1.PlaintextScale) != 0)
			require.True(t, ciphertext0.PlaintextScale.Cmp(ciphertext2.PlaintextScale) != 0)

			require.NoError(t, tc.evaluator.MulRelinThenAdd(ciphertext0, ciphertext1, ciphertext2))
			tc.ringT.MulCoeffsBarrettThenAdd(values0, values1, values2)

			verifyTestVectors(tc, tc.decryptor, values2, ciphertext2, t)
		})
	}

	t.Run("Evaluator/PolyEval", func(t *testing.T) {

		t.Run("Single", func(t *testing.T) {

			if tc.params.MaxLevel() < 4 {
				t.Skip("MaxLevel() to low")
			}

			values, _, ciphertext := newTestVectorsLvl(tc.params.MaxLevel(), tc.params.NewScale(1), tc, tc.encryptorSk)

			coeffs := []uint64{0, 0, 1}

			T := tc.params.T()
			for i := range values.Coeffs[0] {
				values.Coeffs[0][i] = ring.EvalPolyModP(values.Coeffs[0][i], coeffs, T)
			}

			poly := bignum.NewPolynomial(bignum.Monomial, coeffs, nil)

			t.Run(GetTestName("Standard", tc.params, tc.params.MaxLevel()), func(t *testing.T) {

				res, err := tc.evaluator.Polynomial(ciphertext, poly, false, tc.params.PlaintextScale())
				require.NoError(t, err)

				require.True(t, res.PlaintextScale.Cmp(tc.params.PlaintextScale()) == 0)

				verifyTestVectors(tc, tc.decryptor, values, res, t)
			})

			t.Run(GetTestName("Invariant", tc.params, tc.params.MaxLevel()), func(t *testing.T) {

				res, err := tc.evaluator.Polynomial(ciphertext, poly, true, tc.params.PlaintextScale())
				require.NoError(t, err)

				require.True(t, res.PlaintextScale.Cmp(tc.params.PlaintextScale()) == 0)

				verifyTestVectors(tc, tc.decryptor, values, res, t)
			})
		})

		t.Run("Vector", func(t *testing.T) {

			if tc.params.MaxLevel() < 4 {
				t.Skip("MaxLevel() to low")
			}

			values, _, ciphertext := newTestVectorsLvl(tc.params.MaxLevel(), tc.params.NewScale(7), tc, tc.encryptorSk)

			coeffs0 := []uint64{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16}
			coeffs1 := []uint64{2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17}

			slots := values.N()

			slotIndex := make(map[int][]int)
			idx0 := make([]int, slots>>1)
			idx1 := make([]int, slots>>1)
			for i := 0; i < slots>>1; i++ {
				idx0[i] = 2 * i
				idx1[i] = 2*i + 1
			}

			slotIndex[0] = idx0
			slotIndex[1] = idx1

			polyVector, err := rlwe.NewPolynomialVector([]rlwe.Polynomial{
				rlwe.NewPolynomial(bignum.NewPolynomial(bignum.Monomial, coeffs0, nil)),
				rlwe.NewPolynomial(bignum.NewPolynomial(bignum.Monomial, coeffs1, nil)),
			}, slotIndex)
			require.NoError(t, err)

			TInt := new(big.Int).SetUint64(tc.params.T())
			for pol, idx := range slotIndex {
				for _, i := range idx {
					values.Coeffs[0][i] = polyVector.Value[pol].EvaluateModP(new(big.Int).SetUint64(values.Coeffs[0][i]), TInt).Uint64()
				}
			}

			t.Run(GetTestName("Standard", tc.params, tc.params.MaxLevel()), func(t *testing.T) {

				res, err := tc.evaluator.Polynomial(ciphertext, polyVector, false, tc.params.PlaintextScale())
				require.NoError(t, err)

				require.True(t, res.PlaintextScale.Cmp(tc.params.PlaintextScale()) == 0)

				verifyTestVectors(tc, tc.decryptor, values, res, t)
			})

			t.Run(GetTestName("Invariant", tc.params, tc.params.MaxLevel()), func(t *testing.T) {

				res, err := tc.evaluator.Polynomial(ciphertext, polyVector, true, tc.params.PlaintextScale())
				require.NoError(t, err)

				require.True(t, res.PlaintextScale.Cmp(tc.params.PlaintextScale()) == 0)

				verifyTestVectors(tc, tc.decryptor, values, res, t)
			})
		})
	})

	for _, lvl := range tc.testLevel[:] {
		t.Run(GetTestName("Evaluator/Rescale", tc.params, lvl), func(t *testing.T) {

			ringT := tc.params.RingT()

			values0, _, ciphertext0 := newTestVectorsLvl(lvl, tc.params.PlaintextScale(), tc, tc.encryptorPk)

			printNoise := func(msg string, values []uint64, ct *rlwe.Ciphertext) {
				pt := NewPlaintext(tc.params, ct.Level())
				pt.MetaData = ciphertext0.MetaData
				require.NoError(t, tc.encoder.Encode(values0.Coeffs[0], pt))
				ct, err := tc.evaluator.SubNew(ct, pt)
				require.NoError(t, err)
				vartmp, _, _ := rlwe.Norm(ct, tc.decryptor)
				t.Logf("STD(noise) %s: %f\n", msg, vartmp)
			}

			if lvl != 0 {

				values1, _, ciphertext1 := newTestVectorsLvl(lvl, tc.params.PlaintextScale(), tc, tc.encryptorSk)

				if *flagPrintNoise {
					printNoise("0x", values0.Coeffs[0], ciphertext0)
				}

				for i := 0; i < lvl; i++ {
					tc.evaluator.MulRelin(ciphertext0, ciphertext1, ciphertext0)

					ringT.MulCoeffsBarrett(values0, values1, values0)

					if *flagPrintNoise {
						printNoise(fmt.Sprintf("%dx", i+1), values0.Coeffs[0], ciphertext0)
					}

				}

				verifyTestVectors(tc, tc.decryptor, values0, ciphertext0, t)

				require.Nil(t, tc.evaluator.Rescale(ciphertext0, ciphertext0))

				verifyTestVectors(tc, tc.decryptor, values0, ciphertext0, t)

			} else {
				require.NotNil(t, tc.evaluator.Rescale(ciphertext0, ciphertext0))
			}
		})
	}
}

func testLinearTransformation(tc *testContext, t *testing.T) {

	level := tc.params.MaxLevel()
	t.Run(GetTestName("Evaluator/LinearTransformationBSGS=true", tc.params, level), func(t *testing.T) {

		params := tc.params

		values, _, ciphertext := newTestVectorsLvl(level, tc.params.PlaintextScale(), tc, tc.encryptorSk)

		diagonals := make(map[int][]uint64)

		totSlots := values.N()

		diagonals[-15] = make([]uint64, totSlots)
		diagonals[-4] = make([]uint64, totSlots)
		diagonals[-1] = make([]uint64, totSlots)
		diagonals[0] = make([]uint64, totSlots)
		diagonals[1] = make([]uint64, totSlots)
		diagonals[2] = make([]uint64, totSlots)
		diagonals[3] = make([]uint64, totSlots)
		diagonals[4] = make([]uint64, totSlots)
		diagonals[15] = make([]uint64, totSlots)

		for i := 0; i < totSlots; i++ {
			diagonals[-15][i] = 1
			diagonals[-4][i] = 1
			diagonals[-1][i] = 1
			diagonals[0][i] = 1
			diagonals[1][i] = 1
			diagonals[2][i] = 1
			diagonals[3][i] = 1
			diagonals[4][i] = 1
			diagonals[15][i] = 1
		}

		ltparams := rlwe.MemLinearTransformationParameters[uint64]{
			Diagonals:                diagonals,
			Level:                    ciphertext.Level(),
			PlaintextScale:           tc.params.PlaintextScale(),
			PlaintextLogDimensions:   ciphertext.PlaintextLogDimensions,
			LogBabyStepGianStepRatio: 1,
		}

		// Allocate the linear transformation
		linTransf := NewLinearTransformation[uint64](params, ltparams)

		// Encode on the linear transformation
		require.NoError(t, EncodeLinearTransformation[uint64](linTransf, ltparams, tc.encoder))

		galEls := rlwe.GaloisElementsForLinearTransformation[uint64](params, ltparams)

		gks, err := tc.kgen.GenGaloisKeysNew(galEls, tc.sk)
		require.NoError(t, err)
		eval := tc.evaluator.WithKey(rlwe.NewMemEvaluationKeySet(nil, gks...))

		require.NoError(t, eval.LinearTransformation(ciphertext, linTransf, []*rlwe.Ciphertext{ciphertext}))

		tmp := make([]uint64, totSlots)
		copy(tmp, values.Coeffs[0])

		subRing := tc.params.RingT().SubRings[0]

		subRing.Add(values.Coeffs[0], utils.RotateSlotsNew(tmp, -15), values.Coeffs[0])
		subRing.Add(values.Coeffs[0], utils.RotateSlotsNew(tmp, -4), values.Coeffs[0])
		subRing.Add(values.Coeffs[0], utils.RotateSlotsNew(tmp, -1), values.Coeffs[0])
		subRing.Add(values.Coeffs[0], utils.RotateSlotsNew(tmp, 1), values.Coeffs[0])
		subRing.Add(values.Coeffs[0], utils.RotateSlotsNew(tmp, 2), values.Coeffs[0])
		subRing.Add(values.Coeffs[0], utils.RotateSlotsNew(tmp, 3), values.Coeffs[0])
		subRing.Add(values.Coeffs[0], utils.RotateSlotsNew(tmp, 4), values.Coeffs[0])
		subRing.Add(values.Coeffs[0], utils.RotateSlotsNew(tmp, 15), values.Coeffs[0])

		verifyTestVectors(tc, tc.decryptor, values, ciphertext, t)
	})

	t.Run(GetTestName("Evaluator/LinearTransformationBSGS=false", tc.params, level), func(t *testing.T) {

		params := tc.params

		values, _, ciphertext := newTestVectorsLvl(level, tc.params.PlaintextScale(), tc, tc.encryptorSk)

		diagonals := make(map[int][]uint64)

		totSlots := values.N()

		diagonals[-15] = make([]uint64, totSlots)
		diagonals[-4] = make([]uint64, totSlots)
		diagonals[-1] = make([]uint64, totSlots)
		diagonals[0] = make([]uint64, totSlots)
		diagonals[1] = make([]uint64, totSlots)
		diagonals[2] = make([]uint64, totSlots)
		diagonals[3] = make([]uint64, totSlots)
		diagonals[4] = make([]uint64, totSlots)
		diagonals[15] = make([]uint64, totSlots)

		for i := 0; i < totSlots; i++ {
			diagonals[-15][i] = 1
			diagonals[-4][i] = 1
			diagonals[-1][i] = 1
			diagonals[0][i] = 1
			diagonals[1][i] = 1
			diagonals[2][i] = 1
			diagonals[3][i] = 1
			diagonals[4][i] = 1
			diagonals[15][i] = 1
		}

		ltparams := rlwe.MemLinearTransformationParameters[uint64]{
			Diagonals:                diagonals,
			Level:                    ciphertext.Level(),
			PlaintextScale:           tc.params.PlaintextScale(),
			PlaintextLogDimensions:   ciphertext.PlaintextLogDimensions,
			LogBabyStepGianStepRatio: -1,
		}

		// Allocate the linear transformation
		linTransf := NewLinearTransformation[uint64](params, ltparams)

		// Encode on the linear transformation
		require.NoError(t, EncodeLinearTransformation[uint64](linTransf, ltparams, tc.encoder))

		galEls := rlwe.GaloisElementsForLinearTransformation[uint64](params, ltparams)

		gks, err := tc.kgen.GenGaloisKeysNew(galEls, tc.sk)
		require.NoError(t, err)
		eval := tc.evaluator.WithKey(rlwe.NewMemEvaluationKeySet(nil, gks...))

		require.NoError(t, eval.LinearTransformation(ciphertext, linTransf, []*rlwe.Ciphertext{ciphertext}))

		tmp := make([]uint64, totSlots)
		copy(tmp, values.Coeffs[0])

		subRing := tc.params.RingT().SubRings[0]

		subRing.Add(values.Coeffs[0], utils.RotateSlotsNew(tmp, -15), values.Coeffs[0])
		subRing.Add(values.Coeffs[0], utils.RotateSlotsNew(tmp, -4), values.Coeffs[0])
		subRing.Add(values.Coeffs[0], utils.RotateSlotsNew(tmp, -1), values.Coeffs[0])
		subRing.Add(values.Coeffs[0], utils.RotateSlotsNew(tmp, 1), values.Coeffs[0])
		subRing.Add(values.Coeffs[0], utils.RotateSlotsNew(tmp, 2), values.Coeffs[0])
		subRing.Add(values.Coeffs[0], utils.RotateSlotsNew(tmp, 3), values.Coeffs[0])
		subRing.Add(values.Coeffs[0], utils.RotateSlotsNew(tmp, 4), values.Coeffs[0])
		subRing.Add(values.Coeffs[0], utils.RotateSlotsNew(tmp, 15), values.Coeffs[0])

		verifyTestVectors(tc, tc.decryptor, values, ciphertext, t)
	})
}
