package bgv

import (
	"encoding/json"
	"flag"
	"fmt"
	"runtime"
	"testing"

	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/rlwe"
	"github.com/tuneinsight/lattigo/v4/utils"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
	"github.com/tuneinsight/lattigo/v4/utils/sampling"
)

var flagPrintNoise = flag.Bool("print-noise", false, "print the residual noise")
var flagParamString = flag.String("params", "", "specify the test cryptographic parameters as a JSON string. Overrides -short.")

var (
	// TESTN13QP218 is a of 128-bit secure test parameters set with a 32-bit plaintext and depth 4.
	TESTN14QP418 = ParametersLiteral{
		LogN: 13,
		Q:    []uint64{0x3fffffa8001, 0x1000090001, 0x10000c8001, 0x10000f0001, 0xffff00001},
		P:    []uint64{0x7fffffd8001},
		T:    0xffc001,
	}

	// TestParams is a set of test parameters for BGV ensuring 128 bit security in the classic setting.
	TestParams = []ParametersLiteral{TESTN14QP418}
)

func GetTestName(opname string, p Parameters, lvl int) string {
	return fmt.Sprintf("%s/LogN=%d/logQ=%d/logP=%d/logT=%d/Qi=%d/Pi=%d/lvl=%d", opname, p.LogN(), p.LogQ(), p.LogP(), p.LogT(), p.QCount(), p.PCount(), lvl)
}

func TestBGV(t *testing.T) {

	var err error

	paramsLiterals := TestParams

	if *flagParamString != "" {
		var jsonParams ParametersLiteral
		if err = json.Unmarshal([]byte(*flagParamString), &jsonParams); err != nil {
			t.Fatal(err)
		}
		paramsLiterals = []ParametersLiteral{jsonParams} // the custom test suite reads the parameters from the -params flag
	}

	for _, p := range paramsLiterals[:] {

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
			testLinearTransform,
			testMarshalling,
		} {
			testSet(tc, t)
			runtime.GC()
		}
	}
}

type testContext struct {
	params      Parameters
	ringQ       *ring.Ring
	ringT       *ring.Ring
	prng        sampling.PRNG
	uSampler    *ring.UniformSampler
	encoder     Encoder
	kgen        *rlwe.KeyGenerator
	sk          *rlwe.SecretKey
	pk          *rlwe.PublicKey
	encryptorPk rlwe.Encryptor
	encryptorSk rlwe.Encryptor
	decryptor   rlwe.Decryptor
	evaluator   Evaluator
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
	tc.encryptorPk = NewEncryptor(tc.params, tc.pk)
	tc.encryptorSk = NewEncryptor(tc.params, tc.sk)
	tc.decryptor = NewDecryptor(tc.params, tc.sk)
	evk := rlwe.NewEvaluationKeySet()
	evk.RelinearizationKey = tc.kgen.GenRelinearizationKeyNew(tc.sk)
	tc.evaluator = NewEvaluator(tc.params, evk)

	tc.testLevel = []int{0, params.MaxLevel()}

	return
}

func newTestVectorsLvl(level int, scale rlwe.Scale, tc *testContext, encryptor rlwe.Encryptor) (coeffs *ring.Poly, plaintext *rlwe.Plaintext, ciphertext *rlwe.Ciphertext) {
	coeffs = tc.uSampler.ReadNew()
	for i := range coeffs.Coeffs[0] {
		coeffs.Coeffs[0][i] = uint64(i)
	}
	plaintext = NewPlaintext(tc.params, level)
	plaintext.Scale = scale
	tc.encoder.Encode(coeffs.Coeffs[0], plaintext)
	if encryptor != nil {
		ciphertext = encryptor.EncryptNew(plaintext)
	}

	return coeffs, plaintext, ciphertext
}

func verifyTestVectors(tc *testContext, decryptor rlwe.Decryptor, coeffs *ring.Poly, element rlwe.Operand, t *testing.T) {

	var coeffsTest []uint64

	switch el := element.(type) {
	case *rlwe.Plaintext:
		coeffsTest = tc.encoder.DecodeUintNew(el)
	case *rlwe.Ciphertext:

		pt := decryptor.DecryptNew(el)

		coeffsTest = tc.encoder.DecodeUintNew(pt)

		if *flagPrintNoise {
			tc.encoder.Encode(coeffsTest, pt)
			vartmp, _, _ := rlwe.Norm(tc.evaluator.SubNew(el, pt), decryptor)
			t.Logf("STD(noise): %f\n", vartmp)
		}

	default:
		t.Error("invalid test object to verify")
	}

	require.True(t, utils.EqualSlice(coeffs.Coeffs[0], coeffsTest))
}

func testParameters(tc *testContext, t *testing.T) {

	t.Run("Parameters/NewParameters", func(t *testing.T) {
		params, err := NewParametersFromLiteral(ParametersLiteral{
			LogN: 4,
			LogQ: []int{60, 60},
			LogP: []int{60},
			T:    0x10001,
		})
		require.NoError(t, err)
		require.Equal(t, ring.Standard, params.RingType())  // Default ring type should be standard
		require.Equal(t, rlwe.DefaultSigma, params.Sigma()) // Default error std should be rlwe.DefaultSigma
	})

	t.Run("Parameters/CopyNew", func(t *testing.T) {
		params1, params2 := tc.params.CopyNew(), tc.params.CopyNew()
		require.True(t, params1.Equal(tc.params) && params2.Equal(tc.params))
		params1.ringT, _ = ring.NewRing(params1.N(), []uint64{0x40002001})
		require.False(t, params1.Equal(tc.params))
		require.True(t, params2.Equal(tc.params))
	})
}

func testEncoder(tc *testContext, t *testing.T) {

	for _, lvl := range tc.testLevel {
		t.Run(GetTestName("Encoder/Encode&Decode/Uint", tc.params, lvl), func(t *testing.T) {
			values, plaintext, _ := newTestVectorsLvl(lvl, tc.params.DefaultScale(), tc, nil)
			verifyTestVectors(tc, nil, values, plaintext, t)
		})
	}

	for _, lvl := range tc.testLevel {
		t.Run(GetTestName("Encoder/Encode&Decode/Int", tc.params, lvl), func(t *testing.T) {

			T := tc.params.T()
			THalf := T >> 1
			coeffs := tc.uSampler.ReadNew()
			coeffsInt := make([]int64, len(coeffs.Coeffs[0]))
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
			require.True(t, utils.EqualSlice(coeffsInt, tc.encoder.DecodeIntNew(plaintext)))
		})
	}
}

func testEvaluator(tc *testContext, t *testing.T) {

	t.Run("Evaluator", func(t *testing.T) {

		for _, lvl := range tc.testLevel {
			t.Run(GetTestName("AddNew/op0=ct/op2=ct", tc.params, lvl), func(t *testing.T) {

				values0, _, ciphertext0 := newTestVectorsLvl(lvl, tc.params.NewScale(3), tc, tc.encryptorSk)
				values1, _, ciphertext1 := newTestVectorsLvl(lvl, tc.params.NewScale(7), tc, tc.encryptorSk)

				require.True(t, ciphertext0.Scale.Cmp(ciphertext1.Scale) != 0)

				ciphertext2 := tc.evaluator.AddNew(ciphertext0, ciphertext1)
				tc.ringT.Add(values0, values1, values0)

				verifyTestVectors(tc, tc.decryptor, values0, ciphertext2, t)

			})
		}

		for _, lvl := range tc.testLevel {
			t.Run(GetTestName("Add/op0=ct/op2=ct", tc.params, lvl), func(t *testing.T) {

				values0, _, ciphertext0 := newTestVectorsLvl(lvl, tc.params.NewScale(3), tc, tc.encryptorSk)
				values1, _, ciphertext1 := newTestVectorsLvl(lvl, tc.params.NewScale(7), tc, tc.encryptorSk)

				require.True(t, ciphertext0.Scale.Cmp(ciphertext1.Scale) != 0)

				tc.evaluator.Add(ciphertext0, ciphertext1, ciphertext0)
				tc.ringT.Add(values0, values1, values0)

				verifyTestVectors(tc, tc.decryptor, values0, ciphertext0, t)

			})
		}

		for _, lvl := range tc.testLevel {
			t.Run(GetTestName("Add/op0=ct/op2=pt", tc.params, lvl), func(t *testing.T) {

				values0, _, ciphertext0 := newTestVectorsLvl(lvl, tc.params.NewScale(3), tc, tc.encryptorSk)
				values1, plaintext, _ := newTestVectorsLvl(lvl, tc.params.NewScale(7), tc, tc.encryptorSk)

				require.True(t, ciphertext0.Scale.Cmp(plaintext.Scale) != 0)

				tc.evaluator.Add(ciphertext0, plaintext, ciphertext0)
				tc.ringT.Add(values0, values1, values0)

				verifyTestVectors(tc, tc.decryptor, values0, ciphertext0, t)

			})
		}

		for _, lvl := range tc.testLevel {
			t.Run(GetTestName("Add/op0=ct/op2=ct", tc.params, lvl), func(t *testing.T) {

				values0, _, ciphertext0 := newTestVectorsLvl(lvl, tc.params.NewScale(3), tc, tc.encryptorSk)
				values1, _, ciphertext1 := newTestVectorsLvl(lvl, tc.params.NewScale(7), tc, tc.encryptorSk)

				require.True(t, ciphertext0.Scale.Cmp(ciphertext1.Scale) != 0)

				tc.evaluator.Add(ciphertext0, ciphertext1, ciphertext0)
				tc.ringT.Add(values0, values1, values0)

				verifyTestVectors(tc, tc.decryptor, values0, ciphertext0, t)

			})
		}

		for _, lvl := range tc.testLevel {
			t.Run(GetTestName("SubNew/op0=ct/op2=ct", tc.params, lvl), func(t *testing.T) {

				values0, _, ciphertext0 := newTestVectorsLvl(lvl, tc.params.NewScale(3), tc, tc.encryptorSk)
				values1, _, ciphertext1 := newTestVectorsLvl(lvl, tc.params.NewScale(7), tc, tc.encryptorSk)

				require.True(t, ciphertext0.Scale.Cmp(ciphertext1.Scale) != 0)

				ciphertext0 = tc.evaluator.SubNew(ciphertext0, ciphertext1)
				tc.ringT.Sub(values0, values1, values0)

				verifyTestVectors(tc, tc.decryptor, values0, ciphertext0, t)

			})
		}

		for _, lvl := range tc.testLevel {
			t.Run(GetTestName("Sub/op0=ct/op2=ct", tc.params, lvl), func(t *testing.T) {

				values0, _, ciphertext0 := newTestVectorsLvl(lvl, tc.params.NewScale(3), tc, tc.encryptorSk)
				values1, _, ciphertext1 := newTestVectorsLvl(lvl, tc.params.NewScale(7), tc, tc.encryptorSk)

				require.True(t, ciphertext0.Scale.Cmp(ciphertext1.Scale) != 0)

				tc.evaluator.Sub(ciphertext0, ciphertext1, ciphertext0)
				tc.ringT.Sub(values0, values1, values0)

				verifyTestVectors(tc, tc.decryptor, values0, ciphertext0, t)

			})
		}

		for _, lvl := range tc.testLevel {
			t.Run(GetTestName("Sub/op0=ct/op2=pt", tc.params, lvl), func(t *testing.T) {

				values0, _, ciphertext0 := newTestVectorsLvl(lvl, tc.params.NewScale(3), tc, tc.encryptorSk)
				values1, plaintext, _ := newTestVectorsLvl(lvl, tc.params.NewScale(7), tc, tc.encryptorSk)

				require.True(t, ciphertext0.Scale.Cmp(plaintext.Scale) != 0)

				tc.evaluator.Sub(ciphertext0, plaintext, ciphertext0)
				tc.ringT.Sub(values0, values1, values0)

				verifyTestVectors(tc, tc.decryptor, values0, ciphertext0, t)

			})
		}

		for _, lvl := range tc.testLevel {
			t.Run(GetTestName("Sub/op0=ct/op2=ct", tc.params, lvl), func(t *testing.T) {

				values0, _, ciphertext0 := newTestVectorsLvl(lvl, tc.params.NewScale(3), tc, tc.encryptorSk)
				values1, _, ciphertext1 := newTestVectorsLvl(lvl, tc.params.NewScale(7), tc, tc.encryptorSk)

				require.True(t, ciphertext0.Scale.Cmp(ciphertext1.Scale) != 0)

				tc.evaluator.Sub(ciphertext0, ciphertext1, ciphertext0)
				tc.ringT.Sub(values0, values1, values0)

				verifyTestVectors(tc, tc.decryptor, values0, ciphertext0, t)

			})
		}

		for _, lvl := range tc.testLevel {
			t.Run(GetTestName("Neg/op0=ct", tc.params, lvl), func(t *testing.T) {

				values, _, ciphertext := newTestVectorsLvl(lvl, tc.params.DefaultScale(), tc, tc.encryptorSk)

				tc.evaluator.Neg(ciphertext, ciphertext)
				tc.ringT.Neg(values, values)
				tc.ringT.Reduce(values, values)

				verifyTestVectors(tc, tc.decryptor, values, ciphertext, t)

			})
		}

		for _, lvl := range tc.testLevel {
			t.Run(GetTestName("NegNew/op0=ct", tc.params, lvl), func(t *testing.T) {

				values, _, ciphertext := newTestVectorsLvl(lvl, tc.params.DefaultScale(), tc, tc.encryptorSk)

				ciphertext = tc.evaluator.NegNew(ciphertext)
				tc.ringT.Neg(values, values)
				tc.ringT.Reduce(values, values)

				verifyTestVectors(tc, tc.decryptor, values, ciphertext, t)

			})
		}

		for _, lvl := range tc.testLevel {
			t.Run(GetTestName("AddScalar/op0=ct", tc.params, lvl), func(t *testing.T) {

				values, _, ciphertext := newTestVectorsLvl(lvl, tc.params.DefaultScale(), tc, tc.encryptorSk)

				scalar := tc.params.T() >> 1

				tc.evaluator.AddScalar(ciphertext, scalar, ciphertext)
				tc.ringT.AddScalar(values, scalar, values)

				verifyTestVectors(tc, tc.decryptor, values, ciphertext, t)

			})
		}

		for _, lvl := range tc.testLevel {
			t.Run(GetTestName("AddScalarNew/op0=ct", tc.params, lvl), func(t *testing.T) {

				values, _, ciphertext := newTestVectorsLvl(lvl, tc.params.NewScale(7), tc, tc.encryptorSk)

				require.True(t, ciphertext.Scale.Cmp(rlwe.NewScale(1)) != 0)

				scalar := tc.params.T() >> 1

				ciphertext = tc.evaluator.AddScalarNew(ciphertext, scalar)
				tc.ringT.AddScalar(values, scalar, values)

				verifyTestVectors(tc, tc.decryptor, values, ciphertext, t)
			})
		}

		for _, lvl := range tc.testLevel {
			t.Run(GetTestName("MulScalar/op0=ct", tc.params, lvl), func(t *testing.T) {

				if lvl == 0 {
					t.Skip("Level = 0")
				}

				values, _, ciphertext := newTestVectorsLvl(lvl, tc.params.DefaultScale(), tc, tc.encryptorSk)

				scalar := tc.params.T() >> 1

				tc.evaluator.MulScalar(ciphertext, scalar, ciphertext)
				tc.ringT.MulScalar(values, scalar, values)

				verifyTestVectors(tc, tc.decryptor, values, ciphertext, t)

			})
		}

		for _, lvl := range tc.testLevel {
			t.Run(GetTestName("MulScalarNew/op0=ct", tc.params, lvl), func(t *testing.T) {

				if lvl == 0 {
					t.Skip("Level = 0")
				}

				values, _, ciphertext := newTestVectorsLvl(lvl, tc.params.DefaultScale(), tc, tc.encryptorSk)

				scalar := tc.params.T() >> 1

				ciphertext = tc.evaluator.MulScalarNew(ciphertext, scalar)
				tc.ringT.MulScalar(values, scalar, values)

				verifyTestVectors(tc, tc.decryptor, values, ciphertext, t)
			})
		}

		for _, lvl := range tc.testLevel {
			t.Run(GetTestName("MulScalarThenAdd/op0=ct", tc.params, lvl), func(t *testing.T) {

				if lvl == 0 {
					t.Skip("Level = 0")
				}

				values0, _, ciphertext0 := newTestVectorsLvl(lvl, tc.params.NewScale(7), tc, tc.encryptorSk)
				values1, _, ciphertext1 := newTestVectorsLvl(lvl, tc.params.NewScale(3), tc, tc.encryptorSk)

				require.True(t, ciphertext0.Scale.Cmp(ciphertext1.Scale) != 0)

				scalar := tc.params.T() >> 1

				tc.evaluator.MulScalarThenAdd(ciphertext0, scalar, ciphertext1)
				tc.ringT.MulScalarThenAdd(values0, scalar, values1)

				verifyTestVectors(tc, tc.decryptor, values1, ciphertext1, t)
			})
		}

		for _, lvl := range tc.testLevel {
			t.Run(GetTestName("Mul/op0=ct/op2=ct", tc.params, lvl), func(t *testing.T) {

				if lvl == 0 {
					t.Skip("Level = 0")
				}

				values0, _, ciphertext0 := newTestVectorsLvl(lvl, tc.params.NewScale(3), tc, tc.encryptorSk)
				values1, _, ciphertext1 := newTestVectorsLvl(lvl, tc.params.NewScale(7), tc, tc.encryptorSk)

				require.True(t, ciphertext0.Scale.Cmp(ciphertext1.Scale) != 0)

				tc.evaluator.Mul(ciphertext0, ciphertext1, ciphertext0)
				tc.ringT.MulCoeffsBarrett(values0, values1, values0)

				verifyTestVectors(tc, tc.decryptor, values0, ciphertext0, t)

			})
		}

		for _, lvl := range tc.testLevel {
			t.Run(GetTestName("Mul/op0=ct/op2=pt", tc.params, lvl), func(t *testing.T) {

				if lvl == 0 {
					t.Skip("Level = 0")
				}

				values0, _, ciphertext0 := newTestVectorsLvl(lvl, tc.params.NewScale(3), tc, tc.encryptorSk)
				values1, plaintext, _ := newTestVectorsLvl(lvl, tc.params.NewScale(7), tc, tc.encryptorSk)

				require.True(t, ciphertext0.Scale.Cmp(plaintext.Scale) != 0)

				tc.evaluator.Mul(ciphertext0, plaintext, ciphertext0)
				tc.ringT.MulCoeffsBarrett(values0, values1, values0)

				verifyTestVectors(tc, tc.decryptor, values0, ciphertext0, t)

			})
		}

		for _, lvl := range tc.testLevel {
			t.Run(GetTestName("Square/op0=ct/op2=ct", tc.params, lvl), func(t *testing.T) {

				if lvl == 0 {
					t.Skip("Level = 0")
				}

				values0, _, ciphertext0 := newTestVectorsLvl(lvl, tc.params.NewScale(3), tc, tc.encryptorSk)

				tc.evaluator.Mul(ciphertext0, ciphertext0, ciphertext0)
				tc.ringT.MulCoeffsBarrett(values0, values0, values0)

				verifyTestVectors(tc, tc.decryptor, values0, ciphertext0, t)

			})
		}

		for _, lvl := range tc.testLevel {
			t.Run(GetTestName("MulRelin/op0=ct/op2=ct", tc.params, lvl), func(t *testing.T) {

				if lvl == 0 {
					t.Skip("Level = 0")
				}

				values0, _, ciphertext0 := newTestVectorsLvl(lvl, tc.params.NewScale(3), tc, tc.encryptorSk)
				values1, _, ciphertext1 := newTestVectorsLvl(lvl, tc.params.NewScale(7), tc, tc.encryptorSk)

				tc.ringT.MulCoeffsBarrett(values0, values1, values0)

				require.True(t, ciphertext0.Scale.Cmp(ciphertext1.Scale) != 0)

				receiver := NewCiphertext(tc.params, 1, lvl)

				tc.evaluator.MulRelin(ciphertext0, ciphertext1, receiver)

				tc.evaluator.Rescale(receiver, receiver)

				verifyTestVectors(tc, tc.decryptor, values0, receiver, t)

			})
		}

		for _, lvl := range tc.testLevel {
			t.Run(GetTestName("MulThenAdd/op0=ct/op2=ct", tc.params, lvl), func(t *testing.T) {

				if lvl == 0 {
					t.Skip("Level = 0")
				}

				values0, _, ciphertext0 := newTestVectorsLvl(lvl, tc.params.DefaultScale(), tc, tc.encryptorSk)
				values1, _, ciphertext1 := newTestVectorsLvl(lvl, rlwe.NewScale(2), tc, tc.encryptorSk)
				values2, _, ciphertext2 := newTestVectorsLvl(lvl, tc.params.NewScale(7), tc, tc.encryptorSk)

				require.True(t, ciphertext0.Scale.Cmp(ciphertext1.Scale) != 0)
				require.True(t, ciphertext0.Scale.Cmp(ciphertext2.Scale) != 0)

				tc.evaluator.MulThenAdd(ciphertext0, ciphertext1, ciphertext2)
				tc.ringT.MulCoeffsBarrettThenAdd(values0, values1, values2)

				verifyTestVectors(tc, tc.decryptor, values2, ciphertext2, t)

			})
		}

		for _, lvl := range tc.testLevel {
			t.Run(GetTestName("MulRelinThenAdd/op0=ct/op2=ct", tc.params, lvl), func(t *testing.T) {

				if lvl == 0 {
					t.Skip("Level = 0")
				}

				values0, _, ciphertext0 := newTestVectorsLvl(lvl, tc.params.DefaultScale(), tc, tc.encryptorSk)
				values1, _, ciphertext1 := newTestVectorsLvl(lvl, rlwe.NewScale(2), tc, tc.encryptorSk)
				values2, _, ciphertext2 := newTestVectorsLvl(lvl, tc.params.NewScale(7), tc, tc.encryptorSk)

				require.True(t, ciphertext0.Scale.Cmp(ciphertext1.Scale) != 0)
				require.True(t, ciphertext0.Scale.Cmp(ciphertext2.Scale) != 0)

				tc.evaluator.MulRelinThenAdd(ciphertext0, ciphertext1, ciphertext2)
				tc.ringT.MulCoeffsBarrettThenAdd(values0, values1, values2)

				verifyTestVectors(tc, tc.decryptor, values2, ciphertext2, t)

			})
		}

		for _, lvl := range tc.testLevel {
			t.Run(GetTestName("MulThenAdd/op0=ct/op1=pt", tc.params, lvl), func(t *testing.T) {

				if lvl == 0 {
					t.Skip("Level = 0")
				}

				values0, _, ciphertext0 := newTestVectorsLvl(lvl, tc.params.DefaultScale(), tc, tc.encryptorSk)
				values1, plaintext1, _ := newTestVectorsLvl(lvl, rlwe.NewScale(2), tc, tc.encryptorSk)
				values2, _, ciphertext2 := newTestVectorsLvl(lvl, tc.params.NewScale(7), tc, tc.encryptorSk)

				require.True(t, ciphertext0.Scale.Cmp(plaintext1.Scale) != 0)
				require.True(t, ciphertext0.Scale.Cmp(ciphertext2.Scale) != 0)

				tc.evaluator.MulThenAdd(ciphertext0, plaintext1, ciphertext2)
				tc.ringT.MulCoeffsBarrettThenAdd(values0, values1, values2)

				verifyTestVectors(tc, tc.decryptor, values2, ciphertext2, t)

			})
		}

		t.Run(GetTestName("PolyEval/Single", tc.params, tc.params.MaxLevel()), func(t *testing.T) {

			if tc.params.MaxLevel() < 4 {
				t.Skip("MaxLevel() to low")
			}

			values, _, ciphertext := newTestVectorsLvl(tc.params.MaxLevel(), tc.params.NewScale(7), tc, tc.encryptorSk)

			coeffs := []uint64{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16}

			T := tc.params.T()
			for i := range values.Coeffs[0] {
				values.Coeffs[0][i] = ring.EvalPolyModP(values.Coeffs[0][i], coeffs, T)
			}

			poly := NewPoly(coeffs)

			var err error
			if ciphertext, err = tc.evaluator.EvaluatePoly(ciphertext, poly, tc.params.DefaultScale()); err != nil {
				t.Log(err)
				t.Fatal()
			}

			require.True(t, ciphertext.Scale.Cmp(tc.params.DefaultScale()) == 0)

			std, min, max := rlwe.Norm(ciphertext, tc.decryptor)
			t.Logf("Noise -> (std: %f, min: %f, max=%f)\n", std, min, max)

			verifyTestVectors(tc, tc.decryptor, values, ciphertext, t)
		})

		t.Run(GetTestName("PolyEval/Vector", tc.params, tc.params.MaxLevel()), func(t *testing.T) {

			if tc.params.MaxLevel() < 4 {
				t.Skip("MaxLevel() to low")
			}

			values, _, ciphertext := newTestVectorsLvl(tc.params.MaxLevel(), tc.params.NewScale(7), tc, tc.encryptorSk)

			coeffs0 := []uint64{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16}
			coeffs1 := []uint64{2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17}

			slotIndex := make(map[int][]int)
			idx0 := make([]int, tc.params.N()>>1)
			idx1 := make([]int, tc.params.N()>>1)
			for i := 0; i < tc.params.N()>>1; i++ {
				idx0[i] = 2 * i
				idx1[i] = 2*i + 1
			}

			polyVec := []*Polynomial{NewPoly(coeffs0), NewPoly(coeffs1)}

			slotIndex[0] = idx0
			slotIndex[1] = idx1

			T := tc.params.T()
			for pol, idx := range slotIndex {
				for _, i := range idx {
					values.Coeffs[0][i] = ring.EvalPolyModP(values.Coeffs[0][i], polyVec[pol].Coeffs, T)
				}
			}

			var err error
			if ciphertext, err = tc.evaluator.EvaluatePolyVector(ciphertext, polyVec, tc.encoder, slotIndex, tc.params.DefaultScale()); err != nil {
				t.Fail()
			}

			require.True(t, ciphertext.Scale.Cmp(tc.params.DefaultScale()) == 0)

			std, min, max := rlwe.Norm(ciphertext, tc.decryptor)
			t.Logf("Noise -> (std: %f, min: %f, max=%f)\n", std, min, max)

			verifyTestVectors(tc, tc.decryptor, values, ciphertext, t)
		})

		for _, lvl := range tc.testLevel[:] {
			t.Run(GetTestName("Rescale", tc.params, lvl), func(t *testing.T) {

				ringT := tc.params.RingT()

				values0, _, ciphertext0 := newTestVectorsLvl(lvl, tc.params.DefaultScale(), tc, tc.encryptorPk)

				printNoise := func(msg string, values []uint64, ct *rlwe.Ciphertext) {
					pt := NewPlaintext(tc.params, ct.Level())
					pt.MetaData = ciphertext0.MetaData
					tc.encoder.Encode(values0.Coeffs[0], pt)
					vartmp, _, _ := rlwe.Norm(tc.evaluator.SubNew(ct, pt), tc.decryptor)
					t.Logf("STD(noise) %s: %f\n", msg, vartmp)
				}

				if lvl != 0 {

					values1, _, ciphertext1 := newTestVectorsLvl(lvl, tc.params.DefaultScale(), tc, tc.encryptorSk)

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
	})
}

func testLinearTransform(tc *testContext, t *testing.T) {

	t.Run(GetTestName("LinearTransform/Naive", tc.params, tc.params.MaxLevel()), func(t *testing.T) {

		params := tc.params

		values, _, ciphertext := newTestVectorsLvl(tc.params.MaxLevel(), tc.params.DefaultScale(), tc, tc.encryptorSk)

		diagMatrix := make(map[int][]uint64)

		N := params.N()

		diagMatrix[-1] = make([]uint64, N)
		diagMatrix[0] = make([]uint64, N)
		diagMatrix[1] = make([]uint64, N)

		for i := 0; i < N; i++ {
			diagMatrix[-1][i] = 1
			diagMatrix[0][i] = 1
			diagMatrix[1][i] = 1
		}

		linTransf := GenLinearTransform(tc.encoder, diagMatrix, params.MaxLevel(), tc.params.DefaultScale())

		rotations := linTransf.Rotations()

		evk := rlwe.NewEvaluationKeySet()
		for _, galEl := range tc.params.GaloisElementsForRotations(rotations) {
			evk.GaloisKeys[galEl] = tc.kgen.GenGaloisKeyNew(galEl, tc.sk)
		}

		eval := tc.evaluator.WithKey(evk)

		eval.LinearTransform(ciphertext, linTransf, []*rlwe.Ciphertext{ciphertext})

		tmp := make([]uint64, params.N())
		copy(tmp, values.Coeffs[0])

		subRing := tc.params.RingT().SubRings[0]

		subRing.Add(values.Coeffs[0], utils.RotateSlotsNew(tmp, -1), values.Coeffs[0])
		subRing.Add(values.Coeffs[0], utils.RotateSlotsNew(tmp, 1), values.Coeffs[0])

		verifyTestVectors(tc, tc.decryptor, values, ciphertext, t)
	})

	t.Run(GetTestName("LinearTransform/BSGS", tc.params, tc.params.MaxLevel()), func(t *testing.T) {

		params := tc.params

		values, _, ciphertext := newTestVectorsLvl(tc.params.MaxLevel(), tc.params.DefaultScale(), tc, tc.encryptorSk)

		diagMatrix := make(map[int][]uint64)

		N := params.N()

		diagMatrix[-15] = make([]uint64, N)
		diagMatrix[-4] = make([]uint64, N)
		diagMatrix[-1] = make([]uint64, N)
		diagMatrix[0] = make([]uint64, N)
		diagMatrix[1] = make([]uint64, N)
		diagMatrix[2] = make([]uint64, N)
		diagMatrix[3] = make([]uint64, N)
		diagMatrix[4] = make([]uint64, N)
		diagMatrix[15] = make([]uint64, N)

		for i := 0; i < N; i++ {
			diagMatrix[-15][i] = 1
			diagMatrix[-4][i] = 1
			diagMatrix[-1][i] = 1
			diagMatrix[0][i] = 1
			diagMatrix[1][i] = 1
			diagMatrix[2][i] = 1
			diagMatrix[3][i] = 1
			diagMatrix[4][i] = 1
			diagMatrix[15][i] = 1
		}

		linTransf := GenLinearTransformBSGS(tc.encoder, diagMatrix, params.MaxLevel(), tc.params.DefaultScale(), 2.0)

		rotations := linTransf.Rotations()

		evk := rlwe.NewEvaluationKeySet()
		for _, galEl := range tc.params.GaloisElementsForRotations(rotations) {
			evk.GaloisKeys[galEl] = tc.kgen.GenGaloisKeyNew(galEl, tc.sk)
		}

		eval := tc.evaluator.WithKey(evk)

		eval.LinearTransform(ciphertext, linTransf, []*rlwe.Ciphertext{ciphertext})

		tmp := make([]uint64, params.N())
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

func testMarshalling(tc *testContext, t *testing.T) {
	t.Run("Marshalling", func(t *testing.T) {

		t.Run("Parameters/Binary", func(t *testing.T) {

			bytes, err := tc.params.MarshalBinary()
			assert.Nil(t, err)
			var p Parameters
			err = p.UnmarshalBinary(bytes)
			assert.Nil(t, err)
			assert.Equal(t, tc.params, p)
			assert.Equal(t, tc.params.RingQ(), p.RingQ())
		})

		t.Run("Parameters/JSON", func(t *testing.T) {
			// checks that parameters can be marshalled without error
			data, err := json.Marshal(tc.params)
			assert.Nil(t, err)
			assert.NotNil(t, data)

			// checks that ckks.Parameters can be unmarshalled without error
			var paramsRec Parameters
			err = json.Unmarshal(data, &paramsRec)
			assert.Nil(t, err)
			assert.True(t, tc.params.Equal(paramsRec))

			// checks that ckks.Parameters can be unmarshalled with log-moduli definition without error
			dataWithLogModuli := []byte(fmt.Sprintf(`{"LogN":%d,"LogQ":[50,50],"LogP":[60], "T":65537}`, tc.params.LogN()))
			var paramsWithLogModuli Parameters
			err = json.Unmarshal(dataWithLogModuli, &paramsWithLogModuli)
			assert.Nil(t, err)
			assert.Equal(t, 2, paramsWithLogModuli.QCount())
			assert.Equal(t, 1, paramsWithLogModuli.PCount())
			assert.Equal(t, rlwe.DefaultSigma, paramsWithLogModuli.Sigma()) // Omitting sigma should result in Default being used

			// checks that one can provide custom parameters for the secret-key and error distributions
			dataWithCustomSecrets := []byte(fmt.Sprintf(`{"LogN":%d,"LogQ":[50,50],"LogP":[60],"H": 192, "Sigma": 6.6, "T":65537}`, tc.params.LogN()))
			var paramsWithCustomSecrets Parameters
			err = json.Unmarshal(dataWithCustomSecrets, &paramsWithCustomSecrets)
			assert.Nil(t, err)
			assert.Equal(t, 6.6, paramsWithCustomSecrets.Sigma())
			assert.Equal(t, 192, paramsWithCustomSecrets.HammingWeight())
		})
	})
}
