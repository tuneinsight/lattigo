package bgv

import (
	"fmt"
	"runtime"
	"testing"

	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/rlwe"
	"github.com/tuneinsight/lattigo/v4/utils"

	"github.com/stretchr/testify/require"
)

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
	return fmt.Sprintf("%s/LogN=%d/logQ=%d/logP=%d/logT=%d/#Q=%d/#P=%d/lvl=%d", opname, p.LogN(), p.LogQ(), p.LogP(), p.LogT(), p.QCount(), p.PCount(), lvl)
}

func TestBGV(t *testing.T) {

	var err error

	for _, p := range TestParams[:] {

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
			testEncryptor,
			testEvaluator,
			testRotate,
			testInnerSum,
			testLinearTransform,
			testMerge,
			testSwitchKeys,
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
	prng        utils.PRNG
	uSampler    *ring.UniformSampler
	encoder     Encoder
	kgen        rlwe.KeyGenerator
	sk          *rlwe.SecretKey
	pk          *rlwe.PublicKey
	rlk         *rlwe.RelinearizationKey
	encryptorPk Encryptor
	encryptorSk Encryptor
	decryptor   Decryptor
	evaluator   Evaluator
	testLevel   []int
}

func genTestParams(params Parameters) (tc *testContext, err error) {

	tc = new(testContext)
	tc.params = params

	if tc.prng, err = utils.NewPRNG(); err != nil {
		return nil, err
	}

	tc.ringQ = params.RingQ()
	tc.ringT = params.RingT()

	tc.uSampler = ring.NewUniformSampler(tc.prng, tc.ringT)
	tc.kgen = NewKeyGenerator(tc.params)
	tc.sk, tc.pk = tc.kgen.GenKeyPair()
	tc.rlk = tc.kgen.GenRelinearizationKey(tc.sk, 1)
	tc.encoder = NewEncoder(tc.params)
	tc.encryptorPk = NewEncryptor(tc.params, tc.pk)
	tc.encryptorSk = NewEncryptor(tc.params, tc.sk)
	tc.decryptor = NewDecryptor(tc.params, tc.sk)
	tc.evaluator = NewEvaluator(tc.params, rlwe.EvaluationKey{Rlk: tc.rlk})

	tc.testLevel = []int{0, params.MaxLevel()}

	return
}

func newTestVectorsLvl(level int, scale rlwe.Scale, tc *testContext, encryptor Encryptor) (coeffs *ring.Poly, plaintext *rlwe.Plaintext, ciphertext *rlwe.Ciphertext) {
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

func verifyTestVectors(tc *testContext, decryptor Decryptor, coeffs *ring.Poly, element rlwe.Operand, t *testing.T) {

	var coeffsTest []uint64

	switch el := element.(type) {
	case *rlwe.Plaintext:
		coeffsTest = tc.encoder.DecodeUintNew(el)
	case *rlwe.Ciphertext:
		coeffsTest = tc.encoder.DecodeUintNew(decryptor.DecryptNew(el))
	default:
		t.Error("invalid test object to verify")
	}

	require.True(t, utils.EqualSliceUint64(coeffs.Coeffs[0], coeffsTest))
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
		require.True(t, params1.Equals(tc.params) && params2.Equals(tc.params))
		params1.ringT, _ = ring.NewRing(params1.N(), []uint64{0x40002001})
		require.False(t, params1.Equals(tc.params))
		require.True(t, params2.Equals(tc.params))
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
			require.True(t, utils.EqualSliceInt64(coeffsInt, tc.encoder.DecodeIntNew(plaintext)))
		})
	}
}

func testEncryptor(tc *testContext, t *testing.T) {

	for _, lvl := range tc.testLevel {
		t.Run(GetTestName("Encoder/EncryptorPk", tc.params, lvl), func(t *testing.T) {
			values, _, ciphertext := newTestVectorsLvl(lvl, tc.params.DefaultScale(), tc, tc.encryptorPk)
			verifyTestVectors(tc, tc.decryptor, values, ciphertext, t)
		})
	}

	for _, lvl := range tc.testLevel {
		t.Run(GetTestName("Encoder/encryptorSk", tc.params, lvl), func(t *testing.T) {
			values, _, ciphertext := newTestVectorsLvl(lvl, tc.params.DefaultScale(), tc, tc.encryptorSk)
			verifyTestVectors(tc, tc.decryptor, values, ciphertext, t)
		})
	}
}

func testEvaluator(tc *testContext, t *testing.T) {

	t.Run("Evaluator", func(t *testing.T) {

		for _, lvl := range tc.testLevel {
			t.Run(GetTestName("AddNew/op0=ct/op2=ct", tc.params, lvl), func(t *testing.T) {

				values0, _, ciphertext0 := newTestVectorsLvl(lvl, rlwe.NewScale(3), tc, tc.encryptorSk)
				values1, _, ciphertext1 := newTestVectorsLvl(lvl, rlwe.NewScale(7), tc, tc.encryptorSk)

				require.True(t, ciphertext0.Scale.Cmp(ciphertext1.Scale) != 0)

				ciphertext2 := tc.evaluator.AddNew(ciphertext0, ciphertext1)
				tc.ringT.Add(values0, values1, values0)

				verifyTestVectors(tc, tc.decryptor, values0, ciphertext2, t)

			})
		}

		for _, lvl := range tc.testLevel {
			t.Run(GetTestName("Add/op0=ct/op2=ct", tc.params, lvl), func(t *testing.T) {

				values0, _, ciphertext0 := newTestVectorsLvl(lvl, rlwe.NewScale(3), tc, tc.encryptorSk)
				values1, _, ciphertext1 := newTestVectorsLvl(lvl, rlwe.NewScale(7), tc, tc.encryptorSk)

				require.True(t, ciphertext0.Scale.Cmp(ciphertext1.Scale) != 0)

				tc.evaluator.Add(ciphertext0, ciphertext1, ciphertext0)
				tc.ringT.Add(values0, values1, values0)

				verifyTestVectors(tc, tc.decryptor, values0, ciphertext0, t)

			})
		}

		for _, lvl := range tc.testLevel {
			t.Run(GetTestName("Add/op0=ct/op2=pt", tc.params, lvl), func(t *testing.T) {

				values0, _, ciphertext0 := newTestVectorsLvl(lvl, rlwe.NewScale(3), tc, tc.encryptorSk)
				values1, plaintext, _ := newTestVectorsLvl(lvl, rlwe.NewScale(7), tc, tc.encryptorSk)

				require.True(t, ciphertext0.Scale.Cmp(plaintext.Scale) != 0)

				tc.evaluator.Add(ciphertext0, plaintext, ciphertext0)
				tc.ringT.Add(values0, values1, values0)

				verifyTestVectors(tc, tc.decryptor, values0, ciphertext0, t)

			})
		}

		for _, lvl := range tc.testLevel {
			t.Run(GetTestName("Add/op0=ct/op2=ct", tc.params, lvl), func(t *testing.T) {

				values0, _, ciphertext0 := newTestVectorsLvl(lvl, rlwe.NewScale(3), tc, tc.encryptorSk)
				values1, _, ciphertext1 := newTestVectorsLvl(lvl, rlwe.NewScale(7), tc, tc.encryptorSk)

				require.True(t, ciphertext0.Scale.Cmp(ciphertext1.Scale) != 0)

				tc.evaluator.Add(ciphertext0, ciphertext1, ciphertext0)
				tc.ringT.Add(values0, values1, values0)

				verifyTestVectors(tc, tc.decryptor, values0, ciphertext0, t)

			})
		}

		for _, lvl := range tc.testLevel {
			t.Run(GetTestName("SubNew/op0=ct/op2=ct", tc.params, lvl), func(t *testing.T) {

				values0, _, ciphertext0 := newTestVectorsLvl(lvl, rlwe.NewScale(3), tc, tc.encryptorSk)
				values1, _, ciphertext1 := newTestVectorsLvl(lvl, rlwe.NewScale(7), tc, tc.encryptorSk)

				require.True(t, ciphertext0.Scale.Cmp(ciphertext1.Scale) != 0)

				ciphertext0 = tc.evaluator.SubNew(ciphertext0, ciphertext1)
				tc.ringT.Sub(values0, values1, values0)

				verifyTestVectors(tc, tc.decryptor, values0, ciphertext0, t)

			})
		}

		for _, lvl := range tc.testLevel {
			t.Run(GetTestName("Sub/op0=ct/op2=ct", tc.params, lvl), func(t *testing.T) {

				values0, _, ciphertext0 := newTestVectorsLvl(lvl, rlwe.NewScale(3), tc, tc.encryptorSk)
				values1, _, ciphertext1 := newTestVectorsLvl(lvl, rlwe.NewScale(7), tc, tc.encryptorSk)

				require.True(t, ciphertext0.Scale.Cmp(ciphertext1.Scale) != 0)

				tc.evaluator.Sub(ciphertext0, ciphertext1, ciphertext0)
				tc.ringT.Sub(values0, values1, values0)

				verifyTestVectors(tc, tc.decryptor, values0, ciphertext0, t)

			})
		}

		for _, lvl := range tc.testLevel {
			t.Run(GetTestName("Sub/op0=ct/op2=pt", tc.params, lvl), func(t *testing.T) {

				values0, _, ciphertext0 := newTestVectorsLvl(lvl, rlwe.NewScale(3), tc, tc.encryptorSk)
				values1, plaintext, _ := newTestVectorsLvl(lvl, rlwe.NewScale(7), tc, tc.encryptorSk)

				require.True(t, ciphertext0.Scale.Cmp(plaintext.Scale) != 0)

				tc.evaluator.Sub(ciphertext0, plaintext, ciphertext0)
				tc.ringT.Sub(values0, values1, values0)

				verifyTestVectors(tc, tc.decryptor, values0, ciphertext0, t)

			})
		}

		for _, lvl := range tc.testLevel {
			t.Run(GetTestName("Sub/op0=ct/op2=ct", tc.params, lvl), func(t *testing.T) {

				values0, _, ciphertext0 := newTestVectorsLvl(lvl, rlwe.NewScale(3), tc, tc.encryptorSk)
				values1, _, ciphertext1 := newTestVectorsLvl(lvl, rlwe.NewScale(7), tc, tc.encryptorSk)

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

				values, _, ciphertext := newTestVectorsLvl(lvl, rlwe.NewScale(7), tc, tc.encryptorSk)

				require.True(t, ciphertext.Scale.Cmp(1) != 0)

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
			t.Run(GetTestName("MulScalarAndAdd/op0=ct", tc.params, lvl), func(t *testing.T) {

				if lvl == 0 {
					t.Skip("Level = 0")
				}

				values0, _, ciphertext0 := newTestVectorsLvl(lvl, rlwe.NewScale(7), tc, tc.encryptorSk)
				values1, _, ciphertext1 := newTestVectorsLvl(lvl, rlwe.NewScale(3), tc, tc.encryptorSk)

				require.True(t, ciphertext0.Scale.Cmp(ciphertext1.Scale) != 0)

				scalar := tc.params.T() >> 1

				tc.evaluator.MulScalarAndAdd(ciphertext0, scalar, ciphertext1)
				tc.ringT.MulScalarAndAdd(values0, scalar, values1)

				verifyTestVectors(tc, tc.decryptor, values1, ciphertext1, t)
			})
		}

		for _, lvl := range tc.testLevel {
			t.Run(GetTestName("Mul/op0=ct/op2=ct", tc.params, lvl), func(t *testing.T) {

				if lvl == 0 {
					t.Skip("Level = 0")
				}

				values0, _, ciphertext0 := newTestVectorsLvl(lvl, rlwe.NewScale(3), tc, tc.encryptorSk)
				values1, _, ciphertext1 := newTestVectorsLvl(lvl, rlwe.NewScale(7), tc, tc.encryptorSk)

				require.True(t, ciphertext0.Scale.Cmp(ciphertext1.Scale) != 0)

				tc.evaluator.Mul(ciphertext0, ciphertext1, ciphertext0)
				tc.ringT.MulCoeffs(values0, values1, values0)

				verifyTestVectors(tc, tc.decryptor, values0, ciphertext0, t)

			})
		}

		for _, lvl := range tc.testLevel {
			t.Run(GetTestName("Square/op0=ct/op2=ct", tc.params, lvl), func(t *testing.T) {

				if lvl == 0 {
					t.Skip("Level = 0")
				}

				values0, _, ciphertext0 := newTestVectorsLvl(lvl, rlwe.NewScale(3), tc, tc.encryptorSk)

				tc.evaluator.Mul(ciphertext0, ciphertext0, ciphertext0)
				tc.ringT.MulCoeffs(values0, values0, values0)

				verifyTestVectors(tc, tc.decryptor, values0, ciphertext0, t)

			})
		}

		for _, lvl := range tc.testLevel {
			t.Run(GetTestName("MulRelin/op0=ct/op2=ct", tc.params, lvl), func(t *testing.T) {

				if lvl == 0 {
					t.Skip("Level = 0")
				}

				values0, _, ciphertext0 := newTestVectorsLvl(lvl, rlwe.NewScale(3), tc, tc.encryptorSk)
				values1, _, ciphertext1 := newTestVectorsLvl(lvl, rlwe.NewScale(7), tc, tc.encryptorSk)
				tc.ringT.MulCoeffs(values0, values1, values0)

				require.True(t, ciphertext0.Scale.Cmp(ciphertext1.Scale) != 0)

				receiver := NewCiphertext(tc.params, 1, lvl)

				tc.evaluator.Mul(ciphertext0, ciphertext1, receiver)

				require.Equal(t, receiver.Degree(), 2)

				verifyTestVectors(tc, tc.decryptor, values0, receiver, t)

				receiver = tc.evaluator.RelinearizeNew(receiver)

				require.Equal(t, receiver.Degree(), 1)

				verifyTestVectors(tc, tc.decryptor, values0, receiver, t)

				tc.evaluator.MulRelin(ciphertext0, ciphertext1, receiver)

				verifyTestVectors(tc, tc.decryptor, values0, receiver, t)

			})
		}

		for _, lvl := range tc.testLevel {
			t.Run(GetTestName("MulAndAdd/op0=ct/op2=ct", tc.params, lvl), func(t *testing.T) {

				if lvl == 0 {
					t.Skip("Level = 0")
				}

				values0, _, ciphertext0 := newTestVectorsLvl(lvl, tc.params.DefaultScale(), tc, tc.encryptorSk)
				values1, _, ciphertext1 := newTestVectorsLvl(lvl, rlwe.NewScale(2), tc, tc.encryptorSk)
				values2, _, ciphertext2 := newTestVectorsLvl(lvl, rlwe.NewScale(7), tc, tc.encryptorSk)

				require.True(t, ciphertext0.Scale.Cmp(ciphertext1.Scale) != 0)
				require.True(t, ciphertext0.Scale.Cmp(ciphertext2.Scale) != 0)

				tc.evaluator.MulAndAdd(ciphertext0, ciphertext1, ciphertext2)
				tc.ringT.MulCoeffsAndAdd(values0, values1, values2)

				verifyTestVectors(tc, tc.decryptor, values2, ciphertext2, t)

			})
		}

		for _, lvl := range tc.testLevel {
			t.Run(GetTestName("MulRelinAndAdd/op0=ct/op2=ct", tc.params, lvl), func(t *testing.T) {

				if lvl == 0 {
					t.Skip("Level = 0")
				}

				values0, _, ciphertext0 := newTestVectorsLvl(lvl, tc.params.DefaultScale(), tc, tc.encryptorSk)
				values1, _, ciphertext1 := newTestVectorsLvl(lvl, rlwe.NewScale(2), tc, tc.encryptorSk)
				values2, _, ciphertext2 := newTestVectorsLvl(lvl, rlwe.NewScale(7), tc, tc.encryptorSk)

				require.True(t, ciphertext0.Scale.Cmp(ciphertext1.Scale) != 0)
				require.True(t, ciphertext0.Scale.Cmp(ciphertext2.Scale) != 0)

				tc.evaluator.MulRelinAndAdd(ciphertext0, ciphertext1, ciphertext2)
				tc.ringT.MulCoeffsAndAdd(values0, values1, values2)

				verifyTestVectors(tc, tc.decryptor, values2, ciphertext2, t)

			})
		}

		t.Run(GetTestName("PolyEval/Single", tc.params, tc.params.MaxLevel()), func(t *testing.T) {

			if tc.params.MaxLevel() < 4 {
				t.Skip("MaxLevel() to low")
			}

			values, _, ciphertext := newTestVectorsLvl(tc.params.MaxLevel(), rlwe.NewScale(7), tc, tc.encryptorSk)

			coeffs := []uint64{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16}

			T := tc.ringT.Modulus[0]
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

			std, min, max := Norm(ciphertext, tc.decryptor)
			t.Logf("Noise -> (std: %f, min: %f, max=%f)\n", std, min, max)

			verifyTestVectors(tc, tc.decryptor, values, ciphertext, t)
		})

		t.Run(GetTestName("PolyEval/Vector", tc.params, tc.params.MaxLevel()), func(t *testing.T) {

			if tc.params.MaxLevel() < 4 {
				t.Skip("MaxLevel() to low")
			}

			values, _, ciphertext := newTestVectorsLvl(tc.params.MaxLevel(), rlwe.NewScale(7), tc, tc.encryptorSk)

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

			T := tc.ringT.Modulus[0]
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

			std, min, max := Norm(ciphertext, tc.decryptor)
			t.Logf("Noise -> (std: %f, min: %f, max=%f)\n", std, min, max)

			verifyTestVectors(tc, tc.decryptor, values, ciphertext, t)
		})

		for _, lvl := range tc.testLevel[:] {
			t.Run(GetTestName("Rescale", tc.params, lvl), func(t *testing.T) {

				ringT := tc.params.RingT()

				values0, _, ciphertext0 := newTestVectorsLvl(lvl, tc.params.DefaultScale(), tc, tc.encryptorPk)

				if lvl != 0 {

					values1, _, ciphertext1 := newTestVectorsLvl(lvl, tc.params.DefaultScale(), tc, tc.encryptorSk)

					stdErrFresh, _, _ := Norm(ciphertext0, tc.decryptor)

					t.Logf("STDErr 0x: %f\n", stdErrFresh)

					for i := 0; i < lvl; i++ {
						tc.evaluator.MulRelin(ciphertext0, ciphertext1, ciphertext0)

						vartmp, _, _ := Norm(ciphertext0, tc.decryptor)

						t.Logf("STDErr %dx: %f\n", i+1, vartmp)

						ringT.MulCoeffs(values0, values1, values0)
					}

					verifyTestVectors(tc, tc.decryptor, values0, ciphertext0, t)

					stdErrMul, _, _ := Norm(ciphertext0, tc.decryptor)

					ciphertext1 = ciphertext0.CopyNew()

					require.Nil(t, tc.evaluator.Rescale(ciphertext0, ciphertext0))
					require.Nil(t, tc.evaluator.(*evaluator).rescaleOriginal(ciphertext1, ciphertext1))

					verifyTestVectors(tc, tc.decryptor, values0, ciphertext0, t)
					verifyTestVectors(tc, tc.decryptor, values0, ciphertext1, t)

					stdErrNewRescale, _, _ := Norm(ciphertext0, tc.decryptor)
					stdErrOldRescale, _, _ := Norm(ciphertext1, tc.decryptor)

					t.Logf("STDErr (mul): %f\n", stdErrMul)
					t.Logf("STDErr (t(div(xt^-1)): %f\n", stdErrNewRescale)
					t.Logf("STDErr (modswitch(x)): %f\n", stdErrOldRescale)

				} else {
					require.NotNil(t, tc.evaluator.Rescale(ciphertext0, ciphertext0))
				}
			})
		}
	})
}

func testRotate(tc *testContext, t *testing.T) {

	rots := []int{1, -1, 4, -4, 63, -63}
	rotkey := tc.kgen.GenRotationKeysForRotations(rots, true, tc.sk)
	evaluator := tc.evaluator.WithKey(rlwe.EvaluationKey{Rlk: tc.rlk, Rtks: rotkey})

	for _, lvl := range tc.testLevel {
		t.Run(GetTestName("Evaluator/RotateRows", tc.params, lvl), func(t *testing.T) {
			values, _, ciphertext := newTestVectorsLvl(lvl, tc.params.DefaultScale(), tc, tc.encryptorPk)
			evaluator.RotateRows(ciphertext, ciphertext)
			values.Coeffs[0] = append(values.Coeffs[0][tc.params.N()>>1:], values.Coeffs[0][:tc.params.N()>>1]...)
			verifyTestVectors(tc, tc.decryptor, values, ciphertext, t)
		})
	}

	for _, lvl := range tc.testLevel {
		t.Run(GetTestName("Evaluator/RotateRowsNew", tc.params, lvl), func(t *testing.T) {
			values, _, ciphertext := newTestVectorsLvl(lvl, tc.params.DefaultScale(), tc, tc.encryptorPk)
			ciphertext = evaluator.RotateRowsNew(ciphertext)
			values.Coeffs[0] = append(values.Coeffs[0][tc.params.N()>>1:], values.Coeffs[0][:tc.params.N()>>1]...)
			verifyTestVectors(tc, tc.decryptor, values, ciphertext, t)
		})
	}

	for _, lvl := range tc.testLevel {
		t.Run(GetTestName("Evaluator/RotateColumns", tc.params, lvl), func(t *testing.T) {

			values, _, ciphertext := newTestVectorsLvl(lvl, tc.params.DefaultScale(), tc, tc.encryptorPk)

			receiver := NewCiphertext(tc.params, 1, lvl)
			for _, n := range rots {

				evaluator.RotateColumns(ciphertext, n, receiver)
				valuesWant := utils.RotateUint64Slots(values.Coeffs[0], n)

				verifyTestVectors(tc, tc.decryptor, &ring.Poly{Coeffs: [][]uint64{valuesWant}}, receiver, t)
			}
		})
	}

	for _, lvl := range tc.testLevel {
		t.Run(GetTestName("Evaluator/RotateColumnsNew", tc.params, lvl), func(t *testing.T) {

			values, _, ciphertext := newTestVectorsLvl(lvl, tc.params.DefaultScale(), tc, tc.encryptorPk)

			for _, n := range rots {

				receiver := evaluator.RotateColumnsNew(ciphertext, n)
				valuesWant := utils.RotateUint64Slots(values.Coeffs[0], n)

				verifyTestVectors(tc, tc.decryptor, &ring.Poly{Coeffs: [][]uint64{valuesWant}}, receiver, t)
			}
		})
	}
}

func testInnerSum(tc *testContext, t *testing.T) {

	t.Run(GetTestName("InnerSum", tc.params, tc.params.MaxLevel()), func(t *testing.T) {

		batch := 128
		n := tc.params.N() / (2 * batch)

		rotKey := tc.kgen.GenRotationKeysForRotations(tc.params.RotationsForInnerSumLog(batch, n), false, tc.sk)
		eval := tc.evaluator.WithKey(rlwe.EvaluationKey{Rlk: tc.rlk, Rtks: rotKey})

		values, _, ciphertext := newTestVectorsLvl(tc.params.MaxLevel(), tc.params.DefaultScale(), tc, tc.encryptorSk)

		eval.InnerSumLog(ciphertext, batch, n, ciphertext)

		tmp := make([]uint64, tc.params.N())
		copy(tmp, values.Coeffs[0])

		T := tc.params.T()
		for i := 1; i < n; i++ {
			ring.AddVec(values.Coeffs[0], utils.RotateUint64Slots(tmp, i*batch), values.Coeffs[0], T)
		}

		verifyTestVectors(tc, tc.decryptor, values, ciphertext, t)
	})
}

func testLinearTransform(tc *testContext, t *testing.T) {

	t.Run(GetTestName("LinearTransform/Naive", tc.params, tc.params.MaxLevel()), func(t *testing.T) {

		params := tc.params

		values, _, ciphertext := newTestVectorsLvl(tc.params.MaxLevel(), tc.params.DefaultScale(), tc, tc.encryptorSk)

		diagMatrix := make(map[int][]uint64)

		diagMatrix[-1] = make([]uint64, params.N())
		diagMatrix[0] = make([]uint64, params.N())
		diagMatrix[1] = make([]uint64, params.N())

		for i := 0; i < params.N(); i++ {
			diagMatrix[-1][i] = 1
			diagMatrix[0][i] = 1
			diagMatrix[1][i] = 1
		}

		linTransf := GenLinearTransform(tc.encoder, diagMatrix, params.MaxLevel(), tc.params.DefaultScale())

		rots := linTransf.Rotations()

		rotKey := tc.kgen.GenRotationKeysForRotations(rots, false, tc.sk)

		eval := tc.evaluator.WithKey(rlwe.EvaluationKey{Rlk: tc.rlk, Rtks: rotKey})

		eval.LinearTransform(ciphertext, linTransf, []*rlwe.Ciphertext{ciphertext})

		tmp := make([]uint64, params.N())
		copy(tmp, values.Coeffs[0])

		T := tc.params.T()

		ring.AddVec(values.Coeffs[0], utils.RotateUint64Slots(tmp, -1), values.Coeffs[0], T)
		ring.AddVec(values.Coeffs[0], utils.RotateUint64Slots(tmp, 1), values.Coeffs[0], T)

		verifyTestVectors(tc, tc.decryptor, values, ciphertext, t)
	})

	t.Run(GetTestName("LinearTransform/BSGS", tc.params, tc.params.MaxLevel()), func(t *testing.T) {

		params := tc.params

		values, _, ciphertext := newTestVectorsLvl(tc.params.MaxLevel(), tc.params.DefaultScale(), tc, tc.encryptorSk)

		diagMatrix := make(map[int][]uint64)

		diagMatrix[-15] = make([]uint64, params.N())
		diagMatrix[-4] = make([]uint64, params.N())
		diagMatrix[-1] = make([]uint64, params.N())
		diagMatrix[0] = make([]uint64, params.N())
		diagMatrix[1] = make([]uint64, params.N())
		diagMatrix[2] = make([]uint64, params.N())
		diagMatrix[3] = make([]uint64, params.N())
		diagMatrix[4] = make([]uint64, params.N())
		diagMatrix[15] = make([]uint64, params.N())

		for i := 0; i < params.N(); i++ {
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

		rots := linTransf.Rotations()

		rotKey := tc.kgen.GenRotationKeysForRotations(rots, false, tc.sk)

		eval := tc.evaluator.WithKey(rlwe.EvaluationKey{Rlk: tc.rlk, Rtks: rotKey})

		eval.LinearTransform(ciphertext, linTransf, []*rlwe.Ciphertext{ciphertext})

		tmp := make([]uint64, params.N())
		copy(tmp, values.Coeffs[0])

		T := tc.params.T()

		ring.AddVec(values.Coeffs[0], utils.RotateUint64Slots(tmp, -15), values.Coeffs[0], T)
		ring.AddVec(values.Coeffs[0], utils.RotateUint64Slots(tmp, -4), values.Coeffs[0], T)
		ring.AddVec(values.Coeffs[0], utils.RotateUint64Slots(tmp, -1), values.Coeffs[0], T)
		ring.AddVec(values.Coeffs[0], utils.RotateUint64Slots(tmp, 1), values.Coeffs[0], T)
		ring.AddVec(values.Coeffs[0], utils.RotateUint64Slots(tmp, 2), values.Coeffs[0], T)
		ring.AddVec(values.Coeffs[0], utils.RotateUint64Slots(tmp, 3), values.Coeffs[0], T)
		ring.AddVec(values.Coeffs[0], utils.RotateUint64Slots(tmp, 4), values.Coeffs[0], T)
		ring.AddVec(values.Coeffs[0], utils.RotateUint64Slots(tmp, 15), values.Coeffs[0], T)

		verifyTestVectors(tc, tc.decryptor, values, ciphertext, t)
	})
}

func testMerge(tc *testContext, t *testing.T) {

	params := tc.params

	t.Run(GetTestName("Merge", params, params.MaxLevel()), func(t *testing.T) {

		values := make([]uint64, params.N())
		for i := range values {
			values[i] = uint64(i)
		}

		n := 16

		ciphertexts := make(map[int]*rlwe.Ciphertext)
		slotIndex := make(map[int]bool)

		pt := NewPlaintext(params, params.MaxLevel())
		for i := 0; i < params.N(); i += params.N() / n {

			tc.encoder.EncodeCoeffs(append(values[i:], values[i:]...), pt)

			ciphertexts[i] = tc.encryptorSk.EncryptNew(pt)
			slotIndex[i] = true
		}

		// Rotation Keys
		galEls := params.GaloisElementsForMergeRLWE()
		rtks := tc.kgen.GenRotationKeys(galEls, tc.sk)

		eval := NewEvaluator(params, rlwe.EvaluationKey{Rtks: rtks})

		ciphertext := eval.Merge(ciphertexts)

		valuesHave := tc.encoder.DecodeCoeffsNew(tc.decryptor.DecryptNew(ciphertext))

		for i := range values {
			if _, ok := slotIndex[i]; ok {
				require.Equal(t, valuesHave[i], values[i])
			} else {
				require.Equal(t, valuesHave[i], uint64(0))
			}
		}
	})
}

func testSwitchKeys(tc *testContext, t *testing.T) {

	sk2 := tc.kgen.GenSecretKey()
	decryptorSk2 := NewDecryptor(tc.params, sk2)
	switchingKey := tc.kgen.GenSwitchingKey(tc.sk, sk2)

	for _, lvl := range tc.testLevel {
		t.Run(GetTestName("SwitchKeys", tc.params, lvl), func(t *testing.T) {
			values, _, ciphertext := newTestVectorsLvl(lvl, tc.params.DefaultScale(), tc, tc.encryptorPk)
			tc.evaluator.SwitchKeys(ciphertext, switchingKey, ciphertext)
			verifyTestVectors(tc, decryptorSk2, values, ciphertext, t)
		})
	}

	for _, lvl := range tc.testLevel {
		t.Run(GetTestName("SwitchKeysNew", tc.params, lvl), func(t *testing.T) {
			values, _, ciphertext := newTestVectorsLvl(lvl, tc.params.DefaultScale(), tc, tc.encryptorPk)
			ciphertext = tc.evaluator.SwitchKeysNew(ciphertext, switchingKey)
			verifyTestVectors(tc, decryptorSk2, values, ciphertext, t)
		})
	}
}

func testMarshalling(tc *testContext, t *testing.T) {
	t.Run("Marshalling", func(t *testing.T) {
		t.Run(GetTestName("PowerBasis", tc.params, tc.params.MaxLevel()), func(t *testing.T) {

			if tc.params.MaxLevel() < 4 {
				t.Skip("not enough levels")
			}

			_, _, ct := newTestVectorsLvl(tc.params.MaxLevel(), tc.params.DefaultScale(), tc, tc.encryptorPk)

			pb := NewPowerBasis(ct)

			for i := 2; i < 4; i++ {
				pb.GenPower(i, true, tc.evaluator)
			}

			pbBytes, err := pb.MarshalBinary()
			require.Nil(t, err)
			pbNew := new(PowerBasis)
			require.Nil(t, pbNew.UnmarshalBinary(pbBytes))

			for i := range pb.Value {
				ctWant := pb.Value[i]
				ctHave := pbNew.Value[i]
				require.NotNil(t, ctHave)
				for j := range ctWant.Value {
					require.True(t, tc.ringQ.Equal(ctWant.Value[j], ctHave.Value[j]))
				}
			}
		})
	})
}
