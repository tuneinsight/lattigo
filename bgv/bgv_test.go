package bgv

import (
	"fmt"
	"runtime"
	"testing"

	"github.com/tuneinsight/lattigo/v3/ring"
	"github.com/tuneinsight/lattigo/v3/rlwe"
	"github.com/tuneinsight/lattigo/v3/utils"

	"github.com/stretchr/testify/require"
)

var (
	// TESTN13QP218 is a of 128-bit secure test parameters set with a 32-bit plaintext and depth 4.
	TESTN14QP418 = ParametersLiteral{
		ParametersLiteral: rlwe.ParametersLiteral{
			LogN:  13,
			Q:     []uint64{0xffffc4001, 0x1000090001, 0x10000c8001, 0x10000f0001, 0xffff00001},
			P:     []uint64{0x1ffffe0001},
			Sigma: rlwe.DefaultSigma,
		},
		T: 0xffc001,
	}

	// TestParams is a set of test parameters for BGV ensuring 128 bit security in the classic setting.
	TestParams = []ParametersLiteral{TESTN14QP418}
)

func testString(opname string, p Parameters, lvl int) string {
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
	if params.PCount() != 0 {
		tc.rlk = tc.kgen.GenRelinearizationKey(tc.sk, 1)
	}

	tc.encoder = NewEncoder(tc.params)
	tc.encryptorPk = NewEncryptor(tc.params, tc.pk)
	tc.encryptorSk = NewEncryptor(tc.params, tc.sk)
	tc.decryptor = NewDecryptor(tc.params, tc.sk)
	tc.evaluator = NewEvaluator(tc.params, rlwe.EvaluationKey{Rlk: tc.rlk})

	tc.testLevel = []int{0, params.MaxLevel()}

	return
}

func newTestVectorsLvl(level int, scale uint64, tc *testContext, encryptor Encryptor) (coeffs *ring.Poly, plaintext *Plaintext, ciphertext *Ciphertext) {
	coeffs = tc.uSampler.ReadNew()

	plaintext = NewPlaintext(tc.params, level, scale)
	tc.encoder.Encode(coeffs.Coeffs[0], plaintext)
	if encryptor != nil {
		ciphertext = encryptor.EncryptNew(plaintext)
	}

	return coeffs, plaintext, ciphertext
}

func verifyTestVectors(tc *testContext, decryptor Decryptor, coeffs *ring.Poly, element Operand, t *testing.T) {

	var coeffsTest []uint64

	switch el := element.(type) {
	case *Plaintext:
		coeffsTest = tc.encoder.DecodeUintNew(el)
	case *Ciphertext:
		coeffsTest = tc.encoder.DecodeUintNew(decryptor.DecryptNew(el))
	default:
		t.Error("invalid test object to verify")
	}

	require.True(t, utils.EqualSliceUint64(coeffs.Coeffs[0], coeffsTest))
}

func testParameters(tc *testContext, t *testing.T) {

	t.Run("Parameters/NewParameters", func(t *testing.T) {
		params, err := NewParametersFromLiteral(ParametersLiteral{
			ParametersLiteral: rlwe.ParametersLiteral{
				LogN: 4,
				LogQ: []int{60, 60},
				LogP: []int{60},
			},
			T: 0x10001,
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
		t.Run(testString("Encoder/Encode&Decode/Uint", tc.params, lvl), func(t *testing.T) {
			values, plaintext, _ := newTestVectorsLvl(lvl, 1, tc, nil)
			verifyTestVectors(tc, nil, values, plaintext, t)
		})
	}

	for _, lvl := range tc.testLevel {
		t.Run(testString("Encoder/Encode&Decode/Int", tc.params, lvl), func(t *testing.T) {

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

			plaintext := NewPlaintext(tc.params, lvl, 1)
			tc.encoder.Encode(coeffsInt, plaintext)
			require.True(t, utils.EqualSliceInt64(coeffsInt, tc.encoder.DecodeIntNew(plaintext)))
		})
	}
}

func testEncryptor(tc *testContext, t *testing.T) {

	for _, lvl := range tc.testLevel {
		t.Run(testString("Encoder/EncryptorPk", tc.params, lvl), func(t *testing.T) {
			values, _, ciphertext := newTestVectorsLvl(lvl, 1, tc, tc.encryptorPk)
			verifyTestVectors(tc, tc.decryptor, values, ciphertext, t)
		})
	}

	for _, lvl := range tc.testLevel {
		t.Run(testString("Encoder/encryptorSk", tc.params, lvl), func(t *testing.T) {
			values, _, ciphertext := newTestVectorsLvl(lvl, 1, tc, tc.encryptorSk)
			verifyTestVectors(tc, tc.decryptor, values, ciphertext, t)
		})
	}
}

func testEvaluator(tc *testContext, t *testing.T) {

	t.Run("Evaluator", func(t *testing.T) {

		for _, lvl := range tc.testLevel {
			t.Run(testString("AddNew/op0=ct/op2=ct", tc.params, lvl), func(t *testing.T) {

				values0, _, ciphertext0 := newTestVectorsLvl(lvl, 3, tc, tc.encryptorSk)
				values1, _, ciphertext1 := newTestVectorsLvl(lvl, 7, tc, tc.encryptorSk)

				require.NotEqual(t, ciphertext0.Scale, ciphertext1.Scale)

				ciphertext2 := tc.evaluator.AddNew(ciphertext0, ciphertext1)
				tc.ringT.Add(values0, values1, values0)

				verifyTestVectors(tc, tc.decryptor, values0, ciphertext2, t)

			})
		}

		for _, lvl := range tc.testLevel {
			t.Run(testString("Add/op0=ct/op2=ct", tc.params, lvl), func(t *testing.T) {

				values0, _, ciphertext0 := newTestVectorsLvl(lvl, 3, tc, tc.encryptorSk)
				values1, _, ciphertext1 := newTestVectorsLvl(lvl, 7, tc, tc.encryptorSk)

				require.NotEqual(t, ciphertext0.Scale, ciphertext1.Scale)

				tc.evaluator.Add(ciphertext0, ciphertext1, ciphertext0)
				tc.ringT.Add(values0, values1, values0)

				verifyTestVectors(tc, tc.decryptor, values0, ciphertext0, t)

			})
		}

		for _, lvl := range tc.testLevel {
			t.Run(testString("Add/op0=ct/op2=pt", tc.params, lvl), func(t *testing.T) {

				values0, _, ciphertext0 := newTestVectorsLvl(lvl, 3, tc, tc.encryptorSk)
				values1, plaintext, _ := newTestVectorsLvl(lvl, 7, tc, tc.encryptorSk)

				require.NotEqual(t, ciphertext0.Scale, plaintext.Scale)

				tc.evaluator.Add(ciphertext0, plaintext, ciphertext0)
				tc.ringT.Add(values0, values1, values0)

				verifyTestVectors(tc, tc.decryptor, values0, ciphertext0, t)

			})
		}

		for _, lvl := range tc.testLevel {
			t.Run(testString("Add/op0=ct/op2=ct", tc.params, lvl), func(t *testing.T) {

				values0, _, ciphertext0 := newTestVectorsLvl(lvl, 3, tc, tc.encryptorSk)
				values1, _, ciphertext1 := newTestVectorsLvl(lvl, 7, tc, tc.encryptorSk)

				require.NotEqual(t, ciphertext0.Scale, ciphertext1.Scale)

				tc.evaluator.Add(ciphertext0, ciphertext1, ciphertext0)
				tc.ringT.Add(values0, values1, values0)

				verifyTestVectors(tc, tc.decryptor, values0, ciphertext0, t)

			})
		}

		for _, lvl := range tc.testLevel {
			t.Run(testString("SubNew/op0=ct/op2=ct", tc.params, lvl), func(t *testing.T) {

				values0, _, ciphertext0 := newTestVectorsLvl(lvl, 3, tc, tc.encryptorSk)
				values1, _, ciphertext1 := newTestVectorsLvl(lvl, 7, tc, tc.encryptorSk)

				require.NotEqual(t, ciphertext0.Scale, ciphertext1.Scale)

				ciphertext0 = tc.evaluator.SubNew(ciphertext0, ciphertext1)
				tc.ringT.Sub(values0, values1, values0)

				verifyTestVectors(tc, tc.decryptor, values0, ciphertext0, t)

			})
		}

		for _, lvl := range tc.testLevel {
			t.Run(testString("Sub/op0=ct/op2=ct", tc.params, lvl), func(t *testing.T) {

				values0, _, ciphertext0 := newTestVectorsLvl(lvl, 3, tc, tc.encryptorSk)
				values1, _, ciphertext1 := newTestVectorsLvl(lvl, 7, tc, tc.encryptorSk)

				require.NotEqual(t, ciphertext0.Scale, ciphertext1.Scale)

				tc.evaluator.Sub(ciphertext0, ciphertext1, ciphertext0)
				tc.ringT.Sub(values0, values1, values0)

				verifyTestVectors(tc, tc.decryptor, values0, ciphertext0, t)

			})
		}

		for _, lvl := range tc.testLevel {
			t.Run(testString("Sub/op0=ct/op2=pt", tc.params, lvl), func(t *testing.T) {

				values0, _, ciphertext0 := newTestVectorsLvl(lvl, 3, tc, tc.encryptorSk)
				values1, plaintext, _ := newTestVectorsLvl(lvl, 7, tc, tc.encryptorSk)

				require.NotEqual(t, ciphertext0.Scale, plaintext.Scale)

				tc.evaluator.Sub(ciphertext0, plaintext, ciphertext0)
				tc.ringT.Sub(values0, values1, values0)

				verifyTestVectors(tc, tc.decryptor, values0, ciphertext0, t)

			})
		}

		for _, lvl := range tc.testLevel {
			t.Run(testString("Sub/op0=ct/op2=ct", tc.params, lvl), func(t *testing.T) {

				values0, _, ciphertext0 := newTestVectorsLvl(lvl, 3, tc, tc.encryptorSk)
				values1, _, ciphertext1 := newTestVectorsLvl(lvl, 7, tc, tc.encryptorSk)

				require.NotEqual(t, ciphertext0.Scale, ciphertext1.Scale)

				tc.evaluator.Sub(ciphertext0, ciphertext1, ciphertext0)
				tc.ringT.Sub(values0, values1, values0)

				verifyTestVectors(tc, tc.decryptor, values0, ciphertext0, t)

			})
		}

		for _, lvl := range tc.testLevel {
			t.Run(testString("Neg/op0=ct", tc.params, lvl), func(t *testing.T) {

				values, _, ciphertext := newTestVectorsLvl(lvl, 1, tc, tc.encryptorSk)

				tc.evaluator.Neg(ciphertext, ciphertext)
				tc.ringT.Neg(values, values)
				tc.ringT.Reduce(values, values)

				verifyTestVectors(tc, tc.decryptor, values, ciphertext, t)

			})
		}

		for _, lvl := range tc.testLevel {
			t.Run(testString("NegNew/op0=ct", tc.params, lvl), func(t *testing.T) {

				values, _, ciphertext := newTestVectorsLvl(lvl, 1, tc, tc.encryptorSk)

				ciphertext = tc.evaluator.NegNew(ciphertext)
				tc.ringT.Neg(values, values)
				tc.ringT.Reduce(values, values)

				verifyTestVectors(tc, tc.decryptor, values, ciphertext, t)

			})
		}

		for _, lvl := range tc.testLevel {
			t.Run(testString("AddScalar/op0=ct", tc.params, lvl), func(t *testing.T) {

				values, _, ciphertext := newTestVectorsLvl(lvl, 1, tc, tc.encryptorSk)

				scalar := tc.params.T() >> 1

				tc.evaluator.AddScalar(ciphertext, scalar, ciphertext)
				tc.ringT.AddScalar(values, scalar, values)

				verifyTestVectors(tc, tc.decryptor, values, ciphertext, t)

			})
		}

		for _, lvl := range tc.testLevel {
			t.Run(testString("AddScalarNew/op0=ct", tc.params, lvl), func(t *testing.T) {

				values, _, ciphertext := newTestVectorsLvl(lvl, 7, tc, tc.encryptorSk)

				require.NotEqual(t, ciphertext.Scale, 1)

				scalar := tc.params.T() >> 1

				ciphertext = tc.evaluator.AddScalarNew(ciphertext, scalar)
				tc.ringT.AddScalar(values, scalar, values)

				verifyTestVectors(tc, tc.decryptor, values, ciphertext, t)
			})
		}

		for _, lvl := range tc.testLevel {
			t.Run(testString("MulScalar/op0=ct", tc.params, lvl), func(t *testing.T) {

				if lvl == 0 {
					t.Skip("Level = 0")
				}

				values, _, ciphertext := newTestVectorsLvl(lvl, 1, tc, tc.encryptorSk)

				scalar := tc.params.T() >> 1

				tc.evaluator.MulScalar(ciphertext, scalar, ciphertext)
				tc.ringT.MulScalar(values, scalar, values)

				verifyTestVectors(tc, tc.decryptor, values, ciphertext, t)

			})
		}

		for _, lvl := range tc.testLevel {
			t.Run(testString("MulScalarNew/op0=ct", tc.params, lvl), func(t *testing.T) {

				if lvl == 0 {
					t.Skip("Level = 0")
				}

				values, _, ciphertext := newTestVectorsLvl(lvl, 1, tc, tc.encryptorSk)

				scalar := tc.params.T() >> 1

				ciphertext = tc.evaluator.MulScalarNew(ciphertext, scalar)
				tc.ringT.MulScalar(values, scalar, values)

				verifyTestVectors(tc, tc.decryptor, values, ciphertext, t)
			})
		}

		for _, lvl := range tc.testLevel {
			t.Run(testString("MulScalarAndAdd/op0=ct", tc.params, lvl), func(t *testing.T) {

				if lvl == 0 {
					t.Skip("Level = 0")
				}

				values0, _, ciphertext0 := newTestVectorsLvl(lvl, 7, tc, tc.encryptorSk)
				values1, _, ciphertext1 := newTestVectorsLvl(lvl, 3, tc, tc.encryptorSk)

				require.NotEqual(t, ciphertext0.Scale, ciphertext1.Scale)

				scalar := tc.params.T() >> 1

				tc.evaluator.MulScalarAndAdd(ciphertext0, scalar, ciphertext1)
				tc.ringT.MulScalarAndAdd(values0, scalar, values1)

				verifyTestVectors(tc, tc.decryptor, values1, ciphertext1, t)
			})
		}

		for _, lvl := range tc.testLevel {
			t.Run(testString("Mul/op0=ct/op2=ct", tc.params, lvl), func(t *testing.T) {

				if lvl == 0 {
					t.Skip("Level = 0")
				}

				values0, _, ciphertext0 := newTestVectorsLvl(lvl, 3, tc, tc.encryptorSk)
				values1, _, ciphertext1 := newTestVectorsLvl(lvl, 7, tc, tc.encryptorSk)

				require.NotEqual(t, ciphertext0.Scale, ciphertext1.Scale)

				tc.evaluator.Mul(ciphertext0, ciphertext1, ciphertext0)
				tc.ringT.MulCoeffs(values0, values1, values0)

				verifyTestVectors(tc, tc.decryptor, values0, ciphertext0, t)

			})
		}

		for _, lvl := range tc.testLevel {
			t.Run(testString("Square/op0=ct/op2=ct", tc.params, lvl), func(t *testing.T) {

				if lvl == 0 {
					t.Skip("Level = 0")
				}

				values0, _, ciphertext0 := newTestVectorsLvl(lvl, 3, tc, tc.encryptorSk)

				tc.evaluator.Mul(ciphertext0, ciphertext0, ciphertext0)
				tc.ringT.MulCoeffs(values0, values0, values0)

				verifyTestVectors(tc, tc.decryptor, values0, ciphertext0, t)

			})
		}

		for _, lvl := range tc.testLevel {
			t.Run(testString("MulRelin/op0=ct/op2=ct", tc.params, lvl), func(t *testing.T) {

				if lvl == 0 {
					t.Skip("Level = 0")
				}

				values0, _, ciphertext0 := newTestVectorsLvl(lvl, 3, tc, tc.encryptorSk)
				values1, _, ciphertext1 := newTestVectorsLvl(lvl, 7, tc, tc.encryptorSk)
				tc.ringT.MulCoeffs(values0, values1, values0)

				require.NotEqual(t, ciphertext0.Scale, ciphertext1.Scale)

				receiver := NewCiphertext(tc.params, 1, lvl, 1)

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
			t.Run(testString("MulAndAdd/op0=ct/op2=ct", tc.params, lvl), func(t *testing.T) {

				if lvl == 0 {
					t.Skip("Level = 0")
				}

				values0, _, ciphertext0 := newTestVectorsLvl(lvl, 1, tc, tc.encryptorSk)
				values1, _, ciphertext1 := newTestVectorsLvl(lvl, 2, tc, tc.encryptorSk)
				values2, _, ciphertext2 := newTestVectorsLvl(lvl, 7, tc, tc.encryptorSk)

				require.NotEqual(t, ciphertext0.Scale, ciphertext1.Scale)
				require.NotEqual(t, ciphertext0.Scale, ciphertext2.Scale)

				tc.evaluator.MulAndAdd(ciphertext0, ciphertext1, ciphertext2)
				tc.ringT.MulCoeffsAndAdd(values0, values1, values2)

				verifyTestVectors(tc, tc.decryptor, values2, ciphertext2, t)

			})
		}

		for _, lvl := range tc.testLevel {
			t.Run(testString("MulRelinAndAdd/op0=ct/op2=ct", tc.params, lvl), func(t *testing.T) {

				if lvl == 0 {
					t.Skip("Level = 0")
				}

				values0, _, ciphertext0 := newTestVectorsLvl(lvl, 1, tc, tc.encryptorSk)
				values1, _, ciphertext1 := newTestVectorsLvl(lvl, 2, tc, tc.encryptorSk)
				values2, _, ciphertext2 := newTestVectorsLvl(lvl, 7, tc, tc.encryptorSk)

				require.NotEqual(t, ciphertext0.Scale, ciphertext1.Scale)
				require.NotEqual(t, ciphertext0.Scale, ciphertext2.Scale)

				tc.evaluator.MulRelinAndAdd(ciphertext0, ciphertext1, ciphertext2)
				tc.ringT.MulCoeffsAndAdd(values0, values1, values2)

				verifyTestVectors(tc, tc.decryptor, values2, ciphertext2, t)

			})
		}

		t.Run(testString("PolyEval/Single", tc.params, tc.params.MaxLevel()), func(t *testing.T) {

			if tc.params.MaxLevel() < 4 {
				t.Skip("MaxLevel() to low")
			}

			values, _, ciphertext := newTestVectorsLvl(tc.params.MaxLevel(), 7, tc, tc.encryptorSk)

			coeffs := []uint64{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16}

			T := tc.ringT.Modulus[0]
			for i := range values.Coeffs[0] {
				values.Coeffs[0][i] = ring.EvalPolyModP(values.Coeffs[0][i], coeffs, T)
			}

			poly := NewPoly(coeffs)

			var targetScale uint64 = 1

			var err error
			if ciphertext, err = tc.evaluator.EvaluatePoly(ciphertext, poly, targetScale); err != nil {
				t.Log(err)
				t.Fatal()
			}

			require.True(t, ciphertext.Scale == targetScale)

			std, min, max := tc.decryptor.Noise(ciphertext)
			t.Logf("Noise -> (std: %f, min: %f, max=%f)\n", std, min, max)

			verifyTestVectors(tc, tc.decryptor, values, ciphertext, t)
		})

		t.Run(testString("PolyEval/Vector", tc.params, tc.params.MaxLevel()), func(t *testing.T) {

			if tc.params.MaxLevel() < 4 {
				t.Skip("MaxLevel() to low")
			}

			values, _, ciphertext := newTestVectorsLvl(tc.params.MaxLevel(), 7, tc, tc.encryptorSk)

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

			var targetScale uint64 = 1

			var err error
			if ciphertext, err = tc.evaluator.EvaluatePolyVector(ciphertext, polyVec, tc.encoder, slotIndex, targetScale); err != nil {
				t.Fail()
			}

			require.True(t, ciphertext.Scale == targetScale)

			std, min, max := tc.decryptor.Noise(ciphertext)
			t.Logf("Noise -> (std: %f, min: %f, max=%f)\n", std, min, max)

			verifyTestVectors(tc, tc.decryptor, values, ciphertext, t)
		})

		for _, lvl := range tc.testLevel[:] {
			t.Run(testString("Rescale", tc.params, lvl), func(t *testing.T) {

				ringT := tc.params.RingT()

				values0, _, ciphertext0 := newTestVectorsLvl(lvl, 1, tc, tc.encryptorPk)

				if lvl != 0 {

					values1, _, ciphertext1 := newTestVectorsLvl(lvl, 1, tc, tc.encryptorSk)

					stdErrFresh, _, _ := tc.decryptor.Noise(ciphertext0)

					t.Logf("STDErr 0x: %f\n", stdErrFresh)

					for i := 0; i < lvl; i++ {
						tc.evaluator.MulRelin(ciphertext0, ciphertext1, ciphertext0)

						vartmp, _, _ := tc.decryptor.Noise(ciphertext0)

						t.Logf("STDErr %dx: %f\n", i+1, vartmp)

						ringT.MulCoeffs(values0, values1, values0)
					}

					verifyTestVectors(tc, tc.decryptor, values0, ciphertext0, t)

					stdErrMul, _, _ := tc.decryptor.Noise(ciphertext0)

					ciphertext1 = ciphertext0.CopyNew()

					require.Nil(t, tc.evaluator.Rescale(ciphertext0, ciphertext0))
					require.Nil(t, tc.evaluator.(*evaluator).rescaleOriginal(ciphertext1, ciphertext1))

					verifyTestVectors(tc, tc.decryptor, values0, ciphertext0, t)
					verifyTestVectors(tc, tc.decryptor, values0, ciphertext1, t)

					stdErrNewRescale, _, _ := tc.decryptor.Noise(ciphertext0)
					stdErrOldRescale, _, _ := tc.decryptor.Noise(ciphertext1)

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
		t.Run(testString("Evaluator/RotateRows", tc.params, lvl), func(t *testing.T) {
			values, _, ciphertext := newTestVectorsLvl(lvl, 1, tc, tc.encryptorPk)
			evaluator.RotateRows(ciphertext, ciphertext)
			values.Coeffs[0] = append(values.Coeffs[0][tc.params.N()>>1:], values.Coeffs[0][:tc.params.N()>>1]...)
			verifyTestVectors(tc, tc.decryptor, values, ciphertext, t)
		})
	}

	for _, lvl := range tc.testLevel {
		t.Run(testString("Evaluator/RotateRowsNew", tc.params, lvl), func(t *testing.T) {
			values, _, ciphertext := newTestVectorsLvl(lvl, 1, tc, tc.encryptorPk)
			ciphertext = evaluator.RotateRowsNew(ciphertext)
			values.Coeffs[0] = append(values.Coeffs[0][tc.params.N()>>1:], values.Coeffs[0][:tc.params.N()>>1]...)
			verifyTestVectors(tc, tc.decryptor, values, ciphertext, t)
		})
	}

	for _, lvl := range tc.testLevel {
		t.Run(testString("Evaluator/RotateColumns", tc.params, lvl), func(t *testing.T) {

			values, _, ciphertext := newTestVectorsLvl(lvl, 1, tc, tc.encryptorPk)

			receiver := NewCiphertext(tc.params, 1, lvl, 1)
			for _, n := range rots {

				evaluator.RotateColumns(ciphertext, n, receiver)
				valuesWant := utils.RotateUint64Slots(values.Coeffs[0], n)

				verifyTestVectors(tc, tc.decryptor, &ring.Poly{Coeffs: [][]uint64{valuesWant}}, receiver, t)
			}
		})
	}

	for _, lvl := range tc.testLevel {
		t.Run(testString("Evaluator/RotateColumnsNew", tc.params, lvl), func(t *testing.T) {

			values, _, ciphertext := newTestVectorsLvl(lvl, 1, tc, tc.encryptorPk)

			for _, n := range rots {

				receiver := evaluator.RotateColumnsNew(ciphertext, n)
				valuesWant := utils.RotateUint64Slots(values.Coeffs[0], n)

				verifyTestVectors(tc, tc.decryptor, &ring.Poly{Coeffs: [][]uint64{valuesWant}}, receiver, t)
			}
		})
	}
}

func testSwitchKeys(tc *testContext, t *testing.T) {

	sk2 := tc.kgen.GenSecretKey()
	decryptorSk2 := NewDecryptor(tc.params, sk2)
	switchingKey := tc.kgen.GenSwitchingKey(tc.sk, sk2)

	for _, lvl := range tc.testLevel {
		t.Run(testString("SwitchKeys", tc.params, lvl), func(t *testing.T) {
			values, _, ciphertext := newTestVectorsLvl(lvl, 1, tc, tc.encryptorPk)
			tc.evaluator.SwitchKeys(ciphertext, switchingKey, ciphertext)
			verifyTestVectors(tc, decryptorSk2, values, ciphertext, t)
		})
	}

	for _, lvl := range tc.testLevel {
		t.Run(testString("SwitchKeysNew", tc.params, lvl), func(t *testing.T) {
			values, _, ciphertext := newTestVectorsLvl(lvl, 1, tc, tc.encryptorPk)
			ciphertext = tc.evaluator.SwitchKeysNew(ciphertext, switchingKey)
			verifyTestVectors(tc, decryptorSk2, values, ciphertext, t)
		})
	}
}

func testMarshalling(tc *testContext, t *testing.T) {
	t.Run("Marshalling", func(t *testing.T) {
		t.Run(testString("PowerBasis", tc.params, tc.params.MaxLevel()), func(t *testing.T) {

			if tc.params.MaxLevel() < 4 {
				t.Skip("not enough levels")
			}

			_, _, ct := newTestVectorsLvl(tc.params.MaxLevel(), 1, tc, tc.encryptorPk)

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
