package bfv

import (
	"encoding/json"
	"flag"
	"fmt"
	"math/big"
	"math/bits"
	"runtime"
	"testing"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/rlwe"
	"github.com/tuneinsight/lattigo/v4/utils"
)

var flagLongTest = flag.Bool("long", false, "run the long test suite (all parameters). Overrides -short and requires -timeout=0.")
var flagParamString = flag.String("params", "", "specify the test cryptographic parameters as a JSON string. Overrides -short and -long.")

func testString(opname string, p Parameters, lvl int) string {
	return fmt.Sprintf("%s/LogN=%d/logQP=%d/logT=%d/TIsQ0=%t/#Q=%d/#P=%d/lvl=%d", opname, p.LogN(), p.LogQP(), p.LogT(), p.T() == p.Q()[0], p.QCount(), p.PCount(), lvl)
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
	encryptorPk rlwe.Encryptor
	encryptorSk rlwe.Encryptor
	decryptor   rlwe.Decryptor
	evaluator   Evaluator
	testLevel   []int
}

var (
	// TESTTDivQN2Q1P is a set of test parameters where T = Q[0].
	TESTTDivQN2Q1P = ParametersLiteral{
		LogN: 14,
		Q:    []uint64{0x10001, 0xffffffffffe8001, 0xffffffffffd8001, 0xffffffffffc0001, 0xffffffffff28001},
		P:    []uint64{0x1fffffffffe10001, 0x1fffffffffe00001},
		T:    0x10001,
	}

	// TESTTCPrimeQN2Q1P is a set of test parameters where T is coprime with Q.
	TESTTCPrimeQN2Q1P = ParametersLiteral{
		LogN: 14,
		Q:    []uint64{0xffffffffffe8001, 0xffffffffffd8001, 0xffffffffffc0001, 0xffffffffff28001},
		P:    []uint64{0x1fffffffffe10001, 0x1fffffffffe00001},
		T:    0x10001,
	}

	// TestParams is a set of test parameters for BFV ensuring 128 bit security in the classic setting.
	TestParams = []ParametersLiteral{TESTTDivQN2Q1P, TESTTCPrimeQN2Q1P}
)

func TestBFV(t *testing.T) {

	var err error

	var paramsLiterals []ParametersLiteral

	paramsLiterals = append(TestParams, DefaultParams...) // the default test runs for ring degree N=2^12, 2^13, 2^14, 2^15

	if testing.Short() {
		paramsLiterals = TestParams
	}

	if *flagLongTest {
		paramsLiterals = append(paramsLiterals, DefaultPostQuantumParams...) // the long test suite runs for all default parameters
	}

	if *flagParamString != "" {
		var jsonParams ParametersLiteral
		if err = json.Unmarshal([]byte(*flagParamString), &jsonParams); err != nil {
			t.Fatal(err)
		}
		paramsLiterals = []ParametersLiteral{jsonParams} // the custom test suite reads the parameters from the -params flag
	}

	for _, p := range paramsLiterals[:1] {

		var params Parameters
		if params, err = NewParametersFromLiteral(p); err != nil {
			t.Fatal(err)
		}

		var tc *testContext
		if tc, err = genTestParams(params); err != nil {
			t.Fatal(err)
		}

		for _, testSet := range []func(tc *testContext, t *testing.T){
			testParameters,
			testScaler,
			testEncoder,
			testEncryptor,
			testEvaluator,
			testPolyEval,
			testEvaluatorRotate,
			testEvaluatorKeySwitch,
			testMarshaller,
		} {
			testSet(tc, t)
			runtime.GC()
		}
	}
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

	tc.testLevel = []int{params.MaxLevel()}
	if params.T() == params.Q()[0] {
		if params.MaxLevel() != 1 {
			tc.testLevel = append(tc.testLevel, 1)
		}
	} else {
		if 2*bits.Len64(params.T())+params.LogN() > bits.Len64(params.Q()[0]) {
			if params.MaxLevel() != 1 {
				tc.testLevel = append(tc.testLevel, 1)
			}
		} else {
			if params.MaxLevel() != 0 {
				tc.testLevel = append(tc.testLevel, 0)
			}
		}
	}

	return
}

func testParameters(tc *testContext, t *testing.T) {

	t.Run("Parameters/InverseGaloisElement", func(t *testing.T) {
		for i := 1; i < int(tc.params.N()/2); i++ {
			galEl := tc.params.GaloisElementForColumnRotationBy(i)
			mod := uint64(2 * tc.params.N())
			inv := tc.params.InverseGaloisElement(galEl)
			res := (inv * galEl) % mod
			assert.Equal(t, uint64(1), res)
		}
	})

	t.Run(testString("Parameters/CopyNew", tc.params, tc.params.MaxLevel()), func(t *testing.T) {
		params1, params2 := tc.params.CopyNew(), tc.params.CopyNew()
		assert.True(t, params1.Equals(tc.params) && params2.Equals(tc.params))
		params1.ringT, _ = ring.NewRing(tc.params.N(), []uint64{7})
		assert.False(t, params1.Equals(tc.params))
		assert.True(t, params2.Equals(tc.params))
	})
}

func newTestVectorsRingQLvl(level int, tc *testContext, encryptor rlwe.Encryptor, t *testing.T) (coeffs *ring.Poly, pt *rlwe.Plaintext, ct *rlwe.Ciphertext) {
	coeffs = tc.uSampler.ReadNew()
	pt = rlwe.NewPlaintext(tc.params.Parameters, level)
	tc.encoder.Encode(coeffs.Coeffs[0], pt)
	if encryptor != nil {
		ct = encryptor.EncryptNew(pt)
	}
	return
}

func newTestVectorsRingT(tc *testContext, t *testing.T) (coeffs *ring.Poly, pt *PlaintextRingT) {
	coeffs = tc.uSampler.ReadNew()
	pt = NewPlaintextRingT(tc.params)
	tc.encoder.EncodeRingT(coeffs.Coeffs[0], pt)
	return
}

func newTestVectorsMulLvl(level int, tc *testContext, t *testing.T) (coeffs *ring.Poly, pt *PlaintextMul) {
	coeffs = tc.uSampler.ReadNew()
	pt = NewPlaintextMul(tc.params, level)
	tc.encoder.EncodeMul(coeffs.Coeffs[0], pt)
	return
}

func verifyTestVectors(tc *testContext, decryptor rlwe.Decryptor, coeffs *ring.Poly, element rlwe.Operand, t *testing.T) {

	var coeffsTest []uint64

	switch el := element.(type) {
	case *rlwe.Plaintext, *PlaintextMul, *PlaintextRingT:
		coeffsTest = tc.encoder.DecodeUintNew(el)
	case *rlwe.Ciphertext:
		coeffsTest = tc.encoder.DecodeUintNew(decryptor.DecryptNew(el))
	default:
		t.Error("invalid test object to verify")
	}

	require.True(t, utils.EqualSliceUint64(coeffs.Coeffs[0], coeffsTest))
}

func testScaler(tc *testContext, t *testing.T) {
	t.Run(testString("Scaler/DivRoundByQOverT", tc.params, tc.params.MaxLevel()), func(t *testing.T) {

		T := tc.ringT.Modulus[0]
		ringQ := tc.ringQ

		scaler := NewRNSScaler(ringQ, T)

		coeffs := make([]*big.Int, ringQ.N)
		bigQ := ringQ.ModulusAtLevel[tc.params.MaxLevel()]
		for i := 0; i < ringQ.N; i++ {
			coeffs[i] = ring.RandInt(bigQ)
		}

		coeffsWant := make([]*big.Int, ringQ.N)
		bigT := ring.NewUint(T)
		for i := range coeffs {
			coeffsWant[i] = new(big.Int).Set(coeffs[i])
			coeffsWant[i].Mul(coeffsWant[i], bigT)
			ring.DivRound(coeffsWant[i], bigQ, coeffsWant[i])
			coeffsWant[i].Mod(coeffsWant[i], bigT)
		}

		polyQ := ringQ.NewPoly()
		polyT := ring.NewPoly(ringQ.N, 1)
		ringQ.SetCoefficientsBigint(coeffs, polyQ)

		scaler.DivByQOverTRoundedLvl(polyQ.Level(), polyQ, polyT)

		for i := 0; i < ringQ.N; i++ {
			require.Equal(t, polyT.Coeffs[0][i], coeffsWant[i].Uint64())
		}
	})
}

func testEncoder(tc *testContext, t *testing.T) {

	t.Run(testString("Encoder/Encode&Decode/RingT/Uint", tc.params, 0), func(t *testing.T) {
		values, plaintext := newTestVectorsRingT(tc, t)
		verifyTestVectors(tc, nil, values, plaintext, t)

		coeffsInt := make([]uint64, len(values.Coeffs[0]))
		for i, v := range values.Coeffs[0] {
			coeffsInt[i] = v + tc.params.T()*uint64(i%10)
		}

		plaintext = NewPlaintextRingT(tc.params)
		tc.encoder.EncodeRingT(coeffsInt, plaintext)
		//	coeffsTest := tc.encoder.DecodeIntNew(plaintext)

		verifyTestVectors(tc, nil, values, plaintext, t)
	})

	t.Run(testString("Encoder/Encode&Decode/RingT/Int", tc.params, 0), func(t *testing.T) {

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

		plaintext := NewPlaintextRingT(tc.params)
		tc.encoder.EncodeRingT(coeffsInt, plaintext)
		coeffsTest := tc.encoder.DecodeIntNew(plaintext)

		require.True(t, utils.EqualSliceInt64(coeffsInt, coeffsTest))
	})

	for _, lvl := range tc.testLevel {
		t.Run(testString("Encoder/Encode&Decode/RingQ/Uint", tc.params, lvl), func(t *testing.T) {
			values, plaintext, _ := newTestVectorsRingQLvl(lvl, tc, nil, t)
			verifyTestVectors(tc, nil, values, plaintext, t)
		})
	}

	for _, lvl := range tc.testLevel {
		t.Run(testString("Encoder/Encode&Decode/RingQ/Int", tc.params, lvl), func(t *testing.T) {

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

	for _, lvl := range tc.testLevel {
		t.Run(testString("Encoder/Encode&Decode/PlaintextMul", tc.params, lvl), func(t *testing.T) {
			values, plaintext := newTestVectorsMulLvl(lvl, tc, t)
			verifyTestVectors(tc, nil, values, plaintext, t)
		})
	}
}

func testEncryptor(tc *testContext, t *testing.T) {
	for _, lvl := range tc.testLevel {
		t.Run(testString("Encryptor/Encrypt/key=pk", tc.params, lvl), func(t *testing.T) {
			values1, _, ciphertext1 := newTestVectorsRingQLvl(lvl, tc, tc.encryptorPk, t)
			verifyTestVectors(tc, tc.decryptor, values1, ciphertext1, t)
		})
	}
	for _, lvl := range tc.testLevel {
		t.Run(testString("Encryptor/Encrypt/key=sk", tc.params, lvl), func(t *testing.T) {
			values1, _, ciphertext1 := newTestVectorsRingQLvl(lvl, tc, tc.encryptorSk, t)
			verifyTestVectors(tc, tc.decryptor, values1, ciphertext1, t)
		})
	}

	zero := tc.ringT.NewPoly()
	for _, lvl := range tc.testLevel {
		t.Run(testString("Encryptor/EncryptZero/key=pk", tc.params, lvl), func(t *testing.T) {
			ct := tc.encryptorPk.EncryptZeroNew(tc.params.MaxLevel())
			verifyTestVectors(tc, tc.decryptor, zero, ct, t)
		})
	}
	for _, lvl := range tc.testLevel {
		t.Run(testString("Encryptor/EncryptZero/key=sk", tc.params, lvl), func(t *testing.T) {
			ct := tc.encryptorSk.EncryptZeroNew(tc.params.MaxLevel())
			verifyTestVectors(tc, tc.decryptor, zero, ct, t)
		})
	}

	for _, lvl := range tc.testLevel {
		t.Run(testString("Encryptor/WithPRNG/Encrypt", tc.params, lvl), func(t *testing.T) {
			enc := NewPRNGEncryptor(tc.params, tc.sk)
			prng1, _ := utils.NewKeyedPRNG([]byte{'l'})
			prng2, _ := utils.NewKeyedPRNG([]byte{'l'})
			sampler := ring.NewUniformSampler(prng2, tc.ringQ)
			values1, pt, _ := newTestVectorsRingQLvl(lvl, tc, nil, t)
			ciphertext := enc.WithPRNG(prng1).EncryptNew(pt)
			c1Want := sampler.ReadLvlNew(lvl)
			tc.params.RingQ().InvNTTLvl(lvl, c1Want, c1Want)
			assert.True(t, c1Want.Equals(ciphertext.Value[1]))
			verifyTestVectors(tc, tc.decryptor, values1, ciphertext, t)
		})
	}

}

func testEvaluator(tc *testContext, t *testing.T) {

	for _, lvl := range tc.testLevel {
		t.Run(testString("Evaluator/Add/op1=Ciphertext/op2=Ciphertext", tc.params, lvl), func(t *testing.T) {
			values1, _, ciphertext1 := newTestVectorsRingQLvl(lvl, tc, tc.encryptorPk, t)
			values2, _, ciphertext2 := newTestVectorsRingQLvl(lvl, tc, tc.encryptorPk, t)
			tc.evaluator.Add(ciphertext1, ciphertext2, ciphertext1)
			tc.ringT.Add(values1, values2, values1)
			verifyTestVectors(tc, tc.decryptor, values1, ciphertext1, t)
		})
	}

	for _, lvl := range tc.testLevel {
		t.Run(testString("Evaluator/AddNoMod/op1=Ciphertext/op2=Ciphertext", tc.params, lvl), func(t *testing.T) {
			values1, _, ciphertext1 := newTestVectorsRingQLvl(lvl, tc, tc.encryptorPk, t)
			values2, _, ciphertext2 := newTestVectorsRingQLvl(lvl, tc, tc.encryptorPk, t)
			tc.evaluator.AddNoMod(ciphertext1, ciphertext2, ciphertext1)
			tc.evaluator.Reduce(ciphertext1, ciphertext1)
			tc.ringT.Add(values1, values2, values1)
			verifyTestVectors(tc, tc.decryptor, values1, ciphertext1, t)
		})
	}

	for _, lvl := range tc.testLevel {
		t.Run(testString("Evaluator/AddNew/op1=Ciphertext/op2=Ciphertext", tc.params, lvl), func(t *testing.T) {
			values1, _, ciphertext1 := newTestVectorsRingQLvl(lvl, tc, tc.encryptorPk, t)
			values2, _, ciphertext2 := newTestVectorsRingQLvl(lvl, tc, tc.encryptorPk, t)
			ciphertext1 = tc.evaluator.AddNew(ciphertext1, ciphertext2)
			tc.ringT.Add(values1, values2, values1)
			verifyTestVectors(tc, tc.decryptor, values1, ciphertext1, t)
		})
	}

	for _, lvl := range tc.testLevel {
		t.Run(testString("Evaluator/AddNoModNew/op1=Ciphertext/op2=Ciphertext", tc.params, lvl), func(t *testing.T) {
			values1, _, ciphertext1 := newTestVectorsRingQLvl(lvl, tc, tc.encryptorPk, t)
			values2, _, ciphertext2 := newTestVectorsRingQLvl(lvl, tc, tc.encryptorPk, t)
			ciphertext1 = tc.evaluator.AddNoModNew(ciphertext1, ciphertext2)
			ciphertext1 = tc.evaluator.ReduceNew(ciphertext1)
			tc.ringT.Add(values1, values2, values1)
			verifyTestVectors(tc, tc.decryptor, values1, ciphertext1, t)
		})
	}

	for _, lvl := range tc.testLevel {
		t.Run(testString("Evaluator/Add/op1=Ciphertext/op2=PlaintextRingT", tc.params, lvl), func(t *testing.T) {
			values1, plaintextRingT := newTestVectorsRingT(tc, t)
			values2, _, ciphertext := newTestVectorsRingQLvl(lvl, tc, tc.encryptorPk, t)
			ciphertextOut := NewCiphertext(tc.params, 1, lvl)
			tc.evaluator.Add(ciphertext, plaintextRingT, ciphertextOut)
			tc.ringT.Add(values1, values2, values2)
			verifyTestVectors(tc, tc.decryptor, values2, ciphertextOut, t)
		})
	}

	for _, lvl := range tc.testLevel {
		t.Run(testString("Evaluator/Add/op1=Ciphertext/op2=Plaintext", tc.params, lvl), func(t *testing.T) {
			values1, _, ciphertext1 := newTestVectorsRingQLvl(lvl, tc, tc.encryptorPk, t)
			values2, plaintext2, ciphertext2 := newTestVectorsRingQLvl(lvl, tc, tc.encryptorPk, t)
			tc.evaluator.Add(ciphertext1, plaintext2, ciphertext2)
			tc.ringT.Add(values1, values2, values2)
			verifyTestVectors(tc, tc.decryptor, values2, ciphertext2, t)
		})
	}

	for _, lvl := range tc.testLevel {
		t.Run(testString("Evaluator/Sub/op1=Ciphertext/op2=Ciphertext", tc.params, lvl), func(t *testing.T) {
			values1, _, ciphertext1 := newTestVectorsRingQLvl(lvl, tc, tc.encryptorPk, t)
			values2, _, ciphertext2 := newTestVectorsRingQLvl(lvl, tc, tc.encryptorPk, t)
			tc.evaluator.Sub(ciphertext1, ciphertext2, ciphertext1)
			tc.ringT.Sub(values1, values2, values1)
			verifyTestVectors(tc, tc.decryptor, values1, ciphertext1, t)
		})
	}

	for _, lvl := range tc.testLevel {
		t.Run(testString("Evaluator/SubNoMod/op1=Ciphertext/op2=Ciphertext", tc.params, lvl), func(t *testing.T) {
			values1, _, ciphertext1 := newTestVectorsRingQLvl(lvl, tc, tc.encryptorPk, t)
			values2, _, ciphertext2 := newTestVectorsRingQLvl(lvl, tc, tc.encryptorPk, t)
			tc.evaluator.SubNoMod(ciphertext1, ciphertext2, ciphertext1)
			tc.evaluator.Reduce(ciphertext1, ciphertext1)
			tc.ringT.Sub(values1, values2, values1)
			verifyTestVectors(tc, tc.decryptor, values1, ciphertext1, t)
		})
	}

	for _, lvl := range tc.testLevel {
		t.Run(testString("Evaluator/SubNew/op1=Ciphertext/op2=Ciphertext", tc.params, lvl), func(t *testing.T) {
			values1, _, ciphertext1 := newTestVectorsRingQLvl(lvl, tc, tc.encryptorPk, t)
			values2, _, ciphertext2 := newTestVectorsRingQLvl(lvl, tc, tc.encryptorPk, t)
			ciphertext1 = tc.evaluator.SubNew(ciphertext1, ciphertext2)
			tc.ringT.Sub(values1, values2, values1)
			verifyTestVectors(tc, tc.decryptor, values1, ciphertext1, t)
		})
	}

	for _, lvl := range tc.testLevel {
		t.Run(testString("Evaluator/Sub/op1=Ciphertext/op2=PlaintextRingT", tc.params, lvl), func(t *testing.T) {
			values1, plaintextRingT := newTestVectorsRingT(tc, t)
			values2, _, ciphertext := newTestVectorsRingQLvl(lvl, tc, tc.encryptorPk, t)
			ciphertextOut := NewCiphertext(tc.params, 1, lvl)
			plaintextWant := NewPlaintextRingT(tc.params)
			tc.evaluator.Sub(ciphertext, plaintextRingT, ciphertextOut)
			tc.ringT.Sub(values2, values1, plaintextWant.Value)
			verifyTestVectors(tc, tc.decryptor, plaintextWant.Value, ciphertextOut, t)
		})
	}

	for _, lvl := range tc.testLevel {
		t.Run(testString("Evaluator/Sub/op1=Ciphertext/op2=Plaintext", tc.params, lvl), func(t *testing.T) {
			values1, _, ciphertext1 := newTestVectorsRingQLvl(lvl, tc, tc.encryptorPk, t)
			values2, plaintext2, ciphertext2 := newTestVectorsRingQLvl(lvl, tc, tc.encryptorPk, t)
			valuesWant := tc.ringT.NewPoly()
			tc.evaluator.Sub(ciphertext1, plaintext2, ciphertext2)
			tc.ringT.Sub(values1, values2, valuesWant)
			verifyTestVectors(tc, tc.decryptor, valuesWant, ciphertext2, t)
		})
	}

	for _, lvl := range tc.testLevel {
		t.Run(testString("Evaluator/SubNoModNew/op1=Ciphertext/op2=Ciphertext", tc.params, lvl), func(t *testing.T) {
			values1, _, ciphertext1 := newTestVectorsRingQLvl(lvl, tc, tc.encryptorPk, t)
			values2, _, ciphertext2 := newTestVectorsRingQLvl(lvl, tc, tc.encryptorPk, t)
			ciphertext1 = tc.evaluator.SubNoModNew(ciphertext1, ciphertext2)
			ciphertext1 = tc.evaluator.ReduceNew(ciphertext1)
			tc.ringT.Sub(values1, values2, values1)
			verifyTestVectors(tc, tc.decryptor, values1, ciphertext1, t)
		})
	}

	for _, lvl := range tc.testLevel {
		t.Run(testString("Evaluator/SubNoMod/op1=Ciphertext/op2=Plaintext", tc.params, lvl), func(t *testing.T) {
			values1, _, ciphertext1 := newTestVectorsRingQLvl(lvl, tc, tc.encryptorPk, t)
			values2, plaintext2, ciphertext2 := newTestVectorsRingQLvl(lvl, tc, tc.encryptorPk, t)
			valuesWant := tc.ringT.NewPoly()
			tc.evaluator.SubNoMod(ciphertext1, plaintext2, ciphertext2)
			tc.evaluator.Reduce(ciphertext2, ciphertext2)
			tc.ringT.Sub(values1, values2, valuesWant)
			verifyTestVectors(tc, tc.decryptor, valuesWant, ciphertext2, t)
		})
	}

	for _, lvl := range tc.testLevel {
		t.Run(testString("Evaluator/Neg", tc.params, lvl), func(t *testing.T) {
			values1, _, ciphertext1 := newTestVectorsRingQLvl(lvl, tc, tc.encryptorPk, t)
			tc.evaluator.Neg(ciphertext1, ciphertext1)
			tc.ringT.Neg(values1, values1)
			tc.ringT.Reduce(values1, values1)
			verifyTestVectors(tc, tc.decryptor, values1, ciphertext1, t)
		})
	}

	for _, lvl := range tc.testLevel {
		t.Run(testString("Evaluator/NegNew", tc.params, lvl), func(t *testing.T) {
			values1, _, ciphertext1 := newTestVectorsRingQLvl(lvl, tc, tc.encryptorPk, t)
			ciphertext1 = tc.evaluator.NegNew(ciphertext1)
			tc.ringT.Neg(values1, values1)
			tc.ringT.Reduce(values1, values1)
			verifyTestVectors(tc, tc.decryptor, values1, ciphertext1, t)
		})
	}

	for _, lvl := range tc.testLevel {
		t.Run(testString("Evaluator/AddScalar", tc.params, lvl), func(t *testing.T) {
			values1, _, ciphertext1 := newTestVectorsRingQLvl(lvl, tc, tc.encryptorPk, t)
			tc.evaluator.AddScalar(ciphertext1, 37, ciphertext1)
			tc.ringT.AddScalar(values1, 37, values1)
			verifyTestVectors(tc, tc.decryptor, values1, ciphertext1, t)
		})
	}

	for _, lvl := range tc.testLevel {
		t.Run(testString("Evaluator/MulScalar", tc.params, lvl), func(t *testing.T) {
			values1, _, ciphertext1 := newTestVectorsRingQLvl(lvl, tc, tc.encryptorPk, t)
			tc.evaluator.MulScalar(ciphertext1, 37, ciphertext1)
			tc.ringT.MulScalar(values1, 37, values1)
			verifyTestVectors(tc, tc.decryptor, values1, ciphertext1, t)
		})
	}

	for _, lvl := range tc.testLevel {
		t.Run(testString("Evaluator/MulScalarNew", tc.params, lvl), func(t *testing.T) {
			values1, _, ciphertext1 := newTestVectorsRingQLvl(lvl, tc, tc.encryptorPk, t)
			ciphertext1 = tc.evaluator.MulScalarNew(ciphertext1, 37)
			tc.ringT.MulScalar(values1, 37, values1)
			verifyTestVectors(tc, tc.decryptor, values1, ciphertext1, t)
		})
	}

	for _, lvl := range tc.testLevel {

		if lvl == 0 && tc.params.MaxLevel() > 0 {
			lvl++
		}

		t.Run(testString("Evaluator/Mul/op1=Ciphertext/op2=Ciphertext", tc.params, lvl), func(t *testing.T) {
			values1, _, ciphertext1 := newTestVectorsRingQLvl(lvl, tc, tc.encryptorPk, t)
			values2, _, ciphertext2 := newTestVectorsRingQLvl(lvl, tc, tc.encryptorPk, t)
			receiver := NewCiphertext(tc.params, ciphertext1.Degree()+ciphertext2.Degree(), lvl)
			tc.evaluator.Mul(ciphertext1, ciphertext2, receiver)
			tc.ringT.MulCoeffs(values1, values2, values1)
			verifyTestVectors(tc, tc.decryptor, values1, receiver, t)
		})
	}

	for _, lvl := range tc.testLevel {

		if lvl == 0 && tc.params.MaxLevel() > 0 {
			lvl++
		}

		t.Run(testString("Evaluator/MulAndAdd/op1=Ciphertext/op2=Plaintext", tc.params, lvl), func(t *testing.T) {
			values1, _, ciphertext1 := newTestVectorsRingQLvl(lvl, tc, tc.encryptorPk, t)
			values2, plaintext2, ciphertext2 := newTestVectorsRingQLvl(lvl, tc, tc.encryptorPk, t)
			tc.evaluator.MulAndAdd(ciphertext1, plaintext2, ciphertext2)
			tmp := tc.ringT.NewPoly()
			tc.ringT.MulCoeffs(values1, values2, tmp)
			tc.ringT.Add(values2, tmp, values2)
			verifyTestVectors(tc, tc.decryptor, values2, ciphertext2, t)
		})
	}

	for _, lvl := range tc.testLevel {

		if lvl == 0 && tc.params.MaxLevel() > 0 {
			lvl++
		}

		t.Run(testString("Evaluator/MulAndAdd/op1=Ciphertext/op2=Ciphetext", tc.params, lvl), func(t *testing.T) {
			values1, _, ciphertext1 := newTestVectorsRingQLvl(lvl, tc, tc.encryptorPk, t)
			values2, _, ciphertext2 := newTestVectorsRingQLvl(lvl, tc, tc.encryptorPk, t)
			values3, _, ciphertext3 := newTestVectorsRingQLvl(lvl, tc, tc.encryptorPk, t)
			ciphertext3.Resize(2, ciphertext3.Level())
			tc.evaluator.MulAndAdd(ciphertext1, ciphertext2, ciphertext3)
			tc.ringT.MulCoeffs(values1, values2, values1)
			tc.ringT.Add(values3, values1, values3)
			verifyTestVectors(tc, tc.decryptor, values3, ciphertext3, t)
		})
	}

	for _, lvl := range tc.testLevel {

		if lvl == 0 && tc.params.MaxLevel() > 0 {
			lvl++
		}

		t.Run(testString("Evaluator/MulNew/op1=Ciphertext/op2=Ciphertext", tc.params, lvl), func(t *testing.T) {
			values1, _, ciphertext1 := newTestVectorsRingQLvl(lvl, tc, tc.encryptorPk, t)
			values2, _, ciphertext2 := newTestVectorsRingQLvl(lvl, tc, tc.encryptorPk, t)
			receiver := tc.evaluator.MulNew(ciphertext1, ciphertext2)
			tc.ringT.MulCoeffs(values1, values2, values1)
			verifyTestVectors(tc, tc.decryptor, values1, receiver, t)
		})
	}

	for _, lvl := range tc.testLevel {

		if lvl == 0 && tc.params.MaxLevel() > 0 {
			lvl++
		}

		t.Run(testString("Evaluator/MulSquare/op1=Ciphertext/op2=Ciphertext", tc.params, lvl), func(t *testing.T) {
			values1, _, ciphertext1 := newTestVectorsRingQLvl(lvl, tc, tc.encryptorPk, t)
			receiver := NewCiphertext(tc.params, ciphertext1.Degree()+ciphertext1.Degree(), lvl)
			tc.evaluator.Mul(ciphertext1, ciphertext1, receiver)
			tc.ringT.MulCoeffs(values1, values1, values1)
			verifyTestVectors(tc, tc.decryptor, values1, receiver, t)
		})
	}

	for _, lvl := range tc.testLevel {

		if lvl == 0 && tc.params.MaxLevel() > 0 {
			lvl++
		}

		t.Run(testString("Evaluator/Mul/op1=Ciphertext/op2=Plaintext", tc.params, lvl), func(t *testing.T) {
			values1, _, ciphertext1 := newTestVectorsRingQLvl(lvl, tc, tc.encryptorPk, t)
			values2, plaintext2, _ := newTestVectorsRingQLvl(lvl, tc, nil, t)
			tc.evaluator.Mul(ciphertext1, plaintext2, ciphertext1)
			tc.ringT.MulCoeffs(values1, values2, values1)
			verifyTestVectors(tc, tc.decryptor, values1, ciphertext1, t)
		})
	}

	for _, lvl := range tc.testLevel {

		if lvl == 0 && tc.params.MaxLevel() > 0 {
			lvl++
		}

		t.Run(testString("Evaluator/Mul/op1=Ciphertext/op2=PlaintextRingT", tc.params, lvl), func(t *testing.T) {
			values1, plaintextRingT := newTestVectorsRingT(tc, t)
			values2, _, ciphertext := newTestVectorsRingQLvl(lvl, tc, tc.encryptorPk, t)
			ciphertextOut := NewCiphertext(tc.params, 1, lvl)
			tc.evaluator.Mul(ciphertext, plaintextRingT, ciphertextOut)
			tc.ringT.MulCoeffs(values1, values2, values1)
			verifyTestVectors(tc, tc.decryptor, values1, ciphertextOut, t)
		})
	}

	for _, lvl := range tc.testLevel {

		if lvl == 0 && tc.params.MaxLevel() > 0 {
			lvl++
		}

		t.Run(testString("Evaluator/Mul/op1=Ciphertext/op2=PlaintextMul", tc.params, lvl), func(t *testing.T) {
			values1, _, ciphertext1 := newTestVectorsRingQLvl(lvl, tc, tc.encryptorPk, t)
			values2, plaintext2 := newTestVectorsMulLvl(lvl, tc, t)
			tc.evaluator.Mul(ciphertext1, plaintext2, ciphertext1)
			tc.ringT.MulCoeffs(values1, values2, values1)
			verifyTestVectors(tc, tc.decryptor, values1, ciphertext1, t)
		})
	}

	for _, lvl := range tc.testLevel {

		if lvl == 0 && tc.params.MaxLevel() > 0 {
			lvl++
		}

		t.Run(testString("Evaluator/Mul/Relinearize", tc.params, lvl), func(t *testing.T) {
			values1, _, ciphertext1 := newTestVectorsRingQLvl(lvl, tc, tc.encryptorPk, t)
			values2, _, ciphertext2 := newTestVectorsRingQLvl(lvl, tc, tc.encryptorPk, t)
			receiver := NewCiphertext(tc.params, ciphertext1.Degree()+ciphertext2.Degree(), lvl)
			tc.evaluator.Mul(ciphertext1, ciphertext2, receiver)
			tc.ringT.MulCoeffs(values1, values2, values1)
			receiver2 := tc.evaluator.RelinearizeNew(receiver)
			verifyTestVectors(tc, tc.decryptor, values1, receiver2, t)
			tc.evaluator.Relinearize(receiver, receiver)
			verifyTestVectors(tc, tc.decryptor, values1, receiver, t)
		})
	}

	t.Run(testString("Evaluator/RescaleTo", tc.params, 1), func(t *testing.T) {
		values1, _, ciphertext1 := newTestVectorsRingQLvl(tc.params.MaxLevel(), tc, tc.encryptorPk, t)
		tc.evaluator.RescaleTo(1, ciphertext1, ciphertext1)
		verifyTestVectors(tc, tc.decryptor, values1, ciphertext1, t)
		if tc.params.T() != tc.params.RingQ().Modulus[0] { // only happens if T divides Q.
			tc.evaluator.RescaleTo(0, ciphertext1, ciphertext1)
			verifyTestVectors(tc, tc.decryptor, values1, ciphertext1, t)
		}
	})

	t.Run(testString("Evaluator/RescaleTo/ThenAdd", tc.params, 1), func(t *testing.T) {
		values1, _, ciphertext1 := newTestVectorsRingQLvl(tc.params.MaxLevel(), tc, tc.encryptorPk, t)
		values2, _, ciphertext2 := newTestVectorsRingQLvl(tc.params.MaxLevel(), tc, tc.encryptorPk, t)
		tc.evaluator.RescaleTo(1, ciphertext1, ciphertext1)
		tc.evaluator.RescaleTo(1, ciphertext2, ciphertext2)
		assert.True(t, ciphertext1.Level() == 1)
		assert.True(t, ciphertext2.Level() == 1)
		tc.evaluator.Add(ciphertext1, ciphertext2, ciphertext1)
		tc.ringT.Add(values1, values2, values1)
		verifyTestVectors(tc, tc.decryptor, values1, ciphertext1, t)
	})

	t.Run(testString("Evaluator/RescaleTo/MulRelin", tc.params, 1), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectorsRingQLvl(tc.params.MaxLevel(), tc, tc.encryptorPk, t)
		values2, _, ciphertext2 := newTestVectorsRingQLvl(tc.params.MaxLevel(), tc, tc.encryptorPk, t)
		tc.evaluator.RescaleTo(1, ciphertext1, ciphertext1)
		tc.evaluator.RescaleTo(1, ciphertext2, ciphertext2)
		assert.True(t, ciphertext1.Level() == 1)
		assert.True(t, ciphertext2.Level() == 1)
		receiver := NewCiphertext(tc.params, ciphertext1.Degree()+ciphertext2.Degree(), ciphertext2.Level())
		tc.evaluator.Mul(ciphertext1, ciphertext2, receiver)
		tc.ringT.MulCoeffs(values1, values2, values1)
		verifyTestVectors(tc, tc.decryptor, values1, receiver, t)
		receiver2 := tc.evaluator.RelinearizeNew(receiver)
		verifyTestVectors(tc, tc.decryptor, values1, receiver2, t)
	})
}

func testPolyEval(tc *testContext, t *testing.T) {

	t.Run(testString("PowerBasis/Marshalling", tc.params, tc.params.MaxLevel()), func(t *testing.T) {

		_, _, ct := newTestVectorsRingQLvl(tc.params.MaxLevel(), tc, tc.encryptorPk, t)

		pb := NewPowerBasis(ct)

		for i := 2; i < 4; i++ {
			pb.GenPower(i, tc.evaluator)
		}

		pbBytes, err := pb.MarshalBinary()
		assert.Nil(t, err)
		pbNew := new(PowerBasis)
		assert.Nil(t, pbNew.UnmarshalBinary(pbBytes))

		for i := range pb.Value {
			ctWant := pb.Value[i]
			ctHave := pbNew.Value[i]
			require.NotNil(t, ctHave)
			for j := range ct.Value {
				require.True(t, tc.ringQ.Equal(ctWant.Value[j], ctHave.Value[j]))
			}
		}
	})

	for _, lvl := range []int{tc.params.MaxLevel(), tc.params.MaxLevel() - 1} {
		t.Run(testString("PolyEval/Single", tc.params, lvl), func(t *testing.T) {

			if (tc.params.LogQ()-tc.params.LogT())/(tc.params.LogT()+tc.params.LogN()) < 5 {
				t.Skip("Homomorphic Capacity Too Low")
			}

			values, _, ciphertext := newTestVectorsRingQLvl(tc.params.MaxLevel()-1, tc, tc.encryptorPk, t)

			coeffs := []uint64{1, 2, 3, 4, 5, 6, 7, 8}

			T := tc.ringT.Modulus[0]
			for i := range values.Coeffs[0] {
				values.Coeffs[0][i] = ring.EvalPolyModP(values.Coeffs[0][i], coeffs, T)
			}

			poly := NewPoly(coeffs)

			var err error
			if ciphertext, err = tc.evaluator.EvaluatePoly(ciphertext, poly); err != nil {
				t.Fatal(err)
			}

			verifyTestVectors(tc, tc.decryptor, values, ciphertext, t)
		})
	}

	for _, lvl := range []int{tc.params.MaxLevel(), tc.params.MaxLevel() - 1} {
		t.Run(testString("PolyEval/Vector", tc.params, lvl), func(t *testing.T) {

			if (tc.params.LogQ()-tc.params.LogT()-tc.params.LogN())/(tc.params.LogT()+tc.params.LogN()) < 5 {

				t.Skip("Homomorphic Capacity Too Low")
			}

			values, _, ciphertext := newTestVectorsRingQLvl(lvl, tc, tc.encryptorPk, t)

			coeffs0 := []uint64{1, 2, 3, 4, 5, 6, 7, 8}
			coeffs1 := []uint64{2, 3, 4, 5, 6, 7, 8, 9}

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
			if ciphertext, err = tc.evaluator.EvaluatePolyVector(ciphertext, polyVec, tc.encoder, slotIndex); err != nil {
				t.Fatal(err)
			}

			verifyTestVectors(tc, tc.decryptor, values, ciphertext, t)
		})
	}
}

func testEvaluatorKeySwitch(tc *testContext, t *testing.T) {

	sk2 := tc.kgen.GenSecretKey()
	decryptorSk2 := NewDecryptor(tc.params, sk2)
	switchKey := tc.kgen.GenSwitchingKey(tc.sk, sk2)

	for _, lvl := range tc.testLevel {
		t.Run(testString("Evaluator/KeySwitch/InPlace", tc.params, lvl), func(t *testing.T) {
			values, _, ciphertext := newTestVectorsRingQLvl(lvl, tc, tc.encryptorPk, t)
			tc.evaluator.SwitchKeys(ciphertext, switchKey, ciphertext)
			verifyTestVectors(tc, decryptorSk2, values, ciphertext, t)
		})
	}

	for _, lvl := range tc.testLevel {
		t.Run(testString("Evaluator/KeySwitch/New", tc.params, lvl), func(t *testing.T) {
			values, _, ciphertext := newTestVectorsRingQLvl(lvl, tc, tc.encryptorPk, t)
			ciphertext = tc.evaluator.SwitchKeysNew(ciphertext, switchKey)
			verifyTestVectors(tc, decryptorSk2, values, ciphertext, t)
		})
	}

}

func testEvaluatorRotate(tc *testContext, t *testing.T) {

	rots := []int{1, -1, 4, -4, 63, -63}
	rotkey := tc.kgen.GenRotationKeysForRotations(rots, true, tc.sk)
	evaluator := tc.evaluator.WithKey(rlwe.EvaluationKey{Rlk: tc.rlk, Rtks: rotkey})

	for _, lvl := range tc.testLevel {
		t.Run(testString("Evaluator/RotateRows", tc.params, lvl), func(t *testing.T) {
			values, _, ciphertext := newTestVectorsRingQLvl(lvl, tc, tc.encryptorPk, t)
			evaluator.RotateRows(ciphertext, ciphertext)
			values.Coeffs[0] = append(values.Coeffs[0][tc.params.N()>>1:], values.Coeffs[0][:tc.params.N()>>1]...)
			verifyTestVectors(tc, tc.decryptor, values, ciphertext, t)
		})
	}

	for _, lvl := range tc.testLevel {
		t.Run(testString("Evaluator/RotateRowsNew", tc.params, lvl), func(t *testing.T) {
			values, _, ciphertext := newTestVectorsRingQLvl(lvl, tc, tc.encryptorPk, t)
			ciphertext = evaluator.RotateRowsNew(ciphertext)
			values.Coeffs[0] = append(values.Coeffs[0][tc.params.N()>>1:], values.Coeffs[0][:tc.params.N()>>1]...)
			verifyTestVectors(tc, tc.decryptor, values, ciphertext, t)
		})
	}

	for _, lvl := range tc.testLevel {
		t.Run(testString("Evaluator/RotateColumns", tc.params, lvl), func(t *testing.T) {

			values, _, ciphertext := newTestVectorsRingQLvl(lvl, tc, tc.encryptorPk, t)

			receiver := NewCiphertext(tc.params, 1, lvl)
			for _, n := range rots {

				evaluator.RotateColumns(ciphertext, n, receiver)
				valuesWant := utils.RotateUint64Slots(values.Coeffs[0], n)

				verifyTestVectors(tc, tc.decryptor, &ring.Poly{Coeffs: [][]uint64{valuesWant}}, receiver, t)
			}
		})
	}

	for _, lvl := range tc.testLevel {
		t.Run(testString("Evaluator/RotateColumnsNew", tc.params, lvl), func(t *testing.T) {

			values, _, ciphertext := newTestVectorsRingQLvl(lvl, tc, tc.encryptorPk, t)

			for _, n := range rots {

				receiver := evaluator.RotateColumnsNew(ciphertext, n)
				valuesWant := utils.RotateUint64Slots(values.Coeffs[0], n)

				verifyTestVectors(tc, tc.decryptor, &ring.Poly{Coeffs: [][]uint64{valuesWant}}, receiver, t)
			}
		})
	}

	rotkey = tc.kgen.GenRotationKeysForInnerSum(tc.sk)
	evaluator = evaluator.WithKey(rlwe.EvaluationKey{Rlk: tc.rlk, Rtks: rotkey})

	for _, lvl := range tc.testLevel {
		t.Run(testString("Evaluator/Rotate/InnerSum", tc.params, lvl), func(t *testing.T) {
			values, _, ciphertext := newTestVectorsRingQLvl(lvl, tc, tc.encryptorPk, t)

			evaluator.InnerSum(ciphertext, ciphertext)

			var sum uint64
			for _, c := range values.Coeffs[0] {
				sum += c
			}

			sum %= tc.params.T()

			for i := range values.Coeffs[0] {
				values.Coeffs[0][i] = sum
			}
			verifyTestVectors(tc, tc.decryptor, values, ciphertext, t)
		})
	}

	t.Run(testString("Evaluator/RescaleTo/Rotate", tc.params, tc.params.MaxLevel()), func(t *testing.T) {
		values1, _, ciphertext1 := newTestVectorsRingQLvl(tc.params.MaxLevel(), tc, tc.encryptorPk, t)
		rotkey := tc.kgen.GenRotationKeysForRotations(nil, true, tc.sk)
		evaluator := tc.evaluator.WithKey(rlwe.EvaluationKey{Rlk: tc.rlk, Rtks: rotkey})
		tc.evaluator.RescaleTo(1, ciphertext1, ciphertext1)
		assert.True(t, ciphertext1.Level() == 1)
		evaluator.RotateRows(ciphertext1, ciphertext1)
		values1.Coeffs[0] = append(values1.Coeffs[0][tc.params.N()>>1:], values1.Coeffs[0][:tc.params.N()>>1]...)
		verifyTestVectors(tc, tc.decryptor, values1, ciphertext1, t)
	})
}

func testMarshaller(tc *testContext, t *testing.T) {

	t.Run(testString("Marshaller/Parameters/Binary", tc.params, tc.params.MaxLevel()), func(t *testing.T) {
		bytes, err := tc.params.MarshalBinary()
		assert.Nil(t, err)
		var p Parameters
		err = p.UnmarshalBinary(bytes)
		assert.Nil(t, err)
		assert.Equal(t, tc.params, p)
		assert.Equal(t, tc.params.MarshalBinarySize(), len(bytes))
	})

	t.Run(testString("Marshaller/Parameters/JSON", tc.params, tc.params.MaxLevel()), func(t *testing.T) {
		// checks that parameters can be marshalled without error
		data, err := json.Marshal(tc.params)
		assert.Nil(t, err)
		assert.NotNil(t, data)

		// checks that bfv.Parameters can be unmarshalled without error
		var paramsRec Parameters
		err = json.Unmarshal(data, &paramsRec)
		assert.Nil(t, err)
		assert.True(t, tc.params.Equals(paramsRec))

		// checks that bfv.Parameters can be unmarshalled with log-moduli definition without error
		dataWithLogModuli := []byte(fmt.Sprintf(`{"LogN":%d,"LogQ":[50,50],"LogP":[60], "T":65537}`, tc.params.LogN()))
		var paramsWithLogModuli Parameters
		err = json.Unmarshal(dataWithLogModuli, &paramsWithLogModuli)
		assert.Nil(t, err)
		assert.Equal(t, 2, paramsWithLogModuli.QCount())
		assert.Equal(t, 1, paramsWithLogModuli.PCount())
		assert.Equal(t, rlwe.DefaultSigma, paramsWithLogModuli.Sigma()) // ommiting sigma should result in Default being used

		// checks that bfv.Parameters can be unmarshalled with log-moduli definition with empty P without error
		dataWithLogModuliNoP := []byte(fmt.Sprintf(`{"LogN":%d,"LogQ":[50,50],"LogP":[],"T":65537}`, tc.params.LogN()))
		var paramsWithLogModuliNoP Parameters
		err = json.Unmarshal(dataWithLogModuliNoP, &paramsWithLogModuliNoP)
		assert.Nil(t, err)
		assert.Equal(t, 2, paramsWithLogModuliNoP.QCount())
		assert.Equal(t, 0, paramsWithLogModuliNoP.PCount())

		// checks that one can provide custom parameters for the secret-key and error distributions
		dataWithCustomSecrets := []byte(fmt.Sprintf(`{"LogN":%d,"LogQ":[50,50],"LogP":[60], "H": 192, "Sigma": 6.6,"T":65537}`, tc.params.LogN()))
		var paramsWithCustomSecrets Parameters
		err = json.Unmarshal(dataWithCustomSecrets, &paramsWithCustomSecrets)
		assert.Nil(t, err)
		assert.Equal(t, 6.6, paramsWithCustomSecrets.Sigma())
		assert.Equal(t, 192, paramsWithCustomSecrets.HammingWeight())

	})
}
