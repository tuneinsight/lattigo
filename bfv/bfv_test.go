package bfv

import (
	"encoding/json"
	"flag"
	"fmt"
	"math/big"
	"runtime"
	"testing"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
	"github.com/tuneinsight/lattigo/v3/ring"
	"github.com/tuneinsight/lattigo/v3/rlwe"
	"github.com/tuneinsight/lattigo/v3/utils"
)

var flagLongTest = flag.Bool("long", false, "run the long test suite (all parameters). Overrides -short and requires -timeout=0.")
var flagParamString = flag.String("params", "", "specify the test cryptographic parameters as a JSON string. Overrides -short and -long.")

func testString(opname string, p Parameters) string {
	return fmt.Sprintf("%s/LogN=%d/logQ=%d/alpha=%d/beta=%d", opname, p.LogN(), p.LogQP(), p.PCount(), p.Beta())
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
}

func TestBFV(t *testing.T) {

	defaultParams := DefaultParams // the default test runs for ring degree N=2^12, 2^13, 2^14, 2^15
	if testing.Short() {
		defaultParams = DefaultParams[:2] // the short test suite runs for ring degree N=2^12, 2^13
	}
	if *flagLongTest {
		defaultParams = append(defaultParams, DefaultPostQuantumParams...) // the long test suite runs for all default parameters
	}
	if *flagParamString != "" {
		var jsonParams ParametersLiteral
		json.Unmarshal([]byte(*flagParamString), &jsonParams)
		defaultParams = []ParametersLiteral{jsonParams} // the custom test suite reads the parameters from the -params flag
	}

	for _, p := range defaultParams[:] {

		params, err := NewParametersFromLiteral(p)
		if err != nil {
			panic(err)
		}
		var tc *testContext
		if tc, err = genTestParams(params); err != nil {
			panic(err)
		}

		for _, testSet := range []func(tc *testContext, t *testing.T){
			testParameters,
			testEncoder,
			testEvaluator,
			testPolyEval,
			testEvaluatorKeySwitch,
			testEvaluatorRotate,
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
	if params.PCount() != 0 {
		tc.rlk = tc.kgen.GenRelinearizationKey(tc.sk, 1)
	}

	tc.encoder = NewEncoder(tc.params)
	tc.encryptorPk = NewEncryptor(tc.params, tc.pk)
	tc.encryptorSk = NewEncryptor(tc.params, tc.sk)
	tc.decryptor = NewDecryptor(tc.params, tc.sk)
	tc.evaluator = NewEvaluator(tc.params, rlwe.EvaluationKey{Rlk: tc.rlk})
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

	t.Run(testString("Parameters/CopyNew", tc.params), func(t *testing.T) {
		params1, params2 := tc.params.CopyNew(), tc.params.CopyNew()
		assert.True(t, params1.Equals(tc.params) && params2.Equals(tc.params))
		params1.ringT, _ = ring.NewRing(tc.params.N(), []uint64{7})
		assert.False(t, params1.Equals(tc.params))
		assert.True(t, params2.Equals(tc.params))
	})
}

func newTestVectorsRingQ(tc *testContext, encryptor Encryptor, t *testing.T) (coeffs *ring.Poly, plaintext *Plaintext, ciphertext *Ciphertext) {

	coeffs = tc.uSampler.ReadNew()

	plaintext = NewPlaintext(tc.params)

	tc.encoder.EncodeUint(coeffs.Coeffs[0], plaintext)

	if encryptor != nil {
		ciphertext = encryptor.EncryptNew(plaintext)
	}

	return coeffs, plaintext, ciphertext
}

func newTestVectorsRingT(tc *testContext, t *testing.T) (coeffs *ring.Poly, plaintext *PlaintextRingT) {

	coeffs = tc.uSampler.ReadNew()

	plaintext = NewPlaintextRingT(tc.params)

	tc.encoder.EncodeUintRingT(coeffs.Coeffs[0], plaintext)

	return coeffs, plaintext
}

func newTestVectorsMul(tc *testContext, t *testing.T) (coeffs *ring.Poly, plaintext *PlaintextMul) {

	coeffs = tc.uSampler.ReadNew()

	plaintext = NewPlaintextMul(tc.params)

	tc.encoder.EncodeUintMul(coeffs.Coeffs[0], plaintext)

	return coeffs, plaintext
}

func verifyTestVectors(tc *testContext, decryptor Decryptor, coeffs *ring.Poly, element Operand, t *testing.T) {

	var coeffsTest []uint64

	switch el := element.(type) {
	case *Plaintext, *PlaintextMul, *PlaintextRingT:
		coeffsTest = tc.encoder.DecodeUintNew(el)
	case *Ciphertext:
		coeffsTest = tc.encoder.DecodeUintNew(decryptor.DecryptNew(el))
	default:
		t.Error("invalid test object to verify")
	}

	require.True(t, utils.EqualSliceUint64(coeffs.Coeffs[0], coeffsTest))
}

func testEncoder(tc *testContext, t *testing.T) {

	t.Run(testString("Encoder/Scaling", tc.params), func(t *testing.T) {

		T := tc.ringT.Modulus[0]
		ringQ := tc.ringQ

		scaler := NewRNSScaler(ringQ, T)

		coeffs := make([]*big.Int, ringQ.N)
		for i := 0; i < ringQ.N; i++ {
			coeffs[i] = ring.RandInt(ringQ.ModulusBigint)
		}

		coeffsWant := make([]*big.Int, ringQ.N)
		bigT := ring.NewUint(T)
		for i := range coeffs {
			coeffsWant[i] = new(big.Int).Set(coeffs[i])
			coeffsWant[i].Mul(coeffsWant[i], bigT)
			ring.DivRound(coeffsWant[i], ringQ.ModulusBigint, coeffsWant[i])
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

	t.Run(testString("Encoder/Encode&Decode/RingT/Uint", tc.params), func(t *testing.T) {
		values, plaintext := newTestVectorsRingT(tc, t)
		verifyTestVectors(tc, nil, values, plaintext, t)
	})

	t.Run(testString("Encoder/Encode&Decode/RingT/Int", tc.params), func(t *testing.T) {

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
		tc.encoder.EncodeIntRingT(coeffsInt, plaintext)
		coeffsTest := tc.encoder.DecodeIntNew(plaintext)

		require.True(t, utils.EqualSliceInt64(coeffsInt, coeffsTest))
	})

	t.Run(testString("Encoder/Encode&Decode/RingQ/Uint", tc.params), func(t *testing.T) {
		values, plaintext, _ := newTestVectorsRingQ(tc, nil, t)
		verifyTestVectors(tc, nil, values, plaintext, t)
	})

	t.Run(testString("Encoder/Encode&Decode/RingQ/Int", tc.params), func(t *testing.T) {

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
		plaintext := NewPlaintext(tc.params)
		tc.encoder.EncodeInt(coeffsInt, plaintext)
		coeffsTest := tc.encoder.DecodeIntNew(plaintext)

		require.True(t, utils.EqualSliceInt64(coeffsInt, coeffsTest))
	})

	t.Run(testString("Encoder/Encode&Decode/PlaintextMul", tc.params), func(t *testing.T) {
		values, plaintext := newTestVectorsMul(tc, t)
		verifyTestVectors(tc, nil, values, plaintext, t)
	})
}

func testEvaluator(tc *testContext, t *testing.T) {

	t.Run(testString("Evaluator/Add/op1=Ciphertext/op2=Ciphertext", tc.params), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectorsRingQ(tc, tc.encryptorPk, t)
		values2, _, ciphertext2 := newTestVectorsRingQ(tc, tc.encryptorPk, t)

		tc.evaluator.Add(ciphertext1, ciphertext2, ciphertext1)
		tc.ringT.Add(values1, values2, values1)

		verifyTestVectors(tc, tc.decryptor, values1, ciphertext1, t)
	})

	t.Run(testString("Evaluator/AddNoMod/op1=Ciphertext/op2=Ciphertext", tc.params), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectorsRingQ(tc, tc.encryptorPk, t)
		values2, _, ciphertext2 := newTestVectorsRingQ(tc, tc.encryptorPk, t)

		tc.evaluator.AddNoMod(ciphertext1, ciphertext2, ciphertext1)
		tc.evaluator.Reduce(ciphertext1, ciphertext1)
		tc.ringT.Add(values1, values2, values1)

		verifyTestVectors(tc, tc.decryptor, values1, ciphertext1, t)
	})

	t.Run(testString("Evaluator/AddNew/op1=Ciphertext/op2=Ciphertext", tc.params), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectorsRingQ(tc, tc.encryptorPk, t)
		values2, _, ciphertext2 := newTestVectorsRingQ(tc, tc.encryptorPk, t)

		ciphertext1 = tc.evaluator.AddNew(ciphertext1, ciphertext2)
		tc.ringT.Add(values1, values2, values1)

		verifyTestVectors(tc, tc.decryptor, values1, ciphertext1, t)
	})

	t.Run(testString("Evaluator/AddNoModNew/op1=Ciphertext/op2=Ciphertext", tc.params), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectorsRingQ(tc, tc.encryptorPk, t)
		values2, _, ciphertext2 := newTestVectorsRingQ(tc, tc.encryptorPk, t)

		ciphertext1 = tc.evaluator.AddNoModNew(ciphertext1, ciphertext2)
		ciphertext1 = tc.evaluator.ReduceNew(ciphertext1)
		tc.ringT.Add(values1, values2, values1)

		verifyTestVectors(tc, tc.decryptor, values1, ciphertext1, t)
	})

	t.Run(testString("Evaluator/Add/op1=Ciphertext/op2=PlaintextRingT", tc.params), func(t *testing.T) {

		values1, plaintextRingT := newTestVectorsRingT(tc, t)
		values2, _, ciphertext := newTestVectorsRingQ(tc, tc.encryptorPk, t)
		ciphertextOut := NewCiphertext(tc.params, 1)

		tc.evaluator.Add(ciphertext, plaintextRingT, ciphertextOut)
		tc.ringT.Add(values1, values2, values2)

		verifyTestVectors(tc, tc.decryptor, values2, ciphertextOut, t)

		tc.evaluator.Add(plaintextRingT, ciphertext, ciphertextOut)

		verifyTestVectors(tc, tc.decryptor, values2, ciphertextOut, t)
	})

	t.Run(testString("Evaluator/Add/op1=Ciphertext/op2=Plaintext", tc.params), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectorsRingQ(tc, tc.encryptorPk, t)
		values2, plaintext2, ciphertext2 := newTestVectorsRingQ(tc, tc.encryptorPk, t)

		tc.evaluator.Add(ciphertext1, plaintext2, ciphertext2)
		tc.ringT.Add(values1, values2, values2)

		verifyTestVectors(tc, tc.decryptor, values2, ciphertext2, t)

		tc.evaluator.Add(plaintext2, ciphertext1, ciphertext2)

		verifyTestVectors(tc, tc.decryptor, values2, ciphertext2, t)
	})

	t.Run(testString("Evaluator/Sub/op1=Ciphertext/op2=Ciphertext", tc.params), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectorsRingQ(tc, tc.encryptorPk, t)
		values2, _, ciphertext2 := newTestVectorsRingQ(tc, tc.encryptorPk, t)

		tc.evaluator.Sub(ciphertext1, ciphertext2, ciphertext1)
		tc.ringT.Sub(values1, values2, values1)

		verifyTestVectors(tc, tc.decryptor, values1, ciphertext1, t)
	})

	t.Run(testString("Evaluator/SubNoMod/op1=Ciphertext/op2=Ciphertext", tc.params), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectorsRingQ(tc, tc.encryptorPk, t)
		values2, _, ciphertext2 := newTestVectorsRingQ(tc, tc.encryptorPk, t)

		tc.evaluator.SubNoMod(ciphertext1, ciphertext2, ciphertext1)
		tc.evaluator.Reduce(ciphertext1, ciphertext1)
		tc.ringT.Sub(values1, values2, values1)

		verifyTestVectors(tc, tc.decryptor, values1, ciphertext1, t)
	})

	t.Run(testString("Evaluator/SubNew/op1=Ciphertext/op2=Ciphertext", tc.params), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectorsRingQ(tc, tc.encryptorPk, t)
		values2, _, ciphertext2 := newTestVectorsRingQ(tc, tc.encryptorPk, t)

		ciphertext1 = tc.evaluator.SubNew(ciphertext1, ciphertext2)
		tc.ringT.Sub(values1, values2, values1)

		verifyTestVectors(tc, tc.decryptor, values1, ciphertext1, t)
	})

	t.Run(testString("Evaluator/Sub/op1=Ciphertext/op2=PlaintextRingT", tc.params), func(t *testing.T) {

		values1, plaintextRingT := newTestVectorsRingT(tc, t)
		values2, _, ciphertext := newTestVectorsRingQ(tc, tc.encryptorPk, t)
		ciphertextOut := NewCiphertext(tc.params, 1)
		plaintextWant := NewPlaintextRingT(tc.params)

		tc.evaluator.Sub(ciphertext, plaintextRingT, ciphertextOut)
		tc.ringT.Sub(values2, values1, plaintextWant.Value)
		verifyTestVectors(tc, tc.decryptor, plaintextWant.Value, ciphertextOut, t)

		tc.evaluator.Sub(plaintextRingT, ciphertext, ciphertextOut)
		tc.ringT.Sub(values1, values2, plaintextWant.Value)
		verifyTestVectors(tc, tc.decryptor, plaintextWant.Value, ciphertextOut, t)
	})

	t.Run(testString("Evaluator/Sub/op1=Ciphertext/op2=Plaintext", tc.params), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectorsRingQ(tc, tc.encryptorPk, t)
		values2, plaintext2, ciphertext2 := newTestVectorsRingQ(tc, tc.encryptorPk, t)

		valuesWant := tc.ringT.NewPoly()

		tc.evaluator.Sub(ciphertext1, plaintext2, ciphertext2)
		tc.ringT.Sub(values1, values2, valuesWant)
		verifyTestVectors(tc, tc.decryptor, valuesWant, ciphertext2, t)

		tc.evaluator.Sub(plaintext2, ciphertext1, ciphertext2)
		tc.ringT.Sub(values2, values1, valuesWant)
		verifyTestVectors(tc, tc.decryptor, valuesWant, ciphertext2, t)
	})

	t.Run(testString("Evaluator/SubNoModNew/op1=Ciphertext/op2=Ciphertext", tc.params), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectorsRingQ(tc, tc.encryptorPk, t)
		values2, _, ciphertext2 := newTestVectorsRingQ(tc, tc.encryptorPk, t)

		ciphertext1 = tc.evaluator.SubNoModNew(ciphertext1, ciphertext2)
		ciphertext1 = tc.evaluator.ReduceNew(ciphertext1)
		tc.ringT.Sub(values1, values2, values1)

		verifyTestVectors(tc, tc.decryptor, values1, ciphertext1, t)
	})

	t.Run(testString("Evaluator/SubNoMod/op1=Ciphertext/op2=Plaintext", tc.params), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectorsRingQ(tc, tc.encryptorPk, t)
		values2, plaintext2, ciphertext2 := newTestVectorsRingQ(tc, tc.encryptorPk, t)

		valuesWant := tc.ringT.NewPoly()

		tc.evaluator.SubNoMod(ciphertext1, plaintext2, ciphertext2)
		tc.evaluator.Reduce(ciphertext2, ciphertext2)
		tc.ringT.Sub(values1, values2, valuesWant)
		verifyTestVectors(tc, tc.decryptor, valuesWant, ciphertext2, t)

		tc.evaluator.SubNoMod(plaintext2, ciphertext1, ciphertext2)
		tc.evaluator.Reduce(ciphertext2, ciphertext2)
		tc.ringT.Sub(values2, values1, valuesWant)
		verifyTestVectors(tc, tc.decryptor, valuesWant, ciphertext2, t)
	})

	t.Run(testString("Evaluator/Neg", tc.params), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectorsRingQ(tc, tc.encryptorPk, t)

		tc.evaluator.Neg(ciphertext1, ciphertext1)
		tc.ringT.Neg(values1, values1)
		tc.ringT.Reduce(values1, values1)

		verifyTestVectors(tc, tc.decryptor, values1, ciphertext1, t)
	})

	t.Run(testString("Evaluator/NegNew", tc.params), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectorsRingQ(tc, tc.encryptorPk, t)

		ciphertext1 = tc.evaluator.NegNew(ciphertext1)
		tc.ringT.Neg(values1, values1)
		tc.ringT.Reduce(values1, values1)

		verifyTestVectors(tc, tc.decryptor, values1, ciphertext1, t)
	})

	t.Run(testString("Evaluator/AddScalar", tc.params), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectorsRingQ(tc, tc.encryptorPk, t)

		tc.evaluator.AddScalar(ciphertext1, 37, ciphertext1)
		tc.ringT.AddScalar(values1, 37, values1)

		verifyTestVectors(tc, tc.decryptor, values1, ciphertext1, t)
	})

	t.Run(testString("Evaluator/MulScalar", tc.params), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectorsRingQ(tc, tc.encryptorPk, t)

		tc.evaluator.MulScalar(ciphertext1, 37, ciphertext1)
		tc.ringT.MulScalar(values1, 37, values1)

		verifyTestVectors(tc, tc.decryptor, values1, ciphertext1, t)
	})

	t.Run(testString("Evaluator/MulScalarNew", tc.params), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectorsRingQ(tc, tc.encryptorPk, t)

		ciphertext1 = tc.evaluator.MulScalarNew(ciphertext1, 37)
		tc.ringT.MulScalar(values1, 37, values1)

		verifyTestVectors(tc, tc.decryptor, values1, ciphertext1, t)
	})

	t.Run(testString("Evaluator/Mul/op1=Ciphertext/op2=Ciphertext", tc.params), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectorsRingQ(tc, tc.encryptorPk, t)
		values2, _, ciphertext2 := newTestVectorsRingQ(tc, tc.encryptorPk, t)

		receiver := NewCiphertext(tc.params, ciphertext1.Degree()+ciphertext2.Degree())
		tc.evaluator.Mul(ciphertext1, ciphertext2, receiver)
		tc.ringT.MulCoeffs(values1, values2, values1)

		verifyTestVectors(tc, tc.decryptor, values1, receiver, t)
	})

	t.Run(testString("Evaluator/MulAndAdd/op1=Ciphertext/op2=Plaintext", tc.params), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectorsRingQ(tc, tc.encryptorPk, t)
		values2, plaintext2, ciphertext2 := newTestVectorsRingQ(tc, tc.encryptorPk, t)

		tc.evaluator.MulAndAdd(ciphertext1, plaintext2, ciphertext2)
		tmp := tc.ringT.NewPoly()
		tc.ringT.MulCoeffs(values1, values2, tmp)
		tc.ringT.Add(values2, tmp, values2)

		verifyTestVectors(tc, tc.decryptor, values2, ciphertext2, t)
	})

	t.Run(testString("Evaluator/MulAndAdd/op1=Ciphertext/op2=Ciphetext", tc.params), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectorsRingQ(tc, tc.encryptorPk, t)
		values2, _, ciphertext2 := newTestVectorsRingQ(tc, tc.encryptorPk, t)
		values3, _, ciphertext3 := newTestVectorsRingQ(tc, tc.encryptorPk, t)

		ciphertext3.Resize(tc.params.Parameters, 2)
		tc.evaluator.MulAndAdd(ciphertext1, ciphertext2, ciphertext3)
		tc.ringT.MulCoeffs(values1, values2, values1)
		tc.ringT.Add(values3, values1, values3)

		verifyTestVectors(tc, tc.decryptor, values3, ciphertext3, t)
	})

	t.Run(testString("Evaluator/MulNew/op1=Ciphertext/op2=Ciphertext", tc.params), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectorsRingQ(tc, tc.encryptorPk, t)
		values2, _, ciphertext2 := newTestVectorsRingQ(tc, tc.encryptorPk, t)

		receiver := tc.evaluator.MulNew(ciphertext1, ciphertext2)
		tc.ringT.MulCoeffs(values1, values2, values1)

		verifyTestVectors(tc, tc.decryptor, values1, receiver, t)
	})

	t.Run(testString("Evaluator/MulSquare/op1=Ciphertext/op2=Ciphertext", tc.params), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectorsRingQ(tc, tc.encryptorPk, t)

		receiver := NewCiphertext(tc.params, ciphertext1.Degree()+ciphertext1.Degree())
		tc.evaluator.Mul(ciphertext1, ciphertext1, receiver)
		tc.ringT.MulCoeffs(values1, values1, values1)

		verifyTestVectors(tc, tc.decryptor, values1, receiver, t)

		if tc.params.LogN() < 13 {
			t.Skip()
		}
		receiver2 := NewCiphertext(tc.params, receiver.Degree()+receiver.Degree())
		tc.evaluator.Mul(receiver, receiver, receiver2)
		tc.ringT.MulCoeffs(values1, values1, values1)

		verifyTestVectors(tc, tc.decryptor, values1, receiver2, t)
	})

	t.Run(testString("Evaluator/Mul/op1=Ciphertext/op2=Plaintext", tc.params), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectorsRingQ(tc, tc.encryptorPk, t)
		values2, plaintext2, _ := newTestVectorsRingQ(tc, nil, t)

		tc.evaluator.Mul(ciphertext1, plaintext2, ciphertext1)
		tc.ringT.MulCoeffs(values1, values2, values1)

		verifyTestVectors(tc, tc.decryptor, values1, ciphertext1, t)
	})

	t.Run(testString("Evaluator/Mul/op1=Ciphertext/op2=PlaintextRingT", tc.params), func(t *testing.T) {

		values1, plaintextRingT := newTestVectorsRingT(tc, t)
		values2, _, ciphertext := newTestVectorsRingQ(tc, tc.encryptorPk, t)
		ciphertextOut := NewCiphertext(tc.params, 1)

		tc.evaluator.Mul(ciphertext, plaintextRingT, ciphertextOut)
		tc.ringT.MulCoeffs(values1, values2, values1)

		verifyTestVectors(tc, tc.decryptor, values1, ciphertextOut, t)
	})

	t.Run(testString("Evaluator/Mul/op1=Ciphertext/op2=PlaintextMul", tc.params), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectorsRingQ(tc, tc.encryptorPk, t)
		values2, plaintext2 := newTestVectorsMul(tc, t)

		tc.evaluator.Mul(ciphertext1, plaintext2, ciphertext1)
		tc.ringT.MulCoeffs(values1, values2, values1)

		verifyTestVectors(tc, tc.decryptor, values1, ciphertext1, t)
	})

	t.Run(testString("Evaluator/Mul/Relinearize", tc.params), func(t *testing.T) {

		if tc.params.PCount() == 0 {
			t.Skip("#Pi is empty")
		}

		values1, _, ciphertext1 := newTestVectorsRingQ(tc, tc.encryptorPk, t)
		values2, _, ciphertext2 := newTestVectorsRingQ(tc, tc.encryptorPk, t)

		receiver := NewCiphertext(tc.params, ciphertext1.Degree()+ciphertext2.Degree())
		tc.evaluator.Mul(ciphertext1, ciphertext2, receiver)
		tc.ringT.MulCoeffs(values1, values2, values1)

		receiver2 := tc.evaluator.RelinearizeNew(receiver)
		verifyTestVectors(tc, tc.decryptor, values1, receiver2, t)

		tc.evaluator.Relinearize(receiver, receiver)
		verifyTestVectors(tc, tc.decryptor, values1, receiver, t)
	})

	t.Run(testString("Evaluator/QuantizeToLvl", tc.params), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectorsRingQ(tc, tc.encryptorPk, t)

		tc.evaluator.QuantizeToLvl(1, ciphertext1)
		verifyTestVectors(tc, tc.decryptor, values1, ciphertext1, t)

		if tc.params.T() != tc.params.RingQ().Modulus[0] { // only happens if T divides Q.
			tc.evaluator.QuantizeToLvl(0, ciphertext1)
			verifyTestVectors(tc, tc.decryptor, values1, ciphertext1, t)
		}
	})

	t.Run(testString("Evaluator/QuantizeToLvl/Add", tc.params), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectorsRingQ(tc, tc.encryptorPk, t)
		values2, _, ciphertext2 := newTestVectorsRingQ(tc, tc.encryptorPk, t)

		tc.evaluator.QuantizeToLvl(1, ciphertext1)
		tc.evaluator.QuantizeToLvl(1, ciphertext2)

		assert.True(t, ciphertext1.Level() == 1)
		assert.True(t, ciphertext2.Level() == 1)

		tc.evaluator.Add(ciphertext1, ciphertext2, ciphertext1)
		tc.ringT.Add(values1, values2, values1)

		verifyTestVectors(tc, tc.decryptor, values1, ciphertext1, t)
	})

	t.Run(testString("Evaluator/QuantizeToLvl/MulRelin", tc.params), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectorsRingQ(tc, tc.encryptorPk, t)
		values2, _, ciphertext2 := newTestVectorsRingQ(tc, tc.encryptorPk, t)

		tc.evaluator.QuantizeToLvl(1, ciphertext1)
		tc.evaluator.QuantizeToLvl(1, ciphertext2)

		assert.True(t, ciphertext1.Level() == 1)
		assert.True(t, ciphertext2.Level() == 1)

		receiver := NewCiphertext(tc.params, ciphertext1.Degree()+ciphertext2.Degree())
		tc.evaluator.Mul(ciphertext1, ciphertext2, receiver)
		tc.ringT.MulCoeffs(values1, values2, values1)

		verifyTestVectors(tc, tc.decryptor, values1, receiver, t)

		receiver2 := tc.evaluator.RelinearizeNew(receiver)
		verifyTestVectors(tc, tc.decryptor, values1, receiver2, t)
	})
}

func testPolyEval(tc *testContext, t *testing.T) {

	t.Run(testString("PowerBasis/Marshalling", tc.params), func(t *testing.T) {
		_, _, ct := newTestVectorsRingQ(tc, tc.encryptorPk, t)

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

	t.Run(testString("PolyEval/Single", tc.params), func(t *testing.T) {
		if tc.params.PCount() == 0 {
			t.Skip("#Pi is empty")
		}

		values, _, ciphertext := newTestVectorsRingQ(tc, tc.encryptorPk, t)

		coeffs := []uint64{1, 2, 3}

		T := tc.ringT.Modulus[0]
		for i := range values.Coeffs[0] {
			values.Coeffs[0][i] = ring.EvalPolyModP(values.Coeffs[0][i], coeffs, T)
		}

		poly := NewPoly(coeffs)

		var err error
		if ciphertext, err = tc.evaluator.EvaluatePoly(ciphertext, poly); err != nil {
			t.Fail()
		}

		verifyTestVectors(tc, tc.decryptor, values, ciphertext, t)

	})

	t.Run(testString("PolyEval/Vector", tc.params), func(t *testing.T) {

		if tc.params.PCount() == 0 {
			t.Skip("#Pi is empty")
		}

		values, _, ciphertext := newTestVectorsRingQ(tc, tc.encryptorPk, t)

		coeffs0 := []uint64{1, 2, 3}
		coeffs1 := []uint64{2, 3, 4}

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
			t.Fail()
		}

		verifyTestVectors(tc, tc.decryptor, values, ciphertext, t)
	})
}

func testEvaluatorKeySwitch(tc *testContext, t *testing.T) {

	if tc.params.PCount() == 0 {
		t.Skip("#Pi is empty")
	}

	sk2 := tc.kgen.GenSecretKey()
	decryptorSk2 := NewDecryptor(tc.params, sk2)
	switchKey := tc.kgen.GenSwitchingKey(tc.sk, sk2)

	t.Run(testString("Evaluator/KeySwitch/InPlace", tc.params), func(t *testing.T) {
		values, _, ciphertext := newTestVectorsRingQ(tc, tc.encryptorPk, t)
		tc.evaluator.SwitchKeys(ciphertext, switchKey, ciphertext)
		verifyTestVectors(tc, decryptorSk2, values, ciphertext, t)
	})

	t.Run(testString("Evaluator/KeySwitch/New", tc.params), func(t *testing.T) {
		values, _, ciphertext := newTestVectorsRingQ(tc, tc.encryptorPk, t)
		ciphertext = tc.evaluator.SwitchKeysNew(ciphertext, switchKey)
		verifyTestVectors(tc, decryptorSk2, values, ciphertext, t)
	})
}

func testEvaluatorRotate(tc *testContext, t *testing.T) {

	if tc.params.PCount() == 0 {
		t.Skip("#Pi is empty")
	}

	rots := []int{1, -1, 4, -4, 63, -63}
	rotkey := tc.kgen.GenRotationKeysForRotations(rots, true, tc.sk)
	evaluator := tc.evaluator.WithKey(rlwe.EvaluationKey{Rlk: tc.rlk, Rtks: rotkey})

	t.Run(testString("Evaluator/RotateRows", tc.params), func(t *testing.T) {
		values, _, ciphertext := newTestVectorsRingQ(tc, tc.encryptorPk, t)
		evaluator.RotateRows(ciphertext, ciphertext)
		values.Coeffs[0] = append(values.Coeffs[0][tc.params.N()>>1:], values.Coeffs[0][:tc.params.N()>>1]...)
		verifyTestVectors(tc, tc.decryptor, values, ciphertext, t)
	})

	t.Run(testString("Evaluator/RotateRowsNew", tc.params), func(t *testing.T) {
		values, _, ciphertext := newTestVectorsRingQ(tc, tc.encryptorPk, t)
		ciphertext = evaluator.RotateRowsNew(ciphertext)
		values.Coeffs[0] = append(values.Coeffs[0][tc.params.N()>>1:], values.Coeffs[0][:tc.params.N()>>1]...)
		verifyTestVectors(tc, tc.decryptor, values, ciphertext, t)
	})

	t.Run(testString("Evaluator/RotateColumns", tc.params), func(t *testing.T) {

		values, _, ciphertext := newTestVectorsRingQ(tc, tc.encryptorPk, t)

		receiver := NewCiphertext(tc.params, 1)
		for _, n := range rots {

			evaluator.RotateColumns(ciphertext, n, receiver)
			valuesWant := utils.RotateUint64Slots(values.Coeffs[0], n)

			verifyTestVectors(tc, tc.decryptor, &ring.Poly{Coeffs: [][]uint64{valuesWant}}, receiver, t)
		}
	})

	t.Run(testString("Evaluator/RotateColumnsNew", tc.params), func(t *testing.T) {

		values, _, ciphertext := newTestVectorsRingQ(tc, tc.encryptorPk, t)

		for _, n := range rots {

			receiver := evaluator.RotateColumnsNew(ciphertext, n)
			valuesWant := utils.RotateUint64Slots(values.Coeffs[0], n)

			verifyTestVectors(tc, tc.decryptor, &ring.Poly{Coeffs: [][]uint64{valuesWant}}, receiver, t)
		}
	})

	rotkey = tc.kgen.GenRotationKeysForInnerSum(tc.sk)
	evaluator = evaluator.WithKey(rlwe.EvaluationKey{Rlk: tc.rlk, Rtks: rotkey})

	t.Run(testString("Evaluator/Rotate/InnerSum", tc.params), func(t *testing.T) {
		values, _, ciphertext := newTestVectorsRingQ(tc, tc.encryptorPk, t)

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

	t.Run(testString("Evaluator/QuantizeToLvl/Rotate", tc.params), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectorsRingQ(tc, tc.encryptorPk, t)

		tc.evaluator.QuantizeToLvl(1, ciphertext1)

		assert.True(t, ciphertext1.Level() == 1)

		rotkey := tc.kgen.GenRotationKeysForRotations(nil, true, tc.sk)
		evaluator := tc.evaluator.WithKey(rlwe.EvaluationKey{Rlk: tc.rlk, Rtks: rotkey})

		tc.evaluator.QuantizeToLvl(1, ciphertext1)

		assert.True(t, ciphertext1.Level() == 1)

		evaluator.RotateRows(ciphertext1, ciphertext1)
		values1.Coeffs[0] = append(values1.Coeffs[0][tc.params.N()>>1:], values1.Coeffs[0][:tc.params.N()>>1]...)
		verifyTestVectors(tc, tc.decryptor, values1, ciphertext1, t)
	})
}

func testMarshaller(tc *testContext, t *testing.T) {

	t.Run(testString("Marshaller/Parameters/Binary", tc.params), func(t *testing.T) {
		bytes, err := tc.params.MarshalBinary()
		assert.Nil(t, err)
		var p Parameters
		err = p.UnmarshalBinary(bytes)
		assert.Nil(t, err)
		assert.Equal(t, tc.params, p)
		assert.Equal(t, tc.params.MarshalBinarySize(), len(bytes))
	})

	t.Run(testString("Marshaller/Parameters/JSON", tc.params), func(t *testing.T) {
		// checks that parameters can be marshalled without error
		data, err := json.Marshal(tc.params)
		assert.Nil(t, err)
		assert.NotNil(t, data)

		// checks that bfv.Parameters can be unmarshalled without error
		var paramsRec Parameters
		err = json.Unmarshal(data, &paramsRec)
		assert.Nil(t, err)
		assert.True(t, tc.params.Equals(paramsRec))

		// checks that bfv.Paramters can be unmarshalled with log-moduli definition without error
		dataWithLogModuli := []byte(fmt.Sprintf(`{"LogN":%d,"LogQ":[50,50],"LogP":[60], "T":65537}`, tc.params.LogN()))
		var paramsWithLogModuli Parameters
		err = json.Unmarshal(dataWithLogModuli, &paramsWithLogModuli)
		assert.Nil(t, err)
		assert.Equal(t, 2, paramsWithLogModuli.QCount())
		assert.Equal(t, 1, paramsWithLogModuli.PCount())
		assert.Equal(t, rlwe.DefaultSigma, paramsWithLogModuli.Sigma()) // ommiting sigma should result in Default being used

		// checks that bfv.Paramters can be unmarshalled with log-moduli definition with empty P without error
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

	t.Run(testString("Marshaller/Ciphertext", tc.params), func(t *testing.T) {

		ciphertextWant := NewCiphertextRandom(tc.prng, tc.params, 2)

		marshalledCiphertext, err := ciphertextWant.MarshalBinary()
		require.NoError(t, err)

		ciphertextTest := new(Ciphertext)
		err = ciphertextTest.UnmarshalBinary(marshalledCiphertext)
		require.NoError(t, err)

		for i := range ciphertextWant.Value {
			require.True(t, tc.ringQ.Equal(ciphertextWant.Value[i], ciphertextTest.Value[i]))
		}
	})
}
