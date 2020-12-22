package bfv

import (
	"flag"
	"fmt"
	"math/rand"
	"testing"
	"time"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"

	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/utils"
)

var flagLongTest = flag.Bool("long", false, "run the long test suite (all parameters). Overrides -short and requires -timeout=0.")

func testString(opname string, p *Parameters) string {
	return fmt.Sprintf("%sLogN=%d/logQ=%d/alpha=%d/beta=%d", opname, p.logN, p.LogQP(), p.Alpha(), p.Beta())
}

type testContext struct {
	params      *Parameters
	ringQ       *ring.Ring
	ringQP      *ring.Ring
	ringT       *ring.Ring
	prng        utils.PRNG
	uSampler    *ring.UniformSampler
	encoder     Encoder
	kgen        KeyGenerator
	sk          *SecretKey
	pk          *PublicKey
	rlk         *EvaluationKey
	encryptorPk Encryptor
	encryptorSk Encryptor
	decryptor   Decryptor
	evaluator   Evaluator
}

func TestBFV(t *testing.T) {

	var err error
	var testctx = new(testContext)

	rand.Seed(time.Now().UnixNano())

	var defaultParams = DefaultParams[PN12QP109 : PN12QP109+4] // the default test runs for ring degree N=2^12, 2^13, 2^14, 2^15
	if testing.Short() {
		defaultParams = DefaultParams[PN12QP109 : PN12QP109+2] // the short test suite runs for ring degree N=2^12, 2^13
	}
	if *flagLongTest {
		defaultParams = DefaultParams // the long test suite runs for all default parameters
		fmt.Println("bfv runing in long mode")
	}

	for _, p := range defaultParams {

		if testctx, err = genTestParams(p); err != nil {
			panic(err)
		}

		testParameters(testctx, t)
		testEncoder(testctx, t)
		testEncryptor(testctx, t)
		testEvaluator(testctx, t)
		testEvaluatorKeySwitch(testctx, t)
		testEvaluatorRotate(testctx, t)
		testMarshaller(testctx, t)
	}

}

func genTestParams(params *Parameters) (testctx *testContext, err error) {

	testctx = new(testContext)
	testctx.params = params.Copy()

	if testctx.prng, err = utils.NewPRNG(); err != nil {
		return nil, err
	}

	if testctx.ringQ, err = ring.NewRing(params.N(), params.qi); err != nil {
		return nil, err
	}

	if testctx.ringQP, err = ring.NewRing(params.N(), append(params.qi, params.pi...)); err != nil {
		return nil, err
	}

	if testctx.ringT, err = ring.NewRing(params.N(), []uint64{params.t}); err != nil {
		return nil, err
	}

	testctx.uSampler = ring.NewUniformSampler(testctx.prng, testctx.ringT)
	testctx.kgen = NewKeyGenerator(testctx.params)
	testctx.sk, testctx.pk = testctx.kgen.GenKeyPair()
	if params.PiCount() != 0 {
		testctx.rlk = testctx.kgen.GenRelinKey(testctx.sk, 1)
	}

	testctx.encoder = NewEncoder(testctx.params)
	testctx.encryptorPk = NewEncryptorFromPk(testctx.params, testctx.pk)
	testctx.encryptorSk = NewEncryptorFromSk(testctx.params, testctx.sk)
	testctx.decryptor = NewDecryptor(testctx.params, testctx.sk)
	testctx.evaluator = NewEvaluator(testctx.params)
	return

}

func testParameters(testctx *testContext, t *testing.T) {
	t.Run("Parameters/NewParametersFromModuli/", func(t *testing.T) {
		p, err := NewParametersFromModuli(testctx.params.logN, testctx.params.Moduli(), testctx.params.t)
		assert.NoError(t, err)
		assert.True(t, p.Equals(testctx.params))
	})

	t.Run("Parameters/NewParametersFromLogModuli/", func(t *testing.T) {
		p, err := NewParametersFromLogModuli(testctx.params.logN, testctx.params.LogModuli(), testctx.params.t)
		assert.NoError(t, err)
		assert.True(t, p.Equals(testctx.params))
	})
}

func newTestVectorsRingQ(testctx *testContext, encryptor Encryptor, t *testing.T) (coeffs *ring.Poly, plaintext *Plaintext, ciphertext *Ciphertext) {

	coeffs = testctx.uSampler.ReadNew()

	plaintext = NewPlaintext(testctx.params)

	testctx.encoder.EncodeUint(coeffs.Coeffs[0], plaintext)

	if encryptor != nil {
		if testctx.params.PiCount() != 0 {
			ciphertext = testctx.encryptorPk.EncryptNew(plaintext)
		} else {
			ciphertext = testctx.encryptorPk.EncryptFastNew(plaintext)
		}

	}

	return coeffs, plaintext, ciphertext
}

func newTestVectorsRingT(testctx *testContext, t *testing.T) (coeffs *ring.Poly, plaintext *PlaintextRingT) {

	coeffs = testctx.uSampler.ReadNew()

	plaintext = NewPlaintextRingT(testctx.params)

	testctx.encoder.EncodeUintRingT(coeffs.Coeffs[0], plaintext)

	return coeffs, plaintext
}

func newTestVectorsMul(testctx *testContext, t *testing.T) (coeffs *ring.Poly, plaintext *PlaintextMul) {

	coeffs = testctx.uSampler.ReadNew()

	plaintext = NewPlaintextMul(testctx.params)

	testctx.encoder.EncodeUintMul(coeffs.Coeffs[0], plaintext)

	return coeffs, plaintext
}

func verifyTestVectors(testctx *testContext, decryptor Decryptor, coeffs *ring.Poly, element Operand, t *testing.T) {

	var coeffsTest []uint64

	switch el := element.(type) {
	case *Plaintext, *PlaintextMul, *PlaintextRingT:
		coeffsTest = testctx.encoder.DecodeUintNew(el)
	case *Ciphertext:
		coeffsTest = testctx.encoder.DecodeUintNew(decryptor.DecryptNew(el))
	default:
		t.Error("invalid test object to verify")
	}

	require.True(t, utils.EqualSliceUint64(coeffs.Coeffs[0], coeffsTest))
}

func testEncoder(testctx *testContext, t *testing.T) {
	t.Run(testString("Encoder/Encode&Decode/RingT/", testctx.params), func(t *testing.T) {
		values, plaintext := newTestVectorsRingT(testctx, t)
		verifyTestVectors(testctx, nil, values, plaintext, t)
	})

	t.Run(testString("Encoder/Encode&Decode/RingQ/", testctx.params), func(t *testing.T) {
		values, plaintext, _ := newTestVectorsRingQ(testctx, nil, t)
		verifyTestVectors(testctx, nil, values, plaintext, t)
	})

	t.Run(testString("Encoder/Encode&Decode/PlaintextMul/", testctx.params), func(t *testing.T) {
		values, plaintext := newTestVectorsMul(testctx, t)
		verifyTestVectors(testctx, nil, values, plaintext, t)
	})
}

func testEncryptor(testctx *testContext, t *testing.T) {

	coeffs := testctx.uSampler.ReadNew()

	plaintextRingT := NewPlaintextRingT(testctx.params)
	plaintext := NewPlaintext(testctx.params)

	testctx.encoder.EncodeUintRingT(coeffs.Coeffs[0], plaintextRingT)
	testctx.encoder.EncodeUint(coeffs.Coeffs[0], plaintext)

	t.Run(testString("Encryptor/EncryptFromPk/", testctx.params), func(t *testing.T) {

		if testctx.params.PiCount() == 0 {
			t.Skip("#Pi is empty")
		}

		verifyTestVectors(testctx, testctx.decryptor, coeffs, testctx.encryptorPk.EncryptNew(plaintext), t)
	})

	t.Run(testString("Encryptor/EncryptFromPkFast/", testctx.params), func(t *testing.T) {
		verifyTestVectors(testctx, testctx.decryptor, coeffs, testctx.encryptorPk.EncryptFastNew(plaintext), t)
	})

	t.Run(testString("Encryptor/EncryptFromSk/", testctx.params), func(t *testing.T) {
		verifyTestVectors(testctx, testctx.decryptor, coeffs, testctx.encryptorSk.EncryptNew(plaintext), t)
	})
}

func testEvaluator(testctx *testContext, t *testing.T) {

	t.Run(testString("Evaluator/Add/op1=Ciphertext/op2=Ciphertext/", testctx.params), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectorsRingQ(testctx, testctx.encryptorPk, t)
		values2, _, ciphertext2 := newTestVectorsRingQ(testctx, testctx.encryptorPk, t)

		testctx.evaluator.Add(ciphertext1, ciphertext2, ciphertext1)
		testctx.ringT.Add(values1, values2, values1)

		verifyTestVectors(testctx, testctx.decryptor, values1, ciphertext1, t)
	})

	t.Run(testString("Evaluator/AddNew/op1=Ciphertext/op2=Ciphertext/", testctx.params), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectorsRingQ(testctx, testctx.encryptorPk, t)
		values2, _, ciphertext2 := newTestVectorsRingQ(testctx, testctx.encryptorPk, t)

		ciphertext1 = testctx.evaluator.AddNew(ciphertext1, ciphertext2)
		testctx.ringT.Add(values1, values2, values1)

		verifyTestVectors(testctx, testctx.decryptor, values1, ciphertext1, t)
	})

	t.Run(testString("Evaluator/Add/op1=Ciphertext/op2=PlaintextRingT/", testctx.params), func(t *testing.T) {

		values1, plaintextRingT := newTestVectorsRingT(testctx, t)
		values2, _, ciphertext := newTestVectorsRingQ(testctx, testctx.encryptorPk, t)
		ciphertextOut := NewCiphertext(testctx.params, 1)

		testctx.evaluator.Add(ciphertext, plaintextRingT, ciphertextOut)
		testctx.ringT.Add(values1, values2, values2)

		verifyTestVectors(testctx, testctx.decryptor, values2, ciphertextOut, t)

		testctx.evaluator.Add(plaintextRingT, ciphertext, ciphertextOut)

		verifyTestVectors(testctx, testctx.decryptor, values2, ciphertextOut, t)
	})

	t.Run(testString("Evaluator/Add/op1=Ciphertext/op2=Plaintext/", testctx.params), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectorsRingQ(testctx, testctx.encryptorPk, t)
		values2, plaintext2, ciphertext2 := newTestVectorsRingQ(testctx, testctx.encryptorPk, t)

		testctx.evaluator.Add(ciphertext1, plaintext2, ciphertext2)
		testctx.ringT.Add(values1, values2, values2)

		verifyTestVectors(testctx, testctx.decryptor, values2, ciphertext2, t)

		testctx.evaluator.Add(plaintext2, ciphertext1, ciphertext2)

		verifyTestVectors(testctx, testctx.decryptor, values2, ciphertext2, t)
	})

	t.Run(testString("Evaluator/Sub/op1=Ciphertext/op2=Ciphertext/", testctx.params), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectorsRingQ(testctx, testctx.encryptorPk, t)
		values2, _, ciphertext2 := newTestVectorsRingQ(testctx, testctx.encryptorPk, t)

		testctx.evaluator.Sub(ciphertext1, ciphertext2, ciphertext1)
		testctx.ringT.Sub(values1, values2, values1)

		verifyTestVectors(testctx, testctx.decryptor, values1, ciphertext1, t)
	})

	t.Run(testString("Evaluator/SubNew/op1=Ciphertext/op2=Ciphertext/", testctx.params), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectorsRingQ(testctx, testctx.encryptorPk, t)
		values2, _, ciphertext2 := newTestVectorsRingQ(testctx, testctx.encryptorPk, t)

		ciphertext1 = testctx.evaluator.SubNew(ciphertext1, ciphertext2)
		testctx.ringT.Sub(values1, values2, values1)

		verifyTestVectors(testctx, testctx.decryptor, values1, ciphertext1, t)
	})

	t.Run(testString("Evaluator/Sub/op1=Ciphertext/op2=PlaintextRingT/", testctx.params), func(t *testing.T) {

		values1, plaintextRingT := newTestVectorsRingT(testctx, t)
		values2, _, ciphertext := newTestVectorsRingQ(testctx, testctx.encryptorPk, t)
		ciphertextOut := NewCiphertext(testctx.params, 1)
		plaintextWant := NewPlaintextRingT(testctx.params)

		testctx.evaluator.Sub(ciphertext, plaintextRingT, ciphertextOut)
		testctx.ringT.Sub(values2, values1, plaintextWant.value)
		verifyTestVectors(testctx, testctx.decryptor, plaintextWant.value, ciphertextOut, t)

		testctx.evaluator.Sub(plaintextRingT, ciphertext, ciphertextOut)
		testctx.ringT.Sub(values1, values2, plaintextWant.value)
		verifyTestVectors(testctx, testctx.decryptor, plaintextWant.value, ciphertextOut, t)
	})

	t.Run(testString("Evaluator/Sub/op1=Ciphertext/op2=Plaintext/", testctx.params), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectorsRingQ(testctx, testctx.encryptorPk, t)
		values2, plaintext2, ciphertext2 := newTestVectorsRingQ(testctx, testctx.encryptorPk, t)

		valuesWant := testctx.ringT.NewPoly()

		testctx.evaluator.Sub(ciphertext1, plaintext2, ciphertext2)
		testctx.ringT.Sub(values1, values2, valuesWant)
		verifyTestVectors(testctx, testctx.decryptor, valuesWant, ciphertext2, t)

		testctx.evaluator.Sub(plaintext2, ciphertext1, ciphertext2)
		testctx.ringT.Sub(values2, values1, valuesWant)
		verifyTestVectors(testctx, testctx.decryptor, valuesWant, ciphertext2, t)
	})

	t.Run(testString("Evaluator/Mul/op1=Ciphertext/op2=Ciphertext/", testctx.params), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectorsRingQ(testctx, testctx.encryptorPk, t)
		values2, _, ciphertext2 := newTestVectorsRingQ(testctx, testctx.encryptorPk, t)

		receiver := NewCiphertext(testctx.params, ciphertext1.Degree()+ciphertext2.Degree())
		testctx.evaluator.Mul(ciphertext1, ciphertext2, receiver)
		testctx.ringT.MulCoeffs(values1, values2, values1)

		verifyTestVectors(testctx, testctx.decryptor, values1, receiver, t)
	})

	t.Run(testString("Evaluator/Mul/op1=Ciphertext/op2=Plaintext/", testctx.params), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectorsRingQ(testctx, testctx.encryptorPk, t)
		values2, plaintext2, _ := newTestVectorsRingQ(testctx, nil, t)

		testctx.evaluator.Mul(ciphertext1, plaintext2, ciphertext1)
		testctx.ringT.MulCoeffs(values1, values2, values1)

		verifyTestVectors(testctx, testctx.decryptor, values1, ciphertext1, t)
	})

	t.Run(testString("Evaluator/Mul/op1=Ciphertext/op2=PlaintextRingT/", testctx.params), func(t *testing.T) {

		values1, plaintextRingT := newTestVectorsRingT(testctx, t)
		values2, _, ciphertext := newTestVectorsRingQ(testctx, testctx.encryptorPk, t)
		ciphertextOut := NewCiphertext(testctx.params, 1)

		testctx.evaluator.Mul(ciphertext, plaintextRingT, ciphertextOut)
		testctx.ringT.MulCoeffs(values1, values2, values1)

		verifyTestVectors(testctx, testctx.decryptor, values1, ciphertextOut, t)
	})

	t.Run(testString("Evaluator/Mul/op1=Ciphertext/op2=PlaintextMul/", testctx.params), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectorsRingQ(testctx, testctx.encryptorPk, t)
		values2, plaintext2 := newTestVectorsMul(testctx, t)

		testctx.evaluator.Mul(ciphertext1, plaintext2, ciphertext1)
		testctx.ringT.MulCoeffs(values1, values2, values1)

		verifyTestVectors(testctx, testctx.decryptor, values1, ciphertext1, t)
	})

	t.Run(testString("Evaluator/Mul/Relinearize/", testctx.params), func(t *testing.T) {

		if testctx.params.PiCount() == 0 {
			t.Skip("#Pi is empty")
		}

		values1, _, ciphertext1 := newTestVectorsRingQ(testctx, testctx.encryptorPk, t)
		values2, _, ciphertext2 := newTestVectorsRingQ(testctx, testctx.encryptorPk, t)

		receiver := NewCiphertext(testctx.params, ciphertext1.Degree()+ciphertext2.Degree())
		testctx.evaluator.Mul(ciphertext1, ciphertext2, receiver)
		testctx.ringT.MulCoeffs(values1, values2, values1)

		receiver2 := testctx.evaluator.RelinearizeNew(receiver, testctx.rlk)
		verifyTestVectors(testctx, testctx.decryptor, values1, receiver2, t)

		testctx.evaluator.Relinearize(receiver, testctx.rlk, receiver)
		verifyTestVectors(testctx, testctx.decryptor, values1, receiver, t)
	})
}

func testEvaluatorKeySwitch(testctx *testContext, t *testing.T) {

	if testctx.params.PiCount() == 0 {
		t.Skip("#Pi is empty")
	}

	sk2 := testctx.kgen.GenSecretKey()
	decryptorSk2 := NewDecryptor(testctx.params, sk2)
	switchKey := testctx.kgen.GenSwitchingKey(testctx.sk, sk2)

	t.Run(testString("Evaluator/KeySwitch/InPlace/", testctx.params), func(t *testing.T) {
		values, _, ciphertext := newTestVectorsRingQ(testctx, testctx.encryptorPk, t)
		testctx.evaluator.SwitchKeys(ciphertext, switchKey, ciphertext)
		verifyTestVectors(testctx, decryptorSk2, values, ciphertext, t)
	})

	t.Run(testString("Evaluator/KeySwitch/New/", testctx.params), func(t *testing.T) {
		values, _, ciphertext := newTestVectorsRingQ(testctx, testctx.encryptorPk, t)
		ciphertext = testctx.evaluator.SwitchKeysNew(ciphertext, switchKey)
		verifyTestVectors(testctx, decryptorSk2, values, ciphertext, t)
	})
}

func testEvaluatorRotate(testctx *testContext, t *testing.T) {

	if testctx.params.PiCount() == 0 {
		t.Skip("#Pi is empty")
	}

	rotkey := testctx.kgen.GenRotationKeysPow2(testctx.sk)

	t.Run(testString("Evaluator/Rotate/Rows/InPlace/", testctx.params), func(t *testing.T) {
		values, _, ciphertext := newTestVectorsRingQ(testctx, testctx.encryptorPk, t)
		testctx.evaluator.RotateRows(ciphertext, rotkey, ciphertext)
		values.Coeffs[0] = append(values.Coeffs[0][testctx.params.N()>>1:], values.Coeffs[0][:testctx.params.N()>>1]...)
		verifyTestVectors(testctx, testctx.decryptor, values, ciphertext, t)
	})

	t.Run(testString("Evaluator/Rotate/Rows/New/", testctx.params), func(t *testing.T) {
		values, _, ciphertext := newTestVectorsRingQ(testctx, testctx.encryptorPk, t)
		ciphertext = testctx.evaluator.RotateRowsNew(ciphertext, rotkey)
		values.Coeffs[0] = append(values.Coeffs[0][testctx.params.N()>>1:], values.Coeffs[0][:testctx.params.N()>>1]...)
		verifyTestVectors(testctx, testctx.decryptor, values, ciphertext, t)
	})

	valuesWant := testctx.ringT.NewPoly()
	mask := (testctx.params.N() >> 1) - 1
	slots := testctx.params.N() >> 1

	t.Run(testString("Evaluator/Rotate/Cols/InPlace/", testctx.params), func(t *testing.T) {

		values, _, ciphertext := newTestVectorsRingQ(testctx, testctx.encryptorPk, t)

		receiver := NewCiphertext(testctx.params, 1)
		for n := uint64(1); n < slots; n <<= 1 {

			testctx.evaluator.RotateColumns(ciphertext, n, rotkey, receiver)

			for i := uint64(0); i < slots; i++ {
				valuesWant.Coeffs[0][i] = values.Coeffs[0][(i+n)&mask]
				valuesWant.Coeffs[0][i+slots] = values.Coeffs[0][((i+n)&mask)+slots]
			}

			verifyTestVectors(testctx, testctx.decryptor, valuesWant, receiver, t)
		}
	})

	t.Run(testString("Evaluator/Rotate/Cols/New/", testctx.params), func(t *testing.T) {

		values, _, ciphertext := newTestVectorsRingQ(testctx, testctx.encryptorPk, t)

		for n := uint64(1); n < slots; n <<= 1 {

			receiver := testctx.evaluator.RotateColumnsNew(ciphertext, n, rotkey)

			for i := uint64(0); i < slots; i++ {
				valuesWant.Coeffs[0][i] = values.Coeffs[0][(i+n)&mask]
				valuesWant.Coeffs[0][i+slots] = values.Coeffs[0][((i+n)&mask)+slots]
			}

			verifyTestVectors(testctx, testctx.decryptor, valuesWant, receiver, t)
		}
	})

	t.Run(testString("Evaluator/Rotate/Cols/Random/", testctx.params), func(t *testing.T) {

		values, _, ciphertext := newTestVectorsRingQ(testctx, testctx.encryptorPk, t)

		receiver := NewCiphertext(testctx.params, 1)
		prng, err := utils.NewPRNG()
		if err != nil {
			panic(err)
		}

		for n := 0; n < 4; n++ {

			rand := ring.RandUniform(prng, slots, mask)

			testctx.evaluator.RotateColumns(ciphertext, rand, rotkey, receiver)

			for i := uint64(0); i < slots; i++ {
				valuesWant.Coeffs[0][i] = values.Coeffs[0][(i+rand)&mask]
				valuesWant.Coeffs[0][i+slots] = values.Coeffs[0][((i+rand)&mask)+slots]
			}

			verifyTestVectors(testctx, testctx.decryptor, valuesWant, receiver, t)
		}
	})
}

func testMarshaller(testctx *testContext, t *testing.T) {

	t.Run("Marshaller/Parameters/ZeroValue", func(t *testing.T) {
		bytes, err := (&Parameters{}).MarshalBinary()
		assert.Nil(t, err)
		assert.Equal(t, []byte{}, bytes)
		p := new(Parameters)
		err = p.UnmarshalBinary(bytes)
		assert.NotNil(t, err)
	})

	t.Run("Marshaller/Parameters/SupportedParams", func(t *testing.T) {
		bytes, err := testctx.params.MarshalBinary()
		assert.Nil(t, err)
		p := new(Parameters)
		err = p.UnmarshalBinary(bytes)
		assert.Nil(t, err)
		assert.Equal(t, testctx.params, p)
	})

	ringQP := testctx.ringQP

	t.Run(testString("Marshaller/Ciphertext/", testctx.params), func(t *testing.T) {

		ciphertextWant := NewCiphertextRandom(testctx.prng, testctx.params, 2)

		marshalledCiphertext, err := ciphertextWant.MarshalBinary()
		require.NoError(t, err)

		ciphertextTest := new(Ciphertext)
		err = ciphertextTest.UnmarshalBinary(marshalledCiphertext)
		require.NoError(t, err)

		for i := range ciphertextWant.value {
			require.True(t, testctx.ringQ.Equal(ciphertextWant.value[i], ciphertextTest.value[i]))
		}
	})

	t.Run(testString("Marshaller/Sk/", testctx.params), func(t *testing.T) {

		marshalledSk, err := testctx.sk.MarshalBinary()
		require.NoError(t, err)

		sk := new(SecretKey)
		err = sk.UnmarshalBinary(marshalledSk)
		require.NoError(t, err)

		require.True(t, ringQP.Equal(sk.sk, testctx.sk.sk))
	})

	t.Run(testString("Marshaller/Pk/", testctx.params), func(t *testing.T) {

		marshalledPk, err := testctx.pk.MarshalBinary()
		require.NoError(t, err)

		pk := new(PublicKey)
		err = pk.UnmarshalBinary(marshalledPk)
		require.NoError(t, err)

		for k := range testctx.pk.pk {
			require.True(t, ringQP.Equal(pk.pk[k], testctx.pk.pk[k]), k)
		}
	})

	t.Run(testString("Marshaller/EvaluationKey/", testctx.params), func(t *testing.T) {

		if testctx.params.PiCount() == 0 {
			t.Skip("#Pi is empty")
		}

		evalkey := testctx.kgen.GenRelinKey(testctx.sk, 2)
		data, err := evalkey.MarshalBinary()
		require.NoError(t, err)

		resEvalKey := new(EvaluationKey)
		err = resEvalKey.UnmarshalBinary(data)
		require.NoError(t, err)

		for deg := range evalkey.evakey {

			evakeyWant := evalkey.evakey[deg].evakey
			evakeyTest := resEvalKey.evakey[deg].evakey

			for j := range evakeyWant {

				for k := range evakeyWant[j] {
					require.Truef(t, ringQP.Equal(evakeyWant[j][k], evakeyTest[j][k]), "deg %d element [%d][%d]", deg, j, k)
				}
			}
		}
	})

	t.Run(testString("Marshaller/SwitchingKey/", testctx.params), func(t *testing.T) {

		if testctx.params.PiCount() == 0 {
			t.Skip("#Pi is empty")
		}

		skOut := testctx.kgen.GenSecretKey()

		switchingKey := testctx.kgen.GenSwitchingKey(testctx.sk, skOut)
		data, err := switchingKey.MarshalBinary()
		require.NoError(t, err)

		resSwitchingKey := new(SwitchingKey)
		err = resSwitchingKey.UnmarshalBinary(data)
		require.NoError(t, err)

		evakeyWant := switchingKey.evakey
		evakeyTest := resSwitchingKey.evakey

		for j := range evakeyWant {

			for k := range evakeyWant[j] {
				require.Truef(t, ringQP.Equal(evakeyWant[j][k], evakeyTest[j][k]), "marshal SwitchingKey element [%d][%d]", j, k)
			}
		}
	})

	t.Run(testString("Marshaller/RotationKey/", testctx.params), func(t *testing.T) {

		if testctx.params.PiCount() == 0 {
			t.Skip("#Pi is empty")
		}

		rotationKey := NewRotationKeys()

		testctx.kgen.GenRot(RotationRow, testctx.sk, 0, rotationKey)
		testctx.kgen.GenRot(RotationLeft, testctx.sk, 1, rotationKey)
		testctx.kgen.GenRot(RotationLeft, testctx.sk, 2, rotationKey)
		testctx.kgen.GenRot(RotationRight, testctx.sk, 3, rotationKey)
		testctx.kgen.GenRot(RotationRight, testctx.sk, 5, rotationKey)

		data, err := rotationKey.MarshalBinary()
		require.NoError(t, err)

		resRotationKey := new(RotationKeys)
		err = resRotationKey.UnmarshalBinary(data)
		require.NoError(t, err)

		for i := uint64(1); i < testctx.params.N()>>1; i++ {

			if rotationKey.evakeyRotColLeft[i] != nil {

				evakeyWant := rotationKey.evakeyRotColLeft[i].evakey
				evakeyTest := resRotationKey.evakeyRotColLeft[i].evakey

				for j := range evakeyWant {

					for k := range evakeyWant[j] {
						require.Truef(t, ringQP.Equal(evakeyWant[j][k], evakeyTest[j][k]), "marshal RotationKey RotateLeft %d element [%d][%d]", i, j, k)
					}
				}
			}

			if rotationKey.evakeyRotColRight[i] != nil {

				evakeyWant := rotationKey.evakeyRotColRight[i].evakey
				evakeyTest := resRotationKey.evakeyRotColRight[i].evakey

				for j := range evakeyWant {

					for k := range evakeyWant[j] {
						require.Truef(t, ringQP.Equal(evakeyWant[j][k], evakeyTest[j][k]), "marshal RotationKey RotateRight %d element [%d][%d]", i, j, k)
					}
				}
			}
		}

		if rotationKey.evakeyRotRow != nil {

			evakeyWant := rotationKey.evakeyRotRow.evakey
			evakeyTest := resRotationKey.evakeyRotRow.evakey

			for j := range evakeyWant {

				for k := range evakeyWant[j] {
					require.Truef(t, ringQP.Equal(evakeyWant[j][k], evakeyTest[j][k]), "marshal RotationKey RotateRow element [%d][%d]", j, k)
				}
			}
		}
	})
}
