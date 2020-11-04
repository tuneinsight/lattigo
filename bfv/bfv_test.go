package bfv

import (
	"fmt"
	"math/rand"
	"testing"
	"time"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"

	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/utils"
)

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

	var defaultParams []*Parameters

	if testing.Short() {
		defaultParams = DefaultParams[PN12QP109 : PN12QP109+3]
	} else {
		defaultParams = DefaultParams
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
	testctx.rlk = testctx.kgen.GenRelinKey(testctx.sk, 1)
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

func newTestVectorsZQ(testctx *testContext, encryptor Encryptor, t *testing.T) (coeffs *ring.Poly, plaintext *Plaintext, ciphertext *Ciphertext) {

	coeffs = testctx.uSampler.ReadNew()

	plaintext = NewPlaintextZQ(testctx.params)

	testctx.encoder.EncodeUintZQ(coeffs.Coeffs[0], plaintext)

	if encryptor != nil {
		ciphertext = testctx.encryptorPk.EncryptNew(plaintext)
	}

	return coeffs, plaintext, ciphertext
}

func newTestVectorsZT(testctx *testContext, encryptor Encryptor, t *testing.T) (coeffs *ring.Poly, plaintext *Plaintext, ciphertext *Ciphertext) {

	coeffs = testctx.uSampler.ReadNew()

	plaintext = NewPlaintextZT(testctx.params)

	testctx.encoder.EncodeUintZT(coeffs.Coeffs[0], plaintext)

	if encryptor != nil {
		ciphertext = testctx.encryptorPk.EncryptNew(plaintext)
	}

	return coeffs, plaintext, ciphertext
}

func verifyTestVectors(testctx *testContext, decryptor Decryptor, coeffs *ring.Poly, element Operand, t *testing.T) {

	var coeffsTest []uint64

	el := element.El()

	if el.Degree() == 0 {

		coeffsTest = testctx.encoder.DecodeUint(el.Plaintext())

	} else {

		coeffsTest = testctx.encoder.DecodeUint(decryptor.DecryptNew(el.Ciphertext()))
	}

	require.True(t, utils.EqualSliceUint64(coeffs.Coeffs[0], coeffsTest))
}

func testEncoder(testctx *testContext, t *testing.T) {
	t.Run(testString("Encoder/Encode&Decode/ZT/", testctx.params), func(t *testing.T) {
		values, plaintext, _ := newTestVectorsZT(testctx, nil, t)
		verifyTestVectors(testctx, testctx.decryptor, values, plaintext, t)
	})

	t.Run(testString("Encoder/Encode&Decode/ZQ/", testctx.params), func(t *testing.T) {
		values, plaintext, _ := newTestVectorsZQ(testctx, nil, t)
		verifyTestVectors(testctx, testctx.decryptor, values, plaintext, t)
	})
}

func testEncryptor(testctx *testContext, t *testing.T) {

	coeffs := testctx.uSampler.ReadNew()

	plaintext := NewPlaintextZT(testctx.params)
	testctx.encoder.EncodeUintZT(coeffs.Coeffs[0], plaintext)

	t.Run(testString("Encryptor/EncryptFromPk/ZT/", testctx.params), func(t *testing.T) {
		verifyTestVectors(testctx, testctx.decryptor, coeffs, testctx.encryptorPk.EncryptNew(plaintext), t)
	})

	t.Run(testString("Encryptor/EncryptFromPkFast/ZT/", testctx.params), func(t *testing.T) {
		verifyTestVectors(testctx, testctx.decryptor, coeffs, testctx.encryptorPk.EncryptFastNew(plaintext), t)
	})

	t.Run(testString("Encryptor/EncryptFromSk/ZT/", testctx.params), func(t *testing.T) {
		verifyTestVectors(testctx, testctx.decryptor, coeffs, testctx.encryptorSk.EncryptNew(plaintext), t)
	})

	plaintext = NewPlaintextZQ(testctx.params)
	testctx.encoder.EncodeUintZQ(coeffs.Coeffs[0], plaintext)

	t.Run(testString("Encryptor/EncryptFromPk/ZQ/", testctx.params), func(t *testing.T) {
		verifyTestVectors(testctx, testctx.decryptor, coeffs, testctx.encryptorPk.EncryptNew(plaintext), t)
	})

	t.Run(testString("Encryptor/EncryptFromPkFast/ZQ/", testctx.params), func(t *testing.T) {
		verifyTestVectors(testctx, testctx.decryptor, coeffs, testctx.encryptorPk.EncryptFastNew(plaintext), t)
	})

	t.Run(testString("Encryptor/EncryptFromSk/ZQ/", testctx.params), func(t *testing.T) {
		verifyTestVectors(testctx, testctx.decryptor, coeffs, testctx.encryptorSk.EncryptNew(plaintext), t)
	})

}

func testEvaluator(testctx *testContext, t *testing.T) {

	t.Run(testString("Evaluator/Add/CtCtInPlace/", testctx.params), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectorsZQ(testctx, testctx.encryptorPk, t)
		values2, _, ciphertext2 := newTestVectorsZQ(testctx, testctx.encryptorPk, t)

		testctx.evaluator.Add(ciphertext1, ciphertext2, ciphertext1)
		testctx.ringT.Add(values1, values2, values1)

		verifyTestVectors(testctx, testctx.decryptor, values1, ciphertext1, t)
	})

	t.Run(testString("Evaluator/Add/CtCtNew/", testctx.params), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectorsZQ(testctx, testctx.encryptorPk, t)
		values2, _, ciphertext2 := newTestVectorsZQ(testctx, testctx.encryptorPk, t)

		ciphertext1 = testctx.evaluator.AddNew(ciphertext1, ciphertext2)
		testctx.ringT.Add(values1, values2, values1)

		verifyTestVectors(testctx, testctx.decryptor, values1, ciphertext1, t)
	})

	t.Run(testString("Evaluator/Add/CtPlain/ZT/", testctx.params), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectorsZT(testctx, testctx.encryptorPk, t)
		values2, plaintext2, ciphertext2 := newTestVectorsZT(testctx, testctx.encryptorPk, t)

		testctx.evaluator.Add(ciphertext1, plaintext2, ciphertext2)
		testctx.ringT.Add(values1, values2, values2)

		verifyTestVectors(testctx, testctx.decryptor, values2, ciphertext2, t)

		testctx.evaluator.Add(plaintext2, ciphertext1, ciphertext2)

		verifyTestVectors(testctx, testctx.decryptor, values2, ciphertext2, t)
	})

	t.Run(testString("Evaluator/Add/CtPlain/ZQ/", testctx.params), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectorsZQ(testctx, testctx.encryptorPk, t)
		values2, plaintext2, ciphertext2 := newTestVectorsZQ(testctx, testctx.encryptorPk, t)

		testctx.evaluator.Add(ciphertext1, plaintext2, ciphertext2)
		testctx.ringT.Add(values1, values2, values2)

		verifyTestVectors(testctx, testctx.decryptor, values2, ciphertext2, t)

		testctx.evaluator.Add(plaintext2, ciphertext1, ciphertext2)

		verifyTestVectors(testctx, testctx.decryptor, values2, ciphertext2, t)
	})

	t.Run(testString("Evaluator/Sub/CtCtInPlace/", testctx.params), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectorsZQ(testctx, testctx.encryptorPk, t)
		values2, _, ciphertext2 := newTestVectorsZQ(testctx, testctx.encryptorPk, t)

		testctx.evaluator.Sub(ciphertext1, ciphertext2, ciphertext1)
		testctx.ringT.Sub(values1, values2, values1)

		verifyTestVectors(testctx, testctx.decryptor, values1, ciphertext1, t)
	})

	t.Run(testString("Evaluator/Sub/CtCtNew/", testctx.params), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectorsZQ(testctx, testctx.encryptorPk, t)
		values2, _, ciphertext2 := newTestVectorsZQ(testctx, testctx.encryptorPk, t)

		ciphertext1 = testctx.evaluator.SubNew(ciphertext1, ciphertext2)
		testctx.ringT.Sub(values1, values2, values1)

		verifyTestVectors(testctx, testctx.decryptor, values1, ciphertext1, t)
	})

	t.Run(testString("Evaluator/Sub/CtPlain/ZT/", testctx.params), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectorsZT(testctx, testctx.encryptorPk, t)
		values2, plaintext2, ciphertext2 := newTestVectorsZT(testctx, testctx.encryptorPk, t)

		valuesWant := testctx.ringT.NewPoly()

		testctx.evaluator.Sub(ciphertext1, plaintext2, ciphertext2)
		testctx.ringT.Sub(values1, values2, valuesWant)
		verifyTestVectors(testctx, testctx.decryptor, valuesWant, ciphertext2, t)

		testctx.evaluator.Sub(plaintext2, ciphertext1, ciphertext2)
		testctx.ringT.Sub(values2, values1, valuesWant)
		verifyTestVectors(testctx, testctx.decryptor, valuesWant, ciphertext2, t)
	})

	t.Run(testString("Evaluator/Sub/CtPlain/ZQ/", testctx.params), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectorsZQ(testctx, testctx.encryptorPk, t)
		values2, plaintext2, ciphertext2 := newTestVectorsZQ(testctx, testctx.encryptorPk, t)

		valuesWant := testctx.ringT.NewPoly()

		testctx.evaluator.Sub(ciphertext1, plaintext2, ciphertext2)
		testctx.ringT.Sub(values1, values2, valuesWant)
		verifyTestVectors(testctx, testctx.decryptor, valuesWant, ciphertext2, t)

		testctx.evaluator.Sub(plaintext2, ciphertext1, ciphertext2)
		testctx.ringT.Sub(values2, values1, valuesWant)
		verifyTestVectors(testctx, testctx.decryptor, valuesWant, ciphertext2, t)
	})

	t.Run(testString("Evaluator/Mul/CtCt/", testctx.params), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectorsZQ(testctx, testctx.encryptorPk, t)
		values2, _, ciphertext2 := newTestVectorsZQ(testctx, testctx.encryptorPk, t)

		receiver := NewCiphertext(testctx.params, ciphertext1.Degree()+ciphertext2.Degree())
		testctx.evaluator.Mul(ciphertext1, ciphertext2, receiver)
		testctx.ringT.MulCoeffs(values1, values2, values1)

		verifyTestVectors(testctx, testctx.decryptor, values1, receiver, t)
	})

	t.Run(testString("Evaluator/Mul/CtPlain/ZT/", testctx.params), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectorsZT(testctx, testctx.encryptorPk, t)
		values2, plaintext2, _ := newTestVectorsZT(testctx, testctx.encryptorPk, t)

		testctx.evaluator.Mul(ciphertext1, plaintext2, ciphertext1)
		testctx.ringT.MulCoeffs(values1, values2, values1)

		verifyTestVectors(testctx, testctx.decryptor, values1, ciphertext1, t)
	})

	t.Run(testString("Evaluator/Mul/CtPlain/ZQ/", testctx.params), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectorsZQ(testctx, testctx.encryptorPk, t)
		values2, plaintext2, _ := newTestVectorsZQ(testctx, testctx.encryptorPk, t)

		testctx.evaluator.Mul(ciphertext1, plaintext2, ciphertext1)
		testctx.ringT.MulCoeffs(values1, values2, values1)

		verifyTestVectors(testctx, testctx.decryptor, values1, ciphertext1, t)
	})

	t.Run(testString("Evaluator/Mul/Relinearize/", testctx.params), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectorsZQ(testctx, testctx.encryptorPk, t)
		values2, _, ciphertext2 := newTestVectorsZQ(testctx, testctx.encryptorPk, t)

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

	sk2 := testctx.kgen.GenSecretKey()
	decryptorSk2 := NewDecryptor(testctx.params, sk2)
	switchKey := testctx.kgen.GenSwitchingKey(testctx.sk, sk2)

	t.Run(testString("Evaluator/KeySwitch/InPlace/", testctx.params), func(t *testing.T) {
		values, _, ciphertext := newTestVectorsZQ(testctx, testctx.encryptorPk, t)
		testctx.evaluator.SwitchKeys(ciphertext, switchKey, ciphertext)
		verifyTestVectors(testctx, decryptorSk2, values, ciphertext, t)
	})

	t.Run(testString("Evaluator/KeySwitch/New/", testctx.params), func(t *testing.T) {
		values, _, ciphertext := newTestVectorsZQ(testctx, testctx.encryptorPk, t)
		ciphertext = testctx.evaluator.SwitchKeysNew(ciphertext, switchKey)
		verifyTestVectors(testctx, decryptorSk2, values, ciphertext, t)
	})
}

func testEvaluatorRotate(testctx *testContext, t *testing.T) {

	rotkey := testctx.kgen.GenRotationKeysPow2(testctx.sk)

	t.Run(testString("Evaluator/Rotate/Rows/InPlace/", testctx.params), func(t *testing.T) {
		values, _, ciphertext := newTestVectorsZQ(testctx, testctx.encryptorPk, t)
		testctx.evaluator.RotateRows(ciphertext, rotkey, ciphertext)
		values.Coeffs[0] = append(values.Coeffs[0][testctx.params.N()>>1:], values.Coeffs[0][:testctx.params.N()>>1]...)
		verifyTestVectors(testctx, testctx.decryptor, values, ciphertext, t)
	})

	t.Run(testString("Evaluator/Rotate/Rows/New/", testctx.params), func(t *testing.T) {
		values, _, ciphertext := newTestVectorsZQ(testctx, testctx.encryptorPk, t)
		ciphertext = testctx.evaluator.RotateRowsNew(ciphertext, rotkey)
		values.Coeffs[0] = append(values.Coeffs[0][testctx.params.N()>>1:], values.Coeffs[0][:testctx.params.N()>>1]...)
		verifyTestVectors(testctx, testctx.decryptor, values, ciphertext, t)
	})

	valuesWant := testctx.ringT.NewPoly()
	mask := (testctx.params.N() >> 1) - 1
	slots := testctx.params.N() >> 1

	t.Run(testString("Evaluator/Rotate/Cols/InPlace/", testctx.params), func(t *testing.T) {

		values, _, ciphertext := newTestVectorsZQ(testctx, testctx.encryptorPk, t)

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

		values, _, ciphertext := newTestVectorsZQ(testctx, testctx.encryptorPk, t)

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

		values, _, ciphertext := newTestVectorsZQ(testctx, testctx.encryptorPk, t)

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
