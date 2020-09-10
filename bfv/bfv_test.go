package bfv

import (
	"fmt"
	"math/rand"
	"testing"
	"time"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"

	"github.com/ldsec/lattigo/ring"
	"github.com/ldsec/lattigo/utils"
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

var err error
var testctx = new(testContext)

func TestBFV(t *testing.T) {

	rand.Seed(time.Now().UnixNano())

	var defaultParams []*Parameters

	if testing.Short() {
		defaultParams = DefaultParams[PN12QP109 : PN12QP109+3]
	} else {
		defaultParams = DefaultParams
	}

	for _, p := range defaultParams {

		if err = genTestParams(p); err != nil {
			panic(err)
		}

		t.Run("Parameters", testParameters)
		t.Run("Encoder", testEncoder)
		t.Run("Encryptor", testEncryptor)
		t.Run("Evaluator/Add", testEvaluatorAdd)
		t.Run("Evaluator/Sub", testEvaluatorSub)
		t.Run("Evaluator/Mul", testEvaluatorMul)
		t.Run("Evaluator/KeySwitch", testKeySwitch)
		t.Run("Evaluator/RotateRows", testRotateRows)
		t.Run("Evaluator/RotateCols", testRotateCols)
		t.Run("Marshalling", testMarshaller)
	}

}

func genTestParams(params *Parameters) (err error) {

	testctx = new(testContext)
	testctx.params = params.Copy()

	if testctx.prng, err = utils.NewPRNG(); err != nil {
		return err
	}

	if testctx.ringQ, err = ring.NewRing(params.N(), params.qi); err != nil {
		return err
	}

	if testctx.ringQP, err = ring.NewRing(params.N(), append(params.qi, params.pi...)); err != nil {
		return err
	}

	if testctx.ringT, err = ring.NewRing(params.N(), []uint64{params.t}); err != nil {
		return err
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

func testParameters(t *testing.T) {
	t.Run("NewParametersFromModuli", func(t *testing.T) {
		p, err := NewParametersFromModuli(testctx.params.logN, testctx.params.Moduli, testctx.params.t)
		assert.NoError(t, err)
		assert.True(t, p.Equals(testctx.params))
	})

	t.Run("NewParametersFromLogModuli", func(t *testing.T) {
		p, err := NewParametersFromLogModuli(testctx.params.logN, testctx.params.LogModuli(), testctx.params.t)
		assert.NoError(t, err)
		assert.True(t, p.Equals(testctx.params))
	})
}

func newTestVectors(encryptor Encryptor, t *testing.T) (coeffs *ring.Poly, plaintext *Plaintext, ciphertext *Ciphertext) {

	coeffs = testctx.uSampler.ReadNew()

	plaintext = NewPlaintext(testctx.params)

	testctx.encoder.EncodeUint(coeffs.Coeffs[0], plaintext)

	if encryptor != nil {
		ciphertext = testctx.encryptorPk.EncryptNew(plaintext)
	}

	return coeffs, plaintext, ciphertext
}

func verifyTestVectors(decryptor Decryptor, coeffs *ring.Poly, element Operand, t *testing.T) {

	var coeffsTest []uint64

	el := element.Element()

	if el.Degree() == 0 {

		coeffsTest = testctx.encoder.DecodeUint(el.Plaintext())

	} else {

		coeffsTest = testctx.encoder.DecodeUint(decryptor.DecryptNew(el.Ciphertext()))
	}

	require.True(t, utils.EqualSliceUint64(coeffs.Coeffs[0], coeffsTest))
}

func testEncoder(t *testing.T) {
	t.Run(testString("Encode&Decode/", testctx.params), func(t *testing.T) {
		values, plaintext, _ := newTestVectors(nil, t)
		verifyTestVectors(testctx.decryptor, values, plaintext, t)
	})
}

func testEncryptor(t *testing.T) {

	t.Run(testString("EncryptFromPk/", testctx.params), func(t *testing.T) {
		values, _, ciphertext := newTestVectors(testctx.encryptorPk, t)
		verifyTestVectors(testctx.decryptor, values, ciphertext, t)
	})

	t.Run(testString("EncryptFromPkFast/", testctx.params), func(t *testing.T) {
		coeffs := testctx.uSampler.ReadNew()
		plaintext := NewPlaintext(testctx.params)
		testctx.encoder.EncodeUint(coeffs.Coeffs[0], plaintext)
		verifyTestVectors(testctx.decryptor, coeffs, testctx.encryptorPk.EncryptFastNew(plaintext), t)
	})

	t.Run(testString("EncryptFromSk/", testctx.params), func(t *testing.T) {
		coeffs := testctx.uSampler.ReadNew()
		plaintext := NewPlaintext(testctx.params)
		testctx.encoder.EncodeUint(coeffs.Coeffs[0], plaintext)
		verifyTestVectors(testctx.decryptor, coeffs, testctx.encryptorSk.EncryptNew(plaintext), t)
	})
}

func testEvaluatorAdd(t *testing.T) {

	t.Run(testString("CtCtInPlace/", testctx.params), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectors(testctx.encryptorPk, t)
		values2, _, ciphertext2 := newTestVectors(testctx.encryptorPk, t)

		testctx.evaluator.Add(ciphertext1, ciphertext2, ciphertext1)
		testctx.ringT.Add(values1, values2, values1)

		verifyTestVectors(testctx.decryptor, values1, ciphertext1, t)
	})

	t.Run(testString("CtCtNew/", testctx.params), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectors(testctx.encryptorPk, t)
		values2, _, ciphertext2 := newTestVectors(testctx.encryptorPk, t)

		ciphertext1 = testctx.evaluator.AddNew(ciphertext1, ciphertext2)
		testctx.ringT.Add(values1, values2, values1)

		verifyTestVectors(testctx.decryptor, values1, ciphertext1, t)
	})

	t.Run(testString("CtPlain/", testctx.params), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectors(testctx.encryptorPk, t)
		values2, plaintext2, ciphertext2 := newTestVectors(testctx.encryptorPk, t)

		testctx.evaluator.Add(ciphertext1, plaintext2, ciphertext2)
		testctx.ringT.Add(values1, values2, values2)

		verifyTestVectors(testctx.decryptor, values2, ciphertext2, t)

		testctx.evaluator.Add(plaintext2, ciphertext1, ciphertext2)

		verifyTestVectors(testctx.decryptor, values2, ciphertext2, t)
	})
}

func testEvaluatorSub(t *testing.T) {

	t.Run(testString("CtCtInPlace/", testctx.params), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectors(testctx.encryptorPk, t)
		values2, _, ciphertext2 := newTestVectors(testctx.encryptorPk, t)

		testctx.evaluator.Sub(ciphertext1, ciphertext2, ciphertext1)
		testctx.ringT.Sub(values1, values2, values1)

		verifyTestVectors(testctx.decryptor, values1, ciphertext1, t)
	})

	t.Run(testString("CtCtNew/", testctx.params), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectors(testctx.encryptorPk, t)
		values2, _, ciphertext2 := newTestVectors(testctx.encryptorPk, t)

		ciphertext1 = testctx.evaluator.SubNew(ciphertext1, ciphertext2)
		testctx.ringT.Sub(values1, values2, values1)

		verifyTestVectors(testctx.decryptor, values1, ciphertext1, t)
	})

	t.Run(testString("CtPlain/", testctx.params), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectors(testctx.encryptorPk, t)
		values2, plaintext2, ciphertext2 := newTestVectors(testctx.encryptorPk, t)

		valuesWant := testctx.ringT.NewPoly()

		testctx.evaluator.Sub(ciphertext1, plaintext2, ciphertext2)
		testctx.ringT.Sub(values1, values2, valuesWant)
		verifyTestVectors(testctx.decryptor, valuesWant, ciphertext2, t)

		testctx.evaluator.Sub(plaintext2, ciphertext1, ciphertext2)
		testctx.ringT.Sub(values2, values1, valuesWant)
		verifyTestVectors(testctx.decryptor, valuesWant, ciphertext2, t)
	})
}

func testEvaluatorMul(t *testing.T) {

	t.Run(testString("CtCt/", testctx.params), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectors(testctx.encryptorPk, t)
		values2, _, ciphertext2 := newTestVectors(testctx.encryptorPk, t)

		receiver := NewCiphertext(testctx.params, ciphertext1.Degree()+ciphertext2.Degree())
		testctx.evaluator.Mul(ciphertext1, ciphertext2, receiver)
		testctx.ringT.MulCoeffs(values1, values2, values1)

		verifyTestVectors(testctx.decryptor, values1, receiver, t)
	})

	t.Run(testString("CtPlain/", testctx.params), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectors(testctx.encryptorPk, t)
		values2, plaintext2, _ := newTestVectors(testctx.encryptorPk, t)

		testctx.evaluator.Mul(ciphertext1, plaintext2, ciphertext1)
		testctx.ringT.MulCoeffs(values1, values2, values1)

		verifyTestVectors(testctx.decryptor, values1, ciphertext1, t)
	})

	t.Run(testString("Relinearize/", testctx.params), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectors(testctx.encryptorPk, t)
		values2, _, ciphertext2 := newTestVectors(testctx.encryptorPk, t)

		receiver := NewCiphertext(testctx.params, ciphertext1.Degree()+ciphertext2.Degree())
		testctx.evaluator.Mul(ciphertext1, ciphertext2, receiver)
		testctx.ringT.MulCoeffs(values1, values2, values1)

		receiver2 := testctx.evaluator.RelinearizeNew(receiver, testctx.rlk)
		verifyTestVectors(testctx.decryptor, values1, receiver2, t)

		testctx.evaluator.Relinearize(receiver, testctx.rlk, receiver)
		verifyTestVectors(testctx.decryptor, values1, receiver, t)
	})
}

func testKeySwitch(t *testing.T) {

	sk2 := testctx.kgen.GenSecretKey()
	decryptorSk2 := NewDecryptor(testctx.params, sk2)
	switchKey := testctx.kgen.GenSwitchingKey(testctx.sk, sk2)

	t.Run(testString("InPlace/", testctx.params), func(t *testing.T) {
		values, _, ciphertext := newTestVectors(testctx.encryptorPk, t)
		testctx.evaluator.SwitchKeys(ciphertext, switchKey, ciphertext)
		verifyTestVectors(decryptorSk2, values, ciphertext, t)
	})

	t.Run(testString("New/", testctx.params), func(t *testing.T) {
		values, _, ciphertext := newTestVectors(testctx.encryptorPk, t)
		ciphertext = testctx.evaluator.SwitchKeysNew(ciphertext, switchKey)
		verifyTestVectors(decryptorSk2, values, ciphertext, t)
	})
}

func testRotateRows(t *testing.T) {

	rotkey := NewRotationKeys()
	testctx.kgen.GenRot(RotationRow, testctx.sk, 0, rotkey)

	t.Run(testString("InPlace/", testctx.params), func(t *testing.T) {
		values, _, ciphertext := newTestVectors(testctx.encryptorPk, t)
		testctx.evaluator.RotateRows(ciphertext, rotkey, ciphertext)
		values.Coeffs[0] = append(values.Coeffs[0][testctx.params.N()>>1:], values.Coeffs[0][:testctx.params.N()>>1]...)
		verifyTestVectors(testctx.decryptor, values, ciphertext, t)
	})

	t.Run(testString("New/", testctx.params), func(t *testing.T) {
		values, _, ciphertext := newTestVectors(testctx.encryptorPk, t)
		ciphertext = testctx.evaluator.RotateRowsNew(ciphertext, rotkey)
		values.Coeffs[0] = append(values.Coeffs[0][testctx.params.N()>>1:], values.Coeffs[0][:testctx.params.N()>>1]...)
		verifyTestVectors(testctx.decryptor, values, ciphertext, t)
	})
}

func testRotateCols(t *testing.T) {

	rotkey := testctx.kgen.GenRotationKeysPow2(testctx.sk)

	valuesWant := testctx.ringT.NewPoly()
	mask := (testctx.params.N() >> 1) - 1
	slots := testctx.params.N() >> 1

	t.Run(testString("InPlace/", testctx.params), func(t *testing.T) {

		values, _, ciphertext := newTestVectors(testctx.encryptorPk, t)

		receiver := NewCiphertext(testctx.params, 1)
		for n := uint64(1); n < slots; n <<= 1 {

			testctx.evaluator.RotateColumns(ciphertext, n, rotkey, receiver)

			for i := uint64(0); i < slots; i++ {
				valuesWant.Coeffs[0][i] = values.Coeffs[0][(i+n)&mask]
				valuesWant.Coeffs[0][i+slots] = values.Coeffs[0][((i+n)&mask)+slots]
			}

			verifyTestVectors(testctx.decryptor, valuesWant, receiver, t)
		}
	})

	t.Run(testString("New/", testctx.params), func(t *testing.T) {

		values, _, ciphertext := newTestVectors(testctx.encryptorPk, t)

		for n := uint64(1); n < slots; n <<= 1 {

			receiver := testctx.evaluator.RotateColumnsNew(ciphertext, n, rotkey)

			for i := uint64(0); i < slots; i++ {
				valuesWant.Coeffs[0][i] = values.Coeffs[0][(i+n)&mask]
				valuesWant.Coeffs[0][i+slots] = values.Coeffs[0][((i+n)&mask)+slots]
			}

			verifyTestVectors(testctx.decryptor, valuesWant, receiver, t)
		}
	})

	t.Run(testString("Random/", testctx.params), func(t *testing.T) {

		values, _, ciphertext := newTestVectors(testctx.encryptorPk, t)

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

			verifyTestVectors(testctx.decryptor, valuesWant, receiver, t)
		}
	})
}

func testMarshaller(t *testing.T) {

	t.Run("Parameters/ZeroValue", func(t *testing.T) {
		bytes, err := (&Parameters{}).MarshalBinary()
		assert.Nil(t, err)
		assert.Equal(t, []byte{}, bytes)
		p := new(Parameters)
		err = p.UnmarshalBinary(bytes)
		assert.NotNil(t, err)
	})

	t.Run("Parameters/SupportedParams", func(t *testing.T) {
		for _, params := range DefaultParams {
			bytes, err := params.MarshalBinary()
			assert.Nil(t, err)
			p := new(Parameters)
			err = p.UnmarshalBinary(bytes)
			assert.Nil(t, err)
			assert.Equal(t, params, p)
		}
	})

	ringQP := testctx.ringQP

	t.Run(testString("Ciphertext/", testctx.params), func(t *testing.T) {

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

	t.Run(testString("Sk/", testctx.params), func(t *testing.T) {

		marshalledSk, err := testctx.sk.MarshalBinary()
		require.NoError(t, err)

		sk := new(SecretKey)
		err = sk.UnmarshalBinary(marshalledSk)
		require.NoError(t, err)

		require.True(t, ringQP.Equal(sk.sk, testctx.sk.sk))
	})

	t.Run(testString("Pk/", testctx.params), func(t *testing.T) {

		marshalledPk, err := testctx.pk.MarshalBinary()
		require.NoError(t, err)

		pk := new(PublicKey)
		err = pk.UnmarshalBinary(marshalledPk)
		require.NoError(t, err)

		for k := range testctx.pk.pk {
			require.True(t, ringQP.Equal(pk.pk[k], testctx.pk.pk[k]), k)
		}
	})

	t.Run(testString("EvaluationKey/", testctx.params), func(t *testing.T) {

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

	t.Run(testString("SwitchingKey/", testctx.params), func(t *testing.T) {

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

	t.Run(testString("RotationKey/", testctx.params), func(t *testing.T) {

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
