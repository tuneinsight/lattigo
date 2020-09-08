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

func testString(opname string, params *Parameters) string {
	return fmt.Sprintf("%sLogN=%d/logQ=%d", opname, params.logN, params.logQP)
}

type testParams struct {
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
var params = new(testParams)
var defaultParams = DefaultParams[PN12QP109 : PN12QP109+3]

func TestBFV(t *testing.T) {

	rand.Seed(time.Now().UnixNano())

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

func genTestParams(contextParameters *Parameters) (err error) {

	params = new(testParams)
	params.params = contextParameters.Copy()

	if params.prng, err = utils.NewPRNG(); err != nil {
		return err
	}

	if params.ringQ, err = ring.NewRing(contextParameters.N(), contextParameters.qi); err != nil {
		return err
	}

	if params.ringQP, err = ring.NewRing(contextParameters.N(), append(contextParameters.qi, contextParameters.pi...)); err != nil {
		return err
	}

	if params.ringT, err = ring.NewRing(contextParameters.N(), []uint64{contextParameters.t}); err != nil {
		return err
	}

	params.uSampler = ring.NewUniformSampler(params.prng, params.ringT)
	params.kgen = NewKeyGenerator(params.params)
	params.sk, params.pk = params.kgen.GenKeyPair()
	params.rlk = params.kgen.GenRelinKey(params.sk, 1)
	params.encoder = NewEncoder(params.params)
	params.encryptorPk = NewEncryptorFromPk(params.params, params.pk)
	params.encryptorSk = NewEncryptorFromSk(params.params, params.sk)
	params.decryptor = NewDecryptor(params.params, params.sk)
	params.evaluator = NewEvaluator(params.params)
	return

}

func testParameters(t *testing.T) {
	t.Run("NewParametersFromModuli", func(t *testing.T) {
		p, err := NewParametersFromModuli(params.params.logN, params.params.Moduli, params.params.t)
		assert.NoError(t, err)
		assert.True(t, p.Equals(params.params))
	})

	t.Run("NewParametersFromLogModuli", func(t *testing.T) {
		p, err := NewParametersFromLogModuli(params.params.logN, params.params.LogModuli(), params.params.t)
		assert.NoError(t, err)
		assert.True(t, p.Equals(params.params))
	})
}

func newTestVectors(encryptor Encryptor, t *testing.T) (coeffs *ring.Poly, plaintext *Plaintext, ciphertext *Ciphertext) {

	coeffs = params.uSampler.ReadNew()

	plaintext = NewPlaintext(params.params)

	params.encoder.EncodeUint(coeffs.Coeffs[0], plaintext)

	if encryptor != nil {
		ciphertext = params.encryptorPk.EncryptNew(plaintext)
	}

	return coeffs, plaintext, ciphertext
}

func verifyTestVectors(decryptor Decryptor, coeffs *ring.Poly, element Operand, t *testing.T) {

	var coeffsTest []uint64

	el := element.Element()

	if el.Degree() == 0 {

		coeffsTest = params.encoder.DecodeUint(el.Plaintext())

	} else {

		coeffsTest = params.encoder.DecodeUint(decryptor.DecryptNew(el.Ciphertext()))
	}

	require.True(t, utils.EqualSliceUint64(coeffs.Coeffs[0], coeffsTest))
}

func testEncoder(t *testing.T) {
	t.Run(testString("Encode&Decode/", params.params), func(t *testing.T) {
		values, plaintext, _ := newTestVectors(nil, t)
		verifyTestVectors(params.decryptor, values, plaintext, t)
	})
}

func testEncryptor(t *testing.T) {

	t.Run(testString("EncryptFromPk/", params.params), func(t *testing.T) {
		values, _, ciphertext := newTestVectors(params.encryptorPk, t)
		verifyTestVectors(params.decryptor, values, ciphertext, t)
	})

	t.Run(testString("EncryptFromPkFast/", params.params), func(t *testing.T) {
		coeffs := params.uSampler.ReadNew()
		plaintext := NewPlaintext(params.params)
		params.encoder.EncodeUint(coeffs.Coeffs[0], plaintext)
		verifyTestVectors(params.decryptor, coeffs, params.encryptorPk.EncryptFastNew(plaintext), t)
	})

	t.Run(testString("EncryptFromSk/", params.params), func(t *testing.T) {
		coeffs := params.uSampler.ReadNew()
		plaintext := NewPlaintext(params.params)
		params.encoder.EncodeUint(coeffs.Coeffs[0], plaintext)
		verifyTestVectors(params.decryptor, coeffs, params.encryptorSk.EncryptNew(plaintext), t)
	})
}

func testEvaluatorAdd(t *testing.T) {

	t.Run(testString("CtCtInPlace/", params.params), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectors(params.encryptorPk, t)
		values2, _, ciphertext2 := newTestVectors(params.encryptorPk, t)

		params.evaluator.Add(ciphertext1, ciphertext2, ciphertext1)
		params.ringT.Add(values1, values2, values1)

		verifyTestVectors(params.decryptor, values1, ciphertext1, t)
	})

	t.Run(testString("CtCtNew/", params.params), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectors(params.encryptorPk, t)
		values2, _, ciphertext2 := newTestVectors(params.encryptorPk, t)

		ciphertext1 = params.evaluator.AddNew(ciphertext1, ciphertext2)
		params.ringT.Add(values1, values2, values1)

		verifyTestVectors(params.decryptor, values1, ciphertext1, t)
	})

	t.Run(testString("CtPlain/", params.params), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectors(params.encryptorPk, t)
		values2, plaintext2, ciphertext2 := newTestVectors(params.encryptorPk, t)

		params.evaluator.Add(ciphertext1, plaintext2, ciphertext2)
		params.ringT.Add(values1, values2, values2)

		verifyTestVectors(params.decryptor, values2, ciphertext2, t)

		params.evaluator.Add(plaintext2, ciphertext1, ciphertext2)

		verifyTestVectors(params.decryptor, values2, ciphertext2, t)
	})
}

func testEvaluatorSub(t *testing.T) {

	t.Run(testString("CtCtInPlace/", params.params), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectors(params.encryptorPk, t)
		values2, _, ciphertext2 := newTestVectors(params.encryptorPk, t)

		params.evaluator.Sub(ciphertext1, ciphertext2, ciphertext1)
		params.ringT.Sub(values1, values2, values1)

		verifyTestVectors(params.decryptor, values1, ciphertext1, t)
	})

	t.Run(testString("CtCtNew/", params.params), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectors(params.encryptorPk, t)
		values2, _, ciphertext2 := newTestVectors(params.encryptorPk, t)

		ciphertext1 = params.evaluator.SubNew(ciphertext1, ciphertext2)
		params.ringT.Sub(values1, values2, values1)

		verifyTestVectors(params.decryptor, values1, ciphertext1, t)
	})

	t.Run(testString("CtPlain/", params.params), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectors(params.encryptorPk, t)
		values2, plaintext2, ciphertext2 := newTestVectors(params.encryptorPk, t)

		valuesWant := params.ringT.NewPoly()

		params.evaluator.Sub(ciphertext1, plaintext2, ciphertext2)
		params.ringT.Sub(values1, values2, valuesWant)
		verifyTestVectors(params.decryptor, valuesWant, ciphertext2, t)

		params.evaluator.Sub(plaintext2, ciphertext1, ciphertext2)
		params.ringT.Sub(values2, values1, valuesWant)
		verifyTestVectors(params.decryptor, valuesWant, ciphertext2, t)
	})
}

func testEvaluatorMul(t *testing.T) {

	t.Run(testString("CtCt/", params.params), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectors(params.encryptorPk, t)
		values2, _, ciphertext2 := newTestVectors(params.encryptorPk, t)

		receiver := NewCiphertext(params.params, ciphertext1.Degree()+ciphertext2.Degree())
		params.evaluator.Mul(ciphertext1, ciphertext2, receiver)
		params.ringT.MulCoeffs(values1, values2, values1)

		verifyTestVectors(params.decryptor, values1, receiver, t)
	})

	t.Run(testString("CtPlain/", params.params), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectors(params.encryptorPk, t)
		values2, plaintext2, _ := newTestVectors(params.encryptorPk, t)

		params.evaluator.Mul(ciphertext1, plaintext2, ciphertext1)
		params.ringT.MulCoeffs(values1, values2, values1)

		verifyTestVectors(params.decryptor, values1, ciphertext1, t)
	})

	t.Run(testString("Relinearize/", params.params), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectors(params.encryptorPk, t)
		values2, _, ciphertext2 := newTestVectors(params.encryptorPk, t)

		receiver := NewCiphertext(params.params, ciphertext1.Degree()+ciphertext2.Degree())
		params.evaluator.Mul(ciphertext1, ciphertext2, receiver)
		params.ringT.MulCoeffs(values1, values2, values1)

		receiver2 := params.evaluator.RelinearizeNew(receiver, params.rlk)
		verifyTestVectors(params.decryptor, values1, receiver2, t)

		params.evaluator.Relinearize(receiver, params.rlk, receiver)
		verifyTestVectors(params.decryptor, values1, receiver, t)
	})
}

func testKeySwitch(t *testing.T) {

	sk2 := params.kgen.GenSecretKey()
	decryptorSk2 := NewDecryptor(params.params, sk2)
	switchKey := params.kgen.GenSwitchingKey(params.sk, sk2)

	t.Run(testString("InPlace/", params.params), func(t *testing.T) {
		values, _, ciphertext := newTestVectors(params.encryptorPk, t)
		params.evaluator.SwitchKeys(ciphertext, switchKey, ciphertext)
		verifyTestVectors(decryptorSk2, values, ciphertext, t)
	})

	t.Run(testString("New/", params.params), func(t *testing.T) {
		values, _, ciphertext := newTestVectors(params.encryptorPk, t)
		ciphertext = params.evaluator.SwitchKeysNew(ciphertext, switchKey)
		verifyTestVectors(decryptorSk2, values, ciphertext, t)
	})
}

func testRotateRows(t *testing.T) {

	rotkey := NewRotationKeys()
	params.kgen.GenRot(RotationRow, params.sk, 0, rotkey)

	t.Run(testString("InPlace/", params.params), func(t *testing.T) {
		values, _, ciphertext := newTestVectors(params.encryptorPk, t)
		params.evaluator.RotateRows(ciphertext, rotkey, ciphertext)
		values.Coeffs[0] = append(values.Coeffs[0][params.params.N()>>1:], values.Coeffs[0][:params.params.N()>>1]...)
		verifyTestVectors(params.decryptor, values, ciphertext, t)
	})

	t.Run(testString("New/", params.params), func(t *testing.T) {
		values, _, ciphertext := newTestVectors(params.encryptorPk, t)
		ciphertext = params.evaluator.RotateRowsNew(ciphertext, rotkey)
		values.Coeffs[0] = append(values.Coeffs[0][params.params.N()>>1:], values.Coeffs[0][:params.params.N()>>1]...)
		verifyTestVectors(params.decryptor, values, ciphertext, t)
	})
}

func testRotateCols(t *testing.T) {

	rotkey := params.kgen.GenRotationKeysPow2(params.sk)

	valuesWant := params.ringT.NewPoly()
	mask := (params.params.N() >> 1) - 1
	slots := params.params.N() >> 1

	t.Run(testString("InPlace/", params.params), func(t *testing.T) {

		values, _, ciphertext := newTestVectors(params.encryptorPk, t)

		receiver := NewCiphertext(params.params, 1)
		for n := uint64(1); n < slots; n <<= 1 {

			params.evaluator.RotateColumns(ciphertext, n, rotkey, receiver)

			for i := uint64(0); i < slots; i++ {
				valuesWant.Coeffs[0][i] = values.Coeffs[0][(i+n)&mask]
				valuesWant.Coeffs[0][i+slots] = values.Coeffs[0][((i+n)&mask)+slots]
			}

			verifyTestVectors(params.decryptor, valuesWant, receiver, t)
		}
	})

	t.Run(testString("New/", params.params), func(t *testing.T) {

		values, _, ciphertext := newTestVectors(params.encryptorPk, t)

		for n := uint64(1); n < slots; n <<= 1 {

			receiver := params.evaluator.RotateColumnsNew(ciphertext, n, rotkey)

			for i := uint64(0); i < slots; i++ {
				valuesWant.Coeffs[0][i] = values.Coeffs[0][(i+n)&mask]
				valuesWant.Coeffs[0][i+slots] = values.Coeffs[0][((i+n)&mask)+slots]
			}

			verifyTestVectors(params.decryptor, valuesWant, receiver, t)
		}
	})

	t.Run(testString("Random/", params.params), func(t *testing.T) {

		values, _, ciphertext := newTestVectors(params.encryptorPk, t)

		receiver := NewCiphertext(params.params, 1)
		prng, err := utils.NewPRNG()
		if err != nil {
			panic(err)
		}

		for n := 0; n < 4; n++ {

			rand := ring.RandUniform(prng, slots, mask)

			params.evaluator.RotateColumns(ciphertext, rand, rotkey, receiver)

			for i := uint64(0); i < slots; i++ {
				valuesWant.Coeffs[0][i] = values.Coeffs[0][(i+rand)&mask]
				valuesWant.Coeffs[0][i+slots] = values.Coeffs[0][((i+rand)&mask)+slots]
			}

			verifyTestVectors(params.decryptor, valuesWant, receiver, t)
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

	ringQP := params.ringQP

	t.Run(testString("Ciphertext/", params.params), func(t *testing.T) {

		ciphertextWant := NewCiphertextRandom(params.prng, params.params, 2)

		marshalledCiphertext, err := ciphertextWant.MarshalBinary()
		require.NoError(t, err)

		ciphertextTest := new(Ciphertext)
		err = ciphertextTest.UnmarshalBinary(marshalledCiphertext)
		require.NoError(t, err)

		for i := range ciphertextWant.value {
			require.True(t, params.ringQ.Equal(ciphertextWant.value[i], ciphertextTest.value[i]))
		}
	})

	t.Run(testString("Sk/", params.params), func(t *testing.T) {

		marshalledSk, err := params.sk.MarshalBinary()
		require.NoError(t, err)

		sk := new(SecretKey)
		err = sk.UnmarshalBinary(marshalledSk)
		require.NoError(t, err)

		require.True(t, ringQP.Equal(sk.sk, params.sk.sk))
	})

	t.Run(testString("Pk/", params.params), func(t *testing.T) {

		marshalledPk, err := params.pk.MarshalBinary()
		require.NoError(t, err)

		pk := new(PublicKey)
		err = pk.UnmarshalBinary(marshalledPk)
		require.NoError(t, err)

		for k := range params.pk.pk {
			require.True(t, ringQP.Equal(pk.pk[k], params.pk.pk[k]), k)
		}
	})

	t.Run(testString("EvaluationKey/", params.params), func(t *testing.T) {

		evalkey := params.kgen.GenRelinKey(params.sk, 2)
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

	t.Run(testString("SwitchingKey/", params.params), func(t *testing.T) {

		skOut := params.kgen.GenSecretKey()

		switchingKey := params.kgen.GenSwitchingKey(params.sk, skOut)
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

	t.Run(testString("RotationKey/", params.params), func(t *testing.T) {

		rotationKey := NewRotationKeys()

		params.kgen.GenRot(RotationRow, params.sk, 0, rotationKey)
		params.kgen.GenRot(RotationLeft, params.sk, 1, rotationKey)
		params.kgen.GenRot(RotationLeft, params.sk, 2, rotationKey)
		params.kgen.GenRot(RotationRight, params.sk, 3, rotationKey)
		params.kgen.GenRot(RotationRight, params.sk, 5, rotationKey)

		data, err := rotationKey.MarshalBinary()
		require.NoError(t, err)

		resRotationKey := new(RotationKeys)
		err = resRotationKey.UnmarshalBinary(data)
		require.NoError(t, err)

		for i := uint64(1); i < params.params.N()>>1; i++ {

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
