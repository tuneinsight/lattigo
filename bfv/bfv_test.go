package bfv

import (
	"fmt"
	"log"
	"math/rand"
	"testing"
	"time"

	"github.com/ldsec/lattigo/ring"
	"github.com/ldsec/lattigo/utils"
)

func check(t *testing.T, err error) {
	if err != nil {
		t.Error(err)
		log.Fatal(err)
	}
}

func testString(opname string, params *Parameters) string {
	return fmt.Sprintf("%sLogN=%d/logQ=%d", opname, params.LogN, params.logQP)
}

type bfvParams struct {
	params      *Parameters
	bfvContext  *bfvContext
	encoder     *Encoder
	kgen        *KeyGenerator
	sk          *SecretKey
	pk          *PublicKey
	encryptorPk Encryptor
	encryptorSk Encryptor
	decryptor   Decryptor
	evaluator   *Evaluator
}

type bfvTestParameters struct {
	bfvParameters []*Parameters
}

var err error
var testParams = new(bfvTestParameters)

func init() {
	rand.Seed(time.Now().UnixNano())

	testParams.bfvParameters = []*Parameters{
		DefaultParams[PN12QP109],
		DefaultParams[PN13QP218],
		DefaultParams[PN14QP438],
		DefaultParams[PN15QP880],
	}
}

func TestBFV(t *testing.T) {
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

func testMarshaller(t *testing.T) {

	for _, parameters := range testParams.bfvParameters {

		params := genBfvParams(parameters)

		contextQP := params.bfvContext.contextQP

		t.Run(testString("Ciphertext/", parameters), func(t *testing.T) {

			ciphertextWant := NewCiphertextRandom(parameters, 2)

			marshalledCiphertext, err := ciphertextWant.MarshalBinary()
			check(t, err)

			ciphertextTest := newBfvElement(params.params, 1).Ciphertext()
			err = ciphertextTest.UnmarshalBinary(marshalledCiphertext)
			check(t, err)

			for i := range ciphertextWant.value {
				if !params.bfvContext.contextQ.Equal(ciphertextWant.value[i], ciphertextTest.value[i]) {
					t.Errorf("marshal Ciphertext")
				}
			}
		})

		t.Run(testString("Sk/", parameters), func(t *testing.T) {

			marshalledSk, err := params.sk.MarshalBinary()
			check(t, err)

			sk := new(SecretKey)
			err = sk.UnmarshalBinary(marshalledSk)
			check(t, err)

			if !contextQP.Equal(sk.sk, params.sk.sk) {
				t.Errorf("marshal SecretKey")
			}

		})

		t.Run(testString("Pk/", parameters), func(t *testing.T) {

			marshalledPk, err := params.pk.MarshalBinary()
			check(t, err)

			pk := new(PublicKey)
			err = pk.UnmarshalBinary(marshalledPk)
			check(t, err)

			for k := range params.pk.pk {
				if !contextQP.Equal(pk.pk[k], params.pk.pk[k]) {
					t.Errorf("marshal PublicKey element [%d]", k)
				}
			}
		})

		t.Run(testString("EvaluationKey/", parameters), func(t *testing.T) {

			evalkey := params.kgen.NewRelinKey(params.sk, 2)
			data, err := evalkey.MarshalBinary()
			check(t, err)

			resEvalKey := new(EvaluationKey)
			err = resEvalKey.UnmarshalBinary(data)
			check(t, err)

			for deg := range evalkey.evakey {

				evakeyWant := evalkey.evakey[deg].evakey
				evakeyTest := resEvalKey.evakey[deg].evakey

				for j := range evakeyWant {

					for k := range evakeyWant[j] {
						if !contextQP.Equal(evakeyWant[j][k], evakeyTest[j][k]) {
							t.Errorf("marshal EvaluationKey deg %d element [%d][%d]", deg, j, k)
						}
					}
				}
			}
		})

		t.Run(testString("SwitchingKey/", parameters), func(t *testing.T) {

			skOut := params.kgen.NewSecretKey()

			switchingKey := params.kgen.NewSwitchingKey(params.sk, skOut)
			data, err := switchingKey.MarshalBinary()
			check(t, err)

			resSwitchingKey := new(SwitchingKey)
			err = resSwitchingKey.UnmarshalBinary(data)
			check(t, err)

			evakeyWant := switchingKey.evakey
			evakeyTest := resSwitchingKey.evakey

			for j := range evakeyWant {

				for k := range evakeyWant[j] {
					if !contextQP.Equal(evakeyWant[j][k], evakeyTest[j][k]) {
						t.Errorf("marshal SwitchingKey element [%d][%d]", j, k)
					}
				}
			}
		})

		t.Run(testString("RotationKey/", parameters), func(t *testing.T) {

			rotationKey := NewRotationKeys()

			params.kgen.GenRot(RotationRow, params.sk, 0, rotationKey)
			params.kgen.GenRot(RotationLeft, params.sk, 1, rotationKey)
			params.kgen.GenRot(RotationLeft, params.sk, 2, rotationKey)
			params.kgen.GenRot(RotationRight, params.sk, 3, rotationKey)
			params.kgen.GenRot(RotationRight, params.sk, 5, rotationKey)

			data, err := rotationKey.MarshalBinary()
			check(t, err)

			resRotationKey := new(RotationKeys)
			err = resRotationKey.UnmarshalBinary(data)
			check(t, err)

			for i := uint64(1); i < params.bfvContext.n>>1; i++ {

				if rotationKey.evakeyRotColLeft[i] != nil {

					evakeyWant := rotationKey.evakeyRotColLeft[i].evakey
					evakeyTest := resRotationKey.evakeyRotColLeft[i].evakey

					for j := range evakeyWant {

						for k := range evakeyWant[j] {
							if !contextQP.Equal(evakeyWant[j][k], evakeyTest[j][k]) {
								t.Errorf("marshal RotationKey RotateLeft %d element [%d][%d]", i, j, k)
							}
						}
					}
				}

				if rotationKey.evakeyRotColRight[i] != nil {

					evakeyWant := rotationKey.evakeyRotColRight[i].evakey
					evakeyTest := resRotationKey.evakeyRotColRight[i].evakey

					for j := range evakeyWant {

						for k := range evakeyWant[j] {
							if !contextQP.Equal(evakeyWant[j][k], evakeyTest[j][k]) {
								t.Errorf("marshal RotationKey RotateRight %d element [%d][%d]", i, j, k)
							}
						}
					}
				}
			}

			if rotationKey.evakeyRotRow != nil {

				evakeyWant := rotationKey.evakeyRotRow.evakey
				evakeyTest := resRotationKey.evakeyRotRow.evakey

				for j := range evakeyWant {

					for k := range evakeyWant[j] {
						if !contextQP.Equal(evakeyWant[j][k], evakeyTest[j][k]) {
							t.Errorf("marshal RotationKey RotateRow element [%d][%d]", j, k)
						}
					}
				}
			}
		})
	}
}

func genBfvParams(contextParameters *Parameters) (params *bfvParams) {

	params = new(bfvParams)

	params.params = contextParameters.Copy()

	params.bfvContext = newBFVContext(contextParameters)

	params.kgen = NewKeyGenerator(contextParameters)

	params.sk, params.pk = params.kgen.NewKeyPair()

	params.encoder = NewEncoder(contextParameters)

	params.encryptorPk = NewEncryptorFromPk(contextParameters, params.pk)
	params.encryptorSk = NewEncryptorFromSk(contextParameters, params.sk)
	params.decryptor = NewDecryptor(contextParameters, params.sk)

	params.evaluator = NewEvaluator(contextParameters)

	return

}

func newTestVectors(params *bfvParams, encryptor Encryptor, t *testing.T) (coeffs *ring.Poly, plaintext *Plaintext, ciphertext *Ciphertext) {

	coeffs = params.bfvContext.contextT.NewUniformPoly()

	plaintext = NewPlaintext(params.params)

	params.encoder.EncodeUint(coeffs.Coeffs[0], plaintext)

	if encryptor != nil {
		ciphertext = params.encryptorPk.EncryptNew(plaintext)
	}

	return coeffs, plaintext, ciphertext
}

func verifyTestVectors(params *bfvParams, decryptor Decryptor, coeffs *ring.Poly, element Operand, t *testing.T) {

	var coeffsTest []uint64

	el := element.Element()

	if el.Degree() == 0 {

		coeffsTest = params.encoder.DecodeUint(el.Plaintext())

	} else {

		coeffsTest = params.encoder.DecodeUint(decryptor.DecryptNew(el.Ciphertext()))
	}

	if utils.EqualSliceUint64(coeffs.Coeffs[0], coeffsTest) != true {
		t.Errorf("decryption error")
	}
}

func testEncoder(t *testing.T) {

	for _, parameters := range testParams.bfvParameters {

		params := genBfvParams(parameters)

		t.Run(testString("Encode&Decode/", parameters), func(t *testing.T) {

			values, plaintext, _ := newTestVectors(params, nil, t)

			verifyTestVectors(params, params.decryptor, values, plaintext, t)
		})
	}
}

func testEncryptor(t *testing.T) {

	for _, parameters := range testParams.bfvParameters {

		params := genBfvParams(parameters)

		t.Run(testString("EncryptFromPk/", parameters), func(t *testing.T) {

			values, _, ciphertext := newTestVectors(params, params.encryptorPk, t)

			verifyTestVectors(params, params.decryptor, values, ciphertext, t)
		})

		t.Run(testString("EncryptFromSk/", parameters), func(t *testing.T) {

			values, _, ciphertext := newTestVectors(params, params.encryptorSk, t)

			verifyTestVectors(params, params.decryptor, values, ciphertext, t)
		})
	}
}

func testEvaluatorAdd(t *testing.T) {

	for _, parameters := range testParams.bfvParameters {

		params := genBfvParams(parameters)

		t.Run(testString("CtCtInPlace/", parameters), func(t *testing.T) {

			values1, _, ciphertext1 := newTestVectors(params, params.encryptorPk, t)
			values2, _, ciphertext2 := newTestVectors(params, params.encryptorPk, t)

			params.evaluator.Add(ciphertext1, ciphertext2, ciphertext1)
			params.bfvContext.contextT.Add(values1, values2, values1)

			verifyTestVectors(params, params.decryptor, values1, ciphertext1, t)
		})

		t.Run(testString("CtCtNew/", parameters), func(t *testing.T) {

			values1, _, ciphertext1 := newTestVectors(params, params.encryptorPk, t)
			values2, _, ciphertext2 := newTestVectors(params, params.encryptorPk, t)

			ciphertext1 = params.evaluator.AddNew(ciphertext1, ciphertext2)
			params.bfvContext.contextT.Add(values1, values2, values1)

			verifyTestVectors(params, params.decryptor, values1, ciphertext1, t)
		})

		t.Run(testString("CtPlain/", parameters), func(t *testing.T) {

			values1, _, ciphertext1 := newTestVectors(params, params.encryptorPk, t)
			values2, plaintext2, ciphertext2 := newTestVectors(params, params.encryptorPk, t)

			params.evaluator.Add(ciphertext1, plaintext2, ciphertext2)
			params.bfvContext.contextT.Add(values1, values2, values2)

			verifyTestVectors(params, params.decryptor, values2, ciphertext2, t)

			params.evaluator.Add(plaintext2, ciphertext1, ciphertext2)

			verifyTestVectors(params, params.decryptor, values2, ciphertext2, t)
		})
	}
}

func testEvaluatorSub(t *testing.T) {

	for _, parameters := range testParams.bfvParameters {

		params := genBfvParams(parameters)

		t.Run(testString("CtCtInPlace/", parameters), func(t *testing.T) {

			values1, _, ciphertext1 := newTestVectors(params, params.encryptorPk, t)
			values2, _, ciphertext2 := newTestVectors(params, params.encryptorPk, t)

			params.evaluator.Sub(ciphertext1, ciphertext2, ciphertext1)
			params.bfvContext.contextT.Sub(values1, values2, values1)

			verifyTestVectors(params, params.decryptor, values1, ciphertext1, t)
		})

		t.Run(testString("CtCtNew/", parameters), func(t *testing.T) {

			values1, _, ciphertext1 := newTestVectors(params, params.encryptorPk, t)
			values2, _, ciphertext2 := newTestVectors(params, params.encryptorPk, t)

			ciphertext1 = params.evaluator.SubNew(ciphertext1, ciphertext2)
			params.bfvContext.contextT.Sub(values1, values2, values1)

			verifyTestVectors(params, params.decryptor, values1, ciphertext1, t)
		})

		t.Run(testString("CtPlain/", parameters), func(t *testing.T) {

			values1, _, ciphertext1 := newTestVectors(params, params.encryptorPk, t)
			values2, plaintext2, ciphertext2 := newTestVectors(params, params.encryptorPk, t)

			valuesWant := params.bfvContext.contextT.NewPoly()

			params.evaluator.Sub(ciphertext1, plaintext2, ciphertext2)
			params.bfvContext.contextT.Sub(values1, values2, valuesWant)
			verifyTestVectors(params, params.decryptor, valuesWant, ciphertext2, t)

			params.evaluator.Sub(plaintext2, ciphertext1, ciphertext2)
			params.bfvContext.contextT.Sub(values2, values1, valuesWant)
			verifyTestVectors(params, params.decryptor, valuesWant, ciphertext2, t)
		})
	}
}

func testEvaluatorMul(t *testing.T) {

	for _, parameters := range testParams.bfvParameters {

		params := genBfvParams(parameters)

		rlk := params.kgen.NewRelinKey(params.sk, 1)

		t.Run(testString("CtCt/", parameters), func(t *testing.T) {

			values1, _, ciphertext1 := newTestVectors(params, params.encryptorPk, t)
			values2, _, ciphertext2 := newTestVectors(params, params.encryptorPk, t)

			receiver := NewCiphertext(parameters, ciphertext1.Degree()+ciphertext2.Degree())
			params.evaluator.Mul(ciphertext1, ciphertext2, receiver)
			params.bfvContext.contextT.MulCoeffs(values1, values2, values1)

			verifyTestVectors(params, params.decryptor, values1, receiver, t)
		})

		t.Run(testString("CtPlain/", parameters), func(t *testing.T) {

			values1, _, ciphertext1 := newTestVectors(params, params.encryptorPk, t)
			values2, plaintext2, _ := newTestVectors(params, params.encryptorPk, t)

			params.evaluator.Mul(ciphertext1, plaintext2, ciphertext1)
			params.bfvContext.contextT.MulCoeffs(values1, values2, values1)

			verifyTestVectors(params, params.decryptor, values1, ciphertext1, t)
		})

		t.Run(testString("Relinearize/", parameters), func(t *testing.T) {

			values1, _, ciphertext1 := newTestVectors(params, params.encryptorPk, t)
			values2, _, ciphertext2 := newTestVectors(params, params.encryptorPk, t)

			receiver := NewCiphertext(parameters, ciphertext1.Degree()+ciphertext2.Degree())
			params.evaluator.Mul(ciphertext1, ciphertext2, receiver)
			params.bfvContext.contextT.MulCoeffs(values1, values2, values1)

			receiver2 := params.evaluator.RelinearizeNew(receiver, rlk)
			verifyTestVectors(params, params.decryptor, values1, receiver2, t)

			params.evaluator.Relinearize(receiver, rlk, receiver)
			verifyTestVectors(params, params.decryptor, values1, receiver, t)
		})
	}
}

func testKeySwitch(t *testing.T) {

	for _, parameters := range testParams.bfvParameters {

		params := genBfvParams(parameters)

		sk2 := params.kgen.NewSecretKey()
		decryptorSk2 := NewDecryptor(parameters, sk2)
		switchKey := params.kgen.NewSwitchingKey(params.sk, sk2)

		t.Run(testString("InPlace/", parameters), func(t *testing.T) {

			values, _, ciphertext := newTestVectors(params, params.encryptorPk, t)

			params.evaluator.SwitchKeys(ciphertext, switchKey, ciphertext)

			verifyTestVectors(params, decryptorSk2, values, ciphertext, t)
		})

		t.Run(testString("New/", parameters), func(t *testing.T) {

			values, _, ciphertext := newTestVectors(params, params.encryptorPk, t)

			ciphertext = params.evaluator.SwitchKeysNew(ciphertext, switchKey)
			verifyTestVectors(params, decryptorSk2, values, ciphertext, t)
		})
	}
}

func testRotateRows(t *testing.T) {

	for _, parameters := range testParams.bfvParameters {

		params := genBfvParams(parameters)

		rotkey := NewRotationKeys()
		params.kgen.GenRot(RotationRow, params.sk, 0, rotkey)

		t.Run(testString("InPlace/", parameters), func(t *testing.T) {

			values, _, ciphertext := newTestVectors(params, params.encryptorPk, t)

			params.evaluator.RotateRows(ciphertext, rotkey, ciphertext)

			values.Coeffs[0] = append(values.Coeffs[0][params.bfvContext.n>>1:], values.Coeffs[0][:params.bfvContext.n>>1]...)

			verifyTestVectors(params, params.decryptor, values, ciphertext, t)
		})

		t.Run(testString("New/", parameters), func(t *testing.T) {

			values, _, ciphertext := newTestVectors(params, params.encryptorPk, t)

			ciphertext = params.evaluator.RotateRowsNew(ciphertext, rotkey)

			values.Coeffs[0] = append(values.Coeffs[0][params.bfvContext.n>>1:], values.Coeffs[0][:params.bfvContext.n>>1]...)

			verifyTestVectors(params, params.decryptor, values, ciphertext, t)
		})
	}
}

func testRotateCols(t *testing.T) {

	for _, parameters := range testParams.bfvParameters {

		params := genBfvParams(parameters)

		rotkey := params.kgen.NewRotationKeysPow2(params.sk)

		valuesWant := params.bfvContext.contextT.NewPoly()
		mask := (params.bfvContext.n >> 1) - 1
		slots := params.bfvContext.n >> 1

		t.Run(testString("InPlace/", parameters), func(t *testing.T) {

			values, _, ciphertext := newTestVectors(params, params.encryptorPk, t)

			receiver := NewCiphertext(parameters, 1)
			for n := uint64(1); n < slots; n <<= 1 {

				params.evaluator.RotateColumns(ciphertext, n, rotkey, receiver)

				for i := uint64(0); i < slots; i++ {
					valuesWant.Coeffs[0][i] = values.Coeffs[0][(i+n)&mask]
					valuesWant.Coeffs[0][i+slots] = values.Coeffs[0][((i+n)&mask)+slots]
				}

				verifyTestVectors(params, params.decryptor, valuesWant, receiver, t)
			}
		})

		t.Run(testString("New/", parameters), func(t *testing.T) {

			values, _, ciphertext := newTestVectors(params, params.encryptorPk, t)

			for n := uint64(1); n < slots; n <<= 1 {

				receiver := params.evaluator.RotateColumnsNew(ciphertext, n, rotkey)

				for i := uint64(0); i < slots; i++ {
					valuesWant.Coeffs[0][i] = values.Coeffs[0][(i+n)&mask]
					valuesWant.Coeffs[0][i+slots] = values.Coeffs[0][((i+n)&mask)+slots]
				}

				verifyTestVectors(params, params.decryptor, valuesWant, receiver, t)
			}
		})

		t.Run(testString("Random/", parameters), func(t *testing.T) {

			values, _, ciphertext := newTestVectors(params, params.encryptorPk, t)

			receiver := NewCiphertext(parameters, 1)
			for n := 0; n < 4; n++ {

				rand := ring.RandUniform(slots, mask)

				params.evaluator.RotateColumns(ciphertext, rand, rotkey, receiver)

				for i := uint64(0); i < slots; i++ {
					valuesWant.Coeffs[0][i] = values.Coeffs[0][(i+rand)&mask]
					valuesWant.Coeffs[0][i+slots] = values.Coeffs[0][((i+rand)&mask)+slots]
				}

				verifyTestVectors(params, params.decryptor, valuesWant, receiver, t)
			}
		})
	}
}
