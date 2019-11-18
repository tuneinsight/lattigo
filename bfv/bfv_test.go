package bfv

import (
	"fmt"
	"github.com/ldsec/lattigo/ring"
	"log"
	"math/rand"
	"testing"
	"time"
)

func check(t *testing.T, err error) {
	if err != nil {
		t.Error(err)
		log.Fatal(err)
	}
}

func testString2(opname string, params *bfvParams) string {
	return fmt.Sprintf("%sparams=%d", opname, params.bfvContext.N())
}

type bfvParams struct {
	bfvContext  *BfvContext
	encoder     *Encoder
	kgen        *KeyGenerator
	sk          *SecretKey
	pk          *PublicKey
	encryptorPk *Encryptor
	encryptorSk *Encryptor
	decryptor   *Decryptor
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
		&DefaultParams[0],
		&DefaultParams[1],
		&DefaultParams[2],
	}
}

func Test_BFV(t *testing.T) {
	t.Run("Encoder", testEncoder)
	t.Run("Encryptor", testEncryptor)
	t.Run("Evaluator/Add", testEvaluatorAdd)
	t.Run("Evaluator/Sub", testEvaluatorSub)
	t.Run("Evaluator/Mul", testEvaluatorMul)
	t.Run("Evaluator/KeySwitch", testKeySwitch)
	t.Run("Evaluator/RotateRows", testRotateRows)
	t.Run("Evaluator/RotateCols", testRotateCols)

	//tests for marshalling.
	t.Run("Marshalling", testMarshaller)
}

func testMarshaller(t *testing.T) {

	p := testParams.bfvParameters
	params := genBfvParams(p[0])
	//check public key
	pk := new(PublicKey)
	marshalledPk, err := params.pk.MarshalBinary()
	check(t, err)
	err = pk.UnmarshalBinary(marshalledPk)
	check(t, err)

	//check secret keys.
	sk := new(SecretKey)
	marshalledSk, err := params.sk.MarshalBinary()
	check(t, err)
	err = sk.UnmarshalBinary(marshalledSk)
	check(t, err)

	keygen := params.bfvContext.NewKeyGenerator()
	sk_out := keygen.NewSecretKey()
	//Check switching key
	switching_key := keygen.NewSwitchingKey(sk, sk_out)
	data, err := switching_key.MarshalBinary()
	check(t, err)
	res_switching_key := new(SwitchingKey)
	err = res_switching_key.UnmarshalBinary(data)
	check(t, err)
	//check eval key
	eval_key := EvaluationKey{evakey: []*SwitchingKey{switching_key}}
	data, err = eval_key.MarshalBinary()
	check(t, err)
	res_eval_key := new(EvaluationKey)
	err = res_eval_key.UnmarshalBinary(data)
	check(t, err)

}

func genBfvParams(contextParameters *Parameters) (params *bfvParams) {

	params = new(bfvParams)

	if params.bfvContext, err = NewBfvContextWithParam(contextParameters); err != nil {
		log.Fatal(err)
	}

	params.kgen = params.bfvContext.NewKeyGenerator()

	params.sk, params.pk = params.kgen.NewKeyPair()

	params.encoder = params.bfvContext.NewEncoder()

	params.encryptorPk = params.bfvContext.NewEncryptorFromPk(params.pk)
	params.encryptorSk = params.bfvContext.NewEncryptorFromSk(params.sk)
	params.decryptor = params.bfvContext.NewDecryptor(params.sk)

	params.evaluator = params.bfvContext.NewEvaluator()

	return

}

func new_test_vectors(params *bfvParams, encryptor *Encryptor, t *testing.T) (coeffs *ring.Poly, plaintext *Plaintext, ciphertext *Ciphertext) {

	coeffs = params.bfvContext.contextT.NewUniformPoly()

	plaintext = params.bfvContext.NewPlaintext()

	params.encoder.EncodeUint(coeffs.Coeffs[0], plaintext)

	if encryptor != nil {
		ciphertext = params.encryptorPk.EncryptNew(plaintext)
	}

	return coeffs, plaintext, ciphertext
}

func verify_test_vectors(params *bfvParams, decryptor *Decryptor, coeffs *ring.Poly, element Operand, t *testing.T) {

	var coeffsTest []uint64

	el := element.Element()

	if el.Degree() == 0 {

		coeffsTest = params.encoder.DecodeUint(el.Plaintext())

	} else {

		coeffsTest = params.encoder.DecodeUint(decryptor.DecryptNew(el.Ciphertext()))
	}

	if EqualSlice(coeffs.Coeffs[0], coeffsTest) != true {
		t.Errorf("decryption error")
	}
}

func testEncoder(t *testing.T) {

	for _, parameters := range testParams.bfvParameters {

		params := genBfvParams(parameters)

		t.Run(testString2("Encode&Decode/", params), func(t *testing.T) {

			values, plaintext, _ := new_test_vectors(params, nil, t)

			verify_test_vectors(params, params.decryptor, values, plaintext, t)
		})
	}
}

func testEncryptor(t *testing.T) {

	for _, parameters := range testParams.bfvParameters {

		params := genBfvParams(parameters)

		t.Run(testString2("EncryptFromPk/", params), func(t *testing.T) {

			values, _, ciphertext := new_test_vectors(params, params.encryptorPk, t)

			verify_test_vectors(params, params.decryptor, values, ciphertext, t)
		})

		t.Run(testString2("EncryptFromSk/", params), func(t *testing.T) {

			values, _, ciphertext := new_test_vectors(params, params.encryptorSk, t)

			verify_test_vectors(params, params.decryptor, values, ciphertext, t)
		})
	}
}

func testEvaluatorAdd(t *testing.T) {

	for _, parameters := range testParams.bfvParameters {

		params := genBfvParams(parameters)

		t.Run(testString2("CtCtInPlace/", params), func(t *testing.T) {

			values1, _, ciphertext1 := new_test_vectors(params, params.encryptorPk, t)
			values2, _, ciphertext2 := new_test_vectors(params, params.encryptorPk, t)

			params.evaluator.Add(ciphertext1, ciphertext2, ciphertext1)
			params.bfvContext.contextT.Add(values1, values2, values1)

			verify_test_vectors(params, params.decryptor, values1, ciphertext1, t)
		})

		t.Run(testString2("CtCtNew/", params), func(t *testing.T) {

			values1, _, ciphertext1 := new_test_vectors(params, params.encryptorPk, t)
			values2, _, ciphertext2 := new_test_vectors(params, params.encryptorPk, t)

			ciphertext1 = params.evaluator.AddNew(ciphertext1, ciphertext2)
			params.bfvContext.contextT.Add(values1, values2, values1)

			verify_test_vectors(params, params.decryptor, values1, ciphertext1, t)
		})

		t.Run(testString2("CtPlain/", params), func(t *testing.T) {

			values1, _, ciphertext1 := new_test_vectors(params, params.encryptorPk, t)
			values2, plaintext2, ciphertext2 := new_test_vectors(params, params.encryptorPk, t)

			params.evaluator.Add(ciphertext1, plaintext2, ciphertext2)
			params.bfvContext.contextT.Add(values1, values2, values2)

			verify_test_vectors(params, params.decryptor, values2, ciphertext2, t)

			params.evaluator.Add(plaintext2, ciphertext1, ciphertext2)

			verify_test_vectors(params, params.decryptor, values2, ciphertext2, t)
		})
	}
}

func testEvaluatorSub(t *testing.T) {

	for _, parameters := range testParams.bfvParameters {

		params := genBfvParams(parameters)

		t.Run(testString2("CtCtInPlace/", params), func(t *testing.T) {

			values1, _, ciphertext1 := new_test_vectors(params, params.encryptorPk, t)
			values2, _, ciphertext2 := new_test_vectors(params, params.encryptorPk, t)

			params.evaluator.Sub(ciphertext1, ciphertext2, ciphertext1)
			params.bfvContext.contextT.Sub(values1, values2, values1)

			verify_test_vectors(params, params.decryptor, values1, ciphertext1, t)
		})

		t.Run(testString2("CtCtNew/", params), func(t *testing.T) {

			values1, _, ciphertext1 := new_test_vectors(params, params.encryptorPk, t)
			values2, _, ciphertext2 := new_test_vectors(params, params.encryptorPk, t)

			ciphertext1 = params.evaluator.SubNew(ciphertext1, ciphertext2)
			params.bfvContext.contextT.Sub(values1, values2, values1)

			verify_test_vectors(params, params.decryptor, values1, ciphertext1, t)
		})

		t.Run(testString2("CtPlain/", params), func(t *testing.T) {

			values1, _, ciphertext1 := new_test_vectors(params, params.encryptorPk, t)
			values2, plaintext2, ciphertext2 := new_test_vectors(params, params.encryptorPk, t)

			valuesWant := params.bfvContext.contextT.NewPoly()

			params.evaluator.Sub(ciphertext1, plaintext2, ciphertext2)
			params.bfvContext.contextT.Sub(values1, values2, valuesWant)
			verify_test_vectors(params, params.decryptor, valuesWant, ciphertext2, t)

			params.evaluator.Sub(plaintext2, ciphertext1, ciphertext2)
			params.bfvContext.contextT.Sub(values2, values1, valuesWant)
			verify_test_vectors(params, params.decryptor, valuesWant, ciphertext2, t)
		})
	}
}

func testEvaluatorMul(t *testing.T) {

	for _, parameters := range testParams.bfvParameters {

		params := genBfvParams(parameters)

		rlk := params.kgen.NewRelinKey(params.sk, 1)

		t.Run(testString2("CtCt/", params), func(t *testing.T) {

			values1, _, ciphertext1 := new_test_vectors(params, params.encryptorPk, t)
			values2, _, ciphertext2 := new_test_vectors(params, params.encryptorPk, t)

			receiver := params.bfvContext.NewCiphertext(ciphertext1.Degree() + ciphertext2.Degree())
			params.evaluator.Mul(ciphertext1, ciphertext2, receiver)
			params.bfvContext.contextT.MulCoeffs(values1, values2, values1)

			verify_test_vectors(params, params.decryptor, values1, receiver, t)
		})

		t.Run(testString2("CtPlain/", params), func(t *testing.T) {

			values1, _, ciphertext1 := new_test_vectors(params, params.encryptorPk, t)
			values2, plaintext2, _ := new_test_vectors(params, params.encryptorPk, t)

			params.evaluator.Mul(ciphertext1, plaintext2, ciphertext1)
			params.bfvContext.contextT.MulCoeffs(values1, values2, values1)

			verify_test_vectors(params, params.decryptor, values1, ciphertext1, t)
		})

		t.Run(testString2("Relinearize/", params), func(t *testing.T) {

			values1, _, ciphertext1 := new_test_vectors(params, params.encryptorPk, t)
			values2, _, ciphertext2 := new_test_vectors(params, params.encryptorPk, t)

			receiver := params.bfvContext.NewCiphertext(ciphertext1.Degree() + ciphertext2.Degree())
			params.evaluator.Mul(ciphertext1, ciphertext2, receiver)
			params.bfvContext.contextT.MulCoeffs(values1, values2, values1)

			receiver2, err := params.evaluator.RelinearizeNew(receiver, rlk)
			check(t, err)
			verify_test_vectors(params, params.decryptor, values1, receiver2, t)

			check(t, params.evaluator.Relinearize(receiver, rlk, receiver))
			verify_test_vectors(params, params.decryptor, values1, receiver, t)
		})
	}
}

func testKeySwitch(t *testing.T) {

	for _, parameters := range testParams.bfvParameters {

		params := genBfvParams(parameters)

		sk2 := params.kgen.NewSecretKey()
		decryptorSk2 := params.bfvContext.NewDecryptor(sk2)
		switchKey := params.kgen.NewSwitchingKey(params.sk, sk2)

		t.Run(testString2("InPlace/", params), func(t *testing.T) {

			values, _, ciphertext := new_test_vectors(params, params.encryptorPk, t)

			params.evaluator.SwitchKeys(ciphertext, switchKey, ciphertext)

			verify_test_vectors(params, decryptorSk2, values, ciphertext, t)
		})

		t.Run(testString2("New/", params), func(t *testing.T) {

			values, _, ciphertext := new_test_vectors(params, params.encryptorPk, t)

			ciphertext, err = params.evaluator.SwitchKeysNew(ciphertext, switchKey)
			check(t, err)

			verify_test_vectors(params, decryptorSk2, values, ciphertext, t)
		})
	}
}

func testRotateRows(t *testing.T) {

	for _, parameters := range testParams.bfvParameters {

		params := genBfvParams(parameters)

		rotkey := params.bfvContext.NewRotationKeys()
		params.kgen.GenRot(RotationRow, params.sk, 0, rotkey)

		t.Run(testString2("InPlace/", params), func(t *testing.T) {

			values, _, ciphertext := new_test_vectors(params, params.encryptorPk, t)

			check(t, params.evaluator.RotateRows(ciphertext, rotkey, ciphertext))

			values.Coeffs[0] = append(values.Coeffs[0][params.bfvContext.n>>1:], values.Coeffs[0][:params.bfvContext.n>>1]...)

			verify_test_vectors(params, params.decryptor, values, ciphertext, t)
		})

		t.Run(testString2("New/", params), func(t *testing.T) {

			values, _, ciphertext := new_test_vectors(params, params.encryptorPk, t)

			ciphertext, err = params.evaluator.RotateRowsNew(ciphertext, rotkey)
			check(t, err)

			values.Coeffs[0] = append(values.Coeffs[0][params.bfvContext.n>>1:], values.Coeffs[0][:params.bfvContext.n>>1]...)

			verify_test_vectors(params, params.decryptor, values, ciphertext, t)
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

		t.Run(testString2("InPlace/", params), func(t *testing.T) {

			values, _, ciphertext := new_test_vectors(params, params.encryptorPk, t)

			receiver := params.bfvContext.NewCiphertext(1)
			for n := uint64(1); n < slots; n <<= 1 {

				check(t, params.evaluator.RotateColumns(ciphertext, n, rotkey, receiver))

				for i := uint64(0); i < slots; i++ {
					valuesWant.Coeffs[0][i] = values.Coeffs[0][(i+n)&mask]
					valuesWant.Coeffs[0][i+slots] = values.Coeffs[0][((i+n)&mask)+slots]
				}

				verify_test_vectors(params, params.decryptor, valuesWant, receiver, t)
			}
		})

		t.Run(testString2("New/", params), func(t *testing.T) {

			values, _, ciphertext := new_test_vectors(params, params.encryptorPk, t)

			for n := uint64(1); n < slots; n <<= 1 {

				receiver, err := params.evaluator.RotateColumnsNew(ciphertext, n, rotkey)
				check(t, err)

				for i := uint64(0); i < slots; i++ {
					valuesWant.Coeffs[0][i] = values.Coeffs[0][(i+n)&mask]
					valuesWant.Coeffs[0][i+slots] = values.Coeffs[0][((i+n)&mask)+slots]
				}

				verify_test_vectors(params, params.decryptor, valuesWant, receiver, t)
			}
		})

		t.Run(testString2("Random/", params), func(t *testing.T) {

			values, _, ciphertext := new_test_vectors(params, params.encryptorPk, t)

			receiver := params.bfvContext.NewCiphertext(1)
			for n := 0; n < 4; n++ {

				rand := ring.RandUniform(slots, mask)

				check(t, params.evaluator.RotateColumns(ciphertext, rand, rotkey, receiver))

				for i := uint64(0); i < slots; i++ {
					valuesWant.Coeffs[0][i] = values.Coeffs[0][(i+rand)&mask]
					valuesWant.Coeffs[0][i+slots] = values.Coeffs[0][((i+rand)&mask)+slots]
				}

				verify_test_vectors(params, params.decryptor, valuesWant, receiver, t)
			}
		})
	}
}
