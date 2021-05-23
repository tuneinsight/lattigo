package mkckks

import (
	"flag"
	"math"
	"sort"
	"testing"

	"github.com/ldsec/lattigo/v2/ckks"
	"github.com/ldsec/lattigo/v2/mkrlwe"
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/utils"
	"github.com/stretchr/testify/require"
)

var printPrecisionStats = flag.Bool("print-precision", false, "print precision stats")
var minPrec = 13.0

func Test_MKCKKS(t *testing.T) {

	//skip parameter 4 due to memory consumption
	for i, paramLit := range ckks.DefaultParams {

		params, err := ckks.NewParametersFromLiteral(paramLit)
		if err != nil {
			panic(err)
		}
		p := &params

		if i != 4 && i != 9 && i != 0 {

			testEncryptionEqualsDecryption(t, p)
			testEncryptionEqualsDecryptionWithSecretKey(t, p)
			testAdd(t, p)
			testAddManyTimeSameCipher(t, p)
			testAddFourParticipants(t, p)
			testAddPlaintext(t, p)
			testAddPlaintextTwoParticipants(t, p)
			testSub(t, p)
			testSubPlaintext(t, p)
			testNeg(t, p)
			testSubPlaintextTwoParticipants(t, p)
			testMulPlaintext(t, p)
			testMulPlaintextTwoParticipants(t, p)
			testTensor2(t, p)

			testTensor(t, p)

			//testTensorTwoParticipants(t, p)

			testCkksMkbfvBridge(t, p)
			testMarshaler(t, p)
			testRotation(t, p)
			testRotationTwoParticipants(t, p)
			testSquare(t, p)

			//testMul(t, p)
			//testMulAfterAdd(t, p)
			if i != 5 && i != 6 {
				//testAddAfterMul(t, p)
				//testMulFourParticipants(t, p)
			}

		}
	}

}

func testEncryptionEqualsDecryption(t *testing.T, params *ckks.Parameters) {

	sigma := 6.0

	participants := setupPeers(1, params, sigma)

	t.Run(testString("Test encryption equals decryption/", 1, params), func(t *testing.T) {

		// get random value
		value := newTestValue(params, complex(-1, -1), complex(1, 1))

		//encrypt
		cipher := participants[0].Encrypt(value)

		// decrypt
		partialDec := participants[0].GetPartialDecryption(cipher)
		decrypted := participants[0].Decrypt(cipher, []*ring.Poly{partialDec})

		// decode and check
		verifyTestVectors(params, value, decrypted, t)
	})

}

func testEncryptionEqualsDecryptionWithSecretKey(t *testing.T, params *ckks.Parameters) {

	sigma := 6.0

	participants := setupPeers(1, params, sigma)

	t.Run(testString("Test encryption equals decryption with secret key/", 1, params), func(t *testing.T) {

		// get random value
		value := newTestValue(params, complex(-1, -1), complex(1, 1))

		//encrypt
		cipher := participants[0].Encrypt(value)

		// decrypt
		decrypted := DecryptSimple(participants[0].GetSecretKey(), cipher, params)

		// decode and check
		verifyTestVectors(params, value, decrypted, t)
	})

}

func testAdd(t *testing.T, params *ckks.Parameters) {

	sigma := 6.0

	participants := setupPeers(2, params, sigma)

	t.Run(testString("Test add/", 2, params), func(t *testing.T) {

		// generate new values
		values1 := newTestValue(params, complex(-1, -1), complex(1, 1))
		values2 := newTestValue(params, complex(-1, -1), complex(1, 1))

		// Encrypt
		cipher1 := participants[0].Encrypt(values1)
		cipher2 := participants[1].Encrypt(values2)

		// pad and add
		evaluator := NewMKEvaluator(params)
		resCipher := evaluator.Add(cipher1, cipher2)

		// decrypt
		partialDec1 := participants[0].GetPartialDecryption(resCipher)
		partialDec2 := participants[1].GetPartialDecryption(resCipher)

		decrypted := participants[0].Decrypt(resCipher, []*ring.Poly{partialDec1, partialDec2})

		// perform the operation in the plaintext space
		for i := 0; i < len(values1); i++ {
			values1[i] += values2[i]
		}

		// check results
		verifyTestVectors(params, values1, decrypted, t)
	})

}

func testAddManyTimeSameCipher(t *testing.T, params *ckks.Parameters) {

	sigma := 6.0

	participants := setupPeers(2, params, sigma)

	t.Run(testString("Test add many time same cipher/", 2, params), func(t *testing.T) {

		// generate new values
		values1 := newTestValue(params, complex(-1, -1), complex(1, 1))
		values2 := newTestValue(params, complex(-1, -1), complex(1, 1))

		// Encrypt
		cipher1 := participants[0].Encrypt(values1)
		cipher2 := participants[1].Encrypt(values2)

		// pad and add
		evaluator := NewMKEvaluator(params)

		resCipher := evaluator.Add(cipher1, cipher2)
		resCipher = evaluator.Add(resCipher, cipher2)
		resCipher = evaluator.Add(resCipher, cipher2)
		resCipher = evaluator.Add(resCipher, cipher2)

		// decrypt
		partialDec1 := participants[0].GetPartialDecryption(resCipher)
		partialDec2 := participants[1].GetPartialDecryption(resCipher)

		decrypted := participants[0].Decrypt(resCipher, []*ring.Poly{partialDec1, partialDec2})

		// perform the operation in the plaintext space
		for i := 0; i < len(values1); i++ {
			values1[i] += 4 * values2[i]
		}

		// check results
		verifyTestVectors(params, values1, decrypted, t)
	})

}

func testSub(t *testing.T, params *ckks.Parameters) {

	sigma := 6.0

	participants := setupPeers(2, params, sigma)

	t.Run(testString("Test Subtraction/", 2, params), func(t *testing.T) {

		value1 := newTestValue(params, complex(-1, -1), complex(1, 1))
		value2 := newTestValue(params, complex(-1, -1), complex(1, 1))

		// encrypt
		cipher1 := participants[0].Encrypt(value1)
		cipher2 := participants[1].Encrypt(value2)

		// pad and add
		evaluator := NewMKEvaluator(params)

		resCipher := evaluator.Sub(cipher1, cipher2)

		// decrypt
		partialDec1 := participants[0].GetPartialDecryption(resCipher)
		partialDec2 := participants[1].GetPartialDecryption(resCipher)

		decrypted := participants[0].Decrypt(resCipher, []*ring.Poly{partialDec1, partialDec2})

		// perform the operation in the plaintext space
		// perform the operation in the plaintext space
		for i := 0; i < len(value1); i++ {
			value1[i] -= value2[i]
		}

		// check results
		verifyTestVectors(params, value1, decrypted, t)

	})

}

func testAddPlaintext(t *testing.T, params *ckks.Parameters) {

	sigma := 6.0

	participants := setupPeers(1, params, sigma)

	t.Run(testString("Test Plaintext Addition/", 1, params), func(t *testing.T) {

		value1 := newTestValue(params, complex(-1, -1), complex(1, 1))

		// encrypt
		cipher := participants[0].Encrypt(value1)

		// add plaintext to ciphertext
		evaluator := NewMKEvaluator(params)
		value2 := newTestValue(params, complex(-1, -1), complex(1, 1))

		pt := evaluator.NewPlaintextFromValue(value2)
		resCipher := evaluator.AddPlaintext(pt, cipher)

		// decrypt
		partialDec1 := participants[0].GetPartialDecryption(resCipher)

		decrypted := participants[0].Decrypt(resCipher, []*ring.Poly{partialDec1})

		// perform the operation in the plaintext space
		for i := 0; i < len(value1); i++ {
			value1[i] += value2[i]
		}

		// check results
		verifyTestVectors(params, value1, decrypted, t)
	})

}

func testAddPlaintextTwoParticipants(t *testing.T, params *ckks.Parameters) {

	sigma := 6.0

	participants := setupPeers(2, params, sigma)

	t.Run(testString("Test Plaintext Addition/", 2, params), func(t *testing.T) {

		value1 := newTestValue(params, complex(-1, -1), complex(1, 1))
		value2 := newTestValue(params, complex(-1, -1), complex(1, 1))

		// encrypt
		cipher1 := participants[0].Encrypt(value1)
		cipher2 := participants[1].Encrypt(value2)

		// add plaintext to one of the ciphertext then add both ciphertexts
		evaluator := NewMKEvaluator(params)
		value3 := newTestValue(params, complex(-1, -1), complex(1, 1))

		pt := evaluator.NewPlaintextFromValue(value3)
		resCipher1 := evaluator.AddPlaintext(pt, cipher1)

		resCipher := evaluator.Add(resCipher1, cipher2)

		// decrypt
		partialDec1 := participants[0].GetPartialDecryption(resCipher)
		partialDec2 := participants[1].GetPartialDecryption(resCipher)

		decrypted := participants[0].Decrypt(resCipher, []*ring.Poly{partialDec1, partialDec2})

		// perform the operation in the plaintext space
		for i := 0; i < len(value1); i++ {
			value1[i] += value2[i] + value3[i]
		}

		// check results
		verifyTestVectors(params, value1, decrypted, t)

	})

}

func testSubPlaintext(t *testing.T, params *ckks.Parameters) {

	sigma := 6.0

	participants := setupPeers(1, params, sigma)

	t.Run(testString("Test Plaintext Subtraction/", 1, params), func(t *testing.T) {

		value1 := newTestValue(params, complex(-1, -1), complex(1, 1))

		// encrypt
		cipher := participants[0].Encrypt(value1)

		// sub plaintext to ciphertext
		evaluator := NewMKEvaluator(params)
		value2 := newTestValue(params, complex(-1, -1), complex(1, 1))

		pt := evaluator.NewPlaintextFromValue(value2)
		resCipher := evaluator.SubPlaintext(pt, cipher)

		// decrypt
		partialDec1 := participants[0].GetPartialDecryption(resCipher)

		decrypted := participants[0].Decrypt(resCipher, []*ring.Poly{partialDec1})

		// perform the operation in the plaintext space
		for i := 0; i < len(value1); i++ {
			value1[i] -= value2[i]
		}

		// check results
		verifyTestVectors(params, value1, decrypted, t)

	})

}

func testSubPlaintextTwoParticipants(t *testing.T, params *ckks.Parameters) {

	sigma := 6.0

	participants := setupPeers(2, params, sigma)

	t.Run(testString("Test Plaintext Subtraction/", 2, params), func(t *testing.T) {

		value1 := newTestValue(params, complex(-1, -1), complex(1, 1))
		value2 := newTestValue(params, complex(-1, -1), complex(1, 1))

		// encrypt
		cipher1 := participants[0].Encrypt(value1)
		cipher2 := participants[1].Encrypt(value2)

		// sub plaintext to one of the ciphertext then sub both ciphertexts
		evaluator := NewMKEvaluator(params)
		value3 := newTestValue(params, complex(-1, -1), complex(1, 1))

		pt := evaluator.NewPlaintextFromValue(value3)
		resCipher1 := evaluator.SubPlaintext(pt, cipher1)

		resCipher := evaluator.Sub(resCipher1, cipher2)

		// decrypt
		partialDec1 := participants[0].GetPartialDecryption(resCipher)
		partialDec2 := participants[1].GetPartialDecryption(resCipher)

		decrypted := participants[0].Decrypt(resCipher, []*ring.Poly{partialDec1, partialDec2})

		// perform the operation in the plaintext space
		for i := 0; i < len(value1); i++ {
			value1[i] -= value2[i] + value3[i]
		}

		// check results
		verifyTestVectors(params, value1, decrypted, t)

	})

}

func testNeg(t *testing.T, params *ckks.Parameters) {

	sigma := 6.0

	participants := setupPeers(1, params, sigma)

	t.Run(testString("Test Negation/", 1, params), func(t *testing.T) {

		value1 := newTestValue(params, complex(-1, -1), complex(1, 1))

		// encrypt
		cipher := participants[0].Encrypt(value1)

		// add with negated ciphertext
		evaluator := NewMKEvaluator(params)

		resCipher := evaluator.Add(evaluator.Neg(cipher), cipher)

		// decrypt
		partialDec1 := participants[0].GetPartialDecryption(resCipher)

		decrypted := participants[0].Decrypt(resCipher, []*ring.Poly{partialDec1})

		// perform the operation in the plaintext space
		for i := 0; i < len(value1); i++ {
			value1[i] -= value1[i]
		}

		// check results
		verifyTestVectors(params, value1, decrypted, t)

	})
}

func testAddFourParticipants(t *testing.T, params *ckks.Parameters) {

	sigma := 6.0

	participants := setupPeers(4, params, sigma)

	t.Run(testString("Test add/", 4, params), func(t *testing.T) {

		// generate test values
		values1 := newTestValue(params, complex(-1, -1), complex(1, 1))
		values2 := newTestValue(params, complex(-1, -1), complex(1, 1))
		values3 := newTestValue(params, complex(-1, -1), complex(1, 1))
		values4 := newTestValue(params, complex(-1, -1), complex(1, 1))

		// Encrypt
		cipher1 := participants[0].Encrypt(values1)
		cipher2 := participants[1].Encrypt(values2)
		cipher3 := participants[2].Encrypt(values3)
		cipher4 := participants[3].Encrypt(values4)

		// pad and add in 2 steps
		evaluator := NewMKEvaluator(params)

		resCipher1 := evaluator.Add(cipher1, cipher2)
		resCipher2 := evaluator.Add(cipher3, cipher4)

		resCipher := evaluator.Add(resCipher1, resCipher2)

		// decrypt
		partialDec1 := participants[0].GetPartialDecryption(resCipher)
		partialDec2 := participants[1].GetPartialDecryption(resCipher)
		partialDec3 := participants[2].GetPartialDecryption(resCipher)
		partialDec4 := participants[3].GetPartialDecryption(resCipher)

		decrypted := participants[0].Decrypt(resCipher, []*ring.Poly{partialDec1, partialDec2, partialDec3, partialDec4})

		// perform the operation in the plaintext space
		for i := 0; i < len(values1); i++ {
			values1[i] += values2[i] + values3[i] + values4[i]
		}

		// check results
		verifyTestVectors(params, values1, decrypted, t)
	})

}

func testMulPlaintext(t *testing.T, params *ckks.Parameters) {

	sigma := 6.0

	participants := setupPeers(1, params, sigma)

	t.Run(testString("Test Plaintext Multiplication/", 1, params), func(t *testing.T) {

		value1 := newTestValue(params, complex(-1, -1), complex(1, 1))

		// encrypt
		cipher := participants[0].Encrypt(value1)

		// multiply plaintext and ciphertext
		evaluator := NewMKEvaluator(params)
		value2 := newTestValue(params, complex(-1, -1), complex(1, 1))

		pt := evaluator.NewPlaintextFromValue(value2)
		resCipher := evaluator.MultPlaintext(pt, cipher)

		// decrypt
		partialDec1 := participants[0].GetPartialDecryption(resCipher)

		decrypted := participants[0].Decrypt(resCipher, []*ring.Poly{partialDec1})

		// perform the operation in the plaintext space
		for i := 0; i < len(value1); i++ {
			value1[i] *= value2[i]
		}

		// check results
		verifyTestVectors(params, value1, decrypted, t)

	})

}

func testMulPlaintextTwoParticipants(t *testing.T, params *ckks.Parameters) {

	sigma := 6.0

	participants := setupPeers(2, params, sigma)

	t.Run(testString("Test Plaintext Multiplication/", 2, params), func(t *testing.T) {

		value1 := newTestValue(params, complex(-1, -1), complex(1, 1))
		value2 := newTestValue(params, complex(-1, -1), complex(1, 1))

		// encrypt
		cipher1 := participants[0].Encrypt(value1)
		cipher2 := participants[1].Encrypt(value2)

		// add both ciphertexts then multiply by plaintext
		evaluator := NewMKEvaluator(params)
		value3 := newTestValue(params, complex(-1, -1), complex(1, 1))

		resCipherTMP := evaluator.Add(cipher1, cipher2)

		pt := evaluator.NewPlaintextFromValue(value3)
		resCipher := evaluator.MultPlaintext(pt, resCipherTMP)

		// decrypt
		partialDec1 := participants[0].GetPartialDecryption(resCipher)
		partialDec2 := participants[1].GetPartialDecryption(resCipher)

		decrypted := participants[0].Decrypt(resCipher, []*ring.Poly{partialDec1, partialDec2})

		// perform the operation in the plaintext space
		for i := 0; i < len(value1); i++ {
			value1[i] = (value1[i] + value2[i]) * value3[i]
		}

		// check results
		verifyTestVectors(params, value1, decrypted, t)

	})

}

func testTensorTwoParticipants(t *testing.T, params *ckks.Parameters) {
	sigma := 6.0

	participants := setupPeers(2, params, sigma)

	t.Run(testString("Test tensor/", 2, params), func(t *testing.T) {

		// generate test values
		values1 := newTestValue(params, complex(-1, -1), complex(1, 1))
		values2 := newTestValue(params, complex(-1, -1), complex(1, 1))

		// Encrypt
		cipher1 := participants[0].Encrypt(values1)
		cipher2 := participants[1].Encrypt(values2)

		// evaluate multiplication
		evaluator := NewMKEvaluator(params)
		resCipher := evaluator.Mul(cipher1, cipher2)

		// decrypt using all secret keys
		sk1 := participants[0].GetSecretKey()
		sk2 := participants[1].GetSecretKey()

		decrypted := DecryptMul([]*mkrlwe.MKSecretKey{sk1, sk2}, resCipher, params)

		// perform the operation in the plaintext space
		for i := 0; i < len(values1); i++ {
			values1[i] *= values2[i]
		}

		// check results
		verifyTestVectors(params, values1, decrypted, t)
	})
}

func testTensor(t *testing.T, params *ckks.Parameters) {
	sigma := 6.0

	participants := setupPeers(1, params, sigma)

	t.Run(testString("Test tensor/", 1, params), func(t *testing.T) {

		// generate test values
		values1 := newTestValue(params, complex(-1, -1), complex(1, 1))

		// Encrypt
		cipher1 := participants[0].Encrypt(values1)

		// evaluate multiplication
		evaluator := NewMKEvaluator(params)

		resCipher := evaluator.Mul(cipher1, cipher1)

		// decrypt using all secret keys
		sk1 := participants[0].GetSecretKey()

		decrypted := DecryptMul([]*mkrlwe.MKSecretKey{sk1}, resCipher, params)
		// perform the operation in the plaintext space
		for i := 0; i < len(values1); i++ {
			values1[i] *= values1[i]
		}

		// check results
		verifyTestVectors(params, values1, decrypted, t)
	})
}

func testTensor2(t *testing.T, params *ckks.Parameters) {
	sigma := 6.0

	participants := setupPeers(1, params, sigma)

	t.Run(testString("Test tensor 2/", 1, params), func(t *testing.T) {

		// generate test values
		values1 := newTestValue(params, complex(-1, -1), complex(1, 1))

		// Encrypt
		cipher1 := participants[0].Encrypt(values1)

		// evaluate multiplication
		evaluator := NewMKEvaluator(params)

		resCipher := evaluator.Mul(cipher1, cipher1)

		// decrypt using all secret keys
		sk1 := participants[0].GetSecretKey()

		innerP1, innerP2 := DecryptAndCompare([]*mkrlwe.MKSecretKey{sk1}, resCipher, params, cipher1, cipher1)

		plaintext := ckks.NewPlaintext(*params, cipher1.Ciphertexts.Level(), cipher1.Ciphertexts.Scale()*cipher1.Ciphertexts.Scale())

		plaintext.SetValue(innerP2)

		encoder := ckks.NewEncoder(*params)

		decrypted := encoder.Decode(plaintext, params.LogSlots())

		// perform the operation in the plaintext space
		for i := 0; i < len(values1); i++ {
			values1[i] *= values1[i]
		}

		require.Equal(t, true, mkrlwe.EqualsPoly(innerP1, innerP2))

		verifyTestVectors(params, values1, decrypted, t)
	})
}

func testSquare(t *testing.T, params *ckks.Parameters) {

	sigma := 6.0

	participants := setupPeers(1, params, sigma)

	t.Run(testString("Test multiplication/", 1, params), func(t *testing.T) {

		// generate test value
		values1 := newTestValue(params, complex(-1, -1), complex(1, 1))

		// Encrypt
		cipher1 := participants[0].Encrypt(values1)

		// evaluate multiplication
		evaluator := NewMKEvaluator(params)
		evalKeys := []*mkrlwe.MKEvaluationKey{participants[0].GetEvaluationKey()}
		publicKeys := []*mkrlwe.MKPublicKey{participants[0].GetPublicKey()}

		resCipher := evaluator.Mul(cipher1, cipher1)

		evaluator.RelinInPlace(resCipher, evalKeys, publicKeys)

		// decrypt
		partialDec1 := participants[0].GetPartialDecryption(resCipher)

		decrypted := participants[0].Decrypt(resCipher, []*ring.Poly{partialDec1})

		// perform the operation in the plaintext space
		for i := 0; i < len(values1); i++ {
			values1[i] *= values1[i]
		}

		// check results
		verifyTestVectors(params, values1, decrypted, t)
	})

}

func testMul(t *testing.T, params *ckks.Parameters) {

	sigma := 6.0

	participants := setupPeers(2, params, sigma)

	t.Run(testString("Test multiplication/", 2, params), func(t *testing.T) {

		// generate test values
		values1 := newTestValue(params, complex(-1, -1), complex(1, 1))
		values2 := newTestValue(params, complex(-1, -1), complex(1, 1))

		// Encrypt
		cipher1 := participants[0].Encrypt(values1)
		cipher2 := participants[1].Encrypt(values2)

		// evaluate multiplication
		evaluator := NewMKEvaluator(params)
		evalKeys := []*mkrlwe.MKEvaluationKey{participants[0].GetEvaluationKey(), participants[1].GetEvaluationKey()}
		publicKeys := []*mkrlwe.MKPublicKey{participants[0].GetPublicKey(), participants[1].GetPublicKey()}

		resCipher := evaluator.Mul(cipher1, cipher2)
		evaluator.RelinInPlace(resCipher, evalKeys, publicKeys)
		// decrypt
		partialDec1 := participants[0].GetPartialDecryption(resCipher)
		partialDec2 := participants[1].GetPartialDecryption(resCipher)

		decrypted := participants[0].Decrypt(resCipher, []*ring.Poly{partialDec1, partialDec2})
		// perform the operation in the plaintext space
		for i := 0; i < len(values1); i++ {
			values1[i] *= values2[i]
		}

		// check results
		verifyTestVectors(params, values1, decrypted, t)
	})

}

func testAddAfterMul(t *testing.T, params *ckks.Parameters) {

	sigma := 6.0

	participants := setupPeers(4, params, sigma)

	t.Run(testString("Test add after multiplication/", 4, params), func(t *testing.T) {

		// generate test values
		values1 := newTestValue(params, complex(-1, -1), complex(1, 1))
		values2 := newTestValue(params, complex(-1, -1), complex(1, 1))
		values3 := newTestValue(params, complex(-1, -1), complex(1, 1))
		values4 := newTestValue(params, complex(-1, -1), complex(1, 1))

		// Encrypt
		cipher1 := participants[0].Encrypt(values1)
		cipher2 := participants[1].Encrypt(values2)
		cipher3 := participants[2].Encrypt(values3)
		cipher4 := participants[3].Encrypt(values4)

		// pad and add in 2 steps
		evaluator := NewMKEvaluator(params)
		evalKeys := []*mkrlwe.MKEvaluationKey{participants[0].GetEvaluationKey(), participants[1].GetEvaluationKey(), participants[2].GetEvaluationKey(), participants[3].GetEvaluationKey()}
		publicKeys := []*mkrlwe.MKPublicKey{participants[0].GetPublicKey(), participants[1].GetPublicKey(), participants[2].GetPublicKey(), participants[3].GetPublicKey()}

		resCipher1 := evaluator.Mul(cipher1, cipher2)
		resCipher2 := evaluator.Mul(cipher3, cipher4)

		evaluator.RelinInPlace(resCipher1, evalKeys[:2], publicKeys[:2])
		evaluator.RelinInPlace(resCipher2, evalKeys[2:], publicKeys[2:])

		evaluator.DropLevel(resCipher1, 1)
		evaluator.Rescale(resCipher1, resCipher1)
		evaluator.DropLevel(resCipher2, 1)
		evaluator.Rescale(resCipher2, resCipher2)

		// verify intermediate results
		intermediate1 := make([]complex128, len(values1))
		intermediate2 := make([]complex128, len(values1))
		intermediateAdd := make([]complex128, len(values1))

		for i := 0; i < len(values1); i++ {
			intermediate1[i] = values1[i] * values2[i]
			intermediate2[i] = values3[i] * values4[i]
			intermediateAdd[i] = intermediate1[i] + values1[i]
		}

		partialDec1 := participants[0].GetPartialDecryption(resCipher1)
		partialDec2 := participants[1].GetPartialDecryption(resCipher1)
		partialDec3 := participants[2].GetPartialDecryption(resCipher2)
		partialDec4 := participants[3].GetPartialDecryption(resCipher2)

		dec1 := participants[0].Decrypt(resCipher1, []*ring.Poly{partialDec1, partialDec2})
		dec2 := participants[0].Decrypt(resCipher2, []*ring.Poly{partialDec3, partialDec4})

		verifyTestVectors(params, intermediate1, dec1, t)
		verifyTestVectors(params, intermediate2, dec2, t)

		resCipher := evaluator.Add(resCipher1, resCipher2)

		// decrypt
		partialDec1 = participants[0].GetPartialDecryption(resCipher)
		partialDec2 = participants[1].GetPartialDecryption(resCipher)
		partialDec3 = participants[2].GetPartialDecryption(resCipher)
		partialDec4 = participants[3].GetPartialDecryption(resCipher)

		decrypted := participants[0].Decrypt(resCipher, []*ring.Poly{partialDec1, partialDec2, partialDec3, partialDec4})

		// perform the operation in the plaintext space
		for i := 0; i < len(values1); i++ {
			values1[i] = (values1[i] * values2[i]) + (values3[i] * values4[i])
		}

		// check results
		verifyTestVectors(params, values1, decrypted, t)
	})

}

func testMulAfterAdd(t *testing.T, params *ckks.Parameters) {

	sigma := 6.0

	participants := setupPeers(4, params, sigma)

	t.Run(testString("Test multiplication after addition/", 4, params), func(t *testing.T) {

		// generate test values
		values1 := newTestValue(params, complex(-1, -1), complex(1, 1))
		values2 := newTestValue(params, complex(-1, -1), complex(1, 1))
		values3 := newTestValue(params, complex(-1, -1), complex(1, 1))
		values4 := newTestValue(params, complex(-1, -1), complex(1, 1))

		// Encrypt
		cipher1 := participants[0].Encrypt(values1)
		cipher2 := participants[1].Encrypt(values2)
		cipher3 := participants[2].Encrypt(values3)
		cipher4 := participants[3].Encrypt(values4)

		// pad and add in 2 steps
		evaluator := NewMKEvaluator(params)
		evalKeys := []*mkrlwe.MKEvaluationKey{participants[0].GetEvaluationKey(), participants[1].GetEvaluationKey(), participants[2].GetEvaluationKey(), participants[3].GetEvaluationKey()}
		publicKeys := []*mkrlwe.MKPublicKey{participants[0].GetPublicKey(), participants[1].GetPublicKey(), participants[2].GetPublicKey(), participants[3].GetPublicKey()}

		resCipher1 := evaluator.Add(cipher1, cipher2)
		resCipher2 := evaluator.Add(cipher3, cipher4)

		resCipher := evaluator.Mul(resCipher1, resCipher2)

		evaluator.RelinInPlace(resCipher, evalKeys, publicKeys)

		// decrypt
		partialDec1 := participants[0].GetPartialDecryption(resCipher)
		partialDec2 := participants[1].GetPartialDecryption(resCipher)
		partialDec3 := participants[2].GetPartialDecryption(resCipher)
		partialDec4 := participants[3].GetPartialDecryption(resCipher)

		decrypted := participants[0].Decrypt(resCipher, []*ring.Poly{partialDec1, partialDec2, partialDec3, partialDec4})

		// perform the operation in the plaintext space
		for i := 0; i < len(values1); i++ {
			values1[i] = (values1[i] + values2[i]) * (values3[i] + values4[i])
		}

		// check results
		verifyTestVectors(params, values1, decrypted, t)
	})

}

func testMulFourParticipants(t *testing.T, params *ckks.Parameters) {

	sigma := 6.0

	participants := setupPeers(4, params, sigma)

	t.Run(testString("Test multiplication/", 4, params), func(t *testing.T) {

		// generate test values
		values1 := newTestValue(params, complex(-1, -1), complex(1, 1))
		values2 := newTestValue(params, complex(-1, -1), complex(1, 1))
		values3 := newTestValue(params, complex(-1, -1), complex(1, 1))
		values4 := newTestValue(params, complex(-1, -1), complex(1, 1))

		// Encrypt
		cipher1 := participants[0].Encrypt(values1)
		cipher2 := participants[1].Encrypt(values2)
		cipher3 := participants[2].Encrypt(values3)
		cipher4 := participants[3].Encrypt(values4)

		// pad and add in 2 steps
		evaluator := NewMKEvaluator(params)
		evalKeys := []*mkrlwe.MKEvaluationKey{participants[0].GetEvaluationKey(), participants[1].GetEvaluationKey(), participants[2].GetEvaluationKey(), participants[3].GetEvaluationKey()}
		publicKeys := []*mkrlwe.MKPublicKey{participants[0].GetPublicKey(), participants[1].GetPublicKey(), participants[2].GetPublicKey(), participants[3].GetPublicKey()}

		resCipher1 := evaluator.Mul(cipher1, cipher2)
		resCipher2 := evaluator.Mul(cipher3, cipher4)

		evaluator.RelinInPlace(resCipher1, evalKeys[:2], publicKeys[:2])
		evaluator.RelinInPlace(resCipher2, evalKeys[2:], publicKeys[2:])

		evaluator.Rescale(resCipher1, resCipher1)
		evaluator.Rescale(resCipher2, resCipher2)

		// verify intermediate results
		intermediate1 := make([]complex128, len(values1))
		intermediate2 := make([]complex128, len(values1))

		for i := 0; i < len(values1); i++ {
			intermediate1[i] = values1[i] * values2[i]
			intermediate2[i] = values3[i] * values4[i]
		}

		partialDec1 := participants[0].GetPartialDecryption(resCipher1)
		partialDec2 := participants[1].GetPartialDecryption(resCipher1)
		partialDec3 := participants[2].GetPartialDecryption(resCipher2)
		partialDec4 := participants[3].GetPartialDecryption(resCipher2)

		dec1 := participants[0].Decrypt(resCipher1, []*ring.Poly{partialDec1, partialDec2})
		dec2 := participants[0].Decrypt(resCipher2, []*ring.Poly{partialDec3, partialDec4})

		verifyTestVectors(params, intermediate1, dec1, t)
		verifyTestVectors(params, intermediate2, dec2, t)

		resCipher := evaluator.Mul(resCipher1, resCipher2)

		evaluator.RelinInPlace(resCipher, evalKeys, publicKeys)
		evaluator.DropLevel(resCipher, 1)
		evaluator.Rescale(resCipher, resCipher)

		// decrypt
		partialDec1 = participants[0].GetPartialDecryption(resCipher)
		partialDec2 = participants[1].GetPartialDecryption(resCipher)
		partialDec3 = participants[2].GetPartialDecryption(resCipher)
		partialDec4 = participants[3].GetPartialDecryption(resCipher)

		decrypted := participants[0].Decrypt(resCipher, []*ring.Poly{partialDec1, partialDec2, partialDec3, partialDec4})

		// perform the operation in the plaintext space
		for i := 0; i < len(values1); i++ {
			values1[i] *= values2[i] * values3[i] * values4[i]
		}

		// check results
		verifyTestVectors(params, values1, decrypted, t)
	})

}

func testRotation(t *testing.T, params *ckks.Parameters) {
	sigma := 6.0

	participants := setupPeers(1, params, sigma)

	rots := []int{1, -1, 4, -4, 63, -63}

	t.Run(testString("Test Rotation/", 1, params), func(t *testing.T) {

		// generate test values
		values1 := newTestValue(params, complex(-1, -1), complex(1, 1))

		// Encrypt
		cipher1 := participants[0].Encrypt(values1)

		evaluator := NewMKEvaluator(params)

		for _, n := range rots {

			rotKey := participants[0].GetRotationKeys(n)

			values2 := utils.RotateComplex128Slice(values1, n)
			resCipher := evaluator.Rotate(cipher1, n, []*mkrlwe.MKEvalGalKey{rotKey})

			partialDec := participants[0].GetPartialDecryption(resCipher)

			decrypted := participants[0].Decrypt(resCipher, []*ring.Poly{partialDec})
			verifyTestVectors(params, values2, decrypted, t)
		}
	})

}

func testRotationTwoParticipants(t *testing.T, params *ckks.Parameters) {
	sigma := 6.0

	participants := setupPeers(2, params, sigma)

	rots := []int{1, -1, 4, -4, 63, -63}

	t.Run(testString("Test Rotation/", 2, params), func(t *testing.T) {

		// generate test values
		values1 := newTestValue(params, complex(-1, -1), complex(1, 1))
		values2 := newTestValue(params, complex(-1, -1), complex(1, 1))

		// Encrypt
		cipher1 := participants[0].Encrypt(values1)
		cipher2 := participants[1].Encrypt(values2)

		evaluator := NewMKEvaluator(params)

		// add both ciphertexts
		added := evaluator.Add(cipher1, cipher2)

		// operate on plaintext space
		for i := int(0); i < len(values1); i++ {
			values1[i] += values2[i]
		}

		for _, n := range rots {

			// generate rotation keys
			rotKey1 := participants[0].GetRotationKeys(n)
			rotKey2 := participants[1].GetRotationKeys(n)

			rotated := utils.RotateComplex128Slice(values1, n)
			resCipher := evaluator.Rotate(added, n, []*mkrlwe.MKEvalGalKey{rotKey1, rotKey2})

			partialDec1 := participants[0].GetPartialDecryption(resCipher)
			partialDec2 := participants[1].GetPartialDecryption(resCipher)

			decrypted := participants[0].Decrypt(resCipher, []*ring.Poly{partialDec1, partialDec2})

			verifyTestVectors(params, rotated, decrypted, t)
		}
	})

}

func testCkksMkbfvBridge(t *testing.T, params *ckks.Parameters) {

	encoder := ckks.NewEncoder(*params)
	keygen := ckks.NewKeyGenerator(*params)
	sk, pk := keygen.GenKeyPair()
	encryptorPK := ckks.NewEncryptorFromPk(*params, pk)

	t.Run(testString("Test Bridge BFV-MKBFV/", 2, params), func(t *testing.T) {

		// setup bfv environment and encrypt values in bfv
		values1 := newTestValue(params, complex(-1, -1), complex(1, 1))
		plaintext := encoder.EncodeNTTAtLvlNew(params.MaxLevel(), values1, params.LogSlots())

		var ciphertext1 *ckks.Ciphertext

		if params.PCount() != 0 {
			ciphertext1 = encryptorPK.EncryptNew(plaintext)
		} else {
			ciphertext1 = encryptorPK.EncryptFastNew(plaintext)
		}

		// switch to multi key setting and operate with other participant
		prng, err := utils.NewKeyedPRNG([]byte{'l', 'a', 't', 't', 'i', 'g', 'o'})

		if err != nil {
			panic(err)
		}

		// setup keys and public parameters
		a := mkrlwe.GenCommonPublicParam(&params.Parameters, prng)
		part1 := NewParticipantFromSecretKey(params, 6.0, a, sk)
		part2 := NewParticipant(params, 6.0, a)

		// perform addition
		values2 := newTestValue(params, complex(-1, -1), complex(1, 1))
		ciphertext2 := part2.Encrypt(values2)

		evaluator := NewMKEvaluator(params)

		res := evaluator.Add(ciphertext2, &MKCiphertext{Ciphertexts: ciphertext1, PeerID: []uint64{part1.GetID()}})

		// decrypt

		partDec1 := part1.GetPartialDecryption(res)
		partDec2 := part2.GetPartialDecryption(res)

		decrypted := part1.Decrypt(res, []*ring.Poly{partDec1, partDec2})

		//verify
		for i := range values1 {
			values1[i] += values2[i]
		}

		verifyTestVectors(params, values1, decrypted, t)
	})
}

func testMarshaler(t *testing.T, params *ckks.Parameters) {

	sigma := 6.0

	participants := setupPeers(2, params, sigma)

	t.Run(testString("Test Marshaler/", 2, params), func(t *testing.T) {

		// get random value
		value1 := newTestValue(params, complex(-1, -1), complex(1, 1))
		value2 := newTestValue(params, complex(-1, -1), complex(1, 1))

		//encrypt
		cipher1 := participants[0].Encrypt(value1)
		cipher2 := participants[1].Encrypt(value2)

		// add
		evaluator := NewMKEvaluator(params)
		res := evaluator.Add(cipher1, cipher2)

		data := res.MarshalBinary()

		unMarshaled := new(MKCiphertext)

		unMarshaled.UnmarshalBinary(data)

		if !mkrlwe.EqualsSlice(unMarshaled.PeerID, res.PeerID) {
			t.Error("Marshaler error with peer IDs")
		}

		if !mkrlwe.EqualsPoly(res.Ciphertexts.Value[0], unMarshaled.Ciphertexts.Value[0]) ||
			!mkrlwe.EqualsPoly(res.Ciphertexts.Value[1], unMarshaled.Ciphertexts.Value[1]) ||
			!mkrlwe.EqualsPoly(res.Ciphertexts.Value[2], unMarshaled.Ciphertexts.Value[2]) {

			t.Error("Marshaler error with ciphertext")
		}

	})

}

func Test_Utils(t *testing.T) {

	s1 := []uint64{0, 2, 1}

	s2 := []uint64{3, 7, 12, 1, 0}

	expected := []uint64{0, 1, 2, 3, 7, 12}

	res := mkrlwe.MergeSlices(s1, s2)

	if !mkrlwe.EqualsSlice(expected, res) {
		t.Errorf("MergeSlices method failed test")
	}

}

// Generates keys for a set of peers identified by their peerID using a certain ckks parameter index
// returns the slice of keys with the ckks parameters
func setupPeers(peerNbr uint64, params *ckks.Parameters, sigmaSmudging float64) []MKParticipant {

	res := make([]MKParticipant, peerNbr)

	prng, err := utils.NewKeyedPRNG([]byte{'l', 'a', 't', 't', 'i', 'g', 'o'})

	if err != nil {
		panic(err)
	}

	a := mkrlwe.GenCommonPublicParam(&params.Parameters, prng)

	for i := 0; i < int(peerNbr); i++ {

		res[i] = NewParticipant(params, sigmaSmudging, a)
	}

	return res
}

func newTestValue(params *ckks.Parameters, a, b complex128) []complex128 {

	logSlots := params.LogSlots()

	values := make([]complex128, 1<<logSlots)

	for i := uint64(0); i < 1<<logSlots; i++ {
		values[i] = complex(utils.RandFloat64(real(a), real(b)), utils.RandFloat64(imag(a), imag(b)))
	}

	values[0] = complex(0.607538, 0)

	return values
}

func verifyTestVectors(params *ckks.Parameters, valuesWant, valuesTest []complex128, t *testing.T) {

	precStats := GetPrecisionStats(params, valuesWant, valuesTest)

	if *printPrecisionStats {
		t.Log(precStats.String())
	}

	require.GreaterOrEqual(t, real(precStats.MeanPrecision), minPrec)
	require.GreaterOrEqual(t, imag(precStats.MeanPrecision), minPrec)
}

// GetPrecisionStats generates a PrecisionStats struct from the reference values and the decrypted values
func GetPrecisionStats(params *ckks.Parameters, valuesWant, valuesTest []complex128) (prec ckks.PrecisionStats) {

	logSlots := params.LogSlots()
	slots := uint64(1 << logSlots)

	var deltaReal, deltaImag float64

	var delta complex128

	diff := make([]complex128, slots)

	prec.MaxDelta = complex(0, 0)
	prec.MinDelta = complex(1, 1)

	prec.MeanDelta = complex(0, 0)

	cdfResol := 500

	prec.RealDist = make([]struct {
		Prec  float64
		Count int
	}, cdfResol)
	prec.ImagDist = make([]struct {
		Prec  float64
		Count int
	}, cdfResol)

	precReal := make([]float64, len(valuesWant))
	precImag := make([]float64, len(valuesWant))

	for i := range valuesWant {

		delta = valuesTest[i] - valuesWant[i]
		deltaReal = math.Abs(real(delta))
		deltaImag = math.Abs(imag(delta))
		precReal[i] = math.Log2(1 / deltaReal)
		precImag[i] = math.Log2(1 / deltaImag)

		diff[i] += complex(deltaReal, deltaImag)

		prec.MeanDelta += diff[i]

		if deltaReal > real(prec.MaxDelta) {
			prec.MaxDelta = complex(deltaReal, imag(prec.MaxDelta))
		}

		if deltaImag > imag(prec.MaxDelta) {
			prec.MaxDelta = complex(real(prec.MaxDelta), deltaImag)
		}

		if deltaReal < real(prec.MinDelta) {
			prec.MinDelta = complex(deltaReal, imag(prec.MinDelta))
		}

		if deltaImag < imag(prec.MinDelta) {
			prec.MinDelta = complex(real(prec.MinDelta), deltaImag)
		}
	}

	calcCDF(precReal, cdfResol, prec.RealDist)
	calcCDF(precImag, cdfResol, prec.ImagDist)

	prec.MinPrecision = deltaToPrecision(prec.MaxDelta)
	prec.MaxPrecision = deltaToPrecision(prec.MinDelta)
	prec.MeanDelta /= complex(float64(slots), 0)
	prec.MeanPrecision = deltaToPrecision(prec.MeanDelta)
	prec.MedianDelta = calcmedian(diff)
	prec.MedianPrecision = deltaToPrecision(prec.MedianDelta)
	return prec
}

func deltaToPrecision(c complex128) complex128 {
	return complex(math.Log2(1/real(c)), math.Log2(1/imag(c)))
}

func calcCDF(precs []float64, cdfResol int, res []struct {
	Prec  float64
	Count int
}) {
	sortedPrecs := make([]float64, len(precs))
	copy(sortedPrecs, precs)
	sort.Float64s(sortedPrecs)
	minPrec := sortedPrecs[0]
	maxPrec := sortedPrecs[len(sortedPrecs)-1]
	for i := 0; i < cdfResol; i++ {
		curPrec := minPrec + float64(i)*(maxPrec-minPrec)/float64(cdfResol)
		for countSmaller, p := range sortedPrecs {
			if p >= curPrec {
				res[i].Prec = curPrec
				res[i].Count = countSmaller
				break
			}
		}
	}
}

func calcmedian(values []complex128) (median complex128) {

	tmp := make([]float64, len(values))

	for i := range values {
		tmp[i] = real(values[i])
	}

	sort.Float64s(tmp)

	for i := range values {
		values[i] = complex(tmp[i], imag(values[i]))
	}

	for i := range values {
		tmp[i] = imag(values[i])
	}

	sort.Float64s(tmp)

	for i := range values {
		values[i] = complex(real(values[i]), tmp[i])
	}

	index := len(values) / 2

	if len(values)&1 == 1 {
		return values[index]
	}

	if index+1 == len(values) {
		return values[index]
	}

	return (values[index] + values[index+1]) / 2
}

// DecryptMul takes secret keys of k participants and a multikey ciphertext ct of dimension (k+1) ** 2 = ct1 * ct2
// Computes < ct, sk*sk > in order to check that it is equal to the product of plaintexts
func DecryptMul(keys []*mkrlwe.MKSecretKey, ct *MKCiphertext, params *ckks.Parameters) []complex128 {

	ringQ := mkrlwe.GetRingQ(&params.Parameters)
	ringQP := mkrlwe.GetRingQP(&params.Parameters)

	// Compute the tensor product of the secret keys : sk * sk
	level := ct.Ciphertexts.Level()
	nbrElements := len(keys) + 1
	tensorDim := nbrElements * nbrElements

	keyTensor := make([]*ring.Poly, tensorDim)
	for i := 0; i < tensorDim; i++ {
		keyTensor[i] = ringQP.NewPoly()
	}
	if len(ct.PeerID) == 1 {
		keyTensor[1] = keys[0].Key.Value
		keyTensor[2] = keys[0].Key.Value
		ringQP.MulCoeffsMontgomery(keys[0].Key.Value, keys[0].Key.Value, keyTensor[3])
	} else if len(ct.PeerID) == 2 {
		// detect if ciphertext (c01 * c02, 0, ..) or (c01 * c02, c1 * c02, ..)
		if mkrlwe.EqualsPoly(ringQ.NewPoly(), ct.Ciphertexts.Value[1]) {
			keyTensor[2] = keys[1].Key.Value
			keyTensor[3] = keys[0].Key.Value
			ringQP.MulCoeffsMontgomery(keys[0].Key.Value, keys[1].Key.Value, keyTensor[5])
		} else {
			keyTensor[1] = keys[1].Key.Value
			keyTensor[6] = keys[0].Key.Value
			ringQP.MulCoeffsMontgomery(keys[0].Key.Value, keys[1].Key.Value, keyTensor[7])
		}
	} else {
		panic("This function was only designed to process ciphertext up to degree 1")
	}

	// Compute inner product of two matrices -> compute tr(ct^T  sk*sk)
	res := ringQ.NewPoly()
	for i := 1; i < tensorDim; i++ {
		if ct.Ciphertexts.Value[i] == nil {
			continue
		}
		ringQ.MulCoeffsMontgomeryAndAddLvl(level, keyTensor[i], ct.Ciphertexts.Value[i], res)
	}
	ringQ.AddLvl(level, res, ct.Ciphertexts.Value[0], res)

	ringQ.ReduceLvl(level, res, res)
	res.Coeffs = res.Coeffs[:level+1]
	plaintext := ckks.NewPlaintext(*params, level, ct.Ciphertexts.Scale())

	plaintext.SetValue(res)

	encoder := ckks.NewEncoder(*params)

	return encoder.Decode(plaintext, params.LogSlots())
}

// DecryptSimple decrypt a ciphertext with only 1 participant involved using a private key
func DecryptSimple(key *mkrlwe.MKSecretKey, ct *MKCiphertext, params *ckks.Parameters) []complex128 {

	ringQ := mkrlwe.GetRingQ(&params.Parameters)

	// Compute the tensor product of the secret keys : sk * sk
	level := ct.Ciphertexts.Level()

	// derypt with secret key
	res := ringQ.NewPoly()

	ringQ.MulCoeffsMontgomeryAndAddLvl(level, key.Key.Value, ct.Ciphertexts.Value[1], res)
	ringQ.AddLvl(level, res, ct.Ciphertexts.Value[0], res)

	ringQ.ReduceLvl(level, res, res)
	res.Coeffs = res.Coeffs[:level+1]

	plaintext := ckks.NewPlaintext(*params, level, ct.Ciphertexts.Scale())

	plaintext.SetValue(res)

	encoder := ckks.NewEncoder(*params)

	return encoder.Decode(plaintext, params.LogSlots())
}

// DecryptAndCompare takes secret keys of k participants and a multikey ciphertext ct of dimension (k+1) ** 2 = ct1 * ct2
// Computes < ct, sk*sk > in order to check that it is equal to <ct1,sk1>.<ct2,sk2>
func DecryptAndCompare(keys []*mkrlwe.MKSecretKey, ct *MKCiphertext, params *ckks.Parameters, ct1 *MKCiphertext, ct2 *MKCiphertext) (*ring.Poly, *ring.Poly) {

	ringQ := mkrlwe.GetRingQ(&params.Parameters)
	ringQP := mkrlwe.GetRingQP(&params.Parameters)

	// Compute the tensor product of the secret keys : sk * sk
	level := ct.Ciphertexts.Level()
	nbrElements := len(keys) + 1
	tensorDim := nbrElements * nbrElements

	keyTensor := make([]*ring.Poly, tensorDim)
	for i := 0; i < tensorDim; i++ {
		keyTensor[i] = ringQP.NewPoly()
	}
	keyTensor[1] = keys[0].Key.Value
	keyTensor[2] = keys[0].Key.Value
	ringQP.MulCoeffsMontgomery(keys[0].Key.Value, keys[0].Key.Value, keyTensor[3])

	// Compute inner product of two matrices -> compute tr(ct^T  sk*sk)
	res := ringQ.NewPoly()

	for i := 1; i < tensorDim; i++ {
		ringQ.MulCoeffsMontgomeryAndAddLvl(level, keyTensor[i], ct.Ciphertexts.Value[i], res)
	}

	ringQ.AddLvl(level, res, ct.Ciphertexts.Value[0], res)

	ringQ.ReduceLvl(level, res, res)
	res.Coeffs = res.Coeffs[:level+1]

	// Compute the two inner products
	innerProduct1 := ringQ.NewPoly()
	innerProduct2 := ringQ.NewPoly()

	ringQ.MulCoeffsMontgomeryAndAddLvl(level, keys[0].Key.Value, ct1.Ciphertexts.Value[1], innerProduct1)
	ringQ.AddLvl(level, innerProduct1, ct1.Ciphertexts.Value[0], innerProduct1)

	ringQ.MulCoeffsMontgomeryAndAddLvl(level, keys[0].Key.Value, ct2.Ciphertexts.Value[1], innerProduct2)
	ringQ.AddLvl(level, innerProduct2, ct2.Ciphertexts.Value[0], innerProduct2)

	ringQ.ReduceLvl(level, innerProduct1, innerProduct1)
	ringQ.ReduceLvl(level, innerProduct2, innerProduct2)
	innerProduct1.Coeffs = innerProduct1.Coeffs[:level+1]
	innerProduct2.Coeffs = innerProduct2.Coeffs[:level+1]

	mulInnerProduct := ringQ.NewPoly()
	ringQ.MulCoeffsAndAdd(innerProduct1, innerProduct2, mulInnerProduct)
	ringQ.ReduceLvl(level, mulInnerProduct, mulInnerProduct)
	mulInnerProduct.Coeffs = mulInnerProduct.Coeffs[:level+1]
	return res, mulInnerProduct

}
