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
			testAdd(t, p)
			testAddManyTimeSameCipher(t, p)
			testAddFourParticipants(t, p)
			testAddPlaintext(t, p)
			testAddPlaintextTwoParticipants(t, p)
			testSub(t, p)
			testSubPlaintext(t, p)
			testNeg(t, p)
			testNegTwoParticipants(t, p)
			testSubPlaintextTwoParticipants(t, p)
			testMulPlaintext(t, p)
			testMulPlaintextTwoParticipants(t, p)
			//	testCkksMkbfvBridge(t, p)
			//	testMarshaler(t, p)
			testRotation(t, p)
			testRotationTwoParticipants(t, p)
			testSquare(t, p)
			testMul(t, p)
			testMulAfterAdd(t, p)

			if i != 5 && i != 6 {
				testAddAfterMul(t, p)
				testMulFourParticipants(t, p)
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
		ciphers := evaluator.ConvertToMKCiphertext([]*ckks.Ciphertext{cipher1, cipher2}, []uint64{1, 2})
		resCipher := evaluator.Add(ciphers[0], ciphers[1])
		ckksCiphers := evaluator.ConvertToCKKSCiphertext(resCipher)

		// decrypt
		partialDec1 := participants[0].GetPartialDecryption(ckksCiphers[0])
		partialDec2 := participants[1].GetPartialDecryption(ckksCiphers[1])

		decrypted := participants[0].Decrypt(ckksCiphers[0], []*ring.Poly{partialDec1, partialDec2})

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
		ciphers := evaluator.ConvertToMKCiphertext([]*ckks.Ciphertext{cipher1, cipher2}, []uint64{1, 2})

		resCipher := evaluator.Add(ciphers[0], ciphers[1])
		resCipher = evaluator.Add(resCipher, ciphers[1])
		resCipher = evaluator.Add(resCipher, ciphers[1])
		resCipher = evaluator.Add(resCipher, ciphers[1])

		ckksRes := evaluator.ConvertToCKKSCiphertext(resCipher)

		// decrypt
		partialDec1 := participants[0].GetPartialDecryption(ckksRes[0])
		partialDec2 := participants[1].GetPartialDecryption(ckksRes[1])

		decrypted := participants[0].Decrypt(ckksRes[0], []*ring.Poly{partialDec1, partialDec2})

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
		ciphers := evaluator.ConvertToMKCiphertext([]*ckks.Ciphertext{cipher1, cipher2}, []uint64{1, 2})
		resCipher := evaluator.Sub(ciphers[0], ciphers[1])

		ckksRes := evaluator.ConvertToCKKSCiphertext(resCipher)

		// decrypt
		partialDec1 := participants[0].GetPartialDecryption(ckksRes[0])
		partialDec2 := participants[1].GetPartialDecryption(ckksRes[1])

		decrypted := participants[0].Decrypt(ckksRes[0], []*ring.Poly{partialDec1, partialDec2})

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
		value2 := newTestValue(params, complex(-1, -1), complex(1, 1))

		// encrypt
		cipher := participants[0].Encrypt(value1)

		// add plaintext to ciphertext
		evaluator := NewMKEvaluator(params)
		ciphers := evaluator.ConvertToMKCiphertext([]*ckks.Ciphertext{cipher}, []uint64{1})
		pt := evaluator.NewPlaintextFromValue(value2)
		resCipher := evaluator.AddPlaintext(pt, ciphers[0])
		resCKKS := evaluator.ConvertToCKKSCiphertext(resCipher)

		// decrypt
		partialDec1 := participants[0].GetPartialDecryption(resCKKS[0])

		decrypted := participants[0].Decrypt(resCKKS[0], []*ring.Poly{partialDec1})

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
		value3 := newTestValue(params, complex(-1, -1), complex(1, 1))

		// encrypt
		cipher1 := participants[0].Encrypt(value1)
		cipher2 := participants[1].Encrypt(value2)

		// add plaintext to one of the ciphertext then add both ciphertexts
		evaluator := NewMKEvaluator(params)
		ciphers := evaluator.ConvertToMKCiphertext([]*ckks.Ciphertext{cipher1, cipher2}, []uint64{1, 2})

		pt := evaluator.NewPlaintextFromValue(value3)
		resCipher1 := evaluator.AddPlaintext(pt, ciphers[0])

		resCipher := evaluator.Add(resCipher1, ciphers[1])

		resCKKS := evaluator.ConvertToCKKSCiphertext(resCipher)

		// decrypt
		partialDec1 := participants[0].GetPartialDecryption(resCKKS[0])
		partialDec2 := participants[1].GetPartialDecryption(resCKKS[1])

		decrypted := participants[0].Decrypt(resCKKS[0], []*ring.Poly{partialDec1, partialDec2})

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
		value2 := newTestValue(params, complex(-1, -1), complex(1, 1))

		// encrypt
		cipher := participants[0].Encrypt(value1)

		// sub plaintext to ciphertext
		evaluator := NewMKEvaluator(params)
		ciphers := evaluator.ConvertToMKCiphertext([]*ckks.Ciphertext{cipher}, []uint64{1})

		pt := evaluator.NewPlaintextFromValue(value2)
		resCipher := evaluator.SubPlaintext(ciphers[0], pt)

		resCKKS := evaluator.ConvertToCKKSCiphertext(resCipher)

		// decrypt
		partialDec1 := participants[0].GetPartialDecryption(resCKKS[0])

		decrypted := participants[0].Decrypt(resCKKS[0], []*ring.Poly{partialDec1})

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
		value3 := newTestValue(params, complex(-1, -1), complex(1, 1))

		// encrypt
		cipher1 := participants[0].Encrypt(value1)
		cipher2 := participants[1].Encrypt(value2)

		// sub plaintext to one of the ciphertext then sub both ciphertexts
		evaluator := NewMKEvaluator(params)
		ciphers := evaluator.ConvertToMKCiphertext([]*ckks.Ciphertext{cipher1, cipher2}, []uint64{1, 2})

		pt := evaluator.NewPlaintextFromValue(value3)
		resCipher1 := evaluator.SubPlaintext(ciphers[0], pt)

		resCipher := evaluator.Sub(resCipher1, ciphers[1])

		resCKKS := evaluator.ConvertToCKKSCiphertext(resCipher)

		// decrypt
		partialDec1 := participants[0].GetPartialDecryption(resCKKS[0])
		partialDec2 := participants[1].GetPartialDecryption(resCKKS[1])

		decrypted := participants[0].Decrypt(resCKKS[0], []*ring.Poly{partialDec1, partialDec2})

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
		ciphers := evaluator.ConvertToMKCiphertext([]*ckks.Ciphertext{cipher}, []uint64{1})
		resCipher := evaluator.Add(evaluator.Neg(ciphers[0]), ciphers[0])

		resCKKS := evaluator.ConvertToCKKSCiphertext(resCipher)

		// decrypt
		partialDec1 := participants[0].GetPartialDecryption(resCKKS[0])

		decrypted := participants[0].Decrypt(resCKKS[0], []*ring.Poly{partialDec1})

		// perform the operation in the plaintext space
		for i := 0; i < len(value1); i++ {
			value1[i] -= value1[i]
		}

		// check results
		verifyTestVectors(params, value1, decrypted, t)

	})
}

func testNegTwoParticipants(t *testing.T, params *ckks.Parameters) {

	sigma := 6.0

	participants := setupPeers(2, params, sigma)

	t.Run(testString("Test Negation/", 2, params), func(t *testing.T) {

		value1 := newTestValue(params, complex(-1, -1), complex(1, 1))
		value2 := newTestValue(params, complex(-1, -1), complex(1, 1))

		// encrypt
		cipher1 := participants[0].Encrypt(value1)
		cipher2 := participants[1].Encrypt(value2)

		// add with negated ciphertext
		evaluator := NewMKEvaluator(params)
		ciphers := evaluator.ConvertToMKCiphertext([]*ckks.Ciphertext{cipher1, cipher2}, []uint64{1, 2})
		added := evaluator.Add(ciphers[1], ciphers[0])

		resCipher := evaluator.Neg(added)

		resCKKS := evaluator.ConvertToCKKSCiphertext(resCipher)

		// decrypt
		partialDec1 := participants[0].GetPartialDecryption(resCKKS[0])
		partialDec2 := participants[1].GetPartialDecryption(resCKKS[1])

		decrypted := participants[0].Decrypt(resCKKS[0], []*ring.Poly{partialDec1, partialDec2})

		// perform the operation in the plaintext space
		for i := 0; i < len(value1); i++ {
			value1[i] = -value1[i] - value2[i]
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
		ciphers := evaluator.ConvertToMKCiphertext([]*ckks.Ciphertext{cipher1, cipher2, cipher3, cipher4}, []uint64{1, 2, 3, 4})
		resCipher1 := evaluator.Add(ciphers[0], ciphers[1])
		resCipher2 := evaluator.Add(ciphers[2], ciphers[3])

		resCipher := evaluator.Add(resCipher1, resCipher2)

		resCKKS := evaluator.ConvertToCKKSCiphertext(resCipher)

		// decrypt
		partialDec1 := participants[0].GetPartialDecryption(resCKKS[0])
		partialDec2 := participants[1].GetPartialDecryption(resCKKS[1])
		partialDec3 := participants[2].GetPartialDecryption(resCKKS[2])
		partialDec4 := participants[3].GetPartialDecryption(resCKKS[3])

		decrypted := participants[0].Decrypt(resCKKS[0], []*ring.Poly{partialDec1, partialDec2, partialDec3, partialDec4})

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
		value2 := newTestValue(params, complex(-1, -1), complex(1, 1))

		// encrypt
		cipher := participants[0].Encrypt(value1)

		// multiply plaintext and ciphertext
		evaluator := NewMKEvaluator(params)
		ciphers := evaluator.ConvertToMKCiphertext([]*ckks.Ciphertext{cipher}, []uint64{1})

		pt := evaluator.NewPlaintextFromValue(value2)
		resCipher := evaluator.MultPlaintext(pt, ciphers[0])

		resCKKS := evaluator.ConvertToCKKSCiphertext(resCipher)

		// decrypt
		partialDec1 := participants[0].GetPartialDecryption(resCKKS[0])

		decrypted := participants[0].Decrypt(resCKKS[0], []*ring.Poly{partialDec1})

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
		value3 := newTestValue(params, complex(-1, -1), complex(1, 1))

		// encrypt
		cipher1 := participants[0].Encrypt(value1)
		cipher2 := participants[1].Encrypt(value2)

		// add both ciphertexts then multiply by plaintext
		evaluator := NewMKEvaluator(params)
		ciphers := evaluator.ConvertToMKCiphertext([]*ckks.Ciphertext{cipher1, cipher2}, []uint64{1, 2})

		resCipherTMP := evaluator.Add(ciphers[0], ciphers[1])

		pt := evaluator.NewPlaintextFromValue(value3)
		resCipher := evaluator.MultPlaintext(pt, resCipherTMP)

		resCKKS := evaluator.ConvertToCKKSCiphertext(resCipher)

		// decrypt
		partialDec1 := participants[0].GetPartialDecryption(resCKKS[0])
		partialDec2 := participants[1].GetPartialDecryption(resCKKS[1])

		decrypted := participants[0].Decrypt(resCKKS[0], []*ring.Poly{partialDec1, partialDec2})

		// perform the operation in the plaintext space
		for i := 0; i < len(value1); i++ {
			value1[i] = (value1[i] + value2[i]) * value3[i]
		}

		// check results
		verifyTestVectors(params, value1, decrypted, t)

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
		ciphers := evaluator.ConvertToMKCiphertext([]*ckks.Ciphertext{cipher1}, []uint64{1})

		evk := participants[0].GetEvaluationKey()
		evk.PeerID = 1
		evalKeys := []*mkrlwe.MKEvaluationKey{evk}

		pk := participants[0].GetPublicKey()
		pk.PeerID = 1
		publicKeys := []*mkrlwe.MKPublicKey{pk}

		resCipher := evaluator.Mul(ciphers[0], ciphers[0])

		evaluator.RelinInPlace(resCipher, evalKeys, publicKeys)

		resCKKS := evaluator.ConvertToCKKSCiphertext(resCipher)

		// decrypt
		partialDec1 := participants[0].GetPartialDecryption(resCKKS[0])

		decrypted := participants[0].Decrypt(resCKKS[0], []*ring.Poly{partialDec1})

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
		ciphers := evaluator.ConvertToMKCiphertext([]*ckks.Ciphertext{cipher1, cipher2}, []uint64{1, 2})

		evk1 := participants[0].GetEvaluationKey()
		evk1.PeerID = 1
		evk2 := participants[1].GetEvaluationKey()
		evk2.PeerID = 2

		evalKeys := []*mkrlwe.MKEvaluationKey{participants[0].GetEvaluationKey(), participants[1].GetEvaluationKey()}

		pk1 := participants[0].GetPublicKey()
		pk2 := participants[1].GetPublicKey()
		pk1.PeerID = 1
		pk2.PeerID = 2

		publicKeys := []*mkrlwe.MKPublicKey{pk1, pk2}

		resCipher := evaluator.Mul(ciphers[0], ciphers[1])
		evaluator.RelinInPlace(resCipher, evalKeys, publicKeys)

		resCKKS := evaluator.ConvertToCKKSCiphertext(resCipher)

		// decrypt
		partialDec1 := participants[0].GetPartialDecryption(resCKKS[0])
		partialDec2 := participants[1].GetPartialDecryption(resCKKS[1])

		decrypted := participants[0].Decrypt(resCKKS[0], []*ring.Poly{partialDec1, partialDec2})
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
		ciphers := evaluator.ConvertToMKCiphertext([]*ckks.Ciphertext{cipher1, cipher2, cipher3, cipher4}, []uint64{1, 2, 3, 4})

		evk1 := participants[0].GetEvaluationKey()
		evk1.PeerID = 1
		evk2 := participants[1].GetEvaluationKey()
		evk2.PeerID = 2
		evk3 := participants[2].GetEvaluationKey()
		evk3.PeerID = 3
		evk4 := participants[3].GetEvaluationKey()
		evk4.PeerID = 4

		evalKeys := []*mkrlwe.MKEvaluationKey{evk1, evk2, evk3, evk4}

		pk1 := participants[0].GetPublicKey()
		pk1.PeerID = 1
		pk2 := participants[1].GetPublicKey()
		pk2.PeerID = 2
		pk3 := participants[2].GetPublicKey()
		pk3.PeerID = 3
		pk4 := participants[3].GetPublicKey()
		pk4.PeerID = 4

		publicKeys := []*mkrlwe.MKPublicKey{pk1, pk2, pk3, pk4}

		resCipher1 := evaluator.Mul(ciphers[0], ciphers[1])
		resCipher2 := evaluator.Mul(ciphers[2], ciphers[3])

		evaluator.RelinInPlace(resCipher1, evalKeys[:2], publicKeys[:2])
		evaluator.RelinInPlace(resCipher2, evalKeys[2:], publicKeys[2:])

		evaluator.DropLevel(resCipher1, 1)
		evaluator.Rescale(resCipher1, resCipher1)
		evaluator.DropLevel(resCipher2, 1)
		evaluator.Rescale(resCipher2, resCipher2)

		resCipher := evaluator.Add(resCipher1, resCipher2)

		resCKKS := evaluator.ConvertToCKKSCiphertext(resCipher)

		// decrypt
		partialDec1 := participants[0].GetPartialDecryption(resCKKS[0])
		partialDec2 := participants[1].GetPartialDecryption(resCKKS[1])
		partialDec3 := participants[2].GetPartialDecryption(resCKKS[2])
		partialDec4 := participants[3].GetPartialDecryption(resCKKS[3])

		decrypted := participants[0].Decrypt(resCKKS[0], []*ring.Poly{partialDec1, partialDec2, partialDec3, partialDec4})

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
		ciphers := evaluator.ConvertToMKCiphertext([]*ckks.Ciphertext{cipher1, cipher2, cipher3, cipher4}, []uint64{1, 2, 3, 4})

		evk1 := participants[0].GetEvaluationKey()
		evk1.PeerID = 1
		evk2 := participants[1].GetEvaluationKey()
		evk2.PeerID = 2
		evk3 := participants[2].GetEvaluationKey()
		evk3.PeerID = 3
		evk4 := participants[3].GetEvaluationKey()
		evk4.PeerID = 4

		evalKeys := []*mkrlwe.MKEvaluationKey{evk1, evk2, evk3, evk4}

		pk1 := participants[0].GetPublicKey()
		pk1.PeerID = 1
		pk2 := participants[1].GetPublicKey()
		pk2.PeerID = 2
		pk3 := participants[2].GetPublicKey()
		pk3.PeerID = 3
		pk4 := participants[3].GetPublicKey()
		pk4.PeerID = 4

		publicKeys := []*mkrlwe.MKPublicKey{pk1, pk2, pk3, pk4}

		resCipher1 := evaluator.Add(ciphers[0], ciphers[1])
		resCipher2 := evaluator.Add(ciphers[2], ciphers[3])

		resCipher := evaluator.Mul(resCipher1, resCipher2)

		evaluator.RelinInPlace(resCipher, evalKeys, publicKeys)

		resCKKS := evaluator.ConvertToCKKSCiphertext(resCipher)

		// decrypt
		partialDec1 := participants[0].GetPartialDecryption(resCKKS[0])
		partialDec2 := participants[1].GetPartialDecryption(resCKKS[1])
		partialDec3 := participants[2].GetPartialDecryption(resCKKS[2])
		partialDec4 := participants[3].GetPartialDecryption(resCKKS[3])

		decrypted := participants[0].Decrypt(resCKKS[0], []*ring.Poly{partialDec1, partialDec2, partialDec3, partialDec4})

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
		ciphers := evaluator.ConvertToMKCiphertext([]*ckks.Ciphertext{cipher1, cipher2, cipher3, cipher4}, []uint64{1, 2, 3, 4})

		evk1 := participants[0].GetEvaluationKey()
		evk1.PeerID = 1
		evk2 := participants[1].GetEvaluationKey()
		evk2.PeerID = 2
		evk3 := participants[2].GetEvaluationKey()
		evk3.PeerID = 3
		evk4 := participants[3].GetEvaluationKey()
		evk4.PeerID = 4

		evalKeys := []*mkrlwe.MKEvaluationKey{evk1, evk2, evk3, evk4}

		pk1 := participants[0].GetPublicKey()
		pk1.PeerID = 1
		pk2 := participants[1].GetPublicKey()
		pk2.PeerID = 2
		pk3 := participants[2].GetPublicKey()
		pk3.PeerID = 3
		pk4 := participants[3].GetPublicKey()
		pk4.PeerID = 4

		publicKeys := []*mkrlwe.MKPublicKey{pk1, pk2, pk3, pk4}

		resCipher1 := evaluator.Mul(ciphers[0], ciphers[1])
		resCipher2 := evaluator.Mul(ciphers[2], ciphers[3])

		evaluator.RelinInPlace(resCipher1, evalKeys[:2], publicKeys[:2])
		evaluator.RelinInPlace(resCipher2, evalKeys[2:], publicKeys[2:])

		evaluator.Rescale(resCipher1, resCipher1)
		evaluator.Rescale(resCipher2, resCipher2)

		resCipher := evaluator.Mul(resCipher1, resCipher2)

		evaluator.RelinInPlace(resCipher, evalKeys, publicKeys)
		evaluator.DropLevel(resCipher, 1)
		evaluator.Rescale(resCipher, resCipher)

		resCKKS := evaluator.ConvertToCKKSCiphertext(resCipher)

		// decrypt
		partialDec1 := participants[0].GetPartialDecryption(resCKKS[0])
		partialDec2 := participants[1].GetPartialDecryption(resCKKS[1])
		partialDec3 := participants[2].GetPartialDecryption(resCKKS[2])
		partialDec4 := participants[3].GetPartialDecryption(resCKKS[3])

		decrypted := participants[0].Decrypt(resCKKS[0], []*ring.Poly{partialDec1, partialDec2, partialDec3, partialDec4})

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

		ciphers := evaluator.ConvertToMKCiphertext([]*ckks.Ciphertext{cipher1}, []uint64{1})

		for _, n := range rots {

			rotKey := participants[0].GetRotationKeys(n)
			rotKey.PeerID = 1

			values2 := utils.RotateComplex128Slice(values1, n)
			resCipher := evaluator.Rotate(ciphers[0], n, []*mkrlwe.MKEvalGalKey{rotKey})

			resCKKS := evaluator.ConvertToCKKSCiphertext(resCipher)

			partialDec := participants[0].GetPartialDecryption(resCKKS[0])

			decrypted := participants[0].Decrypt(resCKKS[0], []*ring.Poly{partialDec})
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
		ciphers := evaluator.ConvertToMKCiphertext([]*ckks.Ciphertext{cipher1, cipher2}, []uint64{1, 2})

		// add both ciphertexts
		added := evaluator.Add(ciphers[0], ciphers[1])

		// operate on plaintext space
		for i := int(0); i < len(values1); i++ {
			values1[i] += values2[i]
		}

		for _, n := range rots {

			// generate rotation keys
			rotKey1 := participants[0].GetRotationKeys(n)
			rotKey1.PeerID = 1
			rotKey2 := participants[1].GetRotationKeys(n)
			rotKey2.PeerID = 2

			rotated := utils.RotateComplex128Slice(values1, n)
			resCipher := evaluator.Rotate(added, n, []*mkrlwe.MKEvalGalKey{rotKey1, rotKey2})

			resCKKS := evaluator.ConvertToCKKSCiphertext(resCipher)

			partialDec1 := participants[0].GetPartialDecryption(resCKKS[0])
			partialDec2 := participants[1].GetPartialDecryption(resCKKS[1])

			decrypted := participants[0].Decrypt(resCKKS[0], []*ring.Poly{partialDec1, partialDec2})

			verifyTestVectors(params, rotated, decrypted, t)
		}
	})

}

/*
func testCkksMkbfvBridge(t *testing.T, params *ckks.Parameters) {

	encoder := ckks.NewEncoder(*params)
	keygen := ckks.NewKeyGenerator(*params)
	sk, pk := keygen.GenKeyPair()
	encryptorPK := ckks.NewEncryptorFromPk(*params, pk)

	t.Run(testString("Test Bridge BFV-MKBFV/", 2, params), func(t *testing.T) {

		// setup bfv environment and encrypt values in bfv
		values1 := newTestValue(params, complex(-1, -1), complex(1, 1))
		values2 := newTestValue(params, complex(-1, -1), complex(1, 1))
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
		part2 := newParticipant(params, a)

		// keygen
		keysPart2 := mkrlwe.KeyGenWithSecretKey(&params.Parameters, mkrlwe.CopyNewDecomposed(a), sk)

		// perform addition

		ciphertext2 := part2.Encrypt(values2)

		evaluator := NewMKEvaluator(params)

		ciphers := evaluator.ConvertToMKCiphertext([]*ckks.Ciphertext{ciphertext1, ciphertext2}, []uint64{1, 2})

		res := evaluator.Add(ciphers[0], ciphers[1])

		resCKKS := evaluator.ConvertToCKKSCiphertext(res)

		// decrypt
		decryptor := mkrlwe.NewMKDecryptor(&params.Parameters)
		partDec1 := decryptor.PartDec(&resCKKS[0].El().Element, resCKKS[0].Level(), keysPart2.SecretKey, 0.6)

		partDec2 := part2.GetPartialDecryption(resCKKS[1])

		decrypted := part2.Decrypt(resCKKS[1], []*ring.Poly{partDec1, partDec2})

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
		ciphers := evaluator.ConvertToMKCiphertext([]*ckks.Ciphertext{cipher1, cipher2}, []uint64{1, 2})
		res := evaluator.Add(ciphers[0], ciphers[1])

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
*/

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
func setupPeers(peerNbr uint64, params *ckks.Parameters, sigmaSmudging float64) []participant {

	res := make([]participant, peerNbr)

	prng, err := utils.NewKeyedPRNG([]byte{'l', 'a', 't', 't', 'i', 'g', 'o'})
	if err != nil {
		panic(err)
	}

	a := mkrlwe.GenCommonPublicParam(&params.Parameters, prng)
	for i := 0; i < int(peerNbr); i++ {
		res[i] = newParticipant(params, sigmaSmudging, a)
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

//--------------------------------Participants interface for tests------------------------------------------

// MKParticipant is a type for participants in a multi key ckks scheme
type participant interface {
	GetEvaluationKey() *mkrlwe.MKEvaluationKey
	GetPublicKey() *mkrlwe.MKPublicKey
	Encrypt(values []complex128) *ckks.Ciphertext
	Decrypt(cipher *ckks.Ciphertext, partialDecryptions []*ring.Poly) []complex128
	GetPartialDecryption(ciphertext *ckks.Ciphertext) *ring.Poly
	GetRotationKeys(rot int) *mkrlwe.MKEvalGalKey
	GetSecretKey() *mkrlwe.MKSecretKey
}

type mkParticipant struct {
	encryptor MKEncryptor
	decryptor mkrlwe.MKDecryptor
	keys      *mkrlwe.MKKeys
	encoder   ckks.Encoder
	ringQ     *ring.Ring
	params    *ckks.Parameters
}

// GetEvaluationKey returns the evaluation key of the participant
func (participant *mkParticipant) GetEvaluationKey() *mkrlwe.MKEvaluationKey {
	return participant.keys.EvalKey
}

// GetPublicKey returns the publik key of the participant
func (participant *mkParticipant) GetPublicKey() *mkrlwe.MKPublicKey {
	return participant.keys.PublicKey
}

// Encrypt constructs a ciphertext from the given values
func (participant *mkParticipant) Encrypt(values []complex128) *ckks.Ciphertext {
	if values == nil || len(values) <= 0 {
		panic("Cannot encrypt uninitialized or empty values")
	}
	return participant.encryptor.Encrypt(participant.encoder.EncodeNTTAtLvlNew(participant.params.MaxLevel(), values, participant.params.LogSlots()))
}

// Decrypt returns the decryption of the ciphertext given the partial decryption
func (participant *mkParticipant) Decrypt(cipher *ckks.Ciphertext, partialDecryptions []*ring.Poly) []complex128 {

	if cipher == nil || len(cipher.Value) < 1 {
		panic("Cannot decrypt uninitialized ciphertext nor ciphertext containing no value")
	}

	if partialDecryptions == nil || len(partialDecryptions) < 1 {
		panic("Decryption necessitates at least one partialy decrypted ciphertext")
	}

	decrypted := participant.decryptor.MergeDec(&cipher.Element.Element, cipher.Level(), partialDecryptions)

	pt := ckks.NewPlaintext(*participant.params, cipher.Level(), cipher.Scale())

	pt.SetValue(decrypted)

	return participant.encoder.Decode(pt, participant.params.LogSlots())
}

// GetPartialDecryption returns the partial decryption of a ckks ciphertext.
// Value[0] contains the part that belongs to the participant
func (participant *mkParticipant) GetPartialDecryption(ct *ckks.Ciphertext) *ring.Poly {

	if ct == nil || len(ct.Value) < 1 {
		panic("Uninitialized ciphertext")
	}
	return participant.decryptor.PartDec(&ct.Element.Element, ct.Level(), participant.keys.SecretKey)
}

// newParticipant creates a participant for the multi key ckks scheme
// the ckks parameters as well as the standard deviation used for partial decryption must be provided
func newParticipant(params *ckks.Parameters, sigmaSmudging float64, crs *mkrlwe.MKDecomposedPoly) participant {

	if crs == nil || params == nil {
		panic("Uninitialized parameters. Cannot create new participant")
	}

	if sigmaSmudging < params.Sigma() {
		panic("Sigma must be at least greater than the standard deviation of the gaussian distribution")
	}

	if len(crs.Poly) != int(params.Beta()) {
		panic("CRS must be the same dimention as returned by the function ckks.Parameters.Beta()")
	}

	keys := mkrlwe.KeyGen(&params.Parameters, mkrlwe.CopyNewDecomposed(crs))
	encryptor := NewMKEncryptor(keys.PublicKey, params)
	decryptor := mkrlwe.NewMKDecryptor(&params.Parameters, sigmaSmudging)
	encoder := ckks.NewEncoder(*params)
	ringQ := mkrlwe.GetRingQ(&params.Parameters)

	return &mkParticipant{
		encryptor: encryptor,
		decryptor: decryptor,
		keys:      keys,
		encoder:   encoder,
		ringQ:     ringQ,
		params:    params,
	}
}

// GetRotationKeys returns the rotation key set associated with the given rotation
func (participant *mkParticipant) GetRotationKeys(rot int) *mkrlwe.MKEvalGalKey {

	galEl := participant.params.GaloisElementForColumnRotationBy(rot)

	evalKey := mkrlwe.GaloisEvaluationKeyGen(galEl, participant.keys.SecretKey, &participant.params.Parameters)

	return evalKey
}

// GetSecretKey returns the secret key of the participant
func (participant *mkParticipant) GetSecretKey() *mkrlwe.MKSecretKey {
	return participant.keys.SecretKey
}
