package mkbfv

import (
	"testing"

	"github.com/ldsec/lattigo/v2/bfv"
	"github.com/ldsec/lattigo/v2/mkrlwe"
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/utils"
)

func Test_MKBFV(t *testing.T) {

	for i, paramLit := range bfv.DefaultParams {

		params, err := bfv.NewParametersFromLiteral(paramLit)
		if err != nil {
			panic(err)
		}
		p := &params

		testEncryptionEqualsDecryption(t, p)
		testAdd(t, p)
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
		//testBfvMkbfvBridge(t, p)

		if i != 0 && i != 4 && i != 6 {
			testMulFourParticipants(t, p)

		}
		if i != 6 {
			testMul(t, p)
			testAddAfterMul(t, p)
		}

		testRotation(t, p)
		testRotationTwoParticipants(t, p)
		//testMarshaler(t, p)
		testSquare(t, p)
	}
}

func testEncryptionEqualsDecryption(t *testing.T, params *bfv.Parameters) {

	sigma := 6.0

	participants := setupPeers(1, params, sigma)

	ringT := getRingT(params)

	t.Run(testString("Test encryption equals decryption/", 1, params), func(t *testing.T) {

		expected := getRandomPlaintextValue(ringT, params)

		// encrypt
		cipher := participants[0].Encrypt(expected)

		// decrypt
		partialDec := participants[0].GetPartialDecryption(cipher)
		decrypted := participants[0].Decrypt(cipher, []*ring.Poly{partialDec})

		// decode and check
		if !equalsSlice(decrypted, expected) {
			t.Error("Decryption of encryption does not equals plaintext")
		}
	})

}

func testAdd(t *testing.T, params *bfv.Parameters) {

	sigma := 6.0

	participants := setupPeers(2, params, sigma)

	ringT := getRingT(params)

	t.Run(testString("Test Addition/", 2, params), func(t *testing.T) {
		expected1 := getRandomPlaintextValue(ringT, params)
		expected2 := getRandomPlaintextValue(ringT, params)

		// encrypt
		cipher1 := participants[0].Encrypt(expected1)
		cipher2 := participants[1].Encrypt(expected2)

		// pad and add
		evaluator := NewMKEvaluator(params)
		ciphers := evaluator.ConvertToMKCiphertext([]*bfv.Ciphertext{cipher1, cipher2}, []uint64{1, 2})
		resCipher := evaluator.Add(ciphers[0], ciphers[1])

		resBFV := evaluator.ConvertToBFVCiphertext(resCipher)

		// decrypt
		partialDec1 := participants[0].GetPartialDecryption(resBFV[0])
		partialDec2 := participants[1].GetPartialDecryption(resBFV[1])

		decrypted := participants[0].Decrypt(resBFV[0], []*ring.Poly{partialDec1, partialDec2})

		// perform the operation in the plaintext space
		expected := ringT.NewPoly()
		p1 := ringT.NewPoly()
		p2 := ringT.NewPoly()
		copy(p1.Coeffs[0], expected1)
		copy(p2.Coeffs[0], expected2)

		ringT.Add(p1, p2, expected)

		if !equalsSlice(decrypted, expected.Coeffs[0]) {
			t.Error("Homomorphic addition error")
		}

	})

}

func testSub(t *testing.T, params *bfv.Parameters) {

	sigma := 6.0

	participants := setupPeers(2, params, sigma)

	ringT := getRingT(params)

	t.Run(testString("Test Subtraction/", 2, params), func(t *testing.T) {
		expected1 := getRandomPlaintextValue(ringT, params)
		expected2 := getRandomPlaintextValue(ringT, params)

		// encrypt
		cipher1 := participants[0].Encrypt(expected1)
		cipher2 := participants[1].Encrypt(expected2)

		// pad and add
		evaluator := NewMKEvaluator(params)
		ciphers := evaluator.ConvertToMKCiphertext([]*bfv.Ciphertext{cipher1, cipher2}, []uint64{1, 2})
		resCipher := evaluator.Sub(ciphers[0], ciphers[1])

		resBFV := evaluator.ConvertToBFVCiphertext(resCipher)

		// decrypt
		partialDec1 := participants[0].GetPartialDecryption(resBFV[0])
		partialDec2 := participants[1].GetPartialDecryption(resBFV[1])

		decrypted := participants[0].Decrypt(resBFV[0], []*ring.Poly{partialDec1, partialDec2})

		// perform the operation in the plaintext space
		expected := ringT.NewPoly()
		p1 := ringT.NewPoly()
		p2 := ringT.NewPoly()
		copy(p1.Coeffs[0], expected1)
		copy(p2.Coeffs[0], expected2)

		ringT.Sub(p1, p2, expected)

		if !equalsSlice(decrypted, expected.Coeffs[0]) {
			t.Error("Homomorphic subtraction error")
		}

	})

}

func testAddPlaintext(t *testing.T, params *bfv.Parameters) {

	sigma := 6.0

	participants := setupPeers(1, params, sigma)

	ringT := getRingT(params)

	t.Run(testString("Test Plaintext Addition/", 1, params), func(t *testing.T) {

		value1 := getRandomPlaintextValue(ringT, params)
		value2 := getRandomPlaintextValue(ringT, params)

		// encrypt
		cipher := participants[0].Encrypt(value1)

		// add plaintext to ciphertext
		evaluator := NewMKEvaluator(params)
		ciphers := evaluator.ConvertToMKCiphertext([]*bfv.Ciphertext{cipher}, []uint64{1})

		pt := evaluator.NewPlaintextFromValue(value2)
		resCipher := evaluator.AddPlaintext(pt, ciphers[0])

		resBFV := evaluator.ConvertToBFVCiphertext(resCipher)

		// decrypt
		partialDec1 := participants[0].GetPartialDecryption(resBFV[0])

		decrypted := participants[0].Decrypt(resBFV[0], []*ring.Poly{partialDec1})

		// perform the operation in the plaintext space
		expected := ringT.NewPoly()
		p1 := ringT.NewPoly()
		p2 := ringT.NewPoly()
		copy(p1.Coeffs[0], value1)
		copy(p2.Coeffs[0], value2)

		ringT.Add(p1, p2, expected)

		if !equalsSlice(decrypted, expected.Coeffs[0]) {
			t.Error("Homomorphic addition error")
		}

	})

}

func testAddPlaintextTwoParticipants(t *testing.T, params *bfv.Parameters) {

	sigma := 6.0

	participants := setupPeers(2, params, sigma)

	ringT := getRingT(params)

	t.Run(testString("Test Plaintext Addition/", 2, params), func(t *testing.T) {

		value1 := getRandomPlaintextValue(ringT, params)
		value2 := getRandomPlaintextValue(ringT, params)
		value3 := getRandomPlaintextValue(ringT, params)

		// encrypt
		cipher1 := participants[0].Encrypt(value1)
		cipher2 := participants[1].Encrypt(value2)

		// add plaintext to one of the ciphertext then add both ciphertexts
		evaluator := NewMKEvaluator(params)
		ciphers := evaluator.ConvertToMKCiphertext([]*bfv.Ciphertext{cipher1, cipher2}, []uint64{1, 2})

		pt := evaluator.NewPlaintextFromValue(value3)
		resCipher1 := evaluator.AddPlaintext(pt, ciphers[0])

		resCipher := evaluator.Add(resCipher1, ciphers[1])

		resBFV := evaluator.ConvertToBFVCiphertext(resCipher)

		// decrypt
		partialDec1 := participants[0].GetPartialDecryption(resBFV[0])
		partialDec2 := participants[1].GetPartialDecryption(resBFV[1])

		decrypted := participants[0].Decrypt(resBFV[0], []*ring.Poly{partialDec1, partialDec2})

		// perform the operation in the plaintext space
		expected := ringT.NewPoly()

		p1 := ringT.NewPoly()
		p2 := ringT.NewPoly()
		p3 := ringT.NewPoly()
		copy(p1.Coeffs[0], value1)
		copy(p2.Coeffs[0], value2)
		copy(p3.Coeffs[0], value3)

		ringT.Add(p1, p2, expected)
		ringT.Add(expected, p3, expected)

		if !equalsSlice(decrypted, expected.Coeffs[0]) {
			t.Error("Homomorphic addition error")
		}

	})

}

func testSubPlaintext(t *testing.T, params *bfv.Parameters) {

	sigma := 6.0

	participants := setupPeers(1, params, sigma)

	ringT := getRingT(params)

	t.Run(testString("Test Plaintext Subtraction/", 1, params), func(t *testing.T) {

		value1 := getRandomPlaintextValue(ringT, params)
		value2 := getRandomPlaintextValue(ringT, params)

		// encrypt
		cipher := participants[0].Encrypt(value1)

		// sub plaintext to ciphertext
		evaluator := NewMKEvaluator(params)

		ciphers := evaluator.ConvertToMKCiphertext([]*bfv.Ciphertext{cipher}, []uint64{1})

		pt := evaluator.NewPlaintextFromValue(value2)
		resCipher := evaluator.SubPlaintext(pt, ciphers[0])

		resBFV := evaluator.ConvertToBFVCiphertext(resCipher)

		// decrypt
		partialDec1 := participants[0].GetPartialDecryption(resBFV[0])

		decrypted := participants[0].Decrypt(resBFV[0], []*ring.Poly{partialDec1})

		// perform the operation in the plaintext space
		expected := ringT.NewPoly()
		p1 := ringT.NewPoly()
		p2 := ringT.NewPoly()
		copy(p1.Coeffs[0], value1)
		copy(p2.Coeffs[0], value2)

		ringT.Sub(p1, p2, expected)

		if !equalsSlice(decrypted, expected.Coeffs[0]) {
			t.Error("Homomorphic subtraction error")
		}

	})

}

func testSubPlaintextTwoParticipants(t *testing.T, params *bfv.Parameters) {

	sigma := 6.0

	participants := setupPeers(2, params, sigma)

	ringT := getRingT(params)

	t.Run(testString("Test Plaintext Subtraction/", 2, params), func(t *testing.T) {

		value1 := getRandomPlaintextValue(ringT, params)
		value2 := getRandomPlaintextValue(ringT, params)
		value3 := getRandomPlaintextValue(ringT, params)

		// encrypt
		cipher1 := participants[0].Encrypt(value1)
		cipher2 := participants[1].Encrypt(value2)

		// sub plaintext to one of the ciphertext then sub both ciphertexts
		evaluator := NewMKEvaluator(params)

		ciphers := evaluator.ConvertToMKCiphertext([]*bfv.Ciphertext{cipher1, cipher2}, []uint64{1, 2})

		pt := evaluator.NewPlaintextFromValue(value3)
		resCipher1 := evaluator.SubPlaintext(pt, ciphers[0])

		resCipher := evaluator.Sub(resCipher1, ciphers[1])

		resBFV := evaluator.ConvertToBFVCiphertext(resCipher)

		// decrypt
		partialDec1 := participants[0].GetPartialDecryption(resBFV[0])
		partialDec2 := participants[1].GetPartialDecryption(resBFV[1])

		decrypted := participants[0].Decrypt(resBFV[0], []*ring.Poly{partialDec1, partialDec2})

		// perform the operation in the plaintext space
		expected := ringT.NewPoly()

		p1 := ringT.NewPoly()
		p2 := ringT.NewPoly()
		p3 := ringT.NewPoly()
		copy(p1.Coeffs[0], value1)
		copy(p2.Coeffs[0], value2)
		copy(p3.Coeffs[0], value3)

		ringT.Sub(p1, p2, expected)
		ringT.Sub(expected, p3, expected)

		if !equalsSlice(decrypted, expected.Coeffs[0]) {
			t.Error("Homomorphic addition error")
		}

	})

}

func testNeg(t *testing.T, params *bfv.Parameters) {

	sigma := 6.0

	participants := setupPeers(1, params, sigma)

	ringT := getRingT(params)

	t.Run(testString("Test Negation/", 1, params), func(t *testing.T) {

		value1 := getRandomPlaintextValue(ringT, params)

		// encrypt
		cipher := participants[0].Encrypt(value1)

		// add with negated ciphertext
		evaluator := NewMKEvaluator(params)

		ciphers := evaluator.ConvertToMKCiphertext([]*bfv.Ciphertext{cipher}, []uint64{1})

		resCipher := evaluator.Add(evaluator.Neg(ciphers[0]), ciphers[0])

		resBFV := evaluator.ConvertToBFVCiphertext(resCipher)

		// decrypt
		partialDec1 := participants[0].GetPartialDecryption(resBFV[0])

		decrypted := participants[0].Decrypt(resBFV[0], []*ring.Poly{partialDec1})

		// should be 0
		expected := ringT.NewPoly()

		if !equalsSlice(decrypted, expected.Coeffs[0]) {
			t.Error("Homomorphic negation error")
		}

	})
}

func testNegTwoParticipants(t *testing.T, params *bfv.Parameters) {

	sigma := 6.0

	participants := setupPeers(2, params, sigma)

	ringT := getRingT(params)

	t.Run(testString("Test Negation/", 2, params), func(t *testing.T) {

		value1 := getRandomPlaintextValue(ringT, params)
		value2 := getRandomPlaintextValue(ringT, params)

		// encrypt
		cipher1 := participants[0].Encrypt(value1)
		cipher2 := participants[1].Encrypt(value2)

		// add with negated ciphertext
		evaluator := NewMKEvaluator(params)

		ciphers := evaluator.ConvertToMKCiphertext([]*bfv.Ciphertext{cipher1, cipher2}, []uint64{1, 2})

		added := evaluator.Add(ciphers[0], ciphers[1])

		resCipher := evaluator.Neg(added)

		resBFV := evaluator.ConvertToBFVCiphertext(resCipher)

		// decrypt
		partialDec1 := participants[0].GetPartialDecryption(resBFV[0])
		partialDec2 := participants[1].GetPartialDecryption(resBFV[1])

		decrypted := participants[0].Decrypt(resBFV[0], []*ring.Poly{partialDec1, partialDec2})

		// perform plaintext operation
		expected := ringT.NewPoly()
		p1 := ringT.NewPoly()
		p2 := ringT.NewPoly()
		copy(p1.Coeffs[0], value1)
		copy(p2.Coeffs[0], value2)

		ringT.Add(p1, p2, expected)
		ringT.Neg(expected, expected)
		ringT.Reduce(expected, expected)

		if !equalsSlice(decrypted, expected.Coeffs[0]) {
			t.Error("Homomorphic negation error")
		}

	})
}

func testAddFourParticipants(t *testing.T, params *bfv.Parameters) {

	sigma := 6.0

	participants := setupPeers(4, params, sigma)

	ringT := getRingT(params)

	t.Run(testString("Test Addition/", 4, params), func(t *testing.T) {
		expected1 := getRandomPlaintextValue(ringT, params)
		expected2 := getRandomPlaintextValue(ringT, params)
		expected3 := getRandomPlaintextValue(ringT, params)
		expected4 := getRandomPlaintextValue(ringT, params)

		// encrypt
		cipher1 := participants[0].Encrypt(expected1)
		cipher2 := participants[1].Encrypt(expected2)
		cipher3 := participants[2].Encrypt(expected3)
		cipher4 := participants[3].Encrypt(expected4)

		// pad and add in 2 steps
		evaluator := NewMKEvaluator(params)
		ciphers := evaluator.ConvertToMKCiphertext([]*bfv.Ciphertext{cipher1, cipher2, cipher3, cipher4}, []uint64{1, 2, 3, 4})

		resCipher1 := evaluator.Add(ciphers[0], ciphers[1])
		resCipher2 := evaluator.Add(ciphers[2], ciphers[3])

		resCipher := evaluator.Add(resCipher1, resCipher2)

		resBFV := evaluator.ConvertToBFVCiphertext(resCipher)

		// decrypt
		partialDec1 := participants[0].GetPartialDecryption(resBFV[0])
		partialDec2 := participants[1].GetPartialDecryption(resBFV[1])
		partialDec3 := participants[2].GetPartialDecryption(resBFV[2])
		partialDec4 := participants[3].GetPartialDecryption(resBFV[3])

		decrypted := participants[0].Decrypt(resBFV[0], []*ring.Poly{partialDec1, partialDec2, partialDec3, partialDec4})

		// perform the operation in the plaintext space
		expected := ringT.NewPoly()
		p1 := ringT.NewPoly()
		p2 := ringT.NewPoly()
		p3 := ringT.NewPoly()
		p4 := ringT.NewPoly()
		copy(p1.Coeffs[0], expected1)
		copy(p2.Coeffs[0], expected2)
		copy(p3.Coeffs[0], expected3)
		copy(p4.Coeffs[0], expected4)

		ringT.Add(p1, p2, expected)
		ringT.Add(p3, expected, expected)
		ringT.Add(p4, expected, expected)

		if !equalsSlice(decrypted, expected.Coeffs[0]) {
			t.Error("Homomorphic addition error")
		}
	})

}

func testMulPlaintext(t *testing.T, params *bfv.Parameters) {

	sigma := 6.0

	participants := setupPeers(1, params, sigma)

	ringT := getRingT(params)

	t.Run(testString("Test Plaintext Multiplication/", 1, params), func(t *testing.T) {

		value1 := getRandomPlaintextValue(ringT, params)
		value2 := getRandomPlaintextValue(ringT, params)

		// encrypt
		cipher := participants[0].Encrypt(value1)

		// multiply plaintext and ciphertext
		evaluator := NewMKEvaluator(params)
		ciphers := evaluator.ConvertToMKCiphertext([]*bfv.Ciphertext{cipher}, []uint64{1})

		pt := evaluator.NewPlaintextMulFromValue(value2)
		resCipher := evaluator.MultPlaintext(pt, ciphers[0])

		resBFV := evaluator.ConvertToBFVCiphertext(resCipher)

		// decrypt
		partialDec1 := participants[0].GetPartialDecryption(resBFV[0])

		decrypted := participants[0].Decrypt(resBFV[0], []*ring.Poly{partialDec1})

		// perform the operation in the plaintext space
		expected := ringT.NewPoly()
		p1 := ringT.NewPoly()
		p2 := ringT.NewPoly()
		copy(p1.Coeffs[0], value1)
		copy(p2.Coeffs[0], value2)

		ringT.MulCoeffs(p1, p2, expected)

		if !equalsSlice(decrypted, expected.Coeffs[0]) {
			t.Error("Homomorphic multiplication error")
		}

	})

}

func testSquare(t *testing.T, params *bfv.Parameters) {

	sigma := 6.0

	participants := setupPeers(1, params, sigma)

	ringT := getRingT(params)

	t.Run(testString("Test Square/", 1, params), func(t *testing.T) {

		value1 := getRandomPlaintextValue(ringT, params)

		evk := participants[0].GetEvaluationKey()
		evk.PeerID = 1
		evalKeys := []*mkrlwe.MKEvaluationKey{evk}

		pk := participants[0].GetPublicKey()
		pk.PeerID = 1
		publicKeys := []*mkrlwe.MKPublicKey{pk}

		// encrypt
		cipher1 := participants[0].Encrypt(value1)

		// square ciphertexts
		evaluator := NewMKEvaluator(params)
		ciphers := evaluator.ConvertToMKCiphertext([]*bfv.Ciphertext{cipher1}, []uint64{1})

		resCipher := evaluator.Mul(ciphers[0], ciphers[0])
		evaluator.RelinInPlace(resCipher, evalKeys, publicKeys)

		resBFV := evaluator.ConvertToBFVCiphertext(resCipher)

		// decrypt
		partialDec1 := participants[0].GetPartialDecryption(resBFV[0])
		decrypted := participants[0].Decrypt(resBFV[0], []*ring.Poly{partialDec1})

		// perform the operation in the plaintext space
		p1 := ringT.NewPoly()
		copy(p1.Coeffs[0], value1)

		ringT.MulCoeffs(p1, p1, p1)

		if !equalsSlice(decrypted, p1.Coeffs[0]) {
			t.Error("Homomorphic multiplication error")
		}

	})

}

func testMulPlaintextTwoParticipants(t *testing.T, params *bfv.Parameters) {

	sigma := 6.0

	participants := setupPeers(2, params, sigma)

	ringT := getRingT(params)

	t.Run(testString("Test Plaintext Multiplication/", 2, params), func(t *testing.T) {

		value1 := getRandomPlaintextValue(ringT, params)
		value2 := getRandomPlaintextValue(ringT, params)
		value3 := getRandomPlaintextValue(ringT, params)

		// encrypt
		cipher1 := participants[0].Encrypt(value1)
		cipher2 := participants[1].Encrypt(value2)

		// add both ciphertexts then multiply by plaintext
		evaluator := NewMKEvaluator(params)
		ciphers := evaluator.ConvertToMKCiphertext([]*bfv.Ciphertext{cipher1, cipher2}, []uint64{1, 2})

		resCipherTMP := evaluator.Add(ciphers[0], ciphers[1])

		pt := evaluator.NewPlaintextMulFromValue(value3)
		resCipher := evaluator.MultPlaintext(pt, resCipherTMP)

		resBFV := evaluator.ConvertToBFVCiphertext(resCipher)

		// decrypt
		partialDec1 := participants[0].GetPartialDecryption(resBFV[0])
		partialDec2 := participants[1].GetPartialDecryption(resBFV[1])

		decrypted := participants[0].Decrypt(resBFV[0], []*ring.Poly{partialDec1, partialDec2})

		// perform the operation in the plaintext space
		expected := ringT.NewPoly()

		p1 := ringT.NewPoly()
		p2 := ringT.NewPoly()
		p3 := ringT.NewPoly()
		copy(p1.Coeffs[0], value1)
		copy(p2.Coeffs[0], value2)
		copy(p3.Coeffs[0], value3)

		ringT.Add(p1, p2, expected)
		ringT.MulCoeffs(expected, p3, expected)

		if !equalsSlice(decrypted, expected.Coeffs[0]) {
			t.Error("Homomorphic multiplication error")
		}

	})

}

func testMul(t *testing.T, params *bfv.Parameters) {

	sigma := 6.0

	participants := setupPeers(2, params, sigma)

	ringT := getRingT(params)

	t.Run(testString("Test Multiplication/", 2, params), func(t *testing.T) {

		expected1 := getRandomPlaintextValue(ringT, params)
		expected2 := getRandomPlaintextValue(ringT, params)

		// encrypt
		cipher1 := participants[0].Encrypt(expected1)
		cipher2 := participants[1].Encrypt(expected2)

		// multiply using evaluation keys and public keys
		evaluator := NewMKEvaluator(params)

		ciphers := evaluator.ConvertToMKCiphertext([]*bfv.Ciphertext{cipher1, cipher2}, []uint64{1, 2})

		evk1 := participants[0].GetEvaluationKey()
		evk2 := participants[1].GetEvaluationKey()
		evk1.PeerID = 1
		evk2.PeerID = 2
		evalKeys := []*mkrlwe.MKEvaluationKey{evk1, evk2}

		pk1 := participants[0].GetPublicKey()
		pk2 := participants[1].GetPublicKey()
		pk1.PeerID = 1
		pk2.PeerID = 2
		publicKeys := []*mkrlwe.MKPublicKey{pk1, pk2}

		resCipher := evaluator.Mul(ciphers[0], ciphers[1])

		evaluator.RelinInPlace(resCipher, evalKeys, publicKeys)

		resBFV := evaluator.ConvertToBFVCiphertext(resCipher)

		// decrypt
		partialDec1 := participants[0].GetPartialDecryption(resBFV[0])
		partialDec2 := participants[1].GetPartialDecryption(resBFV[1])

		decrypted := participants[0].Decrypt(resBFV[0], []*ring.Poly{partialDec1, partialDec2})

		// perform the operation in the plaintext space
		expected := ringT.NewPoly()
		p1 := ringT.NewPoly()
		p2 := ringT.NewPoly()
		copy(p1.Coeffs[0], expected1)
		copy(p2.Coeffs[0], expected2)
		ringT.MulCoeffs(p1, p2, expected)

		if !equalsSlice(decrypted, expected.Coeffs[0]) {
			t.Error("Homomorphic multiplication error")
		}
	})

}

func testAddAfterMul(t *testing.T, params *bfv.Parameters) {

	sigma := 6.0

	participants := setupPeers(2, params, sigma)

	ringT := getRingT(params)

	t.Run(testString("Test Addition after Multiplication/", 2, params), func(t *testing.T) {

		expected1 := getRandomPlaintextValue(ringT, params)
		expected2 := getRandomPlaintextValue(ringT, params)

		// encrypt
		cipher1 := participants[0].Encrypt(expected1)
		cipher2 := participants[1].Encrypt(expected2)

		// pad
		evaluator := NewMKEvaluator(params)

		ciphers := evaluator.ConvertToMKCiphertext([]*bfv.Ciphertext{cipher1, cipher2}, []uint64{1, 2})

		evk1 := participants[0].GetEvaluationKey()
		evk2 := participants[1].GetEvaluationKey()
		evk1.PeerID = 1
		evk2.PeerID = 2
		evalKeys := []*mkrlwe.MKEvaluationKey{evk1, evk2}

		pk1 := participants[0].GetPublicKey()
		pk2 := participants[1].GetPublicKey()
		pk1.PeerID = 1
		pk2.PeerID = 2
		publicKeys := []*mkrlwe.MKPublicKey{pk1, pk2}

		resCipher := evaluator.Mul(ciphers[0], ciphers[1])

		evaluator.RelinInPlace(resCipher, evalKeys, publicKeys)

		resCipher = evaluator.Add(resCipher, ciphers[0])

		resBFV := evaluator.ConvertToBFVCiphertext(resCipher)

		// decrypt
		partialDec1 := participants[0].GetPartialDecryption(resBFV[0])
		partialDec2 := participants[1].GetPartialDecryption(resBFV[1])

		decrypted := participants[0].Decrypt(resBFV[0], []*ring.Poly{partialDec1, partialDec2})

		// perform the operation in the plaintext space
		expected := ringT.NewPoly()
		p1 := ringT.NewPoly()
		p2 := ringT.NewPoly()
		copy(p1.Coeffs[0], expected1)
		copy(p2.Coeffs[0], expected2)
		ringT.MulCoeffs(p1, p2, expected)
		ringT.Add(expected, p1, expected)
		if !equalsSlice(decrypted, expected.Coeffs[0]) {
			t.Error("Homomorphic multiplication error")
		}
	})

}

func testMulFourParticipants(t *testing.T, params *bfv.Parameters) {

	sigma := 6.0

	participants := setupPeers(4, params, sigma)

	ringT := getRingT(params)

	t.Run(testString("Test Multiplication/", 4, params), func(t *testing.T) {

		expected1 := getRandomPlaintextValue(ringT, params)
		expected2 := getRandomPlaintextValue(ringT, params)
		expected3 := getRandomPlaintextValue(ringT, params)
		expected4 := getRandomPlaintextValue(ringT, params)

		// encrypt
		cipher1 := participants[0].Encrypt(expected1)
		cipher2 := participants[1].Encrypt(expected2)
		cipher3 := participants[2].Encrypt(expected3)
		cipher4 := participants[3].Encrypt(expected4)

		// pad and multiply in 2 steps
		evaluator := NewMKEvaluator(params)

		ciphers := evaluator.ConvertToMKCiphertext([]*bfv.Ciphertext{cipher1, cipher2, cipher3, cipher4}, []uint64{1, 2, 3, 4})

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

		resCipher := evaluator.Mul(resCipher1, resCipher2)

		evaluator.RelinInPlace(resCipher, evalKeys, publicKeys)

		resBFV := evaluator.ConvertToBFVCiphertext(resCipher)

		// decrypt
		partialDec1 := participants[0].GetPartialDecryption(resBFV[0])
		partialDec2 := participants[1].GetPartialDecryption(resBFV[1])
		partialDec3 := participants[2].GetPartialDecryption(resBFV[2])
		partialDec4 := participants[3].GetPartialDecryption(resBFV[3])

		decrypted := participants[0].Decrypt(resBFV[0], []*ring.Poly{partialDec1, partialDec2, partialDec3, partialDec4})

		// perform the operation in the plaintext space
		expected := ringT.NewPoly()
		p1 := ringT.NewPoly()
		p2 := ringT.NewPoly()
		p3 := ringT.NewPoly()
		p4 := ringT.NewPoly()
		copy(p1.Coeffs[0], expected1)
		copy(p2.Coeffs[0], expected2)
		copy(p3.Coeffs[0], expected3)
		copy(p4.Coeffs[0], expected4)

		ringT.MulCoeffs(p1, p2, expected)
		ringT.MulCoeffs(p3, expected, expected)
		ringT.MulCoeffs(p4, expected, expected)

		if !equalsSlice(decrypted, expected.Coeffs[0]) {
			t.Error("Homomorphic multiplication error")
		}
	})

}

func testRotation(t *testing.T, params *bfv.Parameters) {
	sigma := 6.0

	participants := setupPeers(1, params, sigma)

	ringT := getRingT(params)

	rots := []int{1, -1, 4, -4, 63, -63}

	t.Run(testString("Test Rotation/", 1, params), func(t *testing.T) {

		// generate test values
		values1 := getRandomPlaintextValue(ringT, params)

		// Encrypt
		cipher1 := participants[0].Encrypt(values1)

		evaluator := NewMKEvaluator(params)
		ciphers := evaluator.ConvertToMKCiphertext([]*bfv.Ciphertext{cipher1}, []uint64{1})

		for _, n := range rots {

			rotKey := participants[0].GetRotationKeys(n)
			rotKey.PeerID = 1

			resCipher := evaluator.Rotate(ciphers[0], n, []*mkrlwe.MKEvalGalKey{rotKey})

			resBFV := evaluator.ConvertToBFVCiphertext(resCipher)

			partialDec := participants[0].GetPartialDecryption(resBFV[0])
			decrypted := participants[0].Decrypt(resBFV[0], []*ring.Poly{partialDec})

			// perform the operation in the plaintext space
			expected := ringT.NewPoly()
			copy(expected.Coeffs[0], values1)

			nColumns := params.N() >> 1
			valuesWant := append(utils.RotateUint64Slice(expected.Coeffs[0][:nColumns], n), utils.RotateUint64Slice(expected.Coeffs[0][nColumns:], n)...)

			if !utils.EqualSliceUint64(valuesWant, decrypted) {
				t.Errorf("Rotation error")
			}
		}
	})

}

func testRotationTwoParticipants(t *testing.T, params *bfv.Parameters) {
	sigma := 6.0

	participants := setupPeers(2, params, sigma)

	ringT := getRingT(params)

	rots := []int{1, -1, 4, -4, 63, -63}

	t.Run(testString("Test Rotation/", 2, params), func(t *testing.T) {

		// generate test values
		values1 := getRandomPlaintextValue(ringT, params)
		values2 := getRandomPlaintextValue(ringT, params)

		// Encrypt
		cipher1 := participants[0].Encrypt(values1)
		cipher2 := participants[1].Encrypt(values2)

		evaluator := NewMKEvaluator(params)

		ciphers := evaluator.ConvertToMKCiphertext([]*bfv.Ciphertext{cipher1, cipher2}, []uint64{1, 2})

		// add both ciphertexts
		added := evaluator.Add(ciphers[0], ciphers[1])

		for _, n := range rots {

			rotKey1 := participants[0].GetRotationKeys(n)
			rotKey1.PeerID = 1
			rotKey2 := participants[1].GetRotationKeys(n)
			rotKey2.PeerID = 2

			resCipher := evaluator.Rotate(added, n, []*mkrlwe.MKEvalGalKey{rotKey1, rotKey2})

			resBFV := evaluator.ConvertToBFVCiphertext(resCipher)

			partialDec1 := participants[0].GetPartialDecryption(resBFV[0])
			partialDec2 := participants[1].GetPartialDecryption(resBFV[1])

			decrypted := participants[0].Decrypt(resBFV[0], []*ring.Poly{partialDec1, partialDec2})

			// perform the operation in the plaintext space
			expected := ringT.NewPoly()
			p1 := ringT.NewPoly()
			copy(p1.Coeffs[0], values1)
			p2 := ringT.NewPoly()
			copy(p2.Coeffs[0], values2)

			ringT.Add(p1, p2, expected)
			nColumns := params.N() >> 1
			valuesWant := append(utils.RotateUint64Slice(expected.Coeffs[0][:nColumns], n), utils.RotateUint64Slice(expected.Coeffs[0][nColumns:], n)...)
			if !utils.EqualSliceUint64(valuesWant, decrypted) {
				t.Errorf("Rotation error")
			}
		}
	})

}

/*
func testBfvMkbfvBridge(t *testing.T, params *bfv.Parameters) {

	ringT := getRingT(params)
	encoder := bfv.NewEncoder(*params)
	keygen := bfv.NewKeyGenerator(*params)
	sk, pk := keygen.GenKeyPair()
	encryptorPK := bfv.NewEncryptorFromPk(*params, pk)

	t.Run(testString("Test Bridge BFV-MKBFV/", 2, params), func(t *testing.T) {

		// setup bfv environment and encrypt values in bfv
		values1 := getRandomPlaintextValue(ringT, params)
		plaintext := bfv.NewPlaintext(*params)
		encoder.EncodeUint(values1, plaintext)

		var ciphertext1 *bfv.Ciphertext

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
		part2 := newParticipant(params, 6.0, a)

		decryptor := mkrlwe.NewMKDecryptor(&params.Parameters, 6.0)
		keys := mkrlwe.KeyGenWithSecretKey(&params.Parameters, a, sk)

		// perform addition
		values2 := getRandomPlaintextValue(ringT, params)
		ciphertext2 := part2.Encrypt(values2)

		evaluator := NewMKEvaluator(params)

		ciphers := evaluator.ConvertToMKCiphertext([]*bfv.Ciphertext{ciphertext1, ciphertext2}, []uint64{1, 2})

		res := evaluator.Add(ciphers[0], ciphers[1])

		resBFV := evaluator.ConvertToBFVCiphertext(res)

		// decrypt
		partDec1 := decryptor.PartDec(resBFV[0].Element, resBFV[0].Level(), keys.SecretKey)
		partDec2 := part2.GetPartialDecryption(resBFV[1])

		decrypted := part2.Decrypt(resBFV[1], []*ring.Poly{partDec1, partDec2})

		//verify
		expected := ringT.NewPoly()
		p1 := ringT.NewPoly()
		p2 := ringT.NewPoly()
		copy(p1.Coeffs[0], values1)
		copy(p2.Coeffs[0], values2)

		ringT.Add(p1, p2, expected)

		if !equalsSlice(decrypted, expected.Coeffs[0]) {
			t.Error("Homomorphic addition error")
		}
	})
}

func testMarshaler(t *testing.T, params *bfv.Parameters) {

	sigma := 6.0

	ringT := getRingT(params)

	participants := setupPeers(2, params, sigma)

	t.Run(testString("Test Marshaler/", 2, params), func(t *testing.T) {

		// get random value
		value1 := getRandomPlaintextValue(ringT, params)
		value2 := getRandomPlaintextValue(ringT, params)

		//encrypt
		cipher1 := participants[0].Encrypt(value1)
		cipher2 := participants[1].Encrypt(value2)

		// add
		evaluator := NewMKEvaluator(params)
		ciphers := evaluator.ConvertToMKCiphertext([]*bfv.Ciphertext{cipher1, cipher2}, []uint64{1, 2})
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

	if !equalsSlice(expected, res) {
		t.Errorf("MergeSlices method failed test")
	}

}

// equalsSlice returns true if both slices are equal
func equalsSlice(s1, s2 []uint64) bool {

	if len(s1) != len(s2) {
		return false
	}

	for i, e := range s1 {
		if e != s2[i] {

			return false
		}
	}

	return true
}

// equalsPoly returns true if both polynomials are equal
func equalsPoly(p1 *ring.Poly, p2 *ring.Poly) bool {

	if len(p1.Coeffs) != len(p2.Coeffs) {
		return false
	}

	for i, e := range p1.Coeffs {

		if !equalsSlice(e, p2.Coeffs[i]) {
			return false
		}
	}

	return true
}

// Generates keys for a set of peers identified by their peerID using a certain bfv parameter index
// returns the slice of keys with the bfv parameters
func setupPeers(peersNbr uint64, params *bfv.Parameters, sigmaSmudging float64) []participant {

	res := make([]participant, peersNbr)

	prng, err := utils.NewKeyedPRNG([]byte{'l', 'a', 't', 't', 'i', 'g', 'o'})

	if err != nil {
		panic(err)
	}

	// setup keys and public parameters
	a := mkrlwe.GenCommonPublicParam(&params.Parameters, prng)

	for i := 0; i < int(peersNbr); i++ {

		res[i] = newParticipant(params, sigmaSmudging, a)

	}

	return res
}

// returns a uniformly random slice of uint64 in RingT
func getRandomPlaintextValue(ringT *ring.Ring, params *bfv.Parameters) []uint64 {

	return mkrlwe.GetRandomPoly(&params.Parameters, ringT).Coeffs[0]
}

// getRingT returns the ring from which the plaintexts will be sampled
func getRingT(params *bfv.Parameters) *ring.Ring {

	if ringT, err := ring.NewRing(params.N(), []uint64{params.T()}); err != nil {
		panic("Couldn't create ringT with given parameters")

	} else {
		return ringT
	}
}

//--------------------------------Participants interface for tests------------------------------------------

// MKParticipant is a type for participants in a multy key bfv scheme
type participant interface {
	GetEvaluationKey() *mkrlwe.MKEvaluationKey
	GetPublicKey() *mkrlwe.MKPublicKey
	Encrypt(values []uint64) *bfv.Ciphertext
	Decrypt(cipher *bfv.Ciphertext, partialDecryptions []*ring.Poly) []uint64
	GetPartialDecryption(ciphertext *bfv.Ciphertext) *ring.Poly
	GetRotationKeys(rot int) *mkrlwe.MKEvalGalKey
}

type mkParticipant struct {
	encryptor MKEncryptor
	decryptor mkrlwe.MKDecryptor
	params    *bfv.Parameters
	keys      *mkrlwe.MKKeys
	encoder   bfv.Encoder
	ringQ     *ring.Ring
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
func (participant *mkParticipant) Encrypt(values []uint64) *bfv.Ciphertext {
	if values == nil || len(values) <= 0 {
		panic("Cannot encrypt uninitialized or empty values")
	}

	pt := newPlaintext(values, participant.encoder, participant.params)

	return participant.encryptor.Encrypt(pt)
}

// Decrypt returns the decryption of the ciphertext given the partial decryption
func (participant *mkParticipant) Decrypt(cipher *bfv.Ciphertext, partialDecryptions []*ring.Poly) []uint64 {

	if cipher == nil || cipher.Degree() != 1 {
		panic("Cannot decrypt uninitialized ciphertext or cipher of degree greater than 1")
	}

	if partialDecryptions == nil || len(partialDecryptions) < 1 {
		panic("Decryption necessitates at least one partialy decrypted ciphertext")
	}

	decrypted := participant.decryptor.MergeDec(cipher.Element, uint64(len(participant.ringQ.Modulus)-1), partialDecryptions)

	pt := bfv.NewPlaintext(*participant.params)
	pt.SetValue(decrypted)

	return participant.encoder.DecodeUintNew(pt)
}

// GetPartialDecryption returns the partial decryption of an element in the ciphertext
// this function should only be used by participants that were involved in the given ciphertext
func (participant *mkParticipant) GetPartialDecryption(cipher *bfv.Ciphertext) *ring.Poly {

	return participant.decryptor.PartDec(cipher.Element, uint64(len(participant.ringQ.Modulus)-1), participant.keys.SecretKey, 6.0)
}

// newParticipant creates a participant for the multi key bfv scheme
// the bfv parameters as well as the standard deviation used for partial decryption must be provided
func newParticipant(params *bfv.Parameters, sigmaSmudging float64, crs *mkrlwe.MKDecomposedPoly) participant {

	if crs == nil || params == nil {
		panic("Uninitialized parameters. Cannot create new participant")
	}

	if sigmaSmudging < params.Sigma() {
		panic("Sigma must be at least greater than the standard deviation of the gaussian distribution")
	}

	if len(crs.Poly) != int(params.Beta()) {
		panic("CRS must be the same dimention as returned by the function bfv.Parameters.Beta()")
	}

	keys := mkrlwe.KeyGen(&params.Parameters, mkrlwe.CopyNewDecomposed(crs))

	encryptor := NewMKEncryptor(keys.PublicKey, params)
	decryptor := mkrlwe.NewMKDecryptor(&params.Parameters)
	encoder := bfv.NewEncoder(*params)
	ringQ := mkrlwe.GetRingQ(&params.Parameters)

	return &mkParticipant{
		encryptor: encryptor,
		decryptor: decryptor,
		params:    params,
		keys:      keys,
		encoder:   encoder,
		ringQ:     ringQ,
	}
}

// newPlaintext initializes a new bfv Plaintext with an encoded slice of uint64
func newPlaintext(value []uint64, encoder bfv.Encoder, params *bfv.Parameters) *bfv.Plaintext {

	plaintext := bfv.NewPlaintext(*params)

	// Encode
	encoder.EncodeUint(value, plaintext)

	return plaintext
}

// GetRotationKeys returns the rotation key set associated with the given rotation
func (participant *mkParticipant) GetRotationKeys(rot int) *mkrlwe.MKEvalGalKey {

	galEl := participant.params.GaloisElementForColumnRotationBy(rot)

	evalKey := mkrlwe.GaloisEvaluationKeyGen(galEl, participant.keys.SecretKey, &participant.params.Parameters)

	return evalKey
}
