package mkbfv

import (
	"testing"

	"github.com/ldsec/lattigo/v2/bfv"
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/utils"
)

func Test_MKBFV(t *testing.T) {

	for i, p := range bfv.DefaultParams {

		testEncryptionEqualsDecryption(t, p)
		testAdd(t, p)
		testAddFourParticipants(t, p)
		testAddPlaintext(t, p)
		testAddPlaintextTwoParticipants(t, p)
		testSub(t, p)
		testSubPlaintext(t, p)
		testNeg(t, p)
		testSubPlaintextTwoParticipants(t, p)
		testMulPlaintext(t, p)
		testMulPlaintextTwoParticipants(t, p)

		if i != 0 && i != 4 && i != 6 {
			testMulFourParticipants(t, p)

		}
		if i != 6 {
			testMul(t, p)
			testAddAfterMul(t, p)
		}

		//testRotation(t, p)
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

		resCipher := evaluator.Add(cipher1, cipher2)

		// decrypt
		partialDec1 := participants[0].GetPartialDecryption(resCipher)
		partialDec2 := participants[1].GetPartialDecryption(resCipher)

		decrypted := participants[0].Decrypt(resCipher, []*ring.Poly{partialDec1, partialDec2})

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

		resCipher := evaluator.Sub(cipher1, cipher2)

		// decrypt
		partialDec1 := participants[0].GetPartialDecryption(resCipher)
		partialDec2 := participants[1].GetPartialDecryption(resCipher)

		decrypted := participants[0].Decrypt(resCipher, []*ring.Poly{partialDec1, partialDec2})

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

		// encrypt
		cipher := participants[0].Encrypt(value1)

		// add plaintext to ciphertext
		evaluator := NewMKEvaluator(params)
		value2 := getRandomPlaintextValue(ringT, params)

		pt := evaluator.NewPlaintextFromValue(value2)
		resCipher := evaluator.AddPlaintext(pt, cipher)

		// decrypt
		partialDec1 := participants[0].GetPartialDecryption(resCipher)

		decrypted := participants[0].Decrypt(resCipher, []*ring.Poly{partialDec1})

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

		// encrypt
		cipher1 := participants[0].Encrypt(value1)
		cipher2 := participants[1].Encrypt(value2)

		// add plaintext to one of the ciphertext then add both ciphertexts
		evaluator := NewMKEvaluator(params)
		value3 := getRandomPlaintextValue(ringT, params)

		pt := evaluator.NewPlaintextFromValue(value3)
		resCipher1 := evaluator.AddPlaintext(pt, cipher1)

		resCipher := evaluator.Add(resCipher1, cipher2)

		// decrypt
		partialDec1 := participants[0].GetPartialDecryption(resCipher)
		partialDec2 := participants[1].GetPartialDecryption(resCipher)

		decrypted := participants[0].Decrypt(resCipher, []*ring.Poly{partialDec1, partialDec2})

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

		// encrypt
		cipher := participants[0].Encrypt(value1)

		// sub plaintext to ciphertext
		evaluator := NewMKEvaluator(params)
		value2 := getRandomPlaintextValue(ringT, params)

		pt := evaluator.NewPlaintextFromValue(value2)
		resCipher := evaluator.SubPlaintext(pt, cipher)

		// decrypt
		partialDec1 := participants[0].GetPartialDecryption(resCipher)

		decrypted := participants[0].Decrypt(resCipher, []*ring.Poly{partialDec1})

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

		// encrypt
		cipher1 := participants[0].Encrypt(value1)
		cipher2 := participants[1].Encrypt(value2)

		// sub plaintext to one of the ciphertext then sub both ciphertexts
		evaluator := NewMKEvaluator(params)
		value3 := getRandomPlaintextValue(ringT, params)

		pt := evaluator.NewPlaintextFromValue(value3)
		resCipher1 := evaluator.SubPlaintext(pt, cipher1)

		resCipher := evaluator.Sub(resCipher1, cipher2)

		// decrypt
		partialDec1 := participants[0].GetPartialDecryption(resCipher)
		partialDec2 := participants[1].GetPartialDecryption(resCipher)

		decrypted := participants[0].Decrypt(resCipher, []*ring.Poly{partialDec1, partialDec2})

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

		resCipher := evaluator.Add(evaluator.Neg(cipher), cipher)

		// decrypt
		partialDec1 := participants[0].GetPartialDecryption(resCipher)

		decrypted := participants[0].Decrypt(resCipher, []*ring.Poly{partialDec1})

		// should be 0
		expected := ringT.NewPoly()

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

func Test_Dot(t *testing.T) {

	params := bfv.DefaultParams[0]

	ringT := getRingT(params)

	t.Run(testString("Test Dot product/", 1, params), func(t *testing.T) {

		expected1 := getRandomPlaintextValue(ringT, params)
		expected2 := getRandomPlaintextValue(ringT, params)
		// perform the operation in the plaintext space
		expected := ringT.NewPoly()
		p1 := ringT.NewPoly()
		p2 := ringT.NewPoly()
		p3 := ringT.NewPoly()
		copy(p1.Coeffs[0], expected1)
		copy(p2.Coeffs[0], expected2)

		ringT.Add(p1, p2, expected)

		// perform dot product in the plaintext space : multiply (pt1,pt2) by (1,1) and compare with Add(pt1,pt2)

		decomposedPoly1 := NewDecomposedPoly(ringT, 2)
		decomposedPoly1.poly[0] = p1
		decomposedPoly1.poly[1] = p2
		const1Int := make([]uint64, len(p1.Coeffs[0]))
		for i := uint64(0); i < uint64(len(const1Int)); i++ {
			const1Int[i] = 1
		}
		ringT.SetCoefficientsUint64(const1Int, p3)

		decomposedPoly2 := NewDecomposedPoly(ringT, 2)
		decomposedPoly2.poly[0] = p3 // set equal to 1
		decomposedPoly2.poly[1] = p3 // set equal to 1

		dotExpected := Dot(decomposedPoly1, decomposedPoly2, ringT)

		if !equalsSlice(dotExpected.Coeffs[0], expected.Coeffs[0]) {
			t.Error("Dot error")
		}

	})

}

/*
func Test_TensorProduct(t *testing.T) {

	params := bfv.DefaultParams[0]

	sigma := 6.0

	participants := setupPeers(2, params, sigma)

	ringT := getRingT(params)
	ringQ := GetRingQ(params)

	t.Run(testString("Test Tensor product/", 2, params), func(t *testing.T) {

		// Create two ciphertexts (c0,c1,0), (c0,0,c1) with c0 = 0, c1 = 1 and check that tensor product without rescaling = (0,0,0,0,0,0,0,1,0)
		expected1 := getRandomPlaintextValue(ringT, params)
		expected2 := getRandomPlaintextValue(ringT, params)

		// encrypt
		cipher1 := participants[0].Encrypt(expected1)
		cipher2 := participants[1].Encrypt(expected2)

		// Create 3 polynomials in ringQ then set coeffs to 1
		p1 := ringQ.NewPoly()
		const1Int := make([]uint64, len(p1.Coeffs[0]))
		for i := uint64(0); i < uint64(len(const1Int)-1); i++ {
			const1Int[i] = 100
		}
		const1Int[len(const1Int)-1] = 1
		ringQ.SetCoefficientsUint64(const1Int, p1)

		p2 := ringQ.NewPoly()
		ringQ.SetCoefficientsUint64(const1Int, p2)

		cipherPoly := make([]*ring.Poly, 2)

		cipherPoly[0] = p1
		cipherPoly[1] = p1

		cipher1.ciphertexts.SetValue(cipherPoly)
		cipher2.ciphertexts.SetValue(cipherPoly)

		out1, out2 := PadCiphers(cipher1, cipher2, params) // (0,1,0), (0,0,1)

		evaluator := NewMKEvaluator(params)
		testPoly := evaluator.TensorAndRescale(out1.ciphertexts.El(), out2.ciphertexts.El())

		if !equalsSlice(testPoly.ciphertexts.Value()[7].Coeffs[0], []uint64{0, 0, 0, 0, 0, 0, 0, 0, 1}) {
			t.Error("tensor error")
		}

	})

}
*/

func testMulPlaintext(t *testing.T, params *bfv.Parameters) {

	sigma := 6.0

	participants := setupPeers(1, params, sigma)

	ringT := getRingT(params)

	t.Run(testString("Test Plaintext Multiplication/", 1, params), func(t *testing.T) {

		value1 := getRandomPlaintextValue(ringT, params)

		// encrypt
		cipher := participants[0].Encrypt(value1)

		// multiply plaintext and ciphertext
		evaluator := NewMKEvaluator(params)
		value2 := getRandomPlaintextValue(ringT, params)

		pt := evaluator.NewPlaintextMulFromValue(value2)
		resCipher := evaluator.MultPlaintext(pt, cipher)

		// decrypt
		partialDec1 := participants[0].GetPartialDecryption(resCipher)

		decrypted := participants[0].Decrypt(resCipher, []*ring.Poly{partialDec1})

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

func testMulPlaintextTwoParticipants(t *testing.T, params *bfv.Parameters) {

	sigma := 6.0

	participants := setupPeers(2, params, sigma)

	ringT := getRingT(params)

	t.Run(testString("Test Plaintext Multiplication/", 2, params), func(t *testing.T) {

		value1 := getRandomPlaintextValue(ringT, params)
		value2 := getRandomPlaintextValue(ringT, params)

		// encrypt
		cipher1 := participants[0].Encrypt(value1)
		cipher2 := participants[1].Encrypt(value2)

		// add both ciphertexts then multiply by plaintext
		evaluator := NewMKEvaluator(params)
		value3 := getRandomPlaintextValue(ringT, params)

		resCipherTMP := evaluator.Add(cipher1, cipher2)

		pt := evaluator.NewPlaintextMulFromValue(value3)
		resCipher := evaluator.MultPlaintext(pt, resCipherTMP)

		// decrypt
		partialDec1 := participants[0].GetPartialDecryption(resCipher)
		partialDec2 := participants[1].GetPartialDecryption(resCipher)

		decrypted := participants[0].Decrypt(resCipher, []*ring.Poly{partialDec1, partialDec2})

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

		// pad
		evaluator := NewMKEvaluator(params)

		// multiply using evaluation keys and public keys
		evalKeys := []*MKEvaluationKey{participants[0].GetEvaluationKey(), participants[1].GetEvaluationKey()}
		publicKeys := []*MKPublicKey{participants[0].GetPublicKey(), participants[1].GetPublicKey()}

		resCipher := evaluator.Mul(cipher1, cipher2)

		evaluator.RelinInPlace(resCipher, evalKeys, publicKeys)

		// decrypt
		partialDec1 := participants[0].GetPartialDecryption(resCipher)
		partialDec2 := participants[1].GetPartialDecryption(resCipher)

		decrypted := participants[0].Decrypt(resCipher, []*ring.Poly{partialDec1, partialDec2})

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

		// multiply using evaluation keys and publick keys
		evalKeys := []*MKEvaluationKey{participants[0].GetEvaluationKey(), participants[1].GetEvaluationKey()}
		publicKeys := []*MKPublicKey{participants[0].GetPublicKey(), participants[1].GetPublicKey()}

		resCipher := evaluator.Mul(cipher1, cipher2)

		evaluator.RelinInPlace(resCipher, evalKeys, publicKeys)

		resCipher = evaluator.Add(resCipher, cipher1)

		// decrypt
		partialDec1 := participants[0].GetPartialDecryption(resCipher)
		partialDec2 := participants[1].GetPartialDecryption(resCipher)

		decrypted := participants[0].Decrypt(resCipher, []*ring.Poly{partialDec1, partialDec2})

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
		evalKeys := []*MKEvaluationKey{participants[0].GetEvaluationKey(), participants[1].GetEvaluationKey(), participants[2].GetEvaluationKey(), participants[3].GetEvaluationKey()}
		publicKeys := []*MKPublicKey{participants[0].GetPublicKey(), participants[1].GetPublicKey(), participants[2].GetPublicKey(), participants[3].GetPublicKey()}

		resCipher1 := evaluator.Mul(cipher1, cipher2)
		resCipher2 := evaluator.Mul(cipher3, cipher4)

		evaluator.RelinInPlace(resCipher1, evalKeys[:2], publicKeys[:2])
		evaluator.RelinInPlace(resCipher2, evalKeys[2:], publicKeys[2:])

		resCipher := evaluator.Mul(resCipher1, resCipher2)

		evaluator.RelinInPlace(resCipher, evalKeys, publicKeys)

		// decrypt
		partialDec1 := participants[0].GetPartialDecryption(resCipher)
		partialDec2 := participants[1].GetPartialDecryption(resCipher)
		partialDec3 := participants[2].GetPartialDecryption(resCipher)
		partialDec4 := participants[3].GetPartialDecryption(resCipher)

		decrypted := participants[0].Decrypt(resCipher, []*ring.Poly{partialDec1, partialDec2, partialDec3, partialDec4})

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

		for _, n := range rots {

			rotKey := participants[0].GetRotationKeys(n)

			values2 := utils.RotateUint64Slice(values1, n)
			resCipher := evaluator.Rotate(cipher1, n, []*MKEvalGalKey{rotKey})

			partialDec := participants[0].GetPartialDecryption(resCipher)

			decrypted := participants[0].Decrypt(resCipher, []*ring.Poly{partialDec})

			if !utils.EqualSliceUint64(values2, decrypted) {
				t.Errorf("Rotation error")
			}
		}
	})

}

func Test_Utils(t *testing.T) {

	s1 := []uint64{0, 2, 1}

	s2 := []uint64{3, 7, 12, 1, 0}

	expected := []uint64{0, 1, 2, 3, 7, 12}

	res := MergeSlices(s1, s2)

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
func setupPeers(peersNbr uint64, params *bfv.Parameters, sigmaSmudging float64) []MKParticipant {

	res := make([]MKParticipant, peersNbr)

	prng, err := utils.NewKeyedPRNG([]byte{'l', 'a', 't', 't', 'i', 'g', 'o'})

	if err != nil {
		panic(err)
	}

	// setup keys and public parameters
	a := GenCommonPublicParam(params, prng)

	for i := 0; i < int(peersNbr); i++ {

		res[i] = NewParticipant(params, sigmaSmudging, a)

	}

	return res
}

// returns a uniformly random slice of uint64 in RingT
func getRandomPlaintextValue(ringT *ring.Ring, params *bfv.Parameters) []uint64 {

	return GetRandomPoly(params, ringT).Coeffs[0]
}

// getRingT returns the ring from which the plaintexts will be sampled
func getRingT(params *bfv.Parameters) *ring.Ring {

	if ringT, err := ring.NewRing(params.N(), []uint64{params.T()}); err != nil {
		panic("Couldn't create ringT with given parameters")

	} else {
		return ringT
	}
}
