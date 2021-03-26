package mkbfv

import (
	"fmt"
	"testing"

	"github.com/ldsec/lattigo/v2/bfv"
	"github.com/ldsec/lattigo/v2/ring"
)

func Test_MKBFV(t *testing.T) {

	for _, p := range bfv.DefaultParams {
		testEncryptionEqualsDecryption(t, p)
		testAdd(t, p)
		testAddFourParticipants(t, p)
		//testMul(t, i)
		//testMulFourParticipants(t, i)
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

		out1, out2 := PadCiphers(cipher1, cipher2, params)

		resCipher := evaluator.Add(out1, out2)

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

		out1, out2 := PadCiphers(cipher1, cipher2, params)
		out3, out4 := PadCiphers(cipher3, cipher4, params)

		resCipher1 := evaluator.Add(out1, out2)
		resCipher2 := evaluator.Add(out3, out4)

		finalOut1, finalOut2 := PadCiphers(resCipher1, resCipher2, params)

		resCipher := evaluator.Add(finalOut1, finalOut2)

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

		dotExpected := ringT.NewPoly()
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

		Dot(decomposedPoly1, decomposedPoly2, dotExpected, ringT)

		if !equalsSlice(dotExpected.Coeffs[0], expected.Coeffs[0]) {
			t.Error("Dot error")
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
		out1, out2 := PadCiphers(cipher1, cipher2, params)

		// multiply using evaluation keys and publick keys
		evalKeys := []*MKEvaluationKey{participants[0].GetEvaluationKey(), participants[1].GetEvaluationKey()}
		publicKeys := []*MKPublicKey{participants[0].GetPublicKey(), participants[1].GetPublicKey()}

		resCipher := evaluator.MultRelinDynamic(out1, out2, evalKeys, publicKeys)

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

		out1, out2 := PadCiphers(cipher1, cipher2, params)
		out3, out4 := PadCiphers(cipher3, cipher4, params)

		resCipher1 := evaluator.MultRelinDynamic(out1, out2, evalKeys[:2], publicKeys[:2])
		resCipher2 := evaluator.MultRelinDynamic(out3, out4, evalKeys[2:], publicKeys[2:])

		finalOut1, finalOut2 := PadCiphers(resCipher1, resCipher2, params)

		resCipher := evaluator.MultRelinDynamic(finalOut1, finalOut2, evalKeys, publicKeys)

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
			fmt.Printf("%d not equal to %d . Error at index %d", e, s2[i], i)
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

	// setup keys and public parameters
	a := GenCommonPublicParam(params)

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
