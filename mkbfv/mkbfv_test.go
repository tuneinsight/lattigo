package mkbfv

import (
	"fmt"
	"testing"

	"github.com/ldsec/lattigo/v2/bfv"
	"github.com/ldsec/lattigo/v2/ring"
)

func Test_MKBFV(t *testing.T) {

	for i := range bfv.DefaultParams {
		testEncryptionEqualsDecryption(t, i)
		testAdd(t, i)
		testAddFourParticipants(t, i)
	}
}

func testEncryptionEqualsDecryption(t *testing.T, paramsIndex int) {

	sigma := 6.0

	participants, params := setupPeers(1, paramsIndex, sigma)

	ringT := getRingT(params)

	expected := getRandomPlaintextValue(ringT, params)

	// encrypt
	cipher := participants[0].Encrypt(expected)

	// decrypt
	partialDec := participants[0].GetPartialDecryption(cipher)
	decrypted := participants[0].Decrypt(cipher.ciphertexts.Value()[0], []*ring.Poly{partialDec})

	// decode and check
	if !equalsSlice(decrypted, expected) {
		t.Error("Decryption of encryption does not equals plaintext")
	}
}

func testAdd(t *testing.T, paramsIndex int) {

	sigma := 6.0

	participants, params := setupPeers(2, paramsIndex, sigma)

	ringT := getRingT(params)

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

	decrypted := participants[0].Decrypt(resCipher.ciphertexts.Value()[0], []*ring.Poly{partialDec1, partialDec2})

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
}

func testAddFourParticipants(t *testing.T, paramsIndex int) {

	sigma := 6.0

	participants, params := setupPeers(4, paramsIndex, sigma)

	ringT := getRingT(params)

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

	decrypted := participants[0].Decrypt(resCipher.ciphertexts.Value()[0], []*ring.Poly{partialDec1, partialDec2, partialDec3, partialDec4})

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
}

func Test_Dot(t *testing.T) {

	sigma := 6.0

	_, params := setupPeers(2, 0, sigma)

	ringT := getRingT(params)

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

}

/*
func Test_Mul(t *testing.T) {

	ids := []uint64{1, 2}
	sigma := 6.0

	keys, params := setupPeers(ids, 0)

	ringQ := GetRingQ(params)
	ringT := getRingT(params)

	encoder := bfv.NewEncoder(params)

	encryptor1 := NewMKEncryptor(keys[0].publicKey, params, ids[0])
	encryptor2 := NewMKEncryptor(keys[1].publicKey, params, ids[1])
	decryptor := NewMKDecryptor(params, sigma)

	expected1 := getRandomPlaintextValue(ringT, params)
	expected2 := getRandomPlaintextValue(ringT, params)
	plaintext1 := newPlaintext(expected1, ringQ, encoder)
	plaintext2 := newPlaintext(expected2, ringQ, encoder)

	// encrypt
	cipher1 := encryptor1.EncryptMK(plaintext1.Plaintext())
	cipher2 := encryptor2.EncryptMK(plaintext2.Plaintext())

	// pad and mul
	evaluator := NewMKEvaluator(params)
	out1, out2 := PadCiphers(cipher1, cipher2, params)
	resCipher := evaluator.MultRelinDynamic(out1, out2, []*MKEvaluationKey{keys[0].evalKey, keys[1].evalKey}, []*MKPublicKey{keys[0].publicKey, keys[1].publicKey})
	// decrypt
	partialDec1 := decryptor.PartDec(resCipher.ciphertexts.Value()[1], keys[0].secretKey)
	partialDec2 := decryptor.PartDec(resCipher.ciphertexts.Value()[2], keys[1].secretKey)

	decrypted := decryptor.MergeDec(resCipher.ciphertexts.Value()[0], []*ring.Poly{partialDec1, partialDec2})
	// perform operation in plaintext space
	expected := ringT.NewPoly()
	p1 := ringT.NewPoly()
	p2 := ringT.NewPoly()
	copy(p1.Coeffs[0], expected1)
	copy(p2.Coeffs[0], expected2)
	ringT.MulCoeffs(p1, p2, expected)
	if !equalsSlice(encoder.DecodeUintNew(decrypted), expected.Coeffs[0]) {
		t.Error("Homomorphic multiplication error")
	}

}
*/
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
func setupPeers(peersNbr uint64, paramsNbr int, sigmaSmudging float64) ([]MKParticipant, *bfv.Parameters) {

	res := make([]MKParticipant, peersNbr)

	params := bfv.DefaultParams[paramsNbr]

	// setup keys and public parameters
	a := GenCommonPublicParam(params)

	for i := 0; i < int(peersNbr); i++ {

		res[i] = NewParticipant(params, sigmaSmudging, a)

	}

	return res, params
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
