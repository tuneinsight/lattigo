package mkbfv

import (
	"fmt"
	"testing"

	"github.com/ldsec/lattigo/v2/bfv"
	"github.com/ldsec/lattigo/v2/ring"
)

func Test_EncryptionEqualsDecryption(t *testing.T) {

	ids := []uint64{1}

	sigma := 6.0

	keys, params := setupPeers(ids, 0)

	ringQ := GetRingQ(params)
	ringT := getRingT(params)

	encryptor := NewMKEncryptor(keys[0].publicKey, params, ids[0])
	decryptor := NewMKDecryptor(params, sigma)

	encoder := bfv.NewEncoder(params)

	expected := getRandomPlaintextValue(ringT, params)
	plaintext := newPlaintext(expected, ringQ, encoder)

	// encrypt
	cipher := encryptor.EncryptMK(plaintext)

	// decrypt
	partialDec := decryptor.PartDec(cipher.ciphertexts.Value()[1], keys[0].secretKey)
	decrypted := decryptor.MergeDec(cipher.ciphertexts.Value()[0], []*ring.Poly{partialDec})

	// decode and check
	if !equalsSlice(encoder.DecodeUintNew(decrypted), expected) {
		t.Error("Decryption of encryption does not equals plaintext")
	}
}

func Test_Add(t *testing.T) {

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

	// pad and add
	evaluator := NewMKEvaluator(params)

	out1, out2 := PadCiphers(cipher1, cipher2, params)

	resCipher := evaluator.Add(out1, out2)

	// decrypt
	partialDec1 := decryptor.PartDec(resCipher.ciphertexts.Value()[1], keys[0].secretKey)
	partialDec2 := decryptor.PartDec(resCipher.ciphertexts.Value()[2], keys[1].secretKey)

	decrypted := decryptor.MergeDec(resCipher.ciphertexts.Value()[0], []*ring.Poly{partialDec1, partialDec2})

	// perform the operation in the plaintext space
	expected := ringT.NewPoly()
	p1 := ringT.NewPoly()
	p2 := ringT.NewPoly()
	copy(p1.Coeffs[0], expected1)
	copy(p2.Coeffs[0], expected2)

	ringT.Add(p1, p2, expected)

	if !equalsSlice(encoder.DecodeUintNew(decrypted), expected.Coeffs[0]) {
		t.Error("Homomorphic addition error")
	}
}

func Test_AddFourParticipants(t *testing.T) {

	ids := []uint64{1, 2, 5, 7}
	sigma := 6.0

	keys, params := setupPeers(ids, 0)

	ringQ := GetRingQ(params)
	ringT := getRingT(params)

	encoder := bfv.NewEncoder(params)

	encryptor1 := NewMKEncryptor(keys[0].publicKey, params, ids[0])
	encryptor2 := NewMKEncryptor(keys[1].publicKey, params, ids[1])
	encryptor3 := NewMKEncryptor(keys[2].publicKey, params, ids[2])
	encryptor4 := NewMKEncryptor(keys[3].publicKey, params, ids[3])
	decryptor := NewMKDecryptor(params, sigma)

	expected1 := getRandomPlaintextValue(ringT, params)
	expected2 := getRandomPlaintextValue(ringT, params)
	expected3 := getRandomPlaintextValue(ringT, params)
	expected4 := getRandomPlaintextValue(ringT, params)
	plaintext1 := newPlaintext(expected1, ringQ, encoder)
	plaintext2 := newPlaintext(expected2, ringQ, encoder)
	plaintext3 := newPlaintext(expected3, ringQ, encoder)
	plaintext4 := newPlaintext(expected4, ringQ, encoder)

	// encrypt
	cipher1 := encryptor1.EncryptMK(plaintext1.Plaintext())
	cipher2 := encryptor2.EncryptMK(plaintext2.Plaintext())
	cipher3 := encryptor3.EncryptMK(plaintext3.Plaintext())
	cipher4 := encryptor4.EncryptMK(plaintext4.Plaintext())

	// pad and add in 2 steps
	evaluator := NewMKEvaluator(params)

	out1, out2 := PadCiphers(cipher1, cipher2, params)
	out3, out4 := PadCiphers(cipher3, cipher4, params)

	resCipher1 := evaluator.Add(out1, out2)
	resCipher2 := evaluator.Add(out3, out4)

	finalOut1, finalOut2 := PadCiphers(resCipher1, resCipher2, params)

	resCipher := evaluator.Add(finalOut1, finalOut2)

	// decrypt
	partialDec1 := decryptor.PartDec(resCipher.ciphertexts.Value()[1], keys[0].secretKey)
	partialDec2 := decryptor.PartDec(resCipher.ciphertexts.Value()[2], keys[1].secretKey)
	partialDec3 := decryptor.PartDec(resCipher.ciphertexts.Value()[3], keys[2].secretKey)
	partialDec4 := decryptor.PartDec(resCipher.ciphertexts.Value()[4], keys[3].secretKey)

	decrypted := decryptor.MergeDec(resCipher.ciphertexts.Value()[0], []*ring.Poly{partialDec1, partialDec2, partialDec3, partialDec4})

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

	if !equalsSlice(encoder.DecodeUintNew(decrypted), expected.Coeffs[0]) {
		t.Error("Homomorphic addition error")
	}
}

func Test_Dot(t *testing.T) {
	ids := []uint64{1, 2}

	_, params := setupPeers(ids, 0)

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

func Test_KeyGenAndPadCiphertext(t *testing.T) {

	// setup keys and public parameters
	ids := []uint64{12, 3}

	keys, params := setupPeers(ids, 0)
	if keys == nil || params == nil {
		t.Errorf("Generation of common public parameter failed !")
	}

	encoder := bfv.NewEncoder(params)

	//create and encrypt 2 plaintext
	encryptor1 := NewMKEncryptor(keys[0].publicKey, params, ids[0])
	encryptor2 := NewMKEncryptor(keys[1].publicKey, params, ids[1])

	val1 := getRandomPlaintextValue(getRingT(params), params)
	plaintext1 := newPlaintext(val1, GetRingQ(params), encoder)
	val2 := getRandomPlaintextValue(getRingT(params), params)
	plaintext2 := newPlaintext(val2, GetRingQ(params), encoder)

	mkCiphertext1 := encryptor1.EncryptMK(plaintext1)
	mkCiphertext2 := encryptor2.EncryptMK(plaintext2)

	// Pad both ciphertexts
	out1, out2 := PadCiphers(mkCiphertext1, mkCiphertext2, params)

	// check peerIDs
	if !equalsSlice(out1.peerIDs, out2.peerIDs) {
		t.Errorf("PadCipher failed. IDs different")
	}

	// length should be 3 since padded ciphertext should be (c0, c3, c12)
	l1 := len(out1.ciphertexts.Element.Value())
	l2 := len(out2.ciphertexts.Element.Value())
	if l1 != 3 || l2 != 3 {
		t.Errorf("PadCipher failed. Length not equal to 3. Length 1 = %d, length 2 = %d", l1, l2)
	}

	// check content
	ec00 := mkCiphertext1.ciphertexts.Element.Value()[0]
	ec01 := mkCiphertext1.ciphertexts.Element.Value()[1]
	ec10 := mkCiphertext2.ciphertexts.Element.Value()[0]
	ec11 := mkCiphertext2.ciphertexts.Element.Value()[1]

	c00 := out1.ciphertexts.Element.Value()[0]
	c01 := out1.ciphertexts.Element.Value()[1]
	c02 := out1.ciphertexts.Element.Value()[2]
	c10 := out2.ciphertexts.Element.Value()[0]
	c11 := out2.ciphertexts.Element.Value()[1]
	c12 := out2.ciphertexts.Element.Value()[2]

	ringQ := GetRingQ(params)
	pass := equalsPoly(ec00, c00) && equalsPoly(ec10, c10) && equalsPoly(ec01, c02) && equalsPoly(ec11, c11) && equalsPoly(c01, ringQ.NewPoly()) && equalsPoly(c12, ringQ.NewPoly())

	if !pass {
		t.Errorf("PadCipher failed. Content of ciphertext not correct")
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
func setupPeers(peerIDs []uint64, paramsNbr int) ([]*MKKeys, *bfv.Parameters) {

	res := make([]*MKKeys, len(peerIDs))

	params := bfv.DefaultParams[paramsNbr]

	for i := 0; i < len(peerIDs); i++ {

		// setup keys and public parameters
		a := GenCommonPublicParam(params)

		res[i] = KeyGen(params, peerIDs[i], a)

	}

	return res, params
}

// newPlaintext initializes a new bfv Plaintext with an encoded slice of uint64
func newPlaintext(value []uint64, ringQ *ring.Ring, encoder bfv.Encoder) *bfv.Plaintext {

	plaintext := new(bfv.Element)

	ptValues := make([]*ring.Poly, 1)

	ptValues[0] = ringQ.NewPoly()
	plaintext.SetValue(ptValues)

	// Encode
	encoder.EncodeUint(value, plaintext.Plaintext())

	return plaintext.Plaintext()
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
