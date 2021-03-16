package mkbfv

import (
	"testing"

	"github.com/ldsec/lattigo/v2/bfv"
	"github.com/ldsec/lattigo/v2/ring"
)

func Test_EncryptionEqualsDecryption(t *testing.T) {

	ids := []uint64{1}

	keys, params := setupPeers(ids, 0)

	ringQ := GetRingQ(params)

	encryptor := NewMKEncryptor(keys[0].publicKey, params)
	decryptor := NewMKDecryptor(params)

	encoder := bfv.NewEncoder(params)

	plaintext := new(bfv.Element)

	ptValues := make([]*ring.Poly, 1)
	expected := GetRandomPoly(params, ringQ)
	ptValues[0] = ringQ.NewPoly()
	plaintext.SetValue(ptValues)

	// Encode
	encoder.EncodeUint(expected.Coeffs[0], plaintext.Plaintext())

	// encrypt
	cipher := encryptor.EncryptMK(plaintext.Plaintext(), ids[0])

	partialDec := ringQ.NewPoly()

	// decrypt
	decryptor.PartDec(cipher.ciphertexts.Value()[1], keys[0].secretKey, partialDec)

	decrypted := decryptor.MergeDec(cipher.ciphertexts.Value()[0], []*ring.Poly{partialDec})

	// decode and check
	if !equalsSlice(encoder.DecodeUintNew(decrypted), expected.Coeffs[0]) {
		t.Error("Decryption of encryption does not equals plaintext")
	}
}

func Test_Add(t *testing.T) {

	ids := []uint64{1, 2}

	keys, params := setupPeers(ids, 0)

	ringQ := GetRingQ(params)

	encoder := bfv.NewEncoder(params)

	encryptor1 := NewMKEncryptor(keys[0].publicKey, params)
	encryptor2 := NewMKEncryptor(keys[1].publicKey, params)
	decryptor := NewMKDecryptor(params)

	plaintext1 := new(bfv.Element)
	plaintext2 := new(bfv.Element)

	ptValues1 := make([]*ring.Poly, 1)
	ptValues2 := make([]*ring.Poly, 1)
	ptValues1[0] = ringQ.NewPoly()
	ptValues2[0] = ringQ.NewPoly()

	plaintext1.SetValue(ptValues1)
	plaintext2.SetValue(ptValues2)

	expected1 := GetRandomPoly(params, ringQ)
	expected2 := GetRandomPoly(params, ringQ)

	// encode
	encoder.EncodeUint(expected1.Coeffs[0], plaintext1.Plaintext())
	encoder.EncodeUint(expected2.Coeffs[0], plaintext2.Plaintext())

	// encrypt
	cipher1 := encryptor1.EncryptMK(plaintext1.Plaintext(), ids[0])
	cipher2 := encryptor2.EncryptMK(plaintext2.Plaintext(), ids[1])

	// pad and add
	evaluator := NewMKEvaluator(params)

	out1, out2 := PadCiphers(cipher1, cipher2, params)

	resCipher := evaluator.Add(out1, out2, params)

	// decrypt
	partialDec1 := ringQ.NewPoly()
	partialDec2 := ringQ.NewPoly()
	decryptor.PartDec(resCipher.ciphertexts.Value()[1], keys[0].secretKey, partialDec1)
	decryptor.PartDec(resCipher.ciphertexts.Value()[2], keys[1].secretKey, partialDec2)

	decrypted := decryptor.MergeDec(resCipher.ciphertexts.Value()[0], []*ring.Poly{partialDec1, partialDec2})

	expected := ringQ.NewPoly()
	ringQ.AddNoMod(expected1, expected2, expected)

	if !equalsSlice(encoder.DecodeUintNew(decrypted), expected.Coeffs[0]) {
		t.Error("Homomorphic addition error")
	}
}

func Test_Mul(t *testing.T) {

	ids := []uint64{1, 2}

	keys, params := setupPeers(ids, 0)

	ringQ := GetRingQ(params)

	encoder := bfv.NewEncoder(params)

	encryptor1 := NewMKEncryptor(keys[0].publicKey, params)
	encryptor2 := NewMKEncryptor(keys[1].publicKey, params)
	decryptor := NewMKDecryptor(params)

	plaintext1 := new(bfv.Element)
	plaintext2 := new(bfv.Element)

	ptValues1 := make([]*ring.Poly, 1)
	ptValues2 := make([]*ring.Poly, 1)
	ptValues1[0] = ringQ.NewPoly()
	ptValues2[0] = ringQ.NewPoly()

	plaintext1.SetValue(ptValues1)
	plaintext2.SetValue(ptValues2)

	expected1 := GetRandomPoly(params, ringQ)
	expected2 := GetRandomPoly(params, ringQ)

	expected := make([]uint64, len(expected1.Coeffs[0]))

	for i := 0; i < len(expected1.Coeffs[0]); i++ {
		expected[i] = expected1.Coeffs[0][i] * expected2.Coeffs[0][i]
	}

	// encode
	encoder.EncodeUint(expected1.Coeffs[0], plaintext1.Plaintext())
	encoder.EncodeUint(expected2.Coeffs[0], plaintext2.Plaintext())

	// encrypt
	cipher1 := encryptor1.EncryptMK(plaintext1.Plaintext(), ids[0])
	cipher2 := encryptor2.EncryptMK(plaintext2.Plaintext(), ids[1])

	// pad and mul
	evaluator := NewMKEvaluator(params)

	out1, out2 := PadCiphers(cipher1, cipher2, params)

	resCipher := evaluator.MultRelinDynamic(out1, out2, []*MKEvaluationKey{keys[0].evalKey, keys[1].evalKey}, []*MKPublicKey{keys[0].publicKey, keys[1].publicKey}, params)

	// decrypt
	partialDec1 := ringQ.NewPoly()
	partialDec2 := ringQ.NewPoly()
	decryptor.PartDec(resCipher.ciphertexts.Value()[1], keys[0].secretKey, partialDec1)
	decryptor.PartDec(resCipher.ciphertexts.Value()[2], keys[1].secretKey, partialDec2)

	decrypted := decryptor.MergeDec(resCipher.ciphertexts.Value()[0], []*ring.Poly{partialDec1, partialDec2})

	if !equalsSlice(encoder.DecodeUintNew(decrypted), expected) {
		t.Error("Homomorphic multiplication error")
	}

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

func Test_KeyGenAndPadCiphertext(t *testing.T) {

	// setup keys and public parameters
	params := bfv.DefaultParams[0]
	a := GenCommonPublicParam(params)

	if a == nil || params == nil {
		t.Errorf("Generation of common public parameter failed !")
	}

	id1 := uint64(12)
	id2 := uint64(3)

	keys := KeyGen(params, id1, a)
	if keys == nil {
		t.Errorf("Keys generation failed")
	}

	pubKeys := make([]*MKPublicKey, 1)
	pubKeys[0] = keys.publicKey

	evalKeys := make([]*MKEvaluationKey, 1)
	evalKeys[0] = keys.evalKey

	keys.relinKey = GenSharedRelinearizationKey(params, pubKeys, evalKeys)
	if keys.relinKey == nil {
		t.Errorf("Shared relin key generation failed")
	}

	//create and encrypt 2 plaintext
	encryptor := NewMKEncryptor(keys.publicKey, params)

	plaintext := new(bfv.Element)

	values1 := make([]*ring.Poly, 1)
	values2 := make([]*ring.Poly, 1)
	value1 := a.poly[0]
	value2 := a.poly[0]

	values1[0] = value1
	values2[0] = value2

	plaintext.SetValue(values1)
	mkCiphertext1 := encryptor.EncryptMK(plaintext.Plaintext(), id1)

	plaintext.SetValue(values2)
	mkCiphertext2 := encryptor.EncryptMK(plaintext.Plaintext(), id2)

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
