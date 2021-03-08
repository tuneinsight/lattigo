package mkbfv

import (
	"testing"

	"github.com/ldsec/lattigo/v2/bfv"
	"github.com/ldsec/lattigo/v2/ring"
)

func Test_MKBFV(t *testing.T) {

}

func Test_Utils(t *testing.T) {

	s1 := make([]uint64, 3)
	s1[0] = 0
	s1[1] = 2
	s1[2] = 1

	s2 := make([]uint64, 5)
	s2[0] = 3
	s2[1] = 7
	s2[2] = 12
	s2[3] = 1
	s2[4] = 0

	expected := make([]uint64, 6)
	expected[0] = 0
	expected[1] = 1
	expected[2] = 2
	expected[3] = 3
	expected[4] = 7
	expected[5] = 12

	res := MergeSlices(s1, s2)

	if !equalsSlice(expected, res) {
		t.Errorf("MergeSlices method failed test")
	}

}

func Test_PadCiphertext(t *testing.T) {

	// setup keys and public parameters
	params := bfv.DefaultParams[0]
	a := GenCommonPublicParam(params)

	id1 := uint64(12)
	id2 := uint64(3)

	keys := KeyGen(params, id1, a)

	pubKeys := make([]*MKPublicKey, 1)
	pubKeys[0] = keys.publicKey

	evalKeys := make([]*MKEvaluationKey, 1)
	evalKeys[0] = keys.evalKey

	keys.relinKey = GenSharedRelinearizationKey(params, pubKeys, evalKeys)

	//create and encrypt 2 plaintext
	encryptor := NewEncryptor(keys.publicKey, params)

	plaintext := new(bfv.Element)

	values1 := make([]*ring.Poly, 1)
	values2 := make([]*ring.Poly, 1)
	value1 := a.poly[0]
	value2 := a.poly[0]

	values1[0] = value1
	values2[0] = value2

	plaintext.SetValue(values1)
	mkCiphertext1 := EncryptMK(encryptor, plaintext.Plaintext(), id1)

	plaintext.SetValue(values2)
	mkCiphertext2 := EncryptMK(encryptor, plaintext.Plaintext(), id2)

	// Pad both ciphertexts
	out1 := new(MKCiphertext)
	out2 := new(MKCiphertext)

	PadCiphers(mkCiphertext1, mkCiphertext2, out1, out2, params)

	// check peerIDs
	if !equalsSlice(out1.peerIDs, out2.peerIDs) {
		t.Errorf("PadCipher failed. IDs different")
	}

	// check length should be 3 since padded ciphertext should be (c0, c3, c12) since ids are 3 and 12
	l1 := len(mkCiphertext1.ciphertexts.Element.Value())
	l2 := len(mkCiphertext2.ciphertexts.Element.Value())
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
