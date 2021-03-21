package mkbfv

import (
	"flag"
	"fmt"
	"math"
	"sort"
	"testing"

	"github.com/ldsec/lattigo/v2/ckks"
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/utils"
	"github.com/stretchr/testify/require"
)

var printPrecisionStats = flag.Bool("print-precision", false, "print precision stats")
var minPrec float64 = 15.0

func Test_EncryptionEqualsDecryption(t *testing.T) {

	ids := []uint64{1}

	keys, params := setupPeers(ids, 0)

	encryptor := NewMKEncryptor(keys[0].publicKey, params, ids[0])
	decryptor := NewMKDecryptor(params)

	encoder := ckks.NewEncoder(params)

	// encrypt
	values, _, cipher := newTestVectors(params, encoder, encryptor, complex(-1, -1), complex(1, 1), t)

	// decrypt
	scale := cipher.ciphertexts.Scale()
	level := cipher.ciphertexts.Level()
	partialDec := decryptor.PartDec(cipher.ciphertexts.Value()[1], level, keys[0].secretKey)
	decrypted := decryptor.MergeDec(cipher.ciphertexts.Value()[0], scale, level, []*ring.Poly{partialDec})

	// decode and check
	verifyTestVectors(params, encoder, values, decrypted, t)
}

func Test_Add(t *testing.T) {

	ids := []uint64{1, 2}

	keys, params := setupPeers(ids, 0)

	encoder := ckks.NewEncoder(params)

	encryptor1 := NewMKEncryptor(keys[0].publicKey, params, ids[0])
	encryptor2 := NewMKEncryptor(keys[1].publicKey, params, ids[1])
	decryptor := NewMKDecryptor(params)

	// encrypt
	values1, _, cipher1 := newTestVectors(params, encoder, encryptor1, complex(-1, -1), complex(1, 1), t)
	values2, _, cipher2 := newTestVectors(params, encoder, encryptor2, complex(-1, -1), complex(1, 1), t)

	// pad and add
	evaluator := NewMKEvaluator(params)

	out1, out2 := PadCiphers(cipher1, cipher2, params)

	resCipher := evaluator.Add(out1, out2)

	// decrypt
	partialDec1 := decryptor.PartDec(resCipher.ciphertexts.Value()[1], resCipher.ciphertexts.Level(), keys[0].secretKey)
	partialDec2 := decryptor.PartDec(resCipher.ciphertexts.Value()[2], resCipher.ciphertexts.Level(), keys[1].secretKey)

	decrypted := decryptor.MergeDec(resCipher.ciphertexts.Value()[0], resCipher.ciphertexts.Scale(), resCipher.ciphertexts.Level(), []*ring.Poly{partialDec1, partialDec2})

	// perform the operation in the plaintext space
	for i := 0; i < len(values1); i++ {
		values1[i] += values2[i]
	}

	// check results
	verifyTestVectors(params, encoder, values1, decrypted, t)
}

func Test_AddFourParticipants(t *testing.T) {

	ids := []uint64{1, 2, 5, 7}

	keys, params := setupPeers(ids, 0)

	encoder := ckks.NewEncoder(params)

	encryptor1 := NewMKEncryptor(keys[0].publicKey, params, ids[0])
	encryptor2 := NewMKEncryptor(keys[1].publicKey, params, ids[1])
	encryptor3 := NewMKEncryptor(keys[2].publicKey, params, ids[2])
	encryptor4 := NewMKEncryptor(keys[3].publicKey, params, ids[3])
	decryptor := NewMKDecryptor(params)

	// encrypt
	values1, _, cipher1 := newTestVectors(params, encoder, encryptor1, complex(-1, -1), complex(1, 1), t)
	values2, _, cipher2 := newTestVectors(params, encoder, encryptor2, complex(-1, -1), complex(1, 1), t)
	values3, _, cipher3 := newTestVectors(params, encoder, encryptor3, complex(-1, -1), complex(1, 1), t)
	values4, _, cipher4 := newTestVectors(params, encoder, encryptor4, complex(-1, -1), complex(1, 1), t)

	// pad and add in 2 steps
	evaluator := NewMKEvaluator(params)

	out1, out2 := PadCiphers(cipher1, cipher2, params)
	out3, out4 := PadCiphers(cipher3, cipher4, params)

	resCipher1 := evaluator.Add(out1, out2)
	resCipher2 := evaluator.Add(out3, out4)

	finalOut1, finalOut2 := PadCiphers(resCipher1, resCipher2, params)

	resCipher := evaluator.Add(finalOut1, finalOut2)

	// decrypt
	partialDec1 := decryptor.PartDec(resCipher.ciphertexts.Value()[1], resCipher.ciphertexts.Level(), keys[0].secretKey)
	partialDec2 := decryptor.PartDec(resCipher.ciphertexts.Value()[2], resCipher.ciphertexts.Level(), keys[1].secretKey)
	partialDec3 := decryptor.PartDec(resCipher.ciphertexts.Value()[3], resCipher.ciphertexts.Level(), keys[2].secretKey)
	partialDec4 := decryptor.PartDec(resCipher.ciphertexts.Value()[4], resCipher.ciphertexts.Level(), keys[3].secretKey)

	decrypted := decryptor.MergeDec(resCipher.ciphertexts.Value()[0], resCipher.ciphertexts.Scale(), resCipher.ciphertexts.Level(), []*ring.Poly{partialDec1, partialDec2, partialDec3, partialDec4})

	// perform the operation in the plaintext space
	for i := 0; i < len(values1); i++ {
		values1[i] += values2[i] + values3[i] + values4[i]
	}

	// check results
	verifyTestVectors(params, encoder, values1, decrypted, t)

}

/*
func Test_Mul(t *testing.T) {

	ids := []uint64{1, 2}

	keys, params := setupPeers(ids, 0)

	ringQ := GetRingQ(params)
	ringT := getRingT(params)

	encoder := ckks.NewEncoder(params)

	encryptor1 := NewMKEncryptor(keys[0].publicKey, params, ids[0])
	encryptor2 := NewMKEncryptor(keys[1].publicKey, params, ids[1])
	decryptor := NewMKDecryptor(params)

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

	encoder := ckks.NewEncoder(params)

	//create and encrypt 2 plaintext
	encryptor1 := NewMKEncryptor(keys[0].publicKey, params, ids[0])
	encryptor2 := NewMKEncryptor(keys[1].publicKey, params, ids[1])

	_, _, mkCiphertext1 := newTestVectors(params, encoder, encryptor1, complex(-1, -1), complex(1, 1), t)
	_, _, mkCiphertext2 := newTestVectors(params, encoder, encryptor2, complex(-1, -1), complex(1, 1), t)

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
func setupPeers(peerIDs []uint64, paramsNbr int) ([]*MKKeys, *ckks.Parameters) {

	res := make([]*MKKeys, len(peerIDs))

	params := ckks.DefaultParams[paramsNbr]

	for i := 0; i < len(peerIDs); i++ {

		// setup keys and public parameters
		a := GenCommonPublicParam(params)

		res[i] = KeyGen(params, peerIDs[i], a)

	}

	return res, params
}

func newTestVectors(params *ckks.Parameters, encoder ckks.Encoder, encryptor MKEncryptor, a, b complex128, t *testing.T) (values []complex128, plaintext *ckks.Plaintext, ciphertext *MKCiphertext) {

	logSlots := params.LogSlots()

	values = make([]complex128, 1<<logSlots)

	for i := uint64(0); i < 1<<logSlots; i++ {
		values[i] = complex(utils.RandFloat64(real(a), real(b)), utils.RandFloat64(imag(a), imag(b)))
	}

	values[0] = complex(0.607538, 0)

	plaintext = encoder.EncodeNTTAtLvlNew(params.MaxLevel(), values, logSlots)

	if encryptor != nil {
		ciphertext = encryptor.EncryptMK(plaintext)
	}

	return values, plaintext, ciphertext
}

func verifyTestVectors(params *ckks.Parameters, encoder ckks.Encoder, valuesWant []complex128, plaintext *ckks.Plaintext, t *testing.T) {

	precStats := GetPrecisionStats(params, encoder, valuesWant, plaintext)

	if *printPrecisionStats {
		t.Log(precStats.String())
	}

	require.GreaterOrEqual(t, real(precStats.MeanPrecision), minPrec)
	require.GreaterOrEqual(t, imag(precStats.MeanPrecision), minPrec)
}

// GetPrecisionStats generates a PrecisionStats struct from the reference values and the decrypted values
func GetPrecisionStats(params *ckks.Parameters, encoder ckks.Encoder, valuesWant []complex128, plaintext *ckks.Plaintext) (prec ckks.PrecisionStats) {

	logSlots := params.LogSlots()
	slots := uint64(1 << logSlots)

	valuesTest := encoder.Decode(plaintext, logSlots)

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
