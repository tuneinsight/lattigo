package bfv

import (
	"fmt"
	"github.com/ldsec/lattigo/ring"
	"testing"
)

type BFVTESTPARAMS struct {
	bfvcontext   *Context
	batchencoder *BatchEncoder
	kgen         *KeyGenerator
	sk           *SecretKey
	pk           *PublicKey
	encryptorSk  *Encryptor
	encryptorPk  *Encryptor
	decryptor    *Decryptor
	evaluator    *Evaluator
}

func TestBFV(t *testing.T) {

	var err error

	paramSets := DefaultParams[1:2]

	bitDecomps := []uint64{60}

	for _, params := range paramSets {

		bfvTest := new(BFVTESTPARAMS)

		bfvTest.bfvcontext = NewContext()
		if err := bfvTest.bfvcontext.SetParameters(&params); err != nil {
			t.Error(err)
		}

		bfvTest.kgen = bfvTest.bfvcontext.NewKeyGenerator()

		if bfvTest.batchencoder, err = bfvTest.bfvcontext.NewBatchEncoder(); err != nil {
			t.Error(err)
		}

		bfvTest.sk, bfvTest.pk = bfvTest.kgen.NewKeyPair()

		if bfvTest.decryptor, err = bfvTest.bfvcontext.NewDecryptor(bfvTest.sk); err != nil {
			t.Error(err)
		}

		if bfvTest.encryptorPk, err = bfvTest.bfvcontext.NewEncryptorFromPk(bfvTest.pk); err != nil {
			t.Error(err)
		}

		if bfvTest.encryptorSk, err = bfvTest.bfvcontext.NewEncryptorFromSk(bfvTest.sk); err != nil {
			t.Error(err)
		}

		bfvTest.evaluator = bfvTest.bfvcontext.NewEvaluator()

		testEncodeDecode(bfvTest, t)
		testPlaintextBatchEncodeDecode(bfvTest, t)
		testEncryptDecrypt(bfvTest, t)
		testHomomorphicAddition(bfvTest, t)
		testAddWithPlaintext(bfvTest, t)
		testHomomorphicSubtraction(bfvTest, t)
		testSubWithPlaintext(bfvTest, t)
		testHomomorphicMultiplication(bfvTest, t)
		testMulWithPlaintext(bfvTest, t)
		testRelinearization(bfvTest, bitDecomps, t)
		testKeySwitching(bfvTest, bitDecomps, t)
		testGaloisEnd(bfvTest, bitDecomps, t)
		testMarshaler(bfvTest, t)

	}
}

//=============================================================
//=== TESTS FOR ENCODE AND DECODE AN INTEGER ON A PLAINTEXT ===
//=============================================================

func testEncodeDecode(bfvTest *BFVTESTPARAMS, t *testing.T) {

	var base int64
	var value int64
	var sign uint32

	T := bfvTest.bfvcontext.t
	bfvcontext := bfvTest.bfvcontext

	t.Run(fmt.Sprintf("N=%d/T=%d/EncodeDecodeInt", bfvTest.bfvcontext.n, bfvTest.bfvcontext.t), func(t *testing.T) {

		for i := 0; i < 0xFF; i++ {
			//Generates a new plaintext with 0 coefficients
			plaintext := bfvcontext.NewPlaintext()

			//Generates a base between 2 and T/2
			base = int64(ring.RandUniform(256, 255) % (T >> 1))
			if base < 2 {
				base = 2
			}

			encoder := bfvcontext.NewIntEncoder(base)

			//Generates a value which can fit in the ring given the base
			value = int64(ring.RandUniform(0xFFFFFFFF+1, 0xFFFFFFFF))

			if T != 2 {
				sign = uint32(ring.RandUniform(2, 3) & 1)
				if sign&1 == 1 {
					value *= -1
				}
			}

			//Encodes the value on the ring
			encoder.Encode(value, plaintext)

			//Tests the decoding
			if value != encoder.Decode(plaintext) {
				t.Errorf("Error in Encoding/Decoding Function, have %v, want %v, base %v", encoder.Decode(plaintext), value, base)
				break
			}
		}
	})
}

func testMarshaler(bfvTest *BFVTESTPARAMS, t *testing.T) {

	state := true

	bfvContext := bfvTest.bfvcontext
	Sk := bfvTest.sk
	Pk := bfvTest.pk

	t.Run(fmt.Sprintf("N=%d/T=%d/Qi=%dlimbs/MarshalSk",
		bfvTest.bfvcontext.n,
		bfvTest.bfvcontext.t,
		len(bfvTest.bfvcontext.contextQ.Modulus),
	), func(t *testing.T) {
		state = true

		SkBytes, err := Sk.MarshalBinary()
		if err != nil {

		}

		SkTest := bfvTest.kgen.NewSecretKeyEmpty()
		SkTest.UnMarshalBinary(SkBytes)

		if bfvContext.contextQ.Equal(Sk.sk, SkTest.sk) != true {
			t.Errorf("error : binarymarshal secretkey")
		}

	})

	t.Run(fmt.Sprintf("N=%d/T=%d/Qi=%dlimbs/MarshalPk", bfvTest.bfvcontext.n,
		bfvTest.bfvcontext.t,
		len(bfvTest.bfvcontext.contextQ.Modulus)), func(t *testing.T) {
		PkBytes, err := Pk.MarshalBinary()
		if err != nil {
			t.Error(err)
		}

		PkTest := bfvTest.kgen.NewPublicKeyEmpty()
		PkTest.UnMarshalBinary(PkBytes)

		for i := range Pk.pk {
			if bfvContext.contextQ.Equal(Pk.pk[i], PkTest.pk[i]) != true {
				t.Errorf("error : binarymarshal publickey")
				break
			}
		}
	})

	t.Run(fmt.Sprintf("N=%d/T=%d/Qi=%dlimbs/MarshalCiphertext", bfvTest.bfvcontext.n,
		bfvTest.bfvcontext.t,
		len(bfvTest.bfvcontext.contextQ.Modulus)), func(t *testing.T) {

		Ctx := bfvContext.NewRandomCiphertext(4)
		CtxBytes, err := Ctx.MarshalBinary()
		if err != nil {
			t.Error(err)
		}

		CtxTest := bfvContext.NewCiphertext(4)
		if err = CtxTest.UnMarshalBinary(CtxBytes); err != nil {
			t.Error(err)
		}

		for i := range Ctx.Value() {
			if bfvContext.contextQ.Equal(CtxTest.Value()[i], Ctx.Value()[i]) != true {
				t.Errorf("error : binarymarshal ciphertext")
				break
			}
		}

	})

	t.Run(fmt.Sprintf("N=%d/T=%d/Qi=%dlimbs/bitDecomp=%d/Marshalrlk", bfvTest.bfvcontext.n,
		bfvTest.bfvcontext.t,
		len(bfvTest.bfvcontext.contextQ.Modulus),
		15), func(t *testing.T) {

		rlk := bfvTest.kgen.NewRelinKey(bfvTest.sk, 5, 15)

		rlkBytes, err := rlk.MarshalBinary()
		if err != nil {
			t.Error(err)
		}

		rlkTest := bfvTest.kgen.NewRelinKeyEmpty(5, 15)
		rlkTest.UnMarshalBinary(rlkBytes)

		state = true

		for i := range rlkTest.evakey {

			if state != true {
				break
			}

			for j := range rlkTest.evakey[i].evakey {

				if state != true {
					break
				}

				for x := range rlkTest.evakey[i].evakey[j] {

					if bfvContext.contextQ.Equal(rlkTest.evakey[i].evakey[j][x][0], rlk.evakey[i].evakey[j][x][0]) != true {
						t.Errorf("error : binarymarshal rlk")
						state = false
						break
					}

					if bfvContext.contextQ.Equal(rlkTest.evakey[i].evakey[j][x][1], rlk.evakey[i].evakey[j][x][1]) != true {
						t.Errorf("error : binarymarshal rlk")
						state = false
						break
					}
				}
			}
		}
	})

	t.Run(fmt.Sprintf("N=%d/T=%d/Qi=%dlimbs/bitDecomp=%d/MarshalRotKey", bfvTest.bfvcontext.n,
		bfvTest.bfvcontext.t,
		len(bfvTest.bfvcontext.contextQ.Modulus),
		15), func(t *testing.T) {

		rotKey := bfvTest.kgen.NewRotationKeysPow2(Sk, 15, true)

		rotKeyBytes, err := rotKey.MarshalBinary()
		if err != nil {
			t.Error(err)
		}

		rotKeyTest := bfvTest.kgen.NewRotationKeysEmpty()
		rotKeyTest.UnMarshalBinary(rotKeyBytes)

		state = true

		if rotKeyTest.evakeyRotRows != nil {
			for j := range rotKeyTest.evakeyRotRows.evakey {
				for x := range rotKeyTest.evakeyRotRows.evakey[j] {
					if bfvContext.contextQ.Equal(rotKeyTest.evakeyRotRows.evakey[j][x][0], rotKey.evakeyRotRows.evakey[j][x][0]) != true {
						t.Errorf("error : binarymarshal rotation key rot row")
						state = false
						break
					}

					if bfvContext.contextQ.Equal(rotKeyTest.evakeyRotRows.evakey[j][x][1], rotKey.evakeyRotRows.evakey[j][x][1]) != true {
						t.Errorf("error : binarymarshal rotation key rot row")
						state = false
						break
					}

				}
			}
		}

		for k, v := range rotKeyTest.evakeyRotColLeft {
			for j := range v.evakey {
				for x := range v.evakey[j] {
					if bfvContext.contextQ.Equal(v.evakey[j][x][0], rotKey.evakeyRotColLeft[k].evakey[j][x][0]) != true {
						t.Errorf("error : binarymarshal rotation key col L")
						state = false
						break
					}

					if bfvContext.contextQ.Equal(v.evakey[j][x][1], rotKey.evakeyRotColLeft[k].evakey[j][x][1]) != true {
						t.Errorf("error : binarymarshal rotation key col L")
						state = false
						break
					}
				}
			}

		}

		for k, v := range rotKeyTest.evakeyRotColRight {
			for j := range v.evakey {
				for x := range v.evakey[j] {
					if bfvContext.contextQ.Equal(v.evakey[j][x][0], rotKey.evakeyRotColRight[k].evakey[j][x][0]) != true {
						t.Errorf("error : binarymarshal rotation key col R")
						state = false
						break
					}

					if bfvContext.contextQ.Equal(v.evakey[j][x][1], rotKey.evakeyRotColRight[k].evakey[j][x][1]) != true {
						t.Errorf("error : binarymarshal rotation key col R")
						state = false
						break
					}
				}
			}

		}
	})
}

//=====================================
//=== TESTS FOR ENCRYPT AND DECRYPT ===
//=====================================

func testPlaintextBatchEncodeDecode(bfvTest *BFVTESTPARAMS, t *testing.T) {

	t.Run(fmt.Sprintf("N=%d/T=%d/logQ=%d/logP=%d/BatchEncodeDecode", bfvTest.bfvcontext.N(),
		bfvTest.bfvcontext.T(),
		bfvTest.bfvcontext.LogQ(),
		bfvTest.bfvcontext.LogP()), func(t *testing.T) {

		bfvContext := bfvTest.bfvcontext
		encoder := bfvTest.batchencoder

		coeffsWant := bfvContext.contextT.NewUniformPoly()

		plaintextWant := bfvContext.NewPlaintext()

		encoder.EncodeUint(coeffsWant.Coeffs[0], plaintextWant)

		if equalslice(coeffsWant.Coeffs[0], encoder.DecodeUint(plaintextWant)) != true {
			t.Errorf("error : plaintext lift/rescale")
		}

	})

}

func newTestVectors(bfvTest *BFVTESTPARAMS) (coeffs *ring.Poly, plaintext *Plaintext, ciphertext *Ciphertext, err error) {

	coeffs = bfvTest.bfvcontext.contextT.NewUniformPoly()

	plaintext = bfvTest.bfvcontext.NewPlaintext()

	if err = bfvTest.batchencoder.EncodeUint(coeffs.Coeffs[0], plaintext); err != nil {
		return nil, nil, nil, err
	}

	if ciphertext, err = bfvTest.encryptorPk.EncryptNew(plaintext); err != nil {
		return nil, nil, nil, err
	}

	return coeffs, plaintext, ciphertext, nil
}

func verifyTestVectors(bfvTest *BFVTESTPARAMS, coeffs *ring.Poly, element Operand, t *testing.T) {

	var coeffsTest []uint64

	el := element.Element()

	if el.Degree() == 0 {

		coeffsTest = bfvTest.batchencoder.DecodeUint(el.Plaintext())

	} else {

		coeffsTest = bfvTest.batchencoder.DecodeUint(bfvTest.decryptor.DecryptNew(el.Ciphertext()))
	}

	if equalslice(coeffs.Coeffs[0], coeffsTest) != true {
		t.Errorf("decryption error")
	}
}

func testEncryptDecrypt(bfvTest *BFVTESTPARAMS, t *testing.T) {

	t.Run(fmt.Sprintf("N=%d/T=%d/logQ=%d/logP=%d/EncryptPk", bfvTest.bfvcontext.N(),
		bfvTest.bfvcontext.T(),
		bfvTest.bfvcontext.LogQ(),
		bfvTest.bfvcontext.LogP()), func(t *testing.T) {

		coeffs, _, ciphertext, _ := newTestVectors(bfvTest)
		verifyTestVectors(bfvTest, coeffs, ciphertext, t)
	})

	t.Run(fmt.Sprintf("N=%d/T=%d/logQ=%d/logP=%d/EncryptSk", bfvTest.bfvcontext.N(),
		bfvTest.bfvcontext.T(),
		bfvTest.bfvcontext.LogQ(),
		bfvTest.bfvcontext.LogP()), func(t *testing.T) {

		var err error
		var coeffs *ring.Poly
		var plaintext *Plaintext
		var ciphertext *Ciphertext

		coeffs = bfvTest.bfvcontext.contextT.NewUniformPoly()

		plaintext = bfvTest.bfvcontext.NewPlaintext()

		if err = bfvTest.batchencoder.EncodeUint(coeffs.Coeffs[0], plaintext); err != nil {
			t.Error(err)
		}

		if ciphertext, err = bfvTest.encryptorSk.EncryptNew(plaintext); err != nil {
			t.Error(err)
		}

		verifyTestVectors(bfvTest, coeffs, ciphertext, t)
	})
}

func testHomomorphicAddition(bfvTest *BFVTESTPARAMS, t *testing.T) {

	bfvContext := bfvTest.bfvcontext
	evaluator := bfvTest.evaluator

	var ciphertextTest *Ciphertext
	var err error

	t.Run(fmt.Sprintf("N=%d/T=%d/logQ=%d/logP=%d/Add", bfvTest.bfvcontext.N(),
		bfvTest.bfvcontext.T(),
		bfvTest.bfvcontext.LogQ(),
		bfvTest.bfvcontext.LogP()), func(t *testing.T) {

		coeffs0, _, ciphertext0, _ := newTestVectors(bfvTest)
		coeffs1, _, ciphertext1, _ := newTestVectors(bfvTest)

		for i := 0; i < 1; i++ {
			if err = evaluator.Add(ciphertext0, ciphertext1, ciphertext0); err != nil {
				t.Error(err)
			}
			bfvContext.contextT.Add(coeffs0, coeffs1, coeffs0)
		}

		verifyTestVectors(bfvTest, coeffs0, ciphertext0, t)
	})

	t.Run(fmt.Sprintf("N=%d/T=%d/logQ=%d/logP=%d/AddNew", bfvTest.bfvcontext.N(),
		bfvTest.bfvcontext.T(),
		bfvTest.bfvcontext.LogQ(),
		bfvTest.bfvcontext.LogP()), func(t *testing.T) {

		coeffs0, _, ciphertext0, _ := newTestVectors(bfvTest)
		coeffs1, _, ciphertext1, _ := newTestVectors(bfvTest)

		for i := 0; i < 1; i++ {
			if ciphertextTest, err = evaluator.AddNew(ciphertext0, ciphertext1); err != nil {
				t.Error(err)
			}

			bfvContext.contextT.Add(coeffs0, coeffs1, coeffs0)
		}

		verifyTestVectors(bfvTest, coeffs0, ciphertextTest, t)

	})

	t.Run(fmt.Sprintf("N=%d/T=%d/logQ=%d/logP=%d/AddNoMod", bfvTest.bfvcontext.N(),
		bfvTest.bfvcontext.T(),
		bfvTest.bfvcontext.LogQ(),
		bfvTest.bfvcontext.LogP()), func(t *testing.T) {

		coeffs0, _, ciphertext0, _ := newTestVectors(bfvTest)
		coeffs1, _, ciphertext1, _ := newTestVectors(bfvTest)

		for i := 0; i < 1; i++ {

			if err = evaluator.AddNoMod(ciphertext0, ciphertext1, ciphertext0); err != nil {
				t.Error(err)
			}

			evaluator.Reduce(ciphertext0, ciphertext0)
			bfvContext.contextT.Add(coeffs0, coeffs1, coeffs0)

		}

		verifyTestVectors(bfvTest, coeffs0, ciphertext0, t)
	})

	t.Run(fmt.Sprintf("N=%d/T=%d/logQ=%d/logP=%d/AddNoModNew", bfvTest.bfvcontext.N(),
		bfvTest.bfvcontext.T(),
		bfvTest.bfvcontext.LogQ(),
		bfvTest.bfvcontext.LogP()), func(t *testing.T) {

		coeffs0, _, ciphertext0, _ := newTestVectors(bfvTest)
		coeffs1, _, ciphertext1, _ := newTestVectors(bfvTest)

		for i := 0; i < 1; i++ {
			if ciphertextTest, err = evaluator.AddNoModNew(ciphertext0, ciphertext1); err != nil {
				t.Error(err)
			}
			evaluator.Reduce(ciphertextTest, ciphertextTest)
			bfvContext.contextT.Add(coeffs0, coeffs1, coeffs0)
		}

		verifyTestVectors(bfvTest, coeffs0, ciphertextTest, t)
	})
}

//HOMOMORPHIC MULTIPLICATION WITH RELINEARIZATION BASE Qi AND w
func testAddWithPlaintext(bfvTest *BFVTESTPARAMS, t *testing.T) {

	bfvContext := bfvTest.bfvcontext
	evaluator := bfvTest.evaluator

	var ciphertextTest *Ciphertext
	var err error

	t.Run(fmt.Sprintf("N=%d/T=%d/logQ=%d/logP=%d/AddWithPlaintext", bfvTest.bfvcontext.N(),
		bfvTest.bfvcontext.T(),
		bfvTest.bfvcontext.LogQ(),
		bfvTest.bfvcontext.LogP()), func(t *testing.T) {

		coeffs0, _, ciphertext0, _ := newTestVectors(bfvTest)
		coeffs1, plaintext1, _, _ := newTestVectors(bfvTest)

		if err := evaluator.Add(ciphertext0, plaintext1, ciphertext0); err != nil {
			t.Error(err)
		}

		bfvContext.contextT.Add(coeffs0, coeffs1, coeffs0)

		verifyTestVectors(bfvTest, coeffs0, ciphertext0, t)

	})

	t.Run(fmt.Sprintf("N=%d/T=%d/logQ=%d/logP=%d/AddWithPlaintextNew", bfvTest.bfvcontext.N(),
		bfvTest.bfvcontext.T(),
		bfvTest.bfvcontext.LogQ(),
		bfvTest.bfvcontext.LogP()), func(t *testing.T) {

		coeffs0, _, ciphertext0, _ := newTestVectors(bfvTest)
		coeffs1, plaintext1, _, _ := newTestVectors(bfvTest)

		if ciphertextTest, err = evaluator.AddNew(ciphertext0, plaintext1); err != nil {
			t.Error(err)
		}
		bfvContext.contextT.Add(coeffs0, coeffs1, coeffs0)

		verifyTestVectors(bfvTest, coeffs0, ciphertextTest, t)

	})

	t.Run(fmt.Sprintf("N=%d/T=%d/logQ=%d/logP=%d/AddWithPlaintextNoMod", bfvTest.bfvcontext.N(),
		bfvTest.bfvcontext.T(),
		bfvTest.bfvcontext.LogQ(),
		bfvTest.bfvcontext.LogP()), func(t *testing.T) {

		coeffs0, _, ciphertext0, _ := newTestVectors(bfvTest)
		coeffs1, plaintext1, _, _ := newTestVectors(bfvTest)

		if err := evaluator.AddNoMod(ciphertext0, plaintext1, ciphertext0); err != nil {
			t.Error(err)
		}
		evaluator.Reduce(ciphertext0, ciphertext0)
		bfvContext.contextT.Add(coeffs0, coeffs1, coeffs0)

		verifyTestVectors(bfvTest, coeffs0, ciphertext0, t)
	})

	t.Run(fmt.Sprintf("N=%d/T=%d/logQ=%d/logP=%d/AddWithPlaintextNoModNew", bfvTest.bfvcontext.N(),
		bfvTest.bfvcontext.T(),
		bfvTest.bfvcontext.LogQ(),
		bfvTest.bfvcontext.LogP()), func(t *testing.T) {

		coeffs0, _, ciphertext0, _ := newTestVectors(bfvTest)
		coeffs1, plaintext1, _, _ := newTestVectors(bfvTest)

		if ciphertextTest, err = evaluator.AddNoModNew(ciphertext0, plaintext1); err != nil {
			t.Error(err)
		}
		evaluator.Reduce(ciphertextTest, ciphertextTest)
		bfvContext.contextT.Add(coeffs0, coeffs1, coeffs0)

		verifyTestVectors(bfvTest, coeffs0, ciphertextTest, t)
	})
}

func testHomomorphicSubtraction(bfvTest *BFVTESTPARAMS, t *testing.T) {

	bfvContext := bfvTest.bfvcontext
	evaluator := bfvTest.evaluator

	var ciphertextTest *Ciphertext
	var err error

	t.Run(fmt.Sprintf("N=%d/T=%d/logQ=%d/logP=%d/Sub", bfvTest.bfvcontext.N(),
		bfvTest.bfvcontext.T(),
		bfvTest.bfvcontext.LogQ(),
		bfvTest.bfvcontext.LogP()), func(t *testing.T) {

		coeffs0, _, ciphertext0, _ := newTestVectors(bfvTest)
		coeffs1, _, ciphertext1, _ := newTestVectors(bfvTest)

		for i := 0; i < 1; i++ {
			if err := evaluator.Sub(ciphertext0, ciphertext1, ciphertext0); err != nil {
				t.Error(err)
			}
			bfvContext.contextT.Sub(coeffs0, coeffs1, coeffs0)
		}

		verifyTestVectors(bfvTest, coeffs0, ciphertext0, t)

	})

	t.Run(fmt.Sprintf("N=%d/T=%d/logQ=%d/logP=%d/SubNew", bfvTest.bfvcontext.N(),
		bfvTest.bfvcontext.T(),
		bfvTest.bfvcontext.LogQ(),
		bfvTest.bfvcontext.LogP()), func(t *testing.T) {

		coeffs0, _, ciphertext0, _ := newTestVectors(bfvTest)
		coeffs1, _, ciphertext1, _ := newTestVectors(bfvTest)

		for i := 0; i < 1; i++ {
			if ciphertextTest, err = evaluator.SubNew(ciphertext0, ciphertext1); err != nil {
				t.Error(err)
			}
			bfvContext.contextT.Sub(coeffs0, coeffs1, coeffs0)
		}

		verifyTestVectors(bfvTest, coeffs0, ciphertextTest, t)
	})

	t.Run(fmt.Sprintf("N=%d/T=%d/logQ=%d/logP=%d/SubNoMod", bfvTest.bfvcontext.N(),
		bfvTest.bfvcontext.T(),
		bfvTest.bfvcontext.LogQ(),
		bfvTest.bfvcontext.LogP()), func(t *testing.T) {

		coeffs0, _, ciphertext0, _ := newTestVectors(bfvTest)
		coeffs1, _, ciphertext1, _ := newTestVectors(bfvTest)

		for i := 0; i < 1; i++ {
			if err := evaluator.SubNoMod(ciphertext0, ciphertext1, ciphertext0); err != nil {
				t.Error(err)
			}
			evaluator.Reduce(ciphertext0, ciphertext0)
			bfvContext.contextT.Sub(coeffs0, coeffs1, coeffs0)
		}

		verifyTestVectors(bfvTest, coeffs0, ciphertext0, t)

	})

	t.Run(fmt.Sprintf("N=%d/T=%d/logQ=%d/logP=%d/SubNoModNew", bfvTest.bfvcontext.N(),
		bfvTest.bfvcontext.T(),
		bfvTest.bfvcontext.LogQ(),
		bfvTest.bfvcontext.LogP()), func(t *testing.T) {

		coeffs0, _, ciphertext0, _ := newTestVectors(bfvTest)
		coeffs1, _, ciphertext1, _ := newTestVectors(bfvTest)

		for i := 0; i < 1; i++ {
			if ciphertextTest, err = evaluator.SubNoModNew(ciphertext0, ciphertext1); err != nil {
				t.Error(err)
			}
			if err != nil {
				t.Error(err)
			}
			evaluator.Reduce(ciphertextTest, ciphertextTest)
			bfvContext.contextT.Sub(coeffs0, coeffs1, coeffs0)
		}

		verifyTestVectors(bfvTest, coeffs0, ciphertextTest, t)
	})
}

//HOMOMORPHIC MULTIPLICATION WITH RELINEARIZATION BASE Qi AND w
func testSubWithPlaintext(bfvTest *BFVTESTPARAMS, t *testing.T) {

	bfvContext := bfvTest.bfvcontext
	evaluator := bfvTest.evaluator

	var ciphertextTest *Ciphertext
	var err error

	t.Run(fmt.Sprintf("N=%d/T=%d/logQ=%d/logP=%d/SubWithPlaintext", bfvTest.bfvcontext.N(),
		bfvTest.bfvcontext.T(),
		bfvTest.bfvcontext.LogQ(),
		bfvTest.bfvcontext.LogP()), func(t *testing.T) {

		coeffs0, _, ciphertext0, _ := newTestVectors(bfvTest)
		coeffs1, plaintext1, _, _ := newTestVectors(bfvTest)

		if err := evaluator.Sub(ciphertext0, plaintext1, ciphertext0); err != nil {
			t.Error(err)
		}
		bfvContext.contextT.Sub(coeffs0, coeffs1, coeffs0)

		verifyTestVectors(bfvTest, coeffs0, ciphertext0, t)
	})

	t.Run(fmt.Sprintf("N=%d/T=%d/logQ=%d/logP=%d/SubWithPlaintextNew", bfvTest.bfvcontext.N(),
		bfvTest.bfvcontext.T(),
		bfvTest.bfvcontext.LogQ(),
		bfvTest.bfvcontext.LogP()), func(t *testing.T) {

		coeffs0, _, ciphertext0, _ := newTestVectors(bfvTest)
		coeffs1, plaintext1, _, _ := newTestVectors(bfvTest)

		if ciphertextTest, err = evaluator.SubNew(ciphertext0, plaintext1); err != nil {
			t.Error(err)
		}
		bfvContext.contextT.Sub(coeffs0, coeffs1, coeffs0)

		verifyTestVectors(bfvTest, coeffs0, ciphertextTest, t)

	})

	t.Run(fmt.Sprintf("N=%d/T=%d/logQ=%d/logP=%d/SubWithPlaintextNoMod", bfvTest.bfvcontext.N(),
		bfvTest.bfvcontext.T(),
		bfvTest.bfvcontext.LogQ(),
		bfvTest.bfvcontext.LogP()), func(t *testing.T) {

		coeffs0, _, ciphertext0, _ := newTestVectors(bfvTest)
		coeffs1, plaintext1, _, _ := newTestVectors(bfvTest)

		if err := evaluator.SubNoMod(ciphertext0, plaintext1, ciphertext0); err != nil {
			t.Error(err)
		}
		evaluator.Reduce(ciphertext0, ciphertext0)
		bfvContext.contextT.Sub(coeffs0, coeffs1, coeffs0)

		verifyTestVectors(bfvTest, coeffs0, ciphertext0, t)
	})

	t.Run(fmt.Sprintf("N=%d/T=%d/logQ=%d/logP=%d/SubWithPlaintextNoModNew", bfvTest.bfvcontext.N(),
		bfvTest.bfvcontext.T(),
		bfvTest.bfvcontext.LogQ(),
		bfvTest.bfvcontext.LogP()), func(t *testing.T) {

		coeffs0, _, ciphertext0, _ := newTestVectors(bfvTest)
		coeffs1, plaintext1, _, _ := newTestVectors(bfvTest)

		if ciphertextTest, err = evaluator.SubNoModNew(ciphertext0, plaintext1); err != nil {
			t.Error(err)
		}
		evaluator.Reduce(ciphertextTest, ciphertextTest)
		bfvContext.contextT.Sub(coeffs0, coeffs1, coeffs0)

		verifyTestVectors(bfvTest, coeffs0, ciphertextTest, t)
	})
}

//HOMOMORPHIC MULTIPLICATION
func testHomomorphicMultiplication(bfvTest *BFVTESTPARAMS, t *testing.T) {

	bfvContext := bfvTest.bfvcontext
	evaluator := bfvTest.evaluator

	t.Run(fmt.Sprintf("N=%d/T=%d/logQ=%d/logP=%d/Mul", bfvTest.bfvcontext.N(),
		bfvTest.bfvcontext.T(),
		bfvTest.bfvcontext.LogQ(),
		bfvTest.bfvcontext.LogP()), func(t *testing.T) {

		coeffs0, _, ciphertext0, _ := newTestVectors(bfvTest)
		coeffs1, _, ciphertext1, _ := newTestVectors(bfvTest)

		receiverCiphertext := bfvContext.NewCiphertextBig(ciphertext0.Degree() + ciphertext1.Degree())

		if err := evaluator.Mul(ciphertext0, ciphertext1, receiverCiphertext); err != nil {
			t.Error(err)
		}

		bfvContext.contextT.MulCoeffs(coeffs0, coeffs1, coeffs0)

		verifyTestVectors(bfvTest, coeffs0, receiverCiphertext, t)
	})

	t.Run(fmt.Sprintf("N=%d/T=%d/logQ=%d/logP=%d/MulNew", bfvTest.bfvcontext.N(),
		bfvTest.bfvcontext.T(),
		bfvTest.bfvcontext.LogQ(),
		bfvTest.bfvcontext.LogP()), func(t *testing.T) {

		coeffs0, _, ciphertext0, _ := newTestVectors(bfvTest)
		coeffs1, _, ciphertext1, _ := newTestVectors(bfvTest)

		receiverCiphertext, _ := evaluator.MulNew(ciphertext0, ciphertext1)

		bfvContext.contextT.MulCoeffs(coeffs0, coeffs1, coeffs0)

		verifyTestVectors(bfvTest, coeffs0, receiverCiphertext, t)

	})

	t.Run(fmt.Sprintf("N=%d/T=%d/logQ=%d/logP=%d/Square", bfvTest.bfvcontext.N(),
		bfvTest.bfvcontext.T(),
		bfvTest.bfvcontext.LogQ(),
		bfvTest.bfvcontext.LogP()), func(t *testing.T) {

		coeffs0, _, ciphertext0, _ := newTestVectors(bfvTest)

		receiverCiphertext, _ := evaluator.MulNew(ciphertext0, ciphertext0)
		receiverCiphertext, _ = evaluator.MulNew(receiverCiphertext, receiverCiphertext)

		bfvContext.contextT.MulCoeffs(coeffs0, coeffs0, coeffs0)
		bfvContext.contextT.MulCoeffs(coeffs0, coeffs0, coeffs0)

		verifyTestVectors(bfvTest, coeffs0, receiverCiphertext, t)

	})

}

//HOMOMORPHIC MULTIPLICATION WITH RELINEARIZATION BASE Qi AND w
func testMulWithPlaintext(bfvTest *BFVTESTPARAMS, t *testing.T) {

	bfvContext := bfvTest.bfvcontext
	evaluator := bfvTest.evaluator

	t.Run(fmt.Sprintf("N=%d/T=%d/logQ=%d/logP=%d/MulWithPlaintext", bfvTest.bfvcontext.N(),
		bfvTest.bfvcontext.T(),
		bfvTest.bfvcontext.LogQ(),
		bfvTest.bfvcontext.LogP()), func(t *testing.T) {

		coeffs0, _, ciphertext0, _ := newTestVectors(bfvTest)
		coeffs1, plaintext1, _, _ := newTestVectors(bfvTest)

		receiverCiphertext := bfvContext.NewCiphertextBig(ciphertext0.Degree())

		if err := evaluator.Mul(ciphertext0, plaintext1, receiverCiphertext); err != nil {
			t.Error(err)
		}

		bfvContext.contextT.MulCoeffs(coeffs0, coeffs1, coeffs0)

		verifyTestVectors(bfvTest, coeffs0, receiverCiphertext, t)
	})

	t.Run(fmt.Sprintf("N=%d/T=%d/logQ=%d/logP=%d/MulWithPlaintextNew", bfvTest.bfvcontext.N(),
		bfvTest.bfvcontext.T(),
		bfvTest.bfvcontext.LogQ(),
		bfvTest.bfvcontext.LogP()), func(t *testing.T) {

		coeffs0, _, ciphertext0, _ := newTestVectors(bfvTest)
		coeffs1, plaintext1, _, _ := newTestVectors(bfvTest)

		receiverCiphertext, _ := evaluator.MulNew(ciphertext0, plaintext1)
		bfvContext.contextT.MulCoeffs(coeffs0, coeffs1, coeffs0)

		verifyTestVectors(bfvTest, coeffs0, receiverCiphertext, t)

	})
}

func testRelinearization(bfvTest *BFVTESTPARAMS, bitDecomps []uint64, t *testing.T) {

	bfvContext := bfvTest.bfvcontext
	kgen := bfvTest.kgen
	evaluator := bfvTest.evaluator

	for _, bitDecomp := range bitDecomps {

		rlk := kgen.NewRelinKey(bfvTest.sk, 2, bitDecomp)

		t.Run(fmt.Sprintf("N=%d/T=%d/logQ=%d/logP=%d/bitDecomp=%d/Relinearize", bfvTest.bfvcontext.N(),
			bfvTest.bfvcontext.T(),
			bfvTest.bfvcontext.LogQ(),
			bfvTest.bfvcontext.LogP(),
			bitDecomp), func(t *testing.T) {

			coeffs0, _, ciphertext0, _ := newTestVectors(bfvTest)
			coeffs1, _, ciphertext1, _ := newTestVectors(bfvTest)

			ciphertext0, _ = evaluator.MulNew(ciphertext0, ciphertext1)
			bfvContext.contextT.MulCoeffs(coeffs0, coeffs1, coeffs0)

			ciphertext0, _ = evaluator.MulNew(ciphertext0, ciphertext1)
			bfvContext.contextT.MulCoeffs(coeffs0, coeffs1, coeffs0)

			if err := evaluator.Relinearize(ciphertext0, rlk, ciphertext0); err != nil {
				t.Error(err)
			}

			if len(ciphertext0.Value()) > 2 {
				t.Errorf("error : Relinearize")
			}

			verifyTestVectors(bfvTest, coeffs0, ciphertext0, t)

		})

		t.Run(fmt.Sprintf("N=%d/T=%d/logQ=%d/logP=%d/bitDecomp=%d/RelinearizeNew", bfvTest.bfvcontext.N(),
			bfvTest.bfvcontext.T(),
			bfvTest.bfvcontext.LogQ(),
			bfvTest.bfvcontext.LogP(),
			bitDecomp), func(t *testing.T) {

			coeffs0, _, ciphertext0, _ := newTestVectors(bfvTest)
			coeffs1, _, ciphertext1, _ := newTestVectors(bfvTest)

			ciphertext0, _ = evaluator.MulNew(ciphertext0, ciphertext1)
			bfvContext.contextT.MulCoeffs(coeffs0, coeffs1, coeffs0)

			ciphertext0, _ = evaluator.MulNew(ciphertext0, ciphertext1)
			bfvContext.contextT.MulCoeffs(coeffs0, coeffs1, coeffs0)

			ciphertext0, err := evaluator.RelinearizeNew(ciphertext0, rlk)
			if err != nil {
				t.Error(err)
			}

			if len(ciphertext0.Value()) > 2 {
				t.Errorf("error : RelinearizeNew")
			}

			verifyTestVectors(bfvTest, coeffs0, ciphertext0, t)

		})
	}
}

func testKeySwitching(bfvTest *BFVTESTPARAMS, bitDecomps []uint64, t *testing.T) {

	bfvContext := bfvTest.bfvcontext
	encoder := bfvTest.batchencoder
	kgen := bfvTest.kgen
	Sk := bfvTest.sk
	evaluator := bfvTest.evaluator

	SkNew := kgen.NewSecretKey()

	decryptorSkNew, err := bfvContext.NewDecryptor(SkNew)
	if err != nil {
		t.Error(err)
	}

	for _, bitDecomp := range bitDecomps {

		switchingKey := kgen.NewSwitchingKey(Sk, SkNew, bitDecomp)

		t.Run(fmt.Sprintf("N=%d/T=%d/logQ=%d/logP=%d/bitDecomp=%d/SwitchKeys", bfvTest.bfvcontext.N(),
			bfvTest.bfvcontext.T(),
			bfvTest.bfvcontext.LogQ(),
			bfvTest.bfvcontext.LogP(),
			bitDecomp), func(t *testing.T) {

			coeffs0, _, ciphertext0, _ := newTestVectors(bfvTest)

			if err := evaluator.SwitchKeys(ciphertext0, switchingKey, ciphertext0); err != nil {
				t.Error(err)
			}

			if equalslice(coeffs0.Coeffs[0], encoder.DecodeUint(decryptorSkNew.DecryptNew(ciphertext0))) != true {
				t.Errorf("error : switchingKey encrypt/decrypt")
			}

		})

		t.Run(fmt.Sprintf("N=%d/T=%d/logQ=%d/logP=%d/bitDecomp=%d/SwitchKeysNew", bfvTest.bfvcontext.N(),
			bfvTest.bfvcontext.T(),
			bfvTest.bfvcontext.LogQ(),
			bfvTest.bfvcontext.LogP(),
			bitDecomp), func(t *testing.T) {

			coeffs0, _, ciphertext0, _ := newTestVectors(bfvTest)

			ciphertextTest, err := evaluator.SwitchKeysNew(ciphertext0, switchingKey)
			if err != nil {
				t.Error(err)
			}

			if equalslice(coeffs0.Coeffs[0], encoder.DecodeUint(decryptorSkNew.DecryptNew(ciphertextTest))) != true {
				t.Errorf("error : switchingKeyNew encrypt/decrypt")
			}
		})
	}
}

func testGaloisEnd(bfvTest *BFVTESTPARAMS, bitDecomps []uint64, t *testing.T) {

	bfvContext := bfvTest.bfvcontext
	kgen := bfvTest.kgen
	Sk := bfvTest.sk
	evaluator := bfvTest.evaluator

	coeffs, _, ciphertext, _ := newTestVectors(bfvTest)

	receiverCiphertext := bfvContext.NewCiphertext(1)

	slots := bfvContext.n >> 1
	mask := slots - 1

	coeffsWantRotateCol := bfvContext.contextT.NewPoly()
	coeffsWantRotateRow := bfvContext.contextT.NewPoly()
	coeffsWantRotateRow.Coeffs[0] = append(coeffs.Coeffs[0][bfvContext.n>>1:], coeffs.Coeffs[0][:bfvContext.n>>1]...)

	for _, bitDecomp := range bitDecomps {

		rotationKey := kgen.NewRotationKeysPow2(Sk, bitDecomp, true)

		for n := uint64(1); n < bfvContext.n>>1; n <<= 1 {

			for i := uint64(0); i < slots; i++ {
				coeffsWantRotateCol.Coeffs[0][i] = coeffs.Coeffs[0][(i+n)&mask]
				coeffsWantRotateCol.Coeffs[0][i+slots] = coeffs.Coeffs[0][((i+n)&mask)+slots]
			}

			t.Run(fmt.Sprintf("N=%d/T=%d/logQ=%d/logP=%d/bitDecomp=%d/RotateColumns/%d", bfvTest.bfvcontext.N(),
				bfvTest.bfvcontext.T(),
				bfvTest.bfvcontext.LogQ(),
				bfvTest.bfvcontext.LogP(),
				bitDecomp, n), func(t *testing.T) {

				if err := evaluator.RotateColumns(ciphertext, n, rotationKey, receiverCiphertext); err != nil {
					t.Error(err)
				}

				verifyTestVectors(bfvTest, coeffsWantRotateCol, receiverCiphertext, t)
			})
		}

		for n := uint64(1); n < bfvContext.n>>1; n <<= 2 {

			rand := ring.RandUniform(mask+1, mask)

			for i := uint64(0); i < slots; i++ {
				coeffsWantRotateCol.Coeffs[0][i] = coeffs.Coeffs[0][(i+rand)&mask]
				coeffsWantRotateCol.Coeffs[0][i+slots] = coeffs.Coeffs[0][((i+rand)&mask)+slots]
			}

			t.Run(fmt.Sprintf("N=%d/T=%d/logQ=%d/logP=%d/bitDecomp=%d/RotateColumns/%d", bfvTest.bfvcontext.N(),
				bfvTest.bfvcontext.T(),
				bfvTest.bfvcontext.LogQ(),
				bfvTest.bfvcontext.LogP(),
				bitDecomp, rand), func(t *testing.T) {

				if err := evaluator.RotateColumns(ciphertext, rand, rotationKey, receiverCiphertext); err != nil {
					t.Error(err)
				}

				verifyTestVectors(bfvTest, coeffsWantRotateCol, receiverCiphertext, t)
			})
		}

		t.Run(fmt.Sprintf("N=%d/T=%d/logQ=%d/logP=%d/bitDecomp=%d/RotateRows", bfvTest.bfvcontext.N(),
			bfvTest.bfvcontext.T(),
			bfvTest.bfvcontext.LogQ(),
			bfvTest.bfvcontext.LogP(),
			bitDecomp), func(t *testing.T) {

			if err := evaluator.RotateRows(ciphertext, rotationKey, receiverCiphertext); err != nil {
				t.Error(err)
			}

			verifyTestVectors(bfvTest, coeffsWantRotateRow, receiverCiphertext, t)
		})
	}
}
