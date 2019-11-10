package bfv

import (
	"fmt"
	"github.com/ldsec/lattigo/ring"
	"testing"
)

type BFVTESTPARAMS struct {
	bfvcontext   *BfvContext
	batchencoder *BatchEncoder
	kgen         *KeyGenerator
	sk           *SecretKey
	pk           *PublicKey
	encryptorSk  *Encryptor
	encryptorPk  *Encryptor
	decryptor    *Decryptor
	evaluator    *Evaluator
}

func Test_BFV(t *testing.T) {

	var err error

	paramSets := DefaultParams

	for _, params := range paramSets {

		bfvTest := new(BFVTESTPARAMS)

		bfvTest.bfvcontext = NewBfvContext()
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

		test_EncodeDecode(bfvTest, t)
		test_PlaintextBatchEncodeDecode(bfvTest, t)
		test_EncryptDecrypt(bfvTest, t)
		test_HomomorphicAddition(bfvTest, t)
		test_AddWithPlaintext(bfvTest, t)
		test_HomomorphicSubtraction(bfvTest, t)
		test_SubWithPlaintext(bfvTest, t)

		test_HomomorphicMultiplication(bfvTest, t)

		test_MulWithPlaintext(bfvTest, t)
		test_Relinearization(bfvTest, t)
		test_KeySwitching(bfvTest, t)
		test_GaloisEnd(bfvTest, t)
	}
}

//=============================================================
//=== TESTS FOR ENCODE AND DECODE AN INTEGER ON A PLAINTEXT ===
//=============================================================

func test_EncodeDecode(bfvTest *BFVTESTPARAMS, t *testing.T) {

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

//=====================================
//=== TESTS FOR ENCRYPT AND DECRYPT ===
//=====================================

func test_PlaintextBatchEncodeDecode(bfvTest *BFVTESTPARAMS, t *testing.T) {

	t.Run(fmt.Sprintf("N=%d/T=%d/logQ=%d/logP=%d/BatchEncodeDecode", bfvTest.bfvcontext.N(),
		bfvTest.bfvcontext.T(),
		bfvTest.bfvcontext.LogQ(),
		bfvTest.bfvcontext.LogP()), func(t *testing.T) {

		bfvContext := bfvTest.bfvcontext
		encoder := bfvTest.batchencoder

		coeffsWant := bfvContext.contextT.NewUniformPoly()

		plaintextWant := bfvContext.NewPlaintext()

		encoder.EncodeUint(coeffsWant.Coeffs[0], plaintextWant)

		if EqualSlice(coeffsWant.Coeffs[0], encoder.DecodeUint(plaintextWant)) != true {
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

	if EqualSlice(coeffs.Coeffs[0], coeffsTest) != true {
		t.Errorf("decryption error")
	}
}

func test_EncryptDecrypt(bfvTest *BFVTESTPARAMS, t *testing.T) {

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

func test_HomomorphicAddition(bfvTest *BFVTESTPARAMS, t *testing.T) {

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
func test_AddWithPlaintext(bfvTest *BFVTESTPARAMS, t *testing.T) {

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

func test_HomomorphicSubtraction(bfvTest *BFVTESTPARAMS, t *testing.T) {

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
func test_SubWithPlaintext(bfvTest *BFVTESTPARAMS, t *testing.T) {

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
func test_HomomorphicMultiplication(bfvTest *BFVTESTPARAMS, t *testing.T) {

	bfvContext := bfvTest.bfvcontext
	evaluator := bfvTest.evaluator

	t.Run(fmt.Sprintf("N=%d/T=%d/logQ=%d/logP=%d/Mul", bfvTest.bfvcontext.N(),
		bfvTest.bfvcontext.T(),
		bfvTest.bfvcontext.LogQ(),
		bfvTest.bfvcontext.LogP()), func(t *testing.T) {

		coeffs0, _, ciphertext0, _ := newTestVectors(bfvTest)
		coeffs1, _, ciphertext1, _ := newTestVectors(bfvTest)

		receiverCiphertext1 := bfvContext.NewCiphertextBig(ciphertext0.Degree() + ciphertext1.Degree())
		if err := evaluator.Mul(ciphertext0, ciphertext1, receiverCiphertext1); err != nil {
			t.Error(err)
		}
		bfvContext.contextT.MulCoeffs(coeffs0, coeffs1, coeffs0)
		verifyTestVectors(bfvTest, coeffs0, receiverCiphertext1, t)

		receiverCiphertext2 := bfvContext.NewCiphertextBig(receiverCiphertext1.Degree() + ciphertext1.Degree())
		if err := evaluator.Mul(receiverCiphertext1, ciphertext1, receiverCiphertext2); err != nil {
			t.Error(err)
		}
		bfvContext.contextT.MulCoeffs(coeffs0, coeffs1, coeffs0)
		verifyTestVectors(bfvTest, coeffs0, receiverCiphertext2, t)
	})

	t.Run(fmt.Sprintf("N=%d/T=%d/logQ=%d/logP=%d/MulNew", bfvTest.bfvcontext.N(),
		bfvTest.bfvcontext.T(),
		bfvTest.bfvcontext.LogQ(),
		bfvTest.bfvcontext.LogP()), func(t *testing.T) {

		coeffs0, _, ciphertext0, _ := newTestVectors(bfvTest)
		coeffs1, _, ciphertext1, _ := newTestVectors(bfvTest)

		ciphertext0, _ = evaluator.MulNew(ciphertext0, ciphertext1)
		bfvContext.contextT.MulCoeffs(coeffs0, coeffs1, coeffs0)
		verifyTestVectors(bfvTest, coeffs0, ciphertext0, t)

		ciphertext0, _ = evaluator.MulNew(ciphertext0, ciphertext1)
		bfvContext.contextT.MulCoeffs(coeffs0, coeffs1, coeffs0)
		verifyTestVectors(bfvTest, coeffs0, ciphertext0, t)

	})

	t.Run(fmt.Sprintf("N=%d/T=%d/logQ=%d/logP=%d/Square", bfvTest.bfvcontext.N(),
		bfvTest.bfvcontext.T(),
		bfvTest.bfvcontext.LogQ(),
		bfvTest.bfvcontext.LogP()), func(t *testing.T) {

		coeffs0, _, ciphertext0, _ := newTestVectors(bfvTest)

		ciphertext0, _ = evaluator.MulNew(ciphertext0, ciphertext0)
		bfvContext.contextT.MulCoeffs(coeffs0, coeffs0, coeffs0)
		verifyTestVectors(bfvTest, coeffs0, ciphertext0, t)

		ciphertext0, _ = evaluator.MulNew(ciphertext0, ciphertext0)
		bfvContext.contextT.MulCoeffs(coeffs0, coeffs0, coeffs0)
		verifyTestVectors(bfvTest, coeffs0, ciphertext0, t)

	})

}

//HOMOMORPHIC MULTIPLICATION WITH RELINEARIZATION BASE Qi AND w
func test_MulWithPlaintext(bfvTest *BFVTESTPARAMS, t *testing.T) {

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

func test_Relinearization(bfvTest *BFVTESTPARAMS, t *testing.T) {

	bfvContext := bfvTest.bfvcontext
	kgen := bfvTest.kgen
	evaluator := bfvTest.evaluator

	rlk := kgen.NewRelinKey(bfvTest.sk, 2)

	t.Run(fmt.Sprintf("N=%d/T=%d/logQ=%d/logP=%d/Relinearize", bfvTest.bfvcontext.N(),
		bfvTest.bfvcontext.T(),
		bfvTest.bfvcontext.LogQ(),
		bfvTest.bfvcontext.LogP()), func(t *testing.T) {

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

	t.Run(fmt.Sprintf("N=%d/T=%d/logQ=%d/logP=%d/RelinearizeNew", bfvTest.bfvcontext.N(),
		bfvTest.bfvcontext.T(),
		bfvTest.bfvcontext.LogQ(),
		bfvTest.bfvcontext.LogP()), func(t *testing.T) {

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

func test_KeySwitching(bfvTest *BFVTESTPARAMS, t *testing.T) {

	bfvContext := bfvTest.bfvcontext
	encoder := bfvTest.batchencoder
	kgen := bfvTest.kgen
	Sk := bfvTest.sk
	evaluator := bfvTest.evaluator

	SkNew := kgen.NewSecretKey()

	decryptor_SkNew, err := bfvContext.NewDecryptor(SkNew)
	if err != nil {
		t.Error(err)
	}

	switching_key := kgen.NewSwitchingKey(Sk, SkNew)

	t.Run(fmt.Sprintf("N=%d/T=%d/logQ=%d/logP=%d/SwitchKeys", bfvTest.bfvcontext.N(),
		bfvTest.bfvcontext.T(),
		bfvTest.bfvcontext.LogQ(),
		bfvTest.bfvcontext.LogP()), func(t *testing.T) {

		coeffs0, _, ciphertext0, _ := newTestVectors(bfvTest)

		if err := evaluator.SwitchKeys(ciphertext0, switching_key, ciphertext0); err != nil {
			t.Error(err)
		}

		if EqualSlice(coeffs0.Coeffs[0], encoder.DecodeUint(decryptor_SkNew.DecryptNew(ciphertext0))) != true {
			t.Errorf("error : switchingKey encrypt/decrypt")
		}

	})

	t.Run(fmt.Sprintf("N=%d/T=%d/logQ=%d/logP=%d/SwitchKeysNew", bfvTest.bfvcontext.N(),
		bfvTest.bfvcontext.T(),
		bfvTest.bfvcontext.LogQ(),
		bfvTest.bfvcontext.LogP()), func(t *testing.T) {

		coeffs0, _, ciphertext0, _ := newTestVectors(bfvTest)

		ciphertextTest, err := evaluator.SwitchKeysNew(ciphertext0, switching_key)
		if err != nil {
			t.Error(err)
		}

		if EqualSlice(coeffs0.Coeffs[0], encoder.DecodeUint(decryptor_SkNew.DecryptNew(ciphertextTest))) != true {
			t.Errorf("error : switchingKeyNew encrypt/decrypt")
		}
	})

}

func test_GaloisEnd(bfvTest *BFVTESTPARAMS, t *testing.T) {

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

	rotation_key := kgen.NewRotationKeysPow2(Sk, true)

	for n := uint64(1); n < bfvContext.n>>1; n <<= 1 {

		for i := uint64(0); i < slots; i++ {
			coeffsWantRotateCol.Coeffs[0][i] = coeffs.Coeffs[0][(i+n)&mask]
			coeffsWantRotateCol.Coeffs[0][i+slots] = coeffs.Coeffs[0][((i+n)&mask)+slots]
		}

		t.Run(fmt.Sprintf("N=%d/T=%d/logQ=%d/logP=%d/RotateColumns/%d", bfvTest.bfvcontext.N(),
			bfvTest.bfvcontext.T(),
			bfvTest.bfvcontext.LogQ(),
			bfvTest.bfvcontext.LogP(), n), func(t *testing.T) {

			if err := evaluator.RotateColumns(ciphertext, n, rotation_key, receiverCiphertext); err != nil {
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

		t.Run(fmt.Sprintf("N=%d/T=%d/logQ=%d/logP=%d/RotateColumns/%d", bfvTest.bfvcontext.N(),
			bfvTest.bfvcontext.T(),
			bfvTest.bfvcontext.LogQ(),
			bfvTest.bfvcontext.LogP(), rand), func(t *testing.T) {

			if err := evaluator.RotateColumns(ciphertext, rand, rotation_key, receiverCiphertext); err != nil {
				t.Error(err)
			}

			verifyTestVectors(bfvTest, coeffsWantRotateCol, receiverCiphertext, t)
		})
	}

	t.Run(fmt.Sprintf("N=%d/T=%d/logQ=%d/logP=%d/RotateRows", bfvTest.bfvcontext.N(),
		bfvTest.bfvcontext.T(),
		bfvTest.bfvcontext.LogQ(),
		bfvTest.bfvcontext.LogP()), func(t *testing.T) {

		if err := evaluator.RotateRows(ciphertext, rotation_key, receiverCiphertext); err != nil {
			t.Error(err)
		}

		verifyTestVectors(bfvTest, coeffsWantRotateRow, receiverCiphertext, t)
	})

}
