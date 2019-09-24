package bfv

import (
	"fmt"
	"github.com/ldsec/lattigo/ring"
	"log"
	"testing"
)

type BFVTESTPARAMS struct {
	bfvcontext   *BfvContext
	batchencoder *BatchEncoder
	kgen         *KeyGenerator
	sk           *SecretKey
	pk           *PublicKey
	encryptor    *Encryptor
	decryptor    *Decryptor
	evaluator    *Evaluator
}

func Test_BFV(t *testing.T) {

	var err error

	paramSets := DefaultParams[0:2]

	bitDecomps := []uint64{60}

	for _, params := range paramSets {

		log.Printf("Generating testing bfvcontext for N=%d/T=%d/Qi=%dx%dbit/Pi=%dx%dbit", params.N,
			params.T,
			len(params.Qi), 60,
			len(params.Pi), 60)

		bfvTest := new(BFVTESTPARAMS)

		bfvTest.bfvcontext = NewBfvContext()
		if err := bfvTest.bfvcontext.SetParameters(params.N, params.T, params.Qi, params.Pi, params.Sigma); err != nil {
			log.Fatal(err)
		}

		bfvTest.kgen = bfvTest.bfvcontext.NewKeyGenerator()

		bfvTest.batchencoder = bfvTest.bfvcontext.NewBatchEncoder()

		bfvTest.sk, bfvTest.pk, err = bfvTest.kgen.NewKeyPair()

		if err != nil {
			log.Fatal(err)
		}

		bfvTest.decryptor, err = bfvTest.bfvcontext.NewDecryptor(bfvTest.sk)

		if err != nil {
			log.Fatal(err)
		}

		bfvTest.encryptor, err = bfvTest.bfvcontext.NewEncryptor(bfvTest.pk)

		if err != nil {
			log.Fatal(err)
		}

		bfvTest.evaluator, err = bfvTest.bfvcontext.NewEvaluator()

		if err != nil {
			log.Fatal(err)
		}

		test_Marshaler(bfvTest, t)
		test_EncodeDecode(bfvTest, t)
		test_PlaintextBatchEncodeDecode(bfvTest, t)
		test_EncryptDecrypt(bfvTest, t)
		test_HomomorphicAddition(bfvTest, t)
		test_AddWithPlaintext(bfvTest, t)
		test_HomomorphicSubtraction(bfvTest, t)
		test_SubWithPlaintext(bfvTest, t)
		test_HomomorphicMultiplication(bfvTest, t)
		test_MulWithPlaintext(bfvTest, t)
		test_Relinearization(bfvTest, bitDecomps, t)
		test_KeySwitching(bfvTest, bitDecomps, t)
		test_GaloisEnd(bfvTest, bitDecomps, t)

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
			base = int64(ring.RandUniform(256) % (T >> 1))
			if base < 2 {
				base = 2
			}

			encoder := bfvcontext.NewIntEncoder(base)

			//Generates a value which can fit in the ring given the base
			value = int64(ring.RandUniform(0xFFFFFFFF))

			if T != 2 {
				sign = uint32(ring.RandUniform(2) & 1)
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

func test_Marshaler(bfvTest *BFVTESTPARAMS, t *testing.T) {

	state := true

	bfvContext := bfvTest.bfvcontext
	encoder := bfvTest.batchencoder
	Sk := bfvTest.sk
	Pk := bfvTest.pk

	t.Run(fmt.Sprintf("N=%d/T=%d/Qi=%dx%dbit/Pi=%dx%dbit/MarshalBfvContext", bfvTest.bfvcontext.n,
		bfvTest.bfvcontext.t,
		len(bfvTest.bfvcontext.contextQ.Modulus), 60,
		len(bfvTest.bfvcontext.contextP.Modulus), 60), func(t *testing.T) {

		BfvBytes, err := bfvContext.MarshalBinary()
		if err != nil {
			log.Fatal(err)
		}

		BfvTest := NewBfvContext()
		BfvTest.UnMarshalBinary(BfvBytes)

		for i := 0; i < 1; i++ {

			if BfvTest.n != bfvContext.n {
				state = false
				break
			}

			if BfvTest.t != bfvContext.t {
				state = false
				break
			}

			if BfvTest.sigma != bfvContext.sigma {
				state = false
				break
			}

			for i := range BfvTest.contextQ.Modulus {
				if BfvTest.contextQ.Modulus[i] != bfvContext.contextQ.Modulus[i] {
					state = false
					break
				}
			}

			for i := range BfvTest.contextP.Modulus {
				if BfvTest.contextP.Modulus[i] != bfvContext.contextP.Modulus[i] {
					state = false
					break
				}
			}
		}

		if state != true {
			t.Errorf("error : binarymarshal bfvcontext")
		}

	})

	t.Run(fmt.Sprintf("N=%d/T=%d/Qi=%dx%dbit/MarshalSk", bfvTest.bfvcontext.n,
		bfvTest.bfvcontext.t,
		len(bfvTest.bfvcontext.contextQ.Modulus), 60), func(t *testing.T) {
		state = true

		SkBytes, err := Sk.MarshalBinary()
		if err != nil {

		}
		SkTest := new(SecretKey)
		SkTest.UnMarshalBinary(SkBytes)

		if bfvContext.contextQ.Equal(Sk.sk, SkTest.sk) != true {
			t.Errorf("error : binarymarshal secretkey")
		}

	})

	t.Run(fmt.Sprintf("N=%d/T=%d/Qi=%dx%dbit/MarshalPk", bfvTest.bfvcontext.n,
		bfvTest.bfvcontext.t,
		len(bfvTest.bfvcontext.contextQ.Modulus), 60), func(t *testing.T) {
		PkBytes, err := Pk.MarshalBinary()
		if err != nil {
			log.Fatal(err)
		}

		PkTest := new(PublicKey)
		PkTest.UnMarshalBinary(PkBytes)

		for i := range Pk.pk {
			if bfvContext.contextQ.Equal(Pk.pk[i], PkTest.pk[i]) != true {
				t.Errorf("error : binarymarshal publickey")
				break
			}
		}
	})

	t.Run(fmt.Sprintf("N=%d/T=%d/Qi=%dx%dbit/MarshalPlaintext", bfvTest.bfvcontext.n,
		bfvTest.bfvcontext.t,
		len(bfvTest.bfvcontext.contextQ.Modulus), 60), func(t *testing.T) {
		coeffs := bfvContext.NewRandomPlaintextCoeffs()
		Px := bfvContext.NewPlaintext()
		encoder.EncodeUint(coeffs, Px)
		PxBytes, err := Px.MarshalBinary()
		if err != nil {
			log.Fatal(err)
		}

		PxTest := new(Plaintext)
		PxTest.UnMarshalBinary(PxBytes)

		if bfvContext.contextT.Equal(PxTest.Value()[0], Px.Value()[0]) != true {
			t.Errorf("error : binarymarshal plaintext")
		}
	})

	t.Run(fmt.Sprintf("N=%d/T=%d/Qi=%dx%dbit/MarshalCiphertext", bfvTest.bfvcontext.n,
		bfvTest.bfvcontext.t,
		len(bfvTest.bfvcontext.contextQ.Modulus), 60), func(t *testing.T) {

		Ctx := bfvContext.NewRandomCiphertext(4)
		CtxBytes, err := Ctx.MarshalBinary()
		if err != nil {
			log.Fatal(err)
		}

		CtxTest := new(Ciphertext)
		if err = CtxTest.UnmarshalBinary(CtxBytes); err != nil {
			log.Fatal(err)
		}

		for i := range Ctx.Value() {
			if bfvContext.contextQ.Equal(CtxTest.Value()[i], Ctx.Value()[i]) != true {
				t.Errorf("error : binarymarshal ciphertext")
				break
			}
		}

	})

	t.Run(fmt.Sprintf("N=%d/T=%d/Qi=%dx%dbit/bitDecomp=%d/Marshalrlk", bfvTest.bfvcontext.n,
		bfvTest.bfvcontext.t,
		len(bfvTest.bfvcontext.contextQ.Modulus), 60,
		15), func(t *testing.T) {
		rlk, err := bfvTest.kgen.NewRelinKey(bfvTest.sk, 5, 15)
		if err != nil {
			log.Fatal(err)
		}

		rlkBytes, err := rlk.MarshalBinary()
		if err != nil {
			log.Fatal(err)
		}

		rlkTest := new(EvaluationKey)
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

	t.Run(fmt.Sprintf("N=%d/T=%d/Qi=%dx%dbit/bitDecomp=%d/MarshalRotKey", bfvTest.bfvcontext.n,
		bfvTest.bfvcontext.t,
		len(bfvTest.bfvcontext.contextQ.Modulus), 60,
		15), func(t *testing.T) {
		rotKey, err := bfvTest.kgen.NewRotationKeysPow2(Sk, 15, true)
		if err != nil {
			log.Fatal(err)
		}

		rotKeyBytes, err := rotKey.MarshalBinary()

		rotKeyTest := new(RotationKeys)

		rotKeyTest.UnMarshalBinary(rotKeyBytes)

		state = true

		if rotKeyTest.evakey_rot_row != nil {
			for j := range rotKeyTest.evakey_rot_row.evakey {
				for x := range rotKeyTest.evakey_rot_row.evakey[j] {
					if bfvContext.contextQ.Equal(rotKeyTest.evakey_rot_row.evakey[j][x][0], rotKey.evakey_rot_row.evakey[j][x][0]) != true {
						t.Errorf("error : binarymarshal rotation key rot row")
						state = false
						break
					}

					if bfvContext.contextQ.Equal(rotKeyTest.evakey_rot_row.evakey[j][x][1], rotKey.evakey_rot_row.evakey[j][x][1]) != true {
						t.Errorf("error : binarymarshal rotation key rot row")
						state = false
						break
					}

				}
			}
		}

		for k, v := range rotKeyTest.evakey_rot_col_L {
			for j := range v.evakey {
				for x := range v.evakey[j] {
					if bfvContext.contextQ.Equal(v.evakey[j][x][0], rotKey.evakey_rot_col_L[k].evakey[j][x][0]) != true {
						t.Errorf("error : binarymarshal rotation key col L")
						state = false
						break
					}

					if bfvContext.contextQ.Equal(v.evakey[j][x][1], rotKey.evakey_rot_col_L[k].evakey[j][x][1]) != true {
						t.Errorf("error : binarymarshal rotation key col L")
						state = false
						break
					}
				}
			}

		}

		for k, v := range rotKeyTest.evakey_rot_col_R {
			for j := range v.evakey {
				for x := range v.evakey[j] {
					if bfvContext.contextQ.Equal(v.evakey[j][x][0], rotKey.evakey_rot_col_R[k].evakey[j][x][0]) != true {
						t.Errorf("error : binarymarshal rotation key col R")
						state = false
						break
					}

					if bfvContext.contextQ.Equal(v.evakey[j][x][1], rotKey.evakey_rot_col_R[k].evakey[j][x][1]) != true {
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

func test_PlaintextBatchEncodeDecode(bfvTest *BFVTESTPARAMS, t *testing.T) {

	t.Run(fmt.Sprintf("N=%d/T=%d/Qi=%dx%dbit/Pi=%dx%dbit/BatchEncodeDecode", bfvTest.bfvcontext.n,
		bfvTest.bfvcontext.t,
		len(bfvTest.bfvcontext.contextQ.Modulus), 60,
		len(bfvTest.bfvcontext.contextP.Modulus), 60), func(t *testing.T) {

		bfvContext := bfvTest.bfvcontext
		encoder := bfvTest.batchencoder

		coeffsWant := bfvContext.contextT.NewUniformPoly()

		plaintextWant := bfvContext.NewPlaintext()

		encoder.EncodeUint(coeffsWant.Coeffs[0], plaintextWant)

		coeffsTest, err := encoder.DecodeUint(plaintextWant)
		if err != nil {
			log.Fatal(err)
		}

		if EqualSlice(coeffsWant.Coeffs[0], coeffsTest) != true {
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

	if ciphertext, err = bfvTest.encryptor.EncryptNew(plaintext); err != nil {
		return nil, nil, nil, err
	}

	return coeffs, plaintext, ciphertext, nil
}

func verifyTestVectors(bfvTest *BFVTESTPARAMS, coeffs *ring.Poly, element BfvElement, t *testing.T) {

	var coeffsTest []uint64
	var err error

	switch element.(type) {
	case *Ciphertext:

		var plaintext *Plaintext
		if plaintext, err = bfvTest.decryptor.DecryptNew(element.(*Ciphertext)); err != nil {
			log.Fatal(err)
		}

		if coeffsTest, err = bfvTest.batchencoder.DecodeUint(plaintext); err != nil {
			log.Fatal(err)
		}

	case *Plaintext:

		if coeffsTest, err = bfvTest.batchencoder.DecodeUint(element.(*Plaintext)); err != nil {
			log.Fatal(err)
		}
	}

	if EqualSlice(coeffs.Coeffs[0], coeffsTest) != true {
		t.Errorf("decryption error")
	}
}

func test_EncryptDecrypt(bfvTest *BFVTESTPARAMS, t *testing.T) {

	t.Run(fmt.Sprintf("N=%d/T=%d/Qi=%dx%dbit/Pi=%dx%dbit/Encrypt/Decrypt", bfvTest.bfvcontext.n,
		bfvTest.bfvcontext.t,
		len(bfvTest.bfvcontext.contextQ.Modulus), 60,
		len(bfvTest.bfvcontext.contextP.Modulus), 60), func(t *testing.T) {

		coeffs, _, ciphertext, _ := newTestVectors(bfvTest)
		verifyTestVectors(bfvTest, coeffs, ciphertext, t)
	})
}

func test_HomomorphicAddition(bfvTest *BFVTESTPARAMS, t *testing.T) {

	bfvContext := bfvTest.bfvcontext
	evaluator := bfvTest.evaluator

	var ciphertextTest BfvElement
	var err error

	t.Run(fmt.Sprintf("N=%d/T=%d/Qi=%dx%dbit/Pi=%dx%dbit/Add", bfvTest.bfvcontext.n,
		bfvTest.bfvcontext.t,
		len(bfvTest.bfvcontext.contextQ.Modulus), 60,
		len(bfvTest.bfvcontext.contextP.Modulus), 60), func(t *testing.T) {

		coeffs0, _, ciphertext0, _ := newTestVectors(bfvTest)
		coeffs1, _, ciphertext1, _ := newTestVectors(bfvTest)

		for i := 0; i < 1; i++ {
			if err = evaluator.Add(ciphertext0, ciphertext1, ciphertext0); err != nil {
				log.Fatal(err)
			}
			bfvContext.contextT.Add(coeffs0, coeffs1, coeffs0)
		}

		verifyTestVectors(bfvTest, coeffs0, ciphertext0, t)
	})

	t.Run(fmt.Sprintf("N=%d/T=%d/Qi=%dx%dbit/Pi=%dx%dbit/AddNew", bfvTest.bfvcontext.n,
		bfvTest.bfvcontext.t,
		len(bfvTest.bfvcontext.contextQ.Modulus), 60,
		len(bfvTest.bfvcontext.contextP.Modulus), 60), func(t *testing.T) {

		coeffs0, _, ciphertext0, _ := newTestVectors(bfvTest)
		coeffs1, _, ciphertext1, _ := newTestVectors(bfvTest)

		for i := 0; i < 1; i++ {
			if ciphertextTest, err = evaluator.AddNew(ciphertext0, ciphertext1); err != nil {
				log.Fatal(err)
			}

			bfvContext.contextT.Add(coeffs0, coeffs1, coeffs0)
		}

		verifyTestVectors(bfvTest, coeffs0, ciphertextTest, t)

	})

	t.Run(fmt.Sprintf("N=%d/T=%d/Qi=%dx%dbit/Pi=%dx%dbit/AddNoMod", bfvTest.bfvcontext.n,
		bfvTest.bfvcontext.t,
		len(bfvTest.bfvcontext.contextQ.Modulus), 60,
		len(bfvTest.bfvcontext.contextP.Modulus), 60), func(t *testing.T) {

		coeffs0, _, ciphertext0, _ := newTestVectors(bfvTest)
		coeffs1, _, ciphertext1, _ := newTestVectors(bfvTest)

		for i := 0; i < 1; i++ {

			if err = evaluator.AddNoMod(ciphertext0, ciphertext1, ciphertext0); err != nil {
				log.Fatal(err)
			}

			evaluator.Reduce(ciphertext0, ciphertext0)
			bfvContext.contextT.Add(coeffs0, coeffs1, coeffs0)

		}

		verifyTestVectors(bfvTest, coeffs0, ciphertext0, t)
	})

	t.Run(fmt.Sprintf("N=%d/T=%d/Qi=%dx%dbit/Pi=%dx%dbit/AddNoModNew", bfvTest.bfvcontext.n,
		bfvTest.bfvcontext.t,
		len(bfvTest.bfvcontext.contextQ.Modulus), 60,
		len(bfvTest.bfvcontext.contextP.Modulus), 60), func(t *testing.T) {

		coeffs0, _, ciphertext0, _ := newTestVectors(bfvTest)
		coeffs1, _, ciphertext1, _ := newTestVectors(bfvTest)

		for i := 0; i < 1; i++ {
			if ciphertextTest, err = evaluator.AddNoModNew(ciphertext0, ciphertext1); err != nil {
				log.Fatal(err)
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

	var ciphertextTest BfvElement
	var err error

	t.Run(fmt.Sprintf("N=%d/T=%d/Qi=%dx%dbit/Pi=%dx%dbit/AddWithPlaintext", bfvTest.bfvcontext.n,
		bfvTest.bfvcontext.t,
		len(bfvTest.bfvcontext.contextQ.Modulus), 60,
		len(bfvTest.bfvcontext.contextP.Modulus), 60), func(t *testing.T) {

		coeffs0, _, ciphertext0, _ := newTestVectors(bfvTest)
		coeffs1, plaintext1, _, _ := newTestVectors(bfvTest)

		if err := evaluator.Add(ciphertext0, plaintext1, ciphertext0); err != nil {
			log.Fatal(err)
		}

		bfvContext.contextT.Add(coeffs0, coeffs1, coeffs0)

		verifyTestVectors(bfvTest, coeffs0, ciphertext0, t)

	})

	t.Run(fmt.Sprintf("N=%d/T=%d/Qi=%dx%dbit/Pi=%dx%dbit/AddWithPlaintextNew", bfvTest.bfvcontext.n,
		bfvTest.bfvcontext.t,
		len(bfvTest.bfvcontext.contextQ.Modulus), 60,
		len(bfvTest.bfvcontext.contextP.Modulus), 60), func(t *testing.T) {

		coeffs0, _, ciphertext0, _ := newTestVectors(bfvTest)
		coeffs1, plaintext1, _, _ := newTestVectors(bfvTest)

		if ciphertextTest, err = evaluator.AddNew(ciphertext0, plaintext1); err != nil {
			log.Fatal(err)
		}
		bfvContext.contextT.Add(coeffs0, coeffs1, coeffs0)

		verifyTestVectors(bfvTest, coeffs0, ciphertextTest, t)

	})

	t.Run(fmt.Sprintf("N=%d/T=%d/Qi=%dx%dbit/Pi=%dx%dbit/AddWithPlaintextNoMod", bfvTest.bfvcontext.n,
		bfvTest.bfvcontext.t,
		len(bfvTest.bfvcontext.contextQ.Modulus), 60,
		len(bfvTest.bfvcontext.contextP.Modulus), 60), func(t *testing.T) {

		coeffs0, _, ciphertext0, _ := newTestVectors(bfvTest)
		coeffs1, plaintext1, _, _ := newTestVectors(bfvTest)

		if err := evaluator.AddNoMod(ciphertext0, plaintext1, ciphertext0); err != nil {
			log.Fatal(err)
		}
		evaluator.Reduce(ciphertext0, ciphertext0)
		bfvContext.contextT.Add(coeffs0, coeffs1, coeffs0)

		verifyTestVectors(bfvTest, coeffs0, ciphertext0, t)
	})

	t.Run(fmt.Sprintf("N=%d/T=%d/Qi=%dx%dbit/Pi=%dx%dbit/AddWithPlaintextNoModNew", bfvTest.bfvcontext.n,
		bfvTest.bfvcontext.t,
		len(bfvTest.bfvcontext.contextQ.Modulus), 60,
		len(bfvTest.bfvcontext.contextP.Modulus), 60), func(t *testing.T) {

		coeffs0, _, ciphertext0, _ := newTestVectors(bfvTest)
		coeffs1, plaintext1, _, _ := newTestVectors(bfvTest)

		if ciphertextTest, err = evaluator.AddNoModNew(ciphertext0, plaintext1); err != nil {
			log.Fatal(err)
		}
		evaluator.Reduce(ciphertextTest, ciphertextTest)
		bfvContext.contextT.Add(coeffs0, coeffs1, coeffs0)

		verifyTestVectors(bfvTest, coeffs0, ciphertextTest, t)
	})
}

func test_HomomorphicSubtraction(bfvTest *BFVTESTPARAMS, t *testing.T) {

	bfvContext := bfvTest.bfvcontext
	evaluator := bfvTest.evaluator

	var ciphertextTest BfvElement
	var err error

	t.Run(fmt.Sprintf("N=%d/T=%d/Qi=%dx%dbit/Pi=%dx%dbit/Sub", bfvTest.bfvcontext.n,
		bfvTest.bfvcontext.t,
		len(bfvTest.bfvcontext.contextQ.Modulus), 60,
		len(bfvTest.bfvcontext.contextP.Modulus), 60), func(t *testing.T) {

		coeffs0, _, ciphertext0, _ := newTestVectors(bfvTest)
		coeffs1, _, ciphertext1, _ := newTestVectors(bfvTest)

		for i := 0; i < 1; i++ {
			if err := evaluator.Sub(ciphertext0, ciphertext1, ciphertext0); err != nil {
				log.Fatal(err)
			}
			bfvContext.contextT.Sub(coeffs0, coeffs1, coeffs0)
		}

		verifyTestVectors(bfvTest, coeffs0, ciphertext0, t)

	})

	t.Run(fmt.Sprintf("N=%d/T=%d/Qi=%dx%dbit/Pi=%dx%dbit/SubNew", bfvTest.bfvcontext.n,
		bfvTest.bfvcontext.t,
		len(bfvTest.bfvcontext.contextQ.Modulus), 60,
		len(bfvTest.bfvcontext.contextP.Modulus), 60), func(t *testing.T) {

		coeffs0, _, ciphertext0, _ := newTestVectors(bfvTest)
		coeffs1, _, ciphertext1, _ := newTestVectors(bfvTest)

		for i := 0; i < 1; i++ {
			if ciphertextTest, err = evaluator.SubNew(ciphertext0, ciphertext1); err != nil {
				log.Fatal(err)
			}
			bfvContext.contextT.Sub(coeffs0, coeffs1, coeffs0)
		}

		verifyTestVectors(bfvTest, coeffs0, ciphertextTest, t)
	})

	t.Run(fmt.Sprintf("N=%d/T=%d/Qi=%dx%dbit/Pi=%dx%dbit/SubNoMod", bfvTest.bfvcontext.n,
		bfvTest.bfvcontext.t,
		len(bfvTest.bfvcontext.contextQ.Modulus), 60,
		len(bfvTest.bfvcontext.contextP.Modulus), 60), func(t *testing.T) {

		coeffs0, _, ciphertext0, _ := newTestVectors(bfvTest)
		coeffs1, _, ciphertext1, _ := newTestVectors(bfvTest)

		for i := 0; i < 1; i++ {
			if err := evaluator.SubNoMod(ciphertext0, ciphertext1, ciphertext0); err != nil {
				log.Fatal(err)
			}
			evaluator.Reduce(ciphertext0, ciphertext0)
			bfvContext.contextT.Sub(coeffs0, coeffs1, coeffs0)
		}

		verifyTestVectors(bfvTest, coeffs0, ciphertext0, t)

	})

	t.Run(fmt.Sprintf("N=%d/T=%d/Qi=%dx%dbit/Pi=%dx%dbit/SubNoModNew", bfvTest.bfvcontext.n,
		bfvTest.bfvcontext.t,
		len(bfvTest.bfvcontext.contextQ.Modulus), 60,
		len(bfvTest.bfvcontext.contextP.Modulus), 60), func(t *testing.T) {

		coeffs0, _, ciphertext0, _ := newTestVectors(bfvTest)
		coeffs1, _, ciphertext1, _ := newTestVectors(bfvTest)

		for i := 0; i < 1; i++ {
			if ciphertextTest, err = evaluator.SubNoModNew(ciphertext0, ciphertext1); err != nil {
				log.Fatal(err)
			}
			if err != nil {
				log.Fatal(err)
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

	var ciphertextTest BfvElement
	var err error

	t.Run(fmt.Sprintf("N=%d/T=%d/Qi=%dx%dbit/Pi=%dx%dbit/SubWithPlaintext", bfvTest.bfvcontext.n,
		bfvTest.bfvcontext.t,
		len(bfvTest.bfvcontext.contextQ.Modulus), 60,
		len(bfvTest.bfvcontext.contextP.Modulus), 60), func(t *testing.T) {

		coeffs0, _, ciphertext0, _ := newTestVectors(bfvTest)
		coeffs1, plaintext1, _, _ := newTestVectors(bfvTest)

		if err := evaluator.Sub(ciphertext0, plaintext1, ciphertext0); err != nil {
			log.Fatal(err)
		}
		bfvContext.contextT.Sub(coeffs0, coeffs1, coeffs0)

		verifyTestVectors(bfvTest, coeffs0, ciphertext0, t)
	})

	t.Run(fmt.Sprintf("N=%d/T=%d/Qi=%dx%dbit/Pi=%dx%dbit/SubWithPlaintextNew", bfvTest.bfvcontext.n,
		bfvTest.bfvcontext.t,
		len(bfvTest.bfvcontext.contextQ.Modulus), 60,
		len(bfvTest.bfvcontext.contextP.Modulus), 60), func(t *testing.T) {

		coeffs0, _, ciphertext0, _ := newTestVectors(bfvTest)
		coeffs1, plaintext1, _, _ := newTestVectors(bfvTest)

		if ciphertextTest, err = evaluator.SubNew(ciphertext0, plaintext1); err != nil {
			log.Fatal(err)
		}
		bfvContext.contextT.Sub(coeffs0, coeffs1, coeffs0)

		verifyTestVectors(bfvTest, coeffs0, ciphertextTest, t)

	})

	t.Run(fmt.Sprintf("N=%d/T=%d/Qi=%dx%dbit/Pi=%dx%dbit/SubWithPlaintextNoMod", bfvTest.bfvcontext.n,
		bfvTest.bfvcontext.t,
		len(bfvTest.bfvcontext.contextQ.Modulus), 60,
		len(bfvTest.bfvcontext.contextP.Modulus), 60), func(t *testing.T) {

		coeffs0, _, ciphertext0, _ := newTestVectors(bfvTest)
		coeffs1, plaintext1, _, _ := newTestVectors(bfvTest)

		if err := evaluator.SubNoMod(ciphertext0, plaintext1, ciphertext0); err != nil {
			log.Fatal(err)
		}
		evaluator.Reduce(ciphertext0, ciphertext0)
		bfvContext.contextT.Sub(coeffs0, coeffs1, coeffs0)

		verifyTestVectors(bfvTest, coeffs0, ciphertext0, t)
	})

	t.Run(fmt.Sprintf("N=%d/T=%d/Qi=%dx%dbit/Pi=%dx%dbit/SubWithPlaintextNoModNew", bfvTest.bfvcontext.n,
		bfvTest.bfvcontext.t,
		len(bfvTest.bfvcontext.contextQ.Modulus), 60,
		len(bfvTest.bfvcontext.contextP.Modulus), 60), func(t *testing.T) {

		coeffs0, _, ciphertext0, _ := newTestVectors(bfvTest)
		coeffs1, plaintext1, _, _ := newTestVectors(bfvTest)

		if ciphertextTest, err = evaluator.SubNoModNew(ciphertext0, plaintext1); err != nil {
			log.Fatal(err)
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

	t.Run(fmt.Sprintf("N=%d/T=%d/Qi=%dx%dbit/Pi=%dx%dbit/Mul", bfvTest.bfvcontext.n,
		bfvTest.bfvcontext.t,
		len(bfvTest.bfvcontext.contextQ.Modulus), 60,
		len(bfvTest.bfvcontext.contextP.Modulus), 60), func(t *testing.T) {

		coeffs0, _, ciphertext0, _ := newTestVectors(bfvTest)
		coeffs1, _, ciphertext1, _ := newTestVectors(bfvTest)

		receiverCiphertext := bfvContext.NewCiphertextBig(ciphertext0.Degree() + ciphertext1.Degree())

		if err := evaluator.Mul(ciphertext0, ciphertext1, receiverCiphertext); err != nil {
			log.Fatal(err)
		}

		bfvContext.contextT.MulCoeffs(coeffs0, coeffs1, coeffs0)

		verifyTestVectors(bfvTest, coeffs0, receiverCiphertext, t)
	})

	t.Run(fmt.Sprintf("N=%d/T=%d/Qi=%dx%dbit/Pi=%dx%dbit/MulNew", bfvTest.bfvcontext.n,
		bfvTest.bfvcontext.t,
		len(bfvTest.bfvcontext.contextQ.Modulus), 60,
		len(bfvTest.bfvcontext.contextP.Modulus), 60), func(t *testing.T) {

		coeffs0, _, ciphertext0, _ := newTestVectors(bfvTest)
		coeffs1, _, ciphertext1, _ := newTestVectors(bfvTest)

		receiverCiphertext := evaluator.MulNew(ciphertext0, ciphertext1)

		bfvContext.contextT.MulCoeffs(coeffs0, coeffs1, coeffs0)

		verifyTestVectors(bfvTest, coeffs0, receiverCiphertext, t)

	})
}

//HOMOMORPHIC MULTIPLICATION WITH RELINEARIZATION BASE Qi AND w
func test_MulWithPlaintext(bfvTest *BFVTESTPARAMS, t *testing.T) {

	bfvContext := bfvTest.bfvcontext
	evaluator := bfvTest.evaluator

	t.Run(fmt.Sprintf("N=%d/T=%d/Qi=%dx%dbit/Pi=%dx%dbit/MulWithPlaintext", bfvTest.bfvcontext.n,
		bfvTest.bfvcontext.t,
		len(bfvTest.bfvcontext.contextQ.Modulus), 60,
		len(bfvTest.bfvcontext.contextP.Modulus), 60), func(t *testing.T) {

		coeffs0, _, ciphertext0, _ := newTestVectors(bfvTest)
		coeffs1, plaintext1, _, _ := newTestVectors(bfvTest)

		receiverCiphertext := bfvContext.NewCiphertextBig(ciphertext0.Degree())

		if err := evaluator.Mul(ciphertext0, plaintext1, receiverCiphertext); err != nil {
			log.Fatal(err)
		}

		bfvContext.contextT.MulCoeffs(coeffs0, coeffs1, coeffs0)

		verifyTestVectors(bfvTest, coeffs0, receiverCiphertext, t)
	})

	t.Run(fmt.Sprintf("N=%d/T=%d/Qi=%dx%dbit/Pi=%dx%dbit/MulWithPlaintextNew", bfvTest.bfvcontext.n,
		bfvTest.bfvcontext.t,
		len(bfvTest.bfvcontext.contextQ.Modulus), 60,
		len(bfvTest.bfvcontext.contextP.Modulus), 60), func(t *testing.T) {

		coeffs0, _, ciphertext0, _ := newTestVectors(bfvTest)
		coeffs1, plaintext1, _, _ := newTestVectors(bfvTest)

		receiverCiphertext := evaluator.MulNew(ciphertext0, plaintext1)
		bfvContext.contextT.MulCoeffs(coeffs0, coeffs1, coeffs0)

		verifyTestVectors(bfvTest, coeffs0, receiverCiphertext, t)

	})
}

func test_Relinearization(bfvTest *BFVTESTPARAMS, bitDecomps []uint64, t *testing.T) {

	bfvContext := bfvTest.bfvcontext
	kgen := bfvTest.kgen
	evaluator := bfvTest.evaluator

	for _, bitDecomp := range bitDecomps {

		rlk, err := kgen.NewRelinKey(bfvTest.sk, 1, bitDecomp)
		if err != nil {
			log.Fatal(err)
		}

		t.Run(fmt.Sprintf("N=%d/T=%d/Qi=%dx%dbit/Pi=%dx%dbit/bitDecomp=%d/Relinearize", bfvTest.bfvcontext.n,
			bfvTest.bfvcontext.t,
			len(bfvTest.bfvcontext.contextQ.Modulus), 60,
			len(bfvTest.bfvcontext.contextP.Modulus), 60,
			bitDecomp), func(t *testing.T) {

			coeffs0, _, ciphertext0, _ := newTestVectors(bfvTest)
			coeffs1, _, ciphertext1, _ := newTestVectors(bfvTest)

			receiverCiphertext := evaluator.MulNew(ciphertext0, ciphertext1)
			bfvContext.contextT.MulCoeffs(coeffs0, coeffs1, coeffs0)

			if err := evaluator.Relinearize(receiverCiphertext.(*Ciphertext), rlk, receiverCiphertext.(*Ciphertext)); err != nil {
				log.Fatal(err)
			}

			if len(receiverCiphertext.Value()) > 2 {
				t.Errorf("error : Relinearize")
			}

			verifyTestVectors(bfvTest, coeffs0, receiverCiphertext, t)

		})

		t.Run(fmt.Sprintf("N=%d/T=%d/Qi=%dx%dbit/Pi=%dx%dbit/bitDecomp=%d/RelinearizeNew", bfvTest.bfvcontext.n,
			bfvTest.bfvcontext.t,
			len(bfvTest.bfvcontext.contextQ.Modulus), 60,
			len(bfvTest.bfvcontext.contextP.Modulus), 60,
			bitDecomp), func(t *testing.T) {

			coeffs0, _, ciphertext0, _ := newTestVectors(bfvTest)
			coeffs1, _, ciphertext1, _ := newTestVectors(bfvTest)

			receiverCiphertext := evaluator.MulNew(ciphertext0, ciphertext1)
			bfvContext.contextT.MulCoeffs(coeffs0, coeffs1, coeffs0)

			receiverCiphertext, err := evaluator.RelinearizeNew(receiverCiphertext.(*Ciphertext), rlk)
			if err != nil {
				log.Fatal(err)
			}

			if len(receiverCiphertext.Value()) > 2 {
				t.Errorf("error : RelinearizeNew")
			}

			verifyTestVectors(bfvTest, coeffs0, receiverCiphertext, t)

		})
	}
}

func test_KeySwitching(bfvTest *BFVTESTPARAMS, bitDecomps []uint64, t *testing.T) {

	bfvContext := bfvTest.bfvcontext
	encoder := bfvTest.batchencoder
	kgen := bfvTest.kgen
	Sk := bfvTest.sk
	evaluator := bfvTest.evaluator

	SkNew := kgen.NewSecretKey()

	decryptor_SkNew, err := bfvContext.NewDecryptor(SkNew)
	if err != nil {
		log.Fatal(err)
	}

	for _, bitDecomp := range bitDecomps {

		switching_key, err := kgen.NewSwitchingKey(Sk, SkNew, bitDecomp)
		if err != nil {
			log.Fatal(err)
		}

		t.Run(fmt.Sprintf("N=%d/T=%d/Qi=%dx%dbit/Pi=%dx%dbit/bitDecomp=%d/SwitchKeys", bfvTest.bfvcontext.n,
			bfvTest.bfvcontext.t,
			len(bfvTest.bfvcontext.contextQ.Modulus), 60,
			len(bfvTest.bfvcontext.contextP.Modulus), 60,
			bitDecomp), func(t *testing.T) {

			coeffs0, _, ciphertext0, _ := newTestVectors(bfvTest)

			if err := evaluator.SwitchKeys(ciphertext0, switching_key, ciphertext0); err != nil {
				log.Fatal(err)
			}

			plaintext, err := decryptor_SkNew.DecryptNew(ciphertext0)
			if err != nil {
				log.Fatal(err)
			}

			coeffsTest, err := encoder.DecodeUint(plaintext)
			if err != nil {
				log.Fatal(err)
			}

			if EqualSlice(coeffs0.Coeffs[0], coeffsTest) != true {
				t.Errorf("error : switchingKey encrypt/decrypt")
			}

		})

		t.Run(fmt.Sprintf("N=%d/T=%d/Qi=%dx%dbit/Pi=%dx%dbit/bitDecomp=%d/SwitchKeysNew", bfvTest.bfvcontext.n,
			bfvTest.bfvcontext.t,
			len(bfvTest.bfvcontext.contextQ.Modulus), 60,
			len(bfvTest.bfvcontext.contextP.Modulus), 60,
			bitDecomp), func(t *testing.T) {

			coeffs0, _, ciphertext0, _ := newTestVectors(bfvTest)

			ciphertextTest, err := evaluator.SwitchKeysNew(ciphertext0, switching_key)
			if err != nil {
				log.Fatal(err)
			}

			plaintext, err := decryptor_SkNew.DecryptNew(ciphertextTest)
			if err != nil {
				log.Fatal(err)
			}

			coeffsTest, err := encoder.DecodeUint(plaintext)
			if err != nil {
				log.Fatal(err)
			}

			if EqualSlice(coeffs0.Coeffs[0], coeffsTest) != true {
				t.Errorf("error : switchingKeyNew encrypt/decrypt")
			}
		})
	}
}

func test_GaloisEnd(bfvTest *BFVTESTPARAMS, bitDecomps []uint64, t *testing.T) {

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

		rotation_key, err := kgen.NewRotationKeysPow2(Sk, bitDecomp, true)
		if err != nil {
			log.Fatal(err)
		}

		for n := uint64(1); n < bfvContext.n>>1; n <<= 1 {

			for i := uint64(0); i < slots; i++ {
				coeffsWantRotateCol.Coeffs[0][i] = coeffs.Coeffs[0][(i+n)&mask]
				coeffsWantRotateCol.Coeffs[0][i+slots] = coeffs.Coeffs[0][((i+n)&mask)+slots]
			}

			t.Run(fmt.Sprintf("N=%d/T=%d/Qi=%dx%dbit/Pi=%dx%dbit/bitDecomp=%d/RotateColumns/%d", bfvTest.bfvcontext.n,
				bfvTest.bfvcontext.t,
				len(bfvTest.bfvcontext.contextQ.Modulus), 60,
				len(bfvTest.bfvcontext.contextP.Modulus), 60,
				bitDecomp, n), func(t *testing.T) {

				if err := evaluator.RotateColumns(ciphertext, n, rotation_key, receiverCiphertext); err != nil {
					log.Fatal(err)
				}

				verifyTestVectors(bfvTest, coeffsWantRotateCol, receiverCiphertext, t)
			})
		}

		for n := uint64(1); n < bfvContext.n>>1; n <<= 2 {

			rand := ring.RandUniform(mask + 1)

			for i := uint64(0); i < slots; i++ {
				coeffsWantRotateCol.Coeffs[0][i] = coeffs.Coeffs[0][(i+rand)&mask]
				coeffsWantRotateCol.Coeffs[0][i+slots] = coeffs.Coeffs[0][((i+rand)&mask)+slots]
			}

			t.Run(fmt.Sprintf("N=%d/T=%d/Qi=%dx%dbit/Pi=%dx%dbit/bitDecomp=%d/RotateColumns/%d", bfvTest.bfvcontext.n,
				bfvTest.bfvcontext.t,
				len(bfvTest.bfvcontext.contextQ.Modulus), 60,
				len(bfvTest.bfvcontext.contextP.Modulus), 60,
				bitDecomp, rand), func(t *testing.T) {

				if err := evaluator.RotateColumns(ciphertext, rand, rotation_key, receiverCiphertext); err != nil {
					log.Fatal(err)
				}

				verifyTestVectors(bfvTest, coeffsWantRotateCol, receiverCiphertext, t)
			})
		}

		t.Run(fmt.Sprintf("N=%d/T=%d/Qi=%dx%dbit/Pi=%dx%dbit/bitDecomp=%d/RotateRows", bfvTest.bfvcontext.n,
			bfvTest.bfvcontext.t,
			len(bfvTest.bfvcontext.contextQ.Modulus), 60,
			len(bfvTest.bfvcontext.contextP.Modulus), 60,
			bitDecomp), func(t *testing.T) {

			if err := evaluator.RotateRows(ciphertext, rotation_key, receiverCiphertext); err != nil {
				log.Fatal(err)
			}

			verifyTestVectors(bfvTest, coeffsWantRotateRow, receiverCiphertext, t)
		})
	}
}
