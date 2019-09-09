package ckks

import (
	"fmt"
	"github.com/lca1/lattigo/ring"
	"log"
	"math/cmplx"
	"math/rand"
	"testing"
	"time"
)

type CKKSTESTPARAMS struct {
	ckkscontext *CkksContext
	encoder     *Encoder
	kgen        *keygenerator
	sk          *SecretKey
	pk          *PublicKey
	rlk         *EvaluationKey
	rotkey      *RotationKey
	encryptor   *Encryptor
	decryptor   *Decryptor
	evaluator   *Evaluator
}

func randomFloat(min, max float64) float64 {
	return min + rand.Float64()*(max-min)
}

func randomComplex(min, max float64) complex128 {
	return complex(randomFloat(min, max), randomFloat(min, max))
}

func Test_CKKS(t *testing.T) {

	rand.Seed(time.Now().UnixNano())

	var err error

	params := Parameters{9, []uint8{55, 49, 49, 49, 49, 49, 49, 49, 49, 49}, 49, 3.2}

	ckksTest := new(CKKSTESTPARAMS)

	if ckksTest.ckkscontext, err = NewCkksContext(&params); err != nil {
		log.Fatal(err)
	}

	log.Printf("Generating CkksContext for logN=%d/logQ=%d/levels=%d/sigma=%f", ckksTest.ckkscontext.LogN(), ckksTest.ckkscontext.LogQ(), ckksTest.ckkscontext.Levels(), ckksTest.ckkscontext.Sigma())

	ckksTest.encoder = ckksTest.ckkscontext.NewEncoder()

	ckksTest.kgen = ckksTest.ckkscontext.NewKeyGenerator()

	if ckksTest.sk, ckksTest.pk, err = ckksTest.kgen.NewKeyPair(1.0 / 3); err != nil {
		log.Fatal(err)
	}

	log.Printf("Generating relinearization keys")
	if ckksTest.rlk, err = ckksTest.kgen.NewRelinKey(ckksTest.sk, ckksTest.ckkscontext.Scale()); err != nil {
		log.Fatal(err)
	}

	log.Printf("Generating rotation keys for conjugate and powers of 2")
	if ckksTest.rotkey, err = ckksTest.kgen.NewRotationKeysPow2(ckksTest.sk, 15, true); err != nil {
		log.Fatal(err)
	}

	if ckksTest.encryptor, err = ckksTest.ckkscontext.NewEncryptor(ckksTest.pk, ckksTest.sk); err != nil {
		log.Fatal(err)
	}

	if ckksTest.decryptor, err = ckksTest.ckkscontext.NewDecryptor(ckksTest.sk); err != nil {
		log.Fatal(err)
	}

	ckksTest.evaluator = ckksTest.ckkscontext.NewEvaluator()

	test_Encoder(ckksTest, t)

	test_EncryptDecrypt(ckksTest, t)

	test_Add(ckksTest, t)
	test_Sub(ckksTest, t)

	test_AddConst(ckksTest, t)
	test_MulConst(ckksTest, t)
	test_MultByConstAndAdd(ckksTest, t)

	test_ComplexOperations(ckksTest, t)

	test_Rescaling(ckksTest, t)
	test_Mul(ckksTest, t)

	if ckksTest.ckkscontext.Levels() > 8 {
		test_Functions(ckksTest, t)
	} else {
		log.Printf("skipping test_Functions (not enough levels)")
	}

	test_SwitchKeys(ckksTest, t)
	test_Conjugate(ckksTest, t)
	test_RotColumns(ckksTest, t)

	test_MarshalCiphertext(ckksTest, t)
	test_MarshalSecretKey(ckksTest, t)
	test_MarshalPublicKey(ckksTest, t)
	test_MarshalEvaluationKey(ckksTest, t)
	test_MarshalSwitchingKey(ckksTest, t)
	test_MarshalRotationKey(ckksTest, t)

}

func new_test_vectors(params *CKKSTESTPARAMS, a, b float64) (values []complex128, plaintext *Plaintext, ciphertext *Ciphertext, err error) {

	slots := 1 << (params.ckkscontext.logN - 1)

	values = make([]complex128, slots)

	for i := 0; i < slots; i++ {
		values[i] = randomComplex(a, b)
	}

	values[0] = complex(0.607538, 0.555668)

	plaintext = params.ckkscontext.NewPlaintext(params.ckkscontext.Levels()-1, params.ckkscontext.Scale())

	if err = params.encoder.EncodeComplex(plaintext, values); err != nil {
		return nil, nil, nil, err
	}

	ciphertext, err = params.encryptor.EncryptFromPkNew(plaintext)
	if err != nil {
		return nil, nil, nil, err
	}

	return values, plaintext, ciphertext, nil
}

func new_test_vectors_reals(params *CKKSTESTPARAMS, a, b float64) (values []complex128, plaintext *Plaintext, ciphertext *Ciphertext, err error) {

	slots := 1 << (params.ckkscontext.logN - 1)

	values = make([]complex128, slots)

	for i := 0; i < slots; i++ {
		values[i] = complex(randomFloat(a, b), 0)
	}

	values[0] = complex(0.607538, 0)

	plaintext = params.ckkscontext.NewPlaintext(params.ckkscontext.Levels()-1, params.ckkscontext.Scale())

	if err = params.encoder.EncodeComplex(plaintext, values); err != nil {
		return nil, nil, nil, err
	}

	ciphertext, err = params.encryptor.EncryptFromPkNew(plaintext)
	if err != nil {
		return nil, nil, nil, err
	}

	return values, plaintext, ciphertext, nil
}

func verify_test_vectors(params *CKKSTESTPARAMS, valuesWant []complex128, element CkksElement, t *testing.T) (err error) {

	var plaintextTest *Plaintext
	var valuesTest []complex128

	if element.Degree() == 0 {

		plaintextTest = element.(*Plaintext)

	} else {

		if plaintextTest, err = params.decryptor.DecryptNew(element.(*Ciphertext)); err != nil {
			return err
		}
	}

	valuesTest = params.encoder.DecodeComplex(plaintextTest)

	var DeltaReal0, DeltaImag0, DeltaReal1, DeltaImag1 float64

	for i := range valuesWant {

		// Test for big values (> 1)
		DeltaReal0 = real(valuesWant[i]) / real(valuesTest[i])
		DeltaImag0 = imag(valuesWant[i]) / imag(valuesTest[i])

		// Test for small values (< 1)
		DeltaReal1 = real(valuesWant[i]) - real(valuesTest[i])
		DeltaImag1 = imag(valuesWant[i]) - imag(valuesTest[i])

		if DeltaReal1 < 0 {
			DeltaReal1 *= -1
		}
		if DeltaImag1 < 0 {
			DeltaImag1 *= -1
		}

		if (DeltaReal0 < 0.999 || DeltaReal0 > 1.001 || DeltaImag0 < 0.999 || DeltaImag0 > 1.001) && (DeltaReal1 > 0.001 || DeltaImag1 > 0.001) {
			t.Errorf("error : coeff %d, want %f have %f", i, valuesWant[i], valuesTest[i])
			break
		}
	}

	return nil
}

func test_Encoder(params *CKKSTESTPARAMS, t *testing.T) {

	t.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/EncodeDecodeFloat64", params.ckkscontext.logN,
		params.ckkscontext.logQ,
		params.ckkscontext.levels), func(t *testing.T) {
		slots := 1 << (params.ckkscontext.logN - 1)

		valuesWant := make([]float64, slots)
		valuesWantCmplx := make([]complex128, slots)
		for i := 0; i < slots; i++ {
			valuesWant[i] = randomFloat(0.000001, 5)
			valuesWantCmplx[i] = complex(valuesWant[i], 0)
		}

		plaintext := params.ckkscontext.NewPlaintext(params.ckkscontext.Levels()-1, params.ckkscontext.Scale())

		if err := params.encoder.EncodeFloat(plaintext, valuesWant); err != nil {
			log.Fatal(err)
		}

		if err := verify_test_vectors(params, valuesWantCmplx, plaintext, t); err != nil {
			log.Fatal(err)
		}
	})

	t.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/EncodeDecodeComplex128", params.ckkscontext.logN,
		params.ckkscontext.logQ,
		params.ckkscontext.levels), func(t *testing.T) {
		slots := 1 << (params.ckkscontext.logN - 1)

		valuesWant := make([]complex128, slots)

		for i := 0; i < slots; i++ {
			valuesWant[i] = randomComplex(0, 5)
		}

		plaintext := params.ckkscontext.NewPlaintext(params.ckkscontext.Levels()-1, params.ckkscontext.Scale())

		if err := params.encoder.EncodeComplex(plaintext, valuesWant); err != nil {
			log.Fatal(err)
		}

		if err := verify_test_vectors(params, valuesWant, plaintext, t); err != nil {
			log.Fatal(err)
		}
	})
}

func test_EncryptDecrypt(params *CKKSTESTPARAMS, t *testing.T) {

	t.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/EncryptFromPk", params.ckkscontext.logN,
		params.ckkscontext.logQ,
		params.ckkscontext.levels), func(t *testing.T) {
		var err error

		slots := 1 << (params.ckkscontext.logN - 1)

		valuesWant := make([]complex128, slots)

		for i := 0; i < slots; i++ {
			valuesWant[i] = randomComplex(0, 5)
		}

		plaintext := params.ckkscontext.NewPlaintext(params.ckkscontext.Levels()-1, params.ckkscontext.Scale())

		if err = params.encoder.EncodeComplex(plaintext, valuesWant); err != nil {
			log.Fatal(err)
		}

		ciphertext, err := params.encryptor.EncryptFromPkNew(plaintext)
		if err != nil {
			log.Fatal(err)
		}

		if err := verify_test_vectors(params, valuesWant, ciphertext, t); err != nil {
			log.Fatal(err)
		}
	})

	t.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/EncryptFromSk", params.ckkscontext.logN,
		params.ckkscontext.logQ,
		params.ckkscontext.levels), func(t *testing.T) {
		var err error

		slots := 1 << (params.ckkscontext.logN - 1)

		valuesWant := make([]complex128, slots)

		for i := 0; i < slots; i++ {
			valuesWant[i] = randomComplex(0, 5)
		}

		plaintext := params.ckkscontext.NewPlaintext(params.ckkscontext.Levels()-1, params.ckkscontext.Scale())

		if err = params.encoder.EncodeComplex(plaintext, valuesWant); err != nil {
			log.Fatal(err)
		}

		ciphertext, err := params.encryptor.EncryptFromSkNew(plaintext)
		if err != nil {
			log.Fatal(err)
		}

		if err := verify_test_vectors(params, valuesWant, ciphertext, t); err != nil {
			log.Fatal(err)
		}
	})
}

func test_Add(params *CKKSTESTPARAMS, t *testing.T) {

	t.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/AddCtCtInPlace", params.ckkscontext.logN,
		params.ckkscontext.logQ,
		params.ckkscontext.levels), func(t *testing.T) {

		values1, _, ciphertext1, err := new_test_vectors(params, -1, 1)
		if err != nil {
			log.Fatal(err)
		}

		values2, _, ciphertext2, err := new_test_vectors(params, -1, 1)
		if err != nil {
			log.Fatal(err)
		}

		valuesWant := make([]complex128, params.ckkscontext.n>>1)
		for i := 0; i < len(valuesWant); i++ {
			valuesWant[i] = values1[i] + values2[i]
		}

		if err := params.evaluator.Add(ciphertext1, ciphertext2, ciphertext1); err != nil {
			log.Fatal(err)
		}

		if err := verify_test_vectors(params, valuesWant, ciphertext1, t); err != nil {
			log.Fatal(err)
		}
	})

	t.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/AddCtCt", params.ckkscontext.logN,
		params.ckkscontext.logQ,
		params.ckkscontext.levels), func(t *testing.T) {

		values1, _, ciphertext1, err := new_test_vectors(params, -1, 1)
		if err != nil {
			log.Fatal(err)
		}

		values2, _, ciphertext2, err := new_test_vectors(params, -1, 1)
		if err != nil {
			log.Fatal(err)
		}

		receiver := params.ckkscontext.NewCiphertext(1, ciphertext1.Level(), ciphertext1.Scale())

		valuesWant := make([]complex128, params.ckkscontext.n>>1)
		for i := 0; i < len(valuesWant); i++ {
			valuesWant[i] = values1[i] + values2[i]
		}

		if err := params.evaluator.Add(ciphertext1, ciphertext2, receiver); err != nil {
			log.Fatal(err)
		}

		if err := verify_test_vectors(params, valuesWant, receiver, t); err != nil {
			log.Fatal(err)
		}
	})

	t.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/Add(Ct,Plain)", params.ckkscontext.logN,
		params.ckkscontext.logQ,
		params.ckkscontext.levels), func(t *testing.T) {

		values1, _, ciphertext1, err := new_test_vectors(params, -1, 1)
		if err != nil {
			log.Fatal(err)
		}

		values2, plaintext2, _, err := new_test_vectors(params, -1, 1)
		if err != nil {
			log.Fatal(err)
		}

		valuesWant := make([]complex128, params.ckkscontext.n>>1)
		for i := 0; i < len(valuesWant); i++ {
			valuesWant[i] = values1[i] + values2[i]
		}

		if err := params.evaluator.Add(ciphertext1, plaintext2, ciphertext1); err != nil {
			log.Fatal(err)
		}

		if err := verify_test_vectors(params, valuesWant, ciphertext1, t); err != nil {
			log.Fatal(err)
		}
	})

	t.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/Add(Plain,Ct)", params.ckkscontext.logN,
		params.ckkscontext.logQ,
		params.ckkscontext.levels), func(t *testing.T) {

		values1, plaintext1, _, err := new_test_vectors(params, -1, 1)
		if err != nil {
			log.Fatal(err)
		}

		values2, _, ciphertext2, err := new_test_vectors(params, -1, 1)
		if err != nil {
			log.Fatal(err)
		}

		valuesWant := make([]complex128, params.ckkscontext.n>>1)
		for i := 0; i < len(valuesWant); i++ {
			valuesWant[i] = values1[i] + values2[i]
		}

		if err := params.evaluator.Add(plaintext1, ciphertext2, ciphertext2); err != nil {
			log.Fatal(err)
		}

		if err := verify_test_vectors(params, valuesWant, ciphertext2, t); err != nil {
			log.Fatal(err)
		}
	})

	t.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/Add(Plain,Plain)", params.ckkscontext.logN,
		params.ckkscontext.logQ,
		params.ckkscontext.levels), func(t *testing.T) {

		values1, plaintext1, _, err := new_test_vectors(params, -1, 1)
		if err != nil {
			log.Fatal(err)
		}

		values2, plaintext2, _, err := new_test_vectors(params, -1, 1)
		if err != nil {
			log.Fatal(err)
		}

		valuesWant := make([]complex128, params.ckkscontext.n>>1)
		for i := 0; i < len(valuesWant); i++ {
			valuesWant[i] = values1[i] + values2[i]
		}

		if err := params.evaluator.Add(plaintext1, plaintext2, plaintext1); err != nil {
			log.Fatal(err)
		}

		if err := verify_test_vectors(params, valuesWant, plaintext1, t); err != nil {
			log.Fatal(err)
		}
	})
}

func test_Sub(params *CKKSTESTPARAMS, t *testing.T) {

	t.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/Sub(Ct,Ct)", params.ckkscontext.logN,
		params.ckkscontext.logQ,
		params.ckkscontext.levels), func(t *testing.T) {

		values1, _, ciphertext1, err := new_test_vectors(params, -1, 1)
		if err != nil {
			log.Fatal(err)
		}

		values2, _, ciphertext2, err := new_test_vectors(params, -1, 1)
		if err != nil {
			log.Fatal(err)
		}

		valuesWant := make([]complex128, params.ckkscontext.n>>1)
		for i := 0; i < len(valuesWant); i++ {
			valuesWant[i] = values1[i] - values2[i]
		}

		if err := params.evaluator.Sub(ciphertext1, ciphertext2, ciphertext1); err != nil {
			log.Fatal(err)
		}

		if err := verify_test_vectors(params, valuesWant, ciphertext1, t); err != nil {
			log.Fatal(err)
		}
	})

	t.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/Sub(Ct,Plain)", params.ckkscontext.logN,
		params.ckkscontext.logQ,
		params.ckkscontext.levels), func(t *testing.T) {

		values1, _, ciphertext1, err := new_test_vectors(params, -1, 1)
		if err != nil {
			log.Fatal(err)
		}

		values2, plaintext2, _, err := new_test_vectors(params, -1, 1)
		if err != nil {
			log.Fatal(err)
		}

		valuesWant := make([]complex128, params.ckkscontext.n>>1)
		for i := 0; i < len(valuesWant); i++ {
			valuesWant[i] = values1[i] - values2[i]
		}

		if err := params.evaluator.Sub(ciphertext1, plaintext2, ciphertext1); err != nil {
			log.Fatal(err)
		}

		if err := verify_test_vectors(params, valuesWant, ciphertext1, t); err != nil {
			log.Fatal(err)
		}
	})

	t.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/Sub(Plain,Ct)", params.ckkscontext.logN,
		params.ckkscontext.logQ,
		params.ckkscontext.levels), func(t *testing.T) {

		values1, plaintext1, _, err := new_test_vectors(params, -1, 1)
		if err != nil {
			log.Fatal(err)
		}

		values2, _, ciphertext2, err := new_test_vectors(params, -1, 1)
		if err != nil {
			log.Fatal(err)
		}

		valuesWant := make([]complex128, params.ckkscontext.n>>1)
		for i := 0; i < len(valuesWant); i++ {
			valuesWant[i] = values1[i] - values2[i]
		}

		if err := params.evaluator.Sub(plaintext1, ciphertext2, ciphertext2); err != nil {
			log.Fatal(err)
		}

		if err := verify_test_vectors(params, valuesWant, ciphertext2, t); err != nil {
			log.Fatal(err)
		}
	})

	t.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/Sub(Plain,Plain)", params.ckkscontext.logN,
		params.ckkscontext.logQ,
		params.ckkscontext.levels), func(t *testing.T) {

		values1, plaintext1, _, err := new_test_vectors(params, -1, 1)
		if err != nil {
			log.Fatal(err)
		}

		values2, plaintext2, _, err := new_test_vectors(params, -1, 1)
		if err != nil {
			log.Fatal(err)
		}

		valuesWant := make([]complex128, params.ckkscontext.n>>1)
		for i := 0; i < len(valuesWant); i++ {
			valuesWant[i] = values1[i] - values2[i]
		}

		if err := params.evaluator.Sub(plaintext1, plaintext2, plaintext1); err != nil {
			log.Fatal(err)
		}

		if err := verify_test_vectors(params, valuesWant, plaintext1, t); err != nil {
			log.Fatal(err)
		}
	})
}

func test_AddConst(params *CKKSTESTPARAMS, t *testing.T) {

	t.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/AddCmplx(Ct,complex128)", params.ckkscontext.logN,
		params.ckkscontext.logQ,
		params.ckkscontext.levels), func(t *testing.T) {

		values1, _, ciphertext1, err := new_test_vectors(params, -1, 1)
		if err != nil {
			log.Fatal(err)
		}

		constant := complex(3.1415, -1.4142)

		valuesWant := make([]complex128, params.ckkscontext.n>>1)
		for i := 0; i < len(valuesWant); i++ {
			valuesWant[i] = values1[i] + constant
		}

		if err := params.evaluator.AddConst(ciphertext1, constant, ciphertext1); err != nil {
			log.Fatal(err)
		}

		if err := verify_test_vectors(params, valuesWant, ciphertext1, t); err != nil {
			log.Fatal(err)
		}
	})

	t.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/AddCmplx(Plain,complex128)", params.ckkscontext.logN,
		params.ckkscontext.logQ,
		params.ckkscontext.levels), func(t *testing.T) {

		values1, plaintext1, _, err := new_test_vectors(params, -1, 1)
		if err != nil {
			log.Fatal(err)
		}

		constant := complex(3.1415, -1.4142)

		valuesWant := make([]complex128, params.ckkscontext.n>>1)
		for i := 0; i < len(valuesWant); i++ {
			valuesWant[i] = values1[i] + constant
		}

		if err := params.evaluator.AddConst(plaintext1, constant, plaintext1); err != nil {
			log.Fatal(err)
		}

		if err := verify_test_vectors(params, valuesWant, plaintext1, t); err != nil {
			log.Fatal(err)
		}
	})
}

func test_MulConst(params *CKKSTESTPARAMS, t *testing.T) {

	t.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/MultCmplx(Ct,complex128)", params.ckkscontext.logN,
		params.ckkscontext.logQ,
		params.ckkscontext.levels), func(t *testing.T) {

		values1, _, ciphertext1, err := new_test_vectors(params, -5, 5)
		if err != nil {
			log.Fatal(err)
		}

		constant := complex(1.4142, -3.1415)

		valuesWant := make([]complex128, params.ckkscontext.n>>1)
		for i := 0; i < len(valuesWant); i++ {
			valuesWant[i] = values1[i] * (1 / constant)
		}

		if err = params.evaluator.MultConst(ciphertext1, 1/constant, ciphertext1); err != nil {
			log.Fatal(err)
		}

		if err := verify_test_vectors(params, valuesWant, ciphertext1, t); err != nil {
			log.Fatal(err)
		}
	})

	t.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/MultCmplx(Plain,complex128)", params.ckkscontext.logN,
		params.ckkscontext.logQ,
		params.ckkscontext.levels), func(t *testing.T) {

		values1, plaintext1, _, err := new_test_vectors(params, -1, 1)
		if err != nil {
			log.Fatal(err)
		}

		constant := complex(3.1415, -1.4142)

		valuesWant := make([]complex128, params.ckkscontext.n>>1)
		for i := 0; i < len(valuesWant); i++ {
			valuesWant[i] = values1[i] * constant
		}

		if err = params.evaluator.MultConst(plaintext1, constant, plaintext1); err != nil {
			log.Fatal(err)
		}

		if err := verify_test_vectors(params, valuesWant, plaintext1, t); err != nil {
			log.Fatal(err)
		}
	})

	t.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/MultByi", params.ckkscontext.logN,
		params.ckkscontext.logQ,
		params.ckkscontext.levels), func(t *testing.T) {

		values, _, ciphertext1, err := new_test_vectors(params, -1, 1)
		if err != nil {
			log.Fatal(err)
		}

		valuesWant := make([]complex128, params.ckkscontext.n>>1)

		for i := 0; i < len(valuesWant); i++ {
			valuesWant[i] = values[i] * complex(0, 1)
		}

		if err = params.evaluator.MultByi(ciphertext1, ciphertext1); err != nil {
			log.Fatal(err)
		}

		if err := verify_test_vectors(params, valuesWant, ciphertext1, t); err != nil {
			log.Fatal(err)
		}
	})

	t.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/DivByi", params.ckkscontext.logN,
		params.ckkscontext.logQ,
		params.ckkscontext.levels), func(t *testing.T) {

		values, _, ciphertext1, err := new_test_vectors(params, -1, 1)
		if err != nil {
			log.Fatal(err)
		}

		valuesWant := make([]complex128, params.ckkscontext.n>>1)

		for i := 0; i < len(valuesWant); i++ {
			valuesWant[i] = values[i] / complex(0, 1)
		}

		if err = params.evaluator.DivByi(ciphertext1, ciphertext1); err != nil {
			log.Fatal(err)
		}

		if err := verify_test_vectors(params, valuesWant, ciphertext1, t); err != nil {
			log.Fatal(err)
		}
	})
}

func test_MultByConstAndAdd(params *CKKSTESTPARAMS, t *testing.T) {

	t.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/MultByCmplxAndAdd(Ct0, complex128, Ct1)", params.ckkscontext.logN,
		params.ckkscontext.logQ,
		params.ckkscontext.levels), func(t *testing.T) {

		values1, _, ciphertext1, err := new_test_vectors(params, -1, 1)
		if err != nil {
			log.Fatal(err)
		}

		values2, _, ciphertext2, err := new_test_vectors(params, -1, 1)
		if err != nil {
			log.Fatal(err)
		}

		constant := complex(3.1415, -1.4142)

		valuesWant := make([]complex128, params.ckkscontext.n>>1)
		for i := 0; i < len(valuesWant); i++ {
			values2[i] += (values1[i] * constant) + (values1[i] * constant)
		}

		if err = params.evaluator.MultByConstAndAdd(ciphertext1, constant, ciphertext2); err != nil {
			log.Fatal(err)
		}

		params.evaluator.Rescale(ciphertext1, ciphertext1)

		if err = params.evaluator.MultByConstAndAdd(ciphertext1, constant, ciphertext2); err != nil {
			log.Fatal(err)
		}

		if err := verify_test_vectors(params, values2, ciphertext2, t); err != nil {
			log.Fatal(err)
		}
	})
}

func test_ComplexOperations(params *CKKSTESTPARAMS, t *testing.T) {

	t.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/ExtractImag", params.ckkscontext.logN,
		params.ckkscontext.logQ,
		params.ckkscontext.levels), func(t *testing.T) {

		values, _, ciphertext, err := new_test_vectors(params, -1, 1)
		if err != nil {
			log.Fatal(err)
		}

		for i := 0; i < len(values); i++ {
			values[i] = complex(imag(values[i]), 0)
		}

		if err = params.evaluator.ExtractImag(ciphertext, params.rotkey, ciphertext); err != nil {
			log.Fatal(err)
		}

		if err := verify_test_vectors(params, values, ciphertext, t); err != nil {
			log.Fatal(err)
		}
	})

	t.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/SwapRealImag", params.ckkscontext.logN,
		params.ckkscontext.logQ,
		params.ckkscontext.levels), func(t *testing.T) {

		values, _, ciphertext, err := new_test_vectors(params, -1, 1)
		if err != nil {
			log.Fatal(err)
		}

		for i := 0; i < len(values); i++ {
			values[i] = complex(imag(values[i]), real(values[i]))
		}

		if err = params.evaluator.SwapRealImag(ciphertext, params.rotkey, ciphertext); err != nil {
			log.Fatal(err)
		}

		if err := verify_test_vectors(params, values, ciphertext, t); err != nil {
			log.Fatal(err)
		}
	})

	t.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/RemoveReal", params.ckkscontext.logN,
		params.ckkscontext.logQ,
		params.ckkscontext.levels), func(t *testing.T) {

		values, _, ciphertext, err := new_test_vectors(params, -1, 1)
		if err != nil {
			log.Fatal(err)
		}

		for i := 0; i < len(values); i++ {
			values[i] = complex(0, imag(values[i]))
		}

		if err = params.evaluator.RemoveReal(ciphertext, params.rotkey, ciphertext); err != nil {
			log.Fatal(err)
		}

		if err := verify_test_vectors(params, values, ciphertext, t); err != nil {
			log.Fatal(err)
		}
	})

	t.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/RemoveImag", params.ckkscontext.logN,
		params.ckkscontext.logQ,
		params.ckkscontext.levels), func(t *testing.T) {

		values, _, ciphertext, err := new_test_vectors(params, -1, 1)
		if err != nil {
			log.Fatal(err)
		}

		for i := 0; i < len(values); i++ {
			values[i] = complex(real(values[i]), 0)
		}

		if err = params.evaluator.RemoveImag(ciphertext, params.rotkey, ciphertext); err != nil {
			log.Fatal(err)
		}

		if err := verify_test_vectors(params, values, ciphertext, t); err != nil {
			log.Fatal(err)
		}
	})
}

func test_Rescaling(params *CKKSTESTPARAMS, t *testing.T) {

	t.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/Rescaling", params.ckkscontext.logN,
		params.ckkscontext.logQ,
		params.ckkscontext.levels), func(t *testing.T) {

		coeffs := make([]*ring.Int, params.ckkscontext.n)
		for i := uint64(0); i < params.ckkscontext.n; i++ {
			coeffs[i] = ring.RandInt(params.ckkscontext.contextLevel[params.ckkscontext.Levels()-1].ModulusBigint)
			coeffs[i].Div(coeffs[i], ring.NewUint(10))
		}

		coeffsWant := make([]*ring.Int, params.ckkscontext.contextLevel[params.ckkscontext.Levels()-1].N)
		for i := range coeffs {
			coeffsWant[i] = ring.Copy(coeffs[i])
			coeffsWant[i].Div(coeffsWant[i], ring.NewUint(params.ckkscontext.moduli[len(params.ckkscontext.moduli)-1]))
		}

		polTest := params.ckkscontext.contextLevel[params.ckkscontext.Levels()-1].NewPoly()
		polWant := params.ckkscontext.contextLevel[params.ckkscontext.Levels()-1].NewPoly()

		params.ckkscontext.contextLevel[params.ckkscontext.Levels()-1].SetCoefficientsBigint(coeffs, polTest)
		params.ckkscontext.contextLevel[params.ckkscontext.Levels()-1].SetCoefficientsBigint(coeffsWant, polWant)

		params.ckkscontext.contextLevel[params.ckkscontext.Levels()-1].NTT(polTest, polTest)
		params.ckkscontext.contextLevel[params.ckkscontext.Levels()-1].NTT(polWant, polWant)

		rescale(params.evaluator, polTest, polTest)

		for i := uint64(0); i < params.ckkscontext.n; i++ {
			for j := 0; j < len(params.ckkscontext.moduli)-1; j++ {
				if polWant.Coeffs[j][i] != polTest.Coeffs[j][i] {
					t.Errorf("error : coeff %v Qi%v = %s, want %v have %v", i, j, coeffs[i].String(), polWant.Coeffs[j][i], polTest.Coeffs[j][i])
					break
				}
			}
		}
	})
}

func test_Mul(params *CKKSTESTPARAMS, t *testing.T) {

	t.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/Mul(Ct,Ct)", params.ckkscontext.logN,
		params.ckkscontext.logQ,
		params.ckkscontext.levels), func(t *testing.T) {

		values1, _, ciphertext1, err := new_test_vectors(params, -1, 1)
		if err != nil {
			log.Fatal(err)
		}

		values2, _, ciphertext2, err := new_test_vectors(params, -1, 1)
		if err != nil {
			log.Fatal(err)
		}

		// up to a level equal to 2 modulus
		valuesWant := make([]complex128, params.ckkscontext.n>>1)

		for i := 0; i < len(valuesWant); i++ {
			valuesWant[i] = values1[i] * values2[i]
		}

		if err = params.evaluator.MulRelin(ciphertext1, ciphertext2, nil, ciphertext1); err != nil {
			log.Fatal(err)
		}

		if err := verify_test_vectors(params, valuesWant, ciphertext1, t); err != nil {
			log.Fatal(err)
		}
	})

	t.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/Relinearize", params.ckkscontext.logN,
		params.ckkscontext.logQ,
		params.ckkscontext.levels), func(t *testing.T) {

		values1, _, ciphertext1, err := new_test_vectors(params, -1, 1)
		if err != nil {
			log.Fatal(err)
		}

		values2, _, ciphertext2, err := new_test_vectors(params, -1, 1)
		if err != nil {
			log.Fatal(err)
		}

		// up to a level equal to 2 modulus
		valuesWant := make([]complex128, params.ckkscontext.n>>1)

		for i := 0; i < len(valuesWant); i++ {
			valuesWant[i] = values1[i] * values2[i]
		}

		if err = params.evaluator.MulRelin(ciphertext1, ciphertext2, nil, ciphertext1); err != nil {
			log.Fatal(err)
		}

		if err = params.evaluator.Relinearize(ciphertext1, params.rlk, ciphertext1); err != nil {
			log.Fatal(err)
		}

		if err := verify_test_vectors(params, valuesWant, ciphertext1, t); err != nil {
			log.Fatal(err)
		}
	})

	t.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/MulRelin(Ct,Ct)->Rescale", params.ckkscontext.logN,
		params.ckkscontext.logQ,
		params.ckkscontext.levels), func(t *testing.T) {

		values1, _, ciphertext1, err := new_test_vectors(params, -1, 1)
		if err != nil {
			log.Fatal(err)
		}

		values2, _, ciphertext2, err := new_test_vectors(params, -1, 1)
		if err != nil {
			log.Fatal(err)
		}

		// up to a level equal to 2 modulus
		valuesWant := make([]complex128, params.ckkscontext.n>>1)

		for i := 0; i < len(valuesWant); i++ {
			valuesWant[i] = values1[i]
		}

		for i := uint64(0); i < params.ckkscontext.Levels()-1; i++ {

			for i := 0; i < len(valuesWant); i++ {
				valuesWant[i] *= values2[i]
			}

			if err = params.evaluator.MulRelin(ciphertext1, ciphertext2, params.rlk, ciphertext1); err != nil {
				log.Fatal(err)
			}

			if err = params.evaluator.Rescale(ciphertext1, ciphertext1); err != nil {
				log.Fatal(err)
			}
		}

		if err := verify_test_vectors(params, valuesWant, ciphertext1, t); err != nil {
			log.Fatal(err)
		}
	})

	t.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/MulRelin(Ct,Plain)->Rescale", params.ckkscontext.logN,
		params.ckkscontext.logQ,
		params.ckkscontext.levels), func(t *testing.T) {

		values1, _, ciphertext1, err := new_test_vectors(params, -1, 1)
		if err != nil {
			log.Fatal(err)
		}

		values2, plaintext2, _, err := new_test_vectors(params, -1, 1)
		if err != nil {
			log.Fatal(err)
		}

		// up to a level equal to 2 modulus
		valuesWant := make([]complex128, params.ckkscontext.n>>1)

		for i := 0; i < len(valuesWant); i++ {
			valuesWant[i] = values1[i]
		}

		for i := uint64(0); i < params.ckkscontext.Levels()-1; i++ {

			for i := 0; i < len(valuesWant); i++ {
				valuesWant[i] *= values2[i]
			}

			if err = params.evaluator.MulRelin(ciphertext1, plaintext2, params.rlk, ciphertext1); err != nil {
				log.Fatal(err)
			}

			if err = params.evaluator.Rescale(ciphertext1, ciphertext1); err != nil {
				log.Fatal(err)
			}
		}

		if err := verify_test_vectors(params, valuesWant, ciphertext1, t); err != nil {
			log.Fatal(err)
		}
	})

	t.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/MulRelin(Plain,Ct)->Rescale", params.ckkscontext.logN,
		params.ckkscontext.logQ,
		params.ckkscontext.levels), func(t *testing.T) {

		values1, plaintext1, _, err := new_test_vectors(params, -1, 1)
		if err != nil {
			log.Fatal(err)
		}

		values2, _, ciphertext2, err := new_test_vectors(params, -1, 1)
		if err != nil {
			log.Fatal(err)
		}

		// up to a level equal to 2 modulus
		valuesWant := make([]complex128, params.ckkscontext.n>>1)

		for i := 0; i < len(valuesWant); i++ {
			valuesWant[i] = values2[i]
		}

		for i := uint64(0); i < params.ckkscontext.Levels()-1; i++ {

			for i := 0; i < len(valuesWant); i++ {
				valuesWant[i] *= values1[i]
			}

			if err = params.evaluator.MulRelin(plaintext1, ciphertext2, params.rlk, ciphertext2); err != nil {
				log.Fatal(err)
			}

			if err = params.evaluator.Rescale(ciphertext2, ciphertext2); err != nil {
				log.Fatal(err)
			}
		}

		if err := verify_test_vectors(params, valuesWant, ciphertext2, t); err != nil {
			log.Fatal(err)
		}
	})

	t.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/MulRelin(Plain,Plain)->Rescale", params.ckkscontext.logN,
		params.ckkscontext.logQ,
		params.ckkscontext.levels), func(t *testing.T) {

		values1, plaintext1, _, err := new_test_vectors(params, -1, 1)
		if err != nil {
			log.Fatal(err)
		}

		values2, plaintext2, _, err := new_test_vectors(params, -1, 1)
		if err != nil {
			log.Fatal(err)
		}

		// up to a level equal to 2 modulus
		valuesWant := make([]complex128, params.ckkscontext.n>>1)

		for i := 0; i < len(valuesWant); i++ {
			valuesWant[i] = values1[i]
		}

		for i := uint64(0); i < params.ckkscontext.Levels()-1; i++ {

			for i := 0; i < len(valuesWant); i++ {
				valuesWant[i] *= values2[i]
			}

			if err = params.evaluator.MulRelin(plaintext1, plaintext2, params.rlk, plaintext1); err != nil {
				log.Fatal(err)
			}

			if err = params.evaluator.Rescale(plaintext1, plaintext1); err != nil {
				log.Fatal(err)
			}
		}

		if err := verify_test_vectors(params, valuesWant, plaintext1, t); err != nil {
			log.Fatal(err)
		}
	})
}

func test_Functions(params *CKKSTESTPARAMS, t *testing.T) {

	t.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/Square", params.ckkscontext.logN,
		params.ckkscontext.logQ,
		params.ckkscontext.levels), func(t *testing.T) {

		values1, _, ciphertext1, err := new_test_vectors(params, -1, 1)
		if err != nil {
			log.Fatal(err)
		}

		// up to a level equal to 2 modulus
		valuesWant := make([]complex128, params.ckkscontext.n>>1)

		for i := 0; i < len(valuesWant); i++ {
			valuesWant[i] = values1[i]
		}

		for i := uint64(0); i < 1; i++ {

			for j := 0; j < len(valuesWant); j++ {
				valuesWant[j] *= valuesWant[j]
			}

			if err = params.evaluator.MulRelin(ciphertext1, ciphertext1, params.rlk, ciphertext1); err != nil {
				log.Fatal(err)
			}

			params.evaluator.Rescale(ciphertext1, ciphertext1)

		}

		if err := verify_test_vectors(params, valuesWant, ciphertext1, t); err != nil {
			log.Fatal(err)
		}
	})

	t.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/PowerOf2", params.ckkscontext.logN,
		params.ckkscontext.logQ,
		params.ckkscontext.levels), func(t *testing.T) {

		values1, _, ciphertext1, err := new_test_vectors(params, -1, 1)
		if err != nil {
			log.Fatal(err)
		}

		// up to a level equal to 2 modulus
		valuesWant := make([]complex128, params.ckkscontext.n>>1)

		for i := 0; i < len(valuesWant); i++ {
			valuesWant[i] = values1[i]
		}

		var n uint64

		n = 2

		for i := uint64(0); i < n; i++ {
			for j := 0; j < len(valuesWant); j++ {
				valuesWant[j] *= valuesWant[j]
			}
		}

		if err = params.evaluator.PowerOf2(ciphertext1, n, params.rlk, ciphertext1); err != nil {
			log.Fatal(err)
		}

		if ciphertext1.Scale() >= 100 {
			params.evaluator.Rescale(ciphertext1, ciphertext1)
		}

		if err := verify_test_vectors(params, valuesWant, ciphertext1, t); err != nil {
			log.Fatal(err)
		}
	})

	t.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/Power", params.ckkscontext.logN,
		params.ckkscontext.logQ,
		params.ckkscontext.levels), func(t *testing.T) {

		values1, _, ciphertext1, err := new_test_vectors(params, -1, 1)
		if err != nil {
			log.Fatal(err)
		}

		valuesWant := make([]complex128, params.ckkscontext.n>>1)
		tmp := make([]complex128, params.ckkscontext.n>>1)

		for i := 0; i < len(valuesWant); i++ {
			valuesWant[i] = complex(1, 0)
			tmp[i] = values1[i]
		}

		var n uint64

		n = 7

		for j := 0; j < len(valuesWant); j++ {
			for i := n; i > 0; i >>= 1 {

				if i&1 == 1 {
					valuesWant[j] *= tmp[j]
				}

				tmp[j] *= tmp[j]
			}
		}

		if err = params.evaluator.Power(ciphertext1, n, params.rlk, ciphertext1); err != nil {
			log.Fatal(err)
		}

		if ciphertext1.Scale() >= 100 {
			params.evaluator.Rescale(ciphertext1, ciphertext1)
		}

		if err := verify_test_vectors(params, valuesWant, ciphertext1, t); err != nil {
			log.Fatal(err)
		}
	})

	t.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/Inverse", params.ckkscontext.logN,
		params.ckkscontext.logQ,
		params.ckkscontext.levels), func(t *testing.T) {

		values, _, ciphertext1, err := new_test_vectors_reals(params, 0.1, 1)
		if err != nil {
			log.Fatal(err)
		}

		valuesWant := make([]complex128, params.ckkscontext.n>>1)

		for i := 0; i < len(valuesWant); i++ {
			valuesWant[i] = complex(1, 0) / values[i]
		}

		if ciphertext1, err = params.evaluator.InverseNew(ciphertext1, 7, params.rlk); err != nil {
			log.Fatal(err)
		}

		if err := verify_test_vectors(params, valuesWant, ciphertext1, t); err != nil {
			log.Fatal(err)
		}
	})

	t.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/sin(x) [-1-1i, 1+1i] deg16", params.ckkscontext.logN,
		params.ckkscontext.logQ,
		params.ckkscontext.levels), func(t *testing.T) {

		values, _, ciphertext1, err := new_test_vectors(params, -1, 1)
		if err != nil {
			log.Fatal(err)
		}

		valuesWant := make([]complex128, params.ckkscontext.n>>1)

		for i := 0; i < len(valuesWant); i++ {
			valuesWant[i] = cmplx.Sin(values[i])
		}

		cheby := Approximate(cmplx.Sin, complex(-1, -1), complex(1, 1), 16)

		if ciphertext1, err = params.evaluator.EvaluateCheby(ciphertext1, cheby, params.rlk); err != nil {
			log.Fatal(err)
		}

		if err := verify_test_vectors(params, valuesWant, ciphertext1, t); err != nil {
			log.Fatal(err)
		}
	})

	t.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/exp(2*pi*i*x) [-1, 1] deg60", params.ckkscontext.logN,
		params.ckkscontext.logQ,
		params.ckkscontext.levels), func(t *testing.T) {

		values, _, ciphertext1, err := new_test_vectors_reals(params, -1, 1)
		if err != nil {
			log.Fatal(err)
		}

		valuesWant := make([]complex128, params.ckkscontext.n>>1)

		for i := 0; i < len(valuesWant); i++ {
			valuesWant[i] = exp2pi(values[i])
		}

		cheby := Approximate(exp2pi, -1, 1, 60)

		if ciphertext1, err = params.evaluator.EvaluateCheby(ciphertext1, cheby, params.rlk); err != nil {
			log.Fatal(err)
		}

		if err := verify_test_vectors(params, valuesWant, ciphertext1, t); err != nil {
			log.Fatal(err)
		}
	})

	t.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/sin(2*pi*x)/(2*pi) [-15, 15] deg128", params.ckkscontext.logN,
		params.ckkscontext.logQ,
		params.ckkscontext.levels), func(t *testing.T) {

		values, _, ciphertext1, err := new_test_vectors_reals(params, -15, 15)
		if err != nil {
			log.Fatal(err)
		}

		valuesWant := make([]complex128, params.ckkscontext.n>>1)

		for i := 0; i < len(valuesWant); i++ {
			valuesWant[i] = sin2pi2pi(values[i])
		}

		cheby := Approximate(sin2pi2pi, -15, 15, 128)

		if ciphertext1, err = params.evaluator.EvaluateCheby(ciphertext1, cheby, params.rlk); err != nil {
			log.Fatal(err)
		}

		if err := verify_test_vectors(params, valuesWant, ciphertext1, t); err != nil {
			log.Fatal(err)
		}
	})
}

func test_SwitchKeys(params *CKKSTESTPARAMS, t *testing.T) {

	t.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/SwitchKeys", params.ckkscontext.logN,
		params.ckkscontext.logQ,
		params.ckkscontext.levels), func(t *testing.T) {

		slots := 1 << (params.ckkscontext.logN - 1)

		valuesWant := make([]complex128, slots)

		for i := 0; i < slots; i++ {
			valuesWant[i] = randomComplex(0, 5)
		}

		plaintext := params.ckkscontext.NewPlaintext(params.ckkscontext.Levels()-1, params.ckkscontext.Scale())

		if err := params.encoder.EncodeComplex(plaintext, valuesWant); err != nil {
			log.Fatal(err)
		}

		ciphertext, err := params.encryptor.EncryptFromPkNew(plaintext)
		if err != nil {
			log.Fatal(err)
		}

		sk2, _ := params.kgen.NewSecretKey(1.0 / 3)

		switchingkeys, err := params.kgen.NewSwitchingKey(params.sk, sk2, 10)
		if err != nil {
			log.Fatal(err)
		}

		if err = params.evaluator.SwitchKeys(ciphertext, switchingkeys, ciphertext); err != nil {
			log.Fatal(err)
		}

		decryptorSk2, err := params.ckkscontext.NewDecryptor(sk2)
		if err != nil {
			log.Fatal(err)
		}

		plaintextTest, err := decryptorSk2.DecryptNew(ciphertext)
		if err != nil {
			log.Fatal(err)
		}

		valuesTest := params.encoder.DecodeComplex(plaintextTest)

		var DeltaReal0, DeltaImag0, DeltaReal1, DeltaImag1 float64

		for i := range valuesWant {

			// Test for big values (> 1)
			DeltaReal0 = real(valuesWant[i]) / real(valuesTest[i])
			DeltaImag0 = imag(valuesWant[i]) / imag(valuesTest[i])

			// Test for small values (< 1)
			DeltaReal1 = real(valuesWant[i]) - real(valuesTest[i])
			DeltaImag1 = imag(valuesWant[i]) - imag(valuesTest[i])

			if DeltaReal1 < 0 {
				DeltaReal1 *= -1
			}
			if DeltaImag1 < 0 {
				DeltaImag1 *= -1
			}

			if (DeltaReal0 < 0.999 || DeltaReal0 > 1.001 || DeltaImag0 < 0.999 || DeltaImag0 > 1.001) && (DeltaReal1 > 0.001 || DeltaImag1 > 0.001) {
				t.Errorf("error : coeff %d, want %f have %f", i, valuesWant[i], valuesTest[i])
				break
			}
		}
	})
}

func test_Conjugate(params *CKKSTESTPARAMS, t *testing.T) {

	values, plaintext, ciphertext, err := new_test_vectors_reals(params, -15, 15)
	if err != nil {
		log.Fatal(err)
	}

	valuesWant := make([]complex128, len(values))
	for i := range values {
		valuesWant[i] = complex(real(values[i]), -imag(values[i]))
	}

	t.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/Conjugate(Ct)", params.ckkscontext.logN,
		params.ckkscontext.logQ,
		params.ckkscontext.levels), func(t *testing.T) {

		if err := params.evaluator.Conjugate(ciphertext, params.rotkey, ciphertext); err != nil {
			log.Fatal(err)
		}

		if err := verify_test_vectors(params, valuesWant, ciphertext, t); err != nil {
			log.Fatal(err)
		}
	})

	t.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/Conjugate(Plain)", params.ckkscontext.logN,
		params.ckkscontext.logQ,
		params.ckkscontext.levels), func(t *testing.T) {

		if err := params.evaluator.Conjugate(plaintext, params.rotkey, plaintext); err != nil {
			log.Fatal(err)
		}

		if err := verify_test_vectors(params, valuesWant, plaintext, t); err != nil {
			log.Fatal(err)
		}
	})
}

func test_RotColumns(params *CKKSTESTPARAMS, t *testing.T) {

	mask := params.ckkscontext.slots - 1

	values, plaintext, ciphertext, err := new_test_vectors_reals(params, 0.1, 1)
	if err != nil {
		log.Fatal(err)
	}

	valuesWant := make([]complex128, params.ckkscontext.slots)

	plaintextTest := params.ckkscontext.NewPlaintext(ciphertext.Level(), ciphertext.Scale())
	ciphertextTest := params.ckkscontext.NewCiphertext(1, ciphertext.Level(), ciphertext.Scale())
	ciphertextTest.SetScale(ciphertext.Scale())

	t.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/RotColumnsPow2(Ct)", params.ckkscontext.logN,
		params.ckkscontext.logQ,
		params.ckkscontext.levels), func(t *testing.T) {

		for n := uint64(1); n < params.ckkscontext.slots; n <<= 1 {

			// Applies the column rotation to the values
			for i := uint64(0); i < params.ckkscontext.slots; i++ {
				valuesWant[i] = values[(i+n)&mask]
			}

			if err := params.evaluator.RotateColumns(ciphertext, n, params.rotkey, ciphertextTest); err != nil {
				log.Fatal(err)
			}

			if err := verify_test_vectors(params, valuesWant, ciphertextTest, t); err != nil {
				log.Fatal(err)
			}
		}
	})

	t.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/RotColumnsRandom(Ct)", params.ckkscontext.logN,
		params.ckkscontext.logQ,
		params.ckkscontext.levels), func(t *testing.T) {

		for n := uint64(1); n < params.ckkscontext.slots; n <<= 1 {

			rand := ring.RandUniform(params.ckkscontext.slots, mask)

			// Applies the column rotation to the values
			for i := uint64(0); i < params.ckkscontext.slots; i++ {
				valuesWant[i] = values[(i+rand)&mask]
			}

			if err := params.evaluator.RotateColumns(ciphertext, rand, params.rotkey, ciphertextTest); err != nil {
				log.Fatal(err)
			}

			if err := verify_test_vectors(params, valuesWant, ciphertextTest, t); err != nil {
				log.Fatal(err)
			}
		}
	})

	t.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/RotColumnsPow2(Plain)", params.ckkscontext.logN,
		params.ckkscontext.logQ,
		params.ckkscontext.levels), func(t *testing.T) {

		for n := uint64(1); n < params.ckkscontext.slots; n <<= 1 {

			// Applies the column rotation to the values
			for i := uint64(0); i < params.ckkscontext.slots; i++ {
				valuesWant[i] = values[(i+n)&mask]
			}

			if err := params.evaluator.RotateColumns(plaintext, n, params.rotkey, plaintextTest); err != nil {
				log.Fatal(err)
			}

			if err := verify_test_vectors(params, valuesWant, plaintextTest, t); err != nil {
				log.Fatal(err)
			}
		}
	})

	t.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/RotColumnsRandom(Plain)", params.ckkscontext.logN,
		params.ckkscontext.logQ,
		params.ckkscontext.levels), func(t *testing.T) {

		for n := uint64(1); n < params.ckkscontext.slots; n <<= 1 {

			rand := ring.RandUniform(params.ckkscontext.slots, params.ckkscontext.slots-1)

			// Applies the column rotation to the values
			for i := uint64(0); i < params.ckkscontext.slots; i++ {
				valuesWant[i] = values[(i+rand)&mask]
			}

			if err := params.evaluator.RotateColumns(plaintext, rand, params.rotkey, plaintextTest); err != nil {
				log.Fatal(err)
			}

			if err := verify_test_vectors(params, valuesWant, plaintextTest, t); err != nil {
				log.Fatal(err)
			}
		}
	})
}

func test_MarshalCiphertext(params *CKKSTESTPARAMS, t *testing.T) {

	t.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/MarshalCiphertext", params.ckkscontext.logN,
		params.ckkscontext.logQ,
		params.ckkscontext.levels), func(t *testing.T) {

		_, _, ciphertext, err := new_test_vectors_reals(params, 0.1, 1)
		if err != nil {
			log.Fatal(err)
		}

		b, err := ciphertext.MarshalBinary()
		if err != nil {
			log.Fatal(err)
		}

		newCT := params.ckkscontext.NewCiphertext(ciphertext.Degree(), ciphertext.Level(), 0)

		err = newCT.UnMarshalBinary(b)
		if err != nil {
			log.Fatal(err)
		}

		if !params.ckkscontext.keyscontext.Equal(ciphertext.Value()[0], newCT.Value()[0]) {
			t.Errorf("marshal binary ciphertext[0]")
		}

		if !params.ckkscontext.keyscontext.Equal(ciphertext.Value()[1], newCT.Value()[1]) {
			t.Errorf("marshal binary ciphertext[1]")
		}
	})

}

func test_MarshalSecretKey(params *CKKSTESTPARAMS, t *testing.T) {

	t.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/MarshalSK", params.ckkscontext.logN,
		params.ckkscontext.logQ,
		params.ckkscontext.levels), func(t *testing.T) {

		b, _ := params.sk.MarshalBinary(params.ckkscontext)

		newSK := params.kgen.NewSecretKeyEmpty()
		if err := newSK.UnMarshalBinary(b, params.ckkscontext); err != nil {
			log.Fatal(err)
		}

		if !params.ckkscontext.keyscontext.Equal(params.sk.sk, newSK.sk) {
			t.Errorf("marshal binary sk")
		}
	})
}

func test_MarshalPublicKey(params *CKKSTESTPARAMS, t *testing.T) {

	t.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/MarshalPK", params.ckkscontext.logN,
		params.ckkscontext.logQ,
		params.ckkscontext.levels), func(t *testing.T) {

		b, _ := params.pk.MarshalBinary()

		newPK := params.kgen.NewPublicKeyEmpty()
		if err := newPK.UnMarshalBinary(b); err != nil {
			log.Fatal(err)
		}

		if !params.ckkscontext.keyscontext.Equal(params.pk.pk[0], newPK.pk[0]) {
			t.Errorf("marshal binary pk[0]")
		}

		if !params.ckkscontext.keyscontext.Equal(params.pk.pk[1], newPK.pk[1]) {
			t.Errorf("marshal binary pk[1]")
		}
	})
}

func test_MarshalEvaluationKey(params *CKKSTESTPARAMS, t *testing.T) {

	t.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/MarshalRlk", params.ckkscontext.logN,
		params.ckkscontext.logQ,
		params.ckkscontext.levels), func(t *testing.T) {

		b, _ := params.rlk.MarshalBinary()

		newRlk := params.kgen.NewRelinKeyEmpty(params.rlk.evakey.bitDecomp)
		if err := newRlk.UnMarshalBinary(b); err != nil {
			log.Fatal(err)
		}

		for x := range newRlk.evakey.evakey {
			for j := range newRlk.evakey.evakey[x] {
				if !params.ckkscontext.keyscontext.Equal(newRlk.evakey.evakey[x][j][0], params.rlk.evakey.evakey[x][j][0]) {
					t.Errorf("marshal binary rlk[0]")
					break
				}

				if !params.ckkscontext.keyscontext.Equal(newRlk.evakey.evakey[x][j][1], params.rlk.evakey.evakey[x][j][1]) {
					t.Errorf("marshal binary rlk[1]")
					break
				}
			}
		}
	})
}

func test_MarshalSwitchingKey(params *CKKSTESTPARAMS, t *testing.T) {

	t.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/MarshalSwitchingKey", params.ckkscontext.logN,
		params.ckkscontext.logQ,
		params.ckkscontext.levels), func(t *testing.T) {

		bitDecomp := uint64(15)

		s1, _ := params.kgen.NewSecretKey(1.0 / 3)

		switchkey, _ := params.kgen.NewSwitchingKey(params.sk, s1, bitDecomp)

		b, _ := switchkey.MarshalBinary()

		newSwKey := params.kgen.NewSwitchingKeyEmpty(bitDecomp)
		if err := newSwKey.UnMarshalBinary(b); err != nil {
			log.Fatal(err)
		}

		for x := range newSwKey.evakey {
			for j := range newSwKey.evakey[x] {
				if !params.ckkscontext.keyscontext.Equal(newSwKey.evakey[x][j][0], switchkey.evakey[x][j][0]) {
					t.Errorf("marshal binary switchingkey[0]")
					break
				}

				if !params.ckkscontext.keyscontext.Equal(newSwKey.evakey[x][j][1], switchkey.evakey[x][j][1]) {
					t.Errorf("marshal binary switchingkey[1]")
					break
				}
			}
		}

	})
}

func test_MarshalRotationKey(params *CKKSTESTPARAMS, t *testing.T) {

	t.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/MarshalRotKey", params.ckkscontext.logN,
		params.ckkscontext.logQ,
		params.ckkscontext.levels), func(t *testing.T) {

		b, _ := params.rotkey.MarshalBinary()

		newRotKey := params.kgen.NewRotationKeysEmpty()
		if err := newRotKey.UnMarshalBinary(b); err != nil {
			log.Fatal(err)
		}

		for i := range params.rotkey.evakey_rot_col_L {
			for x := range params.rotkey.evakey_rot_col_L[i].evakey {
				for j := range params.rotkey.evakey_rot_col_L[i].evakey[x] {
					if !params.ckkscontext.keyscontext.Equal(params.rotkey.evakey_rot_col_L[i].evakey[x][j][0], newRotKey.evakey_rot_col_L[i].evakey[x][j][0]) {
						t.Errorf("marshal binary rotkeyLeft[0]")
						break
					}

					if !params.ckkscontext.keyscontext.Equal(params.rotkey.evakey_rot_col_L[i].evakey[x][j][1], newRotKey.evakey_rot_col_L[i].evakey[x][j][1]) {
						t.Errorf("marshal binary rotkeyLeft[1]")
						break
					}
				}
			}
		}

		for i := range params.rotkey.evakey_rot_col_R {
			for x := range params.rotkey.evakey_rot_col_R[i].evakey {
				for j := range params.rotkey.evakey_rot_col_R[i].evakey[x] {
					if !params.ckkscontext.keyscontext.Equal(params.rotkey.evakey_rot_col_R[i].evakey[x][j][0], newRotKey.evakey_rot_col_R[i].evakey[x][j][0]) {
						t.Errorf("marshal binary rotkeyRight[0]")
						break
					}

					if !params.ckkscontext.keyscontext.Equal(params.rotkey.evakey_rot_col_R[i].evakey[x][j][1], newRotKey.evakey_rot_col_R[i].evakey[x][j][1]) {
						t.Errorf("marshal binary rotkeyRight[1]")
						break
					}
				}
			}
		}

		if params.rotkey.evakey_rot_row != nil {
			for x := range params.rotkey.evakey_rot_row.evakey {
				for j := range params.rotkey.evakey_rot_row.evakey[x] {
					if !params.ckkscontext.keyscontext.Equal(params.rotkey.evakey_rot_row.evakey[x][j][0], newRotKey.evakey_rot_row.evakey[x][j][0]) {
						t.Errorf("marshal binary rotkeyConjugate[0]")
						break
					}

					if !params.ckkscontext.keyscontext.Equal(params.rotkey.evakey_rot_row.evakey[x][j][1], newRotKey.evakey_rot_row.evakey[x][j][1]) {
						t.Errorf("marshal binary rotkeyConjugate[1]")
						break
					}
				}
			}
		}
	})
}
