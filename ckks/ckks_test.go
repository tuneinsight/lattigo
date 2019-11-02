package ckks

import (
	"fmt"
	"github.com/ldsec/lattigo/ring"
	"log"
	"math"
	"math/cmplx"
	"math/rand"
	"sort"
	"testing"
	"time"
)

type CKKSTESTPARAMS struct {
	medianprec  float64
	ckkscontext *CkksContext
	slots       uint64
	encoder     *Encoder
	levels      uint64
	scale       float64
	kgen        *KeyGenerator
	sk          *SecretKey
	pk          *PublicKey
	rlk         *EvaluationKey
	rotkey      *RotationKey
	encryptorPk *Encryptor
	encryptorSk *Encryptor
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

	medianprec := float64(20) // target median precision in log2 among all the coeffs, determines the success/failure of a test

	ctsDepth := uint64(3) // maximum depth for the bootstrapp LT
	stcDepth := uint64(2)
	showbootintermediateresults := true

	params := Parameters{16, []uint8{55, 40, 40, 40, 40, 40, 40, 40, 40, 45, 45, 45, 55, 55, 55, 55, 55, 55, 55, 55, 55, 45, 45, 45, 45}, []uint8{55, 55, 55, 55, 55}, 1 << 40, 3.2}

	ckksTest := new(CKKSTESTPARAMS)

	ckksTest.medianprec = medianprec

	ckksTest.slots = 1 << 10

	ckksTest.levels = uint64(len(params.Modulichain))
	ckksTest.scale = params.Scale

	if ckksTest.ckkscontext, err = NewCkksContext(&params); err != nil {
		t.Error(err)
	}

	log.Printf("Generated CkksContext for logN=%d/logQ=%d/levels=%d/a=%d/b=%d/sigma=%f",
		ckksTest.ckkscontext.LogN(),
		ckksTest.ckkscontext.LogQ(),
		ckksTest.ckkscontext.Levels(),
		ckksTest.ckkscontext.alpha,
		ckksTest.ckkscontext.beta,
		ckksTest.ckkscontext.Sigma())

	for i, qi := range ckksTest.ckkscontext.ContextKeys().Modulus {
		fmt.Println(i, qi)
	}

	ckksTest.kgen = ckksTest.ckkscontext.NewKeyGenerator()

	ckksTest.sk, ckksTest.pk = ckksTest.kgen.NewKeyPairSparse(64)

	ckksTest.encoder = ckksTest.ckkscontext.NewEncoder()

	if ckksTest.encryptorPk, err = ckksTest.ckkscontext.NewEncryptorFromPk(ckksTest.pk); err != nil {
		t.Error(err)
	}

	if ckksTest.encryptorSk, err = ckksTest.ckkscontext.NewEncryptorFromSk(ckksTest.sk); err != nil {
		t.Error(err)
	}

	if ckksTest.decryptor, err = ckksTest.ckkscontext.NewDecryptor(ckksTest.sk); err != nil {
		t.Error(err)
	}

	ckksTest.evaluator = ckksTest.ckkscontext.NewEvaluator()
/*
	log.Printf("Generating relinearization keys")
	ckksTest.rlk = ckksTest.kgen.NewRelinKey(ckksTest.sk)

	log.Printf("Generating rotation keys for conjugate and powers of 2")
	ckksTest.rotkey = ckksTest.kgen.NewRotationKeysPow2(ckksTest.sk, true)

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

	if len(params.Modulichain) > 9 {
		test_sin2pi2pi(ckksTest, t)
		test_Functions(ckksTest, t)
	}
	test_SwitchKeys(ckksTest, t)
	test_Conjugate(ckksTest, t)
	test_RotColumns(ckksTest, t)
*/
	if len(params.Modulichain) > 10 {
		test_Bootstrapp(ckksTest, ctsDepth, stcDepth, showbootintermediateresults, t)
	}

}

func test_sin2pi2pi(params *CKKSTESTPARAMS, t *testing.T) {

	t.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/a=%d/b=%d/sin(2*pi*x)/(2*pi) [-12, 12] deg128", params.ckkscontext.logN,
		params.ckkscontext.logQ,
		params.ckkscontext.levels,
		params.ckkscontext.alpha,
		params.ckkscontext.beta), func(t *testing.T) {

		values, _, ciphertext1, err := new_test_vectors_reals(params, -0.01, 0.01)
		if err != nil {
			t.Error(err)
		}

		valuesWant := make([]complex128, params.slots)

		for i := 0; i < len(valuesWant); i++ {
			valuesWant[i] = sin2pi2pi(values[i])
		}

		cheby := Approximate(sin2pi2pi, -12, 12, 118)

		if ciphertext1, err = params.evaluator.EvaluateCheby(ciphertext1, cheby, params.rlk); err != nil {
			t.Error(err)
		}

		if err := verify_test_vectors(params, valuesWant, ciphertext1.Element(), t); err != nil {
			t.Error(err)
		}
	})
}

func new_test_vectors(params *CKKSTESTPARAMS, a, b float64) (values []complex128, plaintext *Plaintext, ciphertext *Ciphertext, err error) {

	values = make([]complex128, params.slots)

	for i := uint64(0); i < params.slots; i++ {
		values[i] = randomComplex(a, b)
	}

	values[0] = complex(0.607538, 0.555668)

	plaintext = params.ckkscontext.NewPlaintext(params.ckkscontext.Levels()-1, params.ckkscontext.Scale())

	if err = params.encoder.Encode(plaintext, values, params.slots); err != nil {
		return nil, nil, nil, err
	}

	ciphertext, err = params.encryptorPk.EncryptNew(plaintext)
	if err != nil {
		return nil, nil, nil, err
	}

	return values, plaintext, ciphertext, nil
}

func new_test_vectors_reals(params *CKKSTESTPARAMS, a, b float64) (values []complex128, plaintext *Plaintext, ciphertext *Ciphertext, err error) {

	values = make([]complex128, params.slots)

	for i := uint64(0); i < params.slots; i++ {
		values[i] = complex(randomFloat(a, b), 0)
	}

	values[0] = complex(0.607538, 0)

	plaintext = params.ckkscontext.NewPlaintext(params.ckkscontext.Levels()-1, params.ckkscontext.Scale())

	if err = params.encoder.Encode(plaintext, values, params.slots); err != nil {
		return nil, nil, nil, err
	}

	ciphertext, err = params.encryptorPk.EncryptNew(plaintext)
	if err != nil {
		return nil, nil, nil, err
	}

	return values, plaintext, ciphertext, nil
}

func verify_test_vectors(params *CKKSTESTPARAMS, valuesWant []complex128, element *ckksElement, t *testing.T) (err error) {

	var plaintextTest *Plaintext
	var valuesTest []complex128

	if element.Degree() == 0 {

		plaintextTest = element.Plaintext()

	} else {

		plaintextTest = params.decryptor.DecryptNew(element.Ciphertext())
	}

	valuesTest = params.encoder.Decode(plaintextTest, params.slots)

	var deltaReal, deltaImag float64

	var minprec, maxprec, meanprec, medianprec complex128

	diff := make([]complex128, params.slots)

	minprec = complex(0, 0)
	maxprec = complex(1, 1)

	meanprec = complex(0, 0)

	distrib_real := make(map[uint64]uint64)
	distrib_imag := make(map[uint64]uint64)

	for i := range valuesWant {

		// Test the ratio for big values (> 1) and difference for small values (< 1)

		deltaReal = math.Abs(real(valuesTest[i]) - real(valuesWant[i]))
		deltaImag = math.Abs(imag(valuesTest[i]) - imag(valuesWant[i]))

		diff[i] += complex(deltaReal, 0)
		diff[i] += complex(0, deltaImag)

		meanprec += diff[i]

		if real(diff[i]) > real(minprec) {
			minprec = complex(real(diff[i]), 0)
		}

		if imag(diff[i]) > imag(minprec) {
			minprec = complex(real(minprec), imag(diff[i]))
		}

		if real(diff[i]) < real(maxprec) {
			maxprec = complex(real(diff[i]), 0)
		}

		if imag(diff[i]) < imag(maxprec) {
			maxprec = complex(real(maxprec), imag(diff[i]))
		}

		distrib_real[uint64(math.Floor(math.Log2(1/real(diff[i]))))] += 1
		distrib_imag[uint64(math.Floor(math.Log2(1/imag(diff[i]))))] += 1
	}

	meanprec /= complex(float64(params.slots), 0)
	medianprec = calcmedian(diff)

	t.Log()
	t.Logf("Minimum precision : (%.2f, %.2f) bits \n", math.Log2(1/real(minprec)), math.Log2(1/imag(minprec)))
	t.Logf("Maximum precision : (%.2f, %.2f) bits \n", math.Log2(1/real(maxprec)), math.Log2(1/imag(maxprec)))
	t.Logf("Mean    precision : (%.2f, %.2f) bits \n", math.Log2(1/real(meanprec)), math.Log2(1/imag(meanprec)))
	t.Logf("Median  precision : (%.2f, %.2f) bits \n", math.Log2(1/real(medianprec)), math.Log2(1/imag(medianprec)))
	t.Log()

	if math.Log2(1/real(medianprec)) < params.medianprec || math.Log2(1/imag(medianprec)) < params.medianprec {
		t.Errorf("Mean precision error : target (%.2f, %.2f) > result (%.2f, %.2f)", params.medianprec, params.medianprec, math.Log2(1/real(medianprec)), math.Log2(1/imag(medianprec)))
	}

	/*
		fmt.Println()
		fmt.Println("Distribution of the precision :")
		keys_real := []int{}
		keys_imag := []int{}
		for i := range distrib_real {
			keys_real = append(keys_real, int(i))
		}
		for i := range distrib_imag {
			keys_imag = append(keys_imag, int(i))
		}
		sort.Ints(keys_real)
		sort.Ints(keys_imag)
		for _, i := range keys_real {
			fmt.Printf("bits %d : %.2f %% \n", i, (float64(distrib_real[uint64(i)])/float64(params.slots))*100)
		}
		for _, i := range keys_imag {
			fmt.Printf("bits %d : %.2f %% \n", i, (float64(distrib_imag[uint64(i)])/float64(params.slots))*100)
		}
	*/
	return nil
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

func test_Encoder(params *CKKSTESTPARAMS, t *testing.T) {

	t.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/a=%d/b=%d/EncodeDecodeComplex128", params.ckkscontext.logN,
		params.ckkscontext.logQ,
		params.ckkscontext.levels, params.ckkscontext.alpha, params.ckkscontext.beta), func(t *testing.T) {

		valuesWant := make([]complex128, params.slots)

		for i := uint64(0); i < params.slots; i++ {
			valuesWant[i] = randomComplex(0, 5)
		}

		plaintext := params.ckkscontext.NewPlaintext(params.ckkscontext.Levels()-1, params.ckkscontext.Scale())

		if err := params.encoder.Encode(plaintext, valuesWant, params.slots); err != nil {
			t.Error(err)
		}

		if err := verify_test_vectors(params, valuesWant, plaintext.Element(), t); err != nil {
			t.Error(err)
		}
	})
}

func test_EncryptDecrypt(params *CKKSTESTPARAMS, t *testing.T) {

	t.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/a=%d/b=%d/EncryptFromPk", params.ckkscontext.logN,
		params.ckkscontext.logQ,
		params.ckkscontext.levels, params.ckkscontext.alpha, params.ckkscontext.beta), func(t *testing.T) {
		var err error

		valuesWant := make([]complex128, params.slots)

		for i := uint64(0); i < params.slots; i++ {
			valuesWant[i] = randomComplex(0, 5)
		}

		plaintext := params.ckkscontext.NewPlaintext(params.ckkscontext.Levels()-1, params.ckkscontext.Scale())

		if err = params.encoder.Encode(plaintext, valuesWant, params.slots); err != nil {
			t.Error(err)
		}

		ciphertext, err := params.encryptorPk.EncryptNew(plaintext)
		if err != nil {
			t.Error(err)
		}

		if err := verify_test_vectors(params, valuesWant, ciphertext.Element(), t); err != nil {
			t.Error(err)
		}
	})

	t.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/a=%d/b=%d/EncryptFromSk", params.ckkscontext.logN,
		params.ckkscontext.logQ,
		params.ckkscontext.levels, params.ckkscontext.alpha, params.ckkscontext.beta), func(t *testing.T) {
		var err error

		valuesWant := make([]complex128, params.slots)

		for i := uint64(0); i < params.slots; i++ {
			valuesWant[i] = randomComplex(0, 5)
		}

		plaintext := params.ckkscontext.NewPlaintext(params.ckkscontext.Levels()-1, params.ckkscontext.Scale())

		if err = params.encoder.Encode(plaintext, valuesWant, params.slots); err != nil {
			t.Error(err)
		}

		ciphertext, err := params.encryptorSk.EncryptNew(plaintext)
		if err != nil {
			t.Error(err)
		}

		if err := verify_test_vectors(params, valuesWant, ciphertext.Element(), t); err != nil {
			t.Error(err)
		}
	})
}

func test_Add(params *CKKSTESTPARAMS, t *testing.T) {

	t.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/a=%d/b=%d/AddCtCtInPlace", params.ckkscontext.logN,
		params.ckkscontext.logQ,
		params.ckkscontext.levels, params.ckkscontext.alpha, params.ckkscontext.beta), func(t *testing.T) {

		values1, _, ciphertext1, err := new_test_vectors(params, -1, 1)
		if err != nil {
			t.Error(err)
		}

		values2, _, ciphertext2, err := new_test_vectors(params, -1, 1)
		if err != nil {
			t.Error(err)
		}

		valuesWant := make([]complex128, params.slots)
		for i := 0; i < len(valuesWant); i++ {
			valuesWant[i] = values1[i] + values2[i]
		}

		if err := params.evaluator.Add(ciphertext1, ciphertext2, ciphertext1); err != nil {
			t.Error(err)
		}

		if err := verify_test_vectors(params, valuesWant, ciphertext1.Element(), t); err != nil {
			t.Error(err)
		}
	})

	t.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/a=%d/b=%d/AddCtCt", params.ckkscontext.logN,
		params.ckkscontext.logQ,
		params.ckkscontext.levels, params.ckkscontext.alpha, params.ckkscontext.beta), func(t *testing.T) {

		values1, _, ciphertext1, err := new_test_vectors(params, -1, 1)
		if err != nil {
			t.Error(err)
		}

		values2, _, ciphertext2, err := new_test_vectors(params, -1, 1)
		if err != nil {
			t.Error(err)
		}

		receiver := params.ckkscontext.NewCiphertext(1, ciphertext1.Level(), ciphertext1.Scale())

		valuesWant := make([]complex128, params.slots)
		for i := 0; i < len(valuesWant); i++ {
			valuesWant[i] = values1[i] + values2[i]
		}

		if err := params.evaluator.Add(ciphertext1, ciphertext2, receiver); err != nil {
			t.Error(err)
		}

		if err := verify_test_vectors(params, valuesWant, receiver.Element(), t); err != nil {
			t.Error(err)
		}
	})

	t.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/a=%d/b=%d/Add(Ct,Plain)", params.ckkscontext.logN,
		params.ckkscontext.logQ,
		params.ckkscontext.levels, params.ckkscontext.alpha, params.ckkscontext.beta), func(t *testing.T) {

		values1, _, ciphertext1, err := new_test_vectors(params, -1, 1)
		if err != nil {
			t.Error(err)
		}

		values2, plaintext2, _, err := new_test_vectors(params, -1, 1)
		if err != nil {
			t.Error(err)
		}

		valuesWant := make([]complex128, params.slots)
		for i := 0; i < len(valuesWant); i++ {
			valuesWant[i] = values1[i] + values2[i]
		}

		if err := params.evaluator.Add(ciphertext1, plaintext2, ciphertext1); err != nil {
			t.Error(err)
		}

		if err := verify_test_vectors(params, valuesWant, ciphertext1.Element(), t); err != nil {
			t.Error(err)
		}
	})

	t.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/a=%d/b=%d/Add(Plain,Ct)", params.ckkscontext.logN,
		params.ckkscontext.logQ,
		params.ckkscontext.levels, params.ckkscontext.alpha, params.ckkscontext.beta), func(t *testing.T) {

		values1, plaintext1, _, err := new_test_vectors(params, -1, 1)
		if err != nil {
			t.Error(err)
		}

		values2, _, ciphertext2, err := new_test_vectors(params, -1, 1)
		if err != nil {
			t.Error(err)
		}

		valuesWant := make([]complex128, params.slots)
		for i := 0; i < len(valuesWant); i++ {
			valuesWant[i] = values1[i] + values2[i]
		}

		if err := params.evaluator.Add(plaintext1, ciphertext2, ciphertext2); err != nil {
			t.Error(err)
		}

		if err := verify_test_vectors(params, valuesWant, ciphertext2.Element(), t); err != nil {
			t.Error(err)
		}
	})
}

func test_Sub(params *CKKSTESTPARAMS, t *testing.T) {

	t.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/a=%d/b=%d/Sub(Ct,Ct)", params.ckkscontext.logN,
		params.ckkscontext.logQ,
		params.ckkscontext.levels, params.ckkscontext.alpha, params.ckkscontext.beta), func(t *testing.T) {

		values1, _, ciphertext1, err := new_test_vectors(params, -1, 1)
		if err != nil {
			t.Error(err)
		}

		values2, _, ciphertext2, err := new_test_vectors(params, -1, 1)
		if err != nil {
			t.Error(err)
		}

		valuesWant := make([]complex128, params.slots)
		for i := 0; i < len(valuesWant); i++ {
			valuesWant[i] = values1[i] - values2[i]
		}

		if err := params.evaluator.Sub(ciphertext1, ciphertext2, ciphertext1); err != nil {
			t.Error(err)
		}

		if err := verify_test_vectors(params, valuesWant, ciphertext1.Element(), t); err != nil {
			t.Error(err)
		}
	})

	t.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/a=%d/b=%d/Sub(Ct,Plain)", params.ckkscontext.logN,
		params.ckkscontext.logQ,
		params.ckkscontext.levels, params.ckkscontext.alpha, params.ckkscontext.beta), func(t *testing.T) {

		values1, _, ciphertext1, err := new_test_vectors(params, -1, 1)
		if err != nil {
			t.Error(err)
		}

		values2, plaintext2, _, err := new_test_vectors(params, -1, 1)
		if err != nil {
			t.Error(err)
		}

		valuesWant := make([]complex128, params.slots)
		for i := 0; i < len(valuesWant); i++ {
			valuesWant[i] = values1[i] - values2[i]
		}

		if err := params.evaluator.Sub(ciphertext1, plaintext2, ciphertext1); err != nil {
			t.Error(err)
		}

		if err := verify_test_vectors(params, valuesWant, ciphertext1.Element(), t); err != nil {
			t.Error(err)
		}
	})

	t.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/a=%d/b=%d/Sub(Plain,Ct)", params.ckkscontext.logN,
		params.ckkscontext.logQ,
		params.ckkscontext.levels, params.ckkscontext.alpha, params.ckkscontext.beta), func(t *testing.T) {

		values1, plaintext1, _, err := new_test_vectors(params, -1, 1)
		if err != nil {
			t.Error(err)
		}

		values2, _, ciphertext2, err := new_test_vectors(params, -1, 1)
		if err != nil {
			t.Error(err)
		}

		valuesWant := make([]complex128, params.slots)
		for i := 0; i < len(valuesWant); i++ {
			valuesWant[i] = values1[i] - values2[i]
		}

		if err := params.evaluator.Sub(plaintext1, ciphertext2, ciphertext2); err != nil {
			t.Error(err)
		}

		if err := verify_test_vectors(params, valuesWant, ciphertext2.Element(), t); err != nil {
			t.Error(err)
		}
	})
}

func test_AddConst(params *CKKSTESTPARAMS, t *testing.T) {

	t.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/a=%d/b=%d/AddCmplx(Ct,complex128)", params.ckkscontext.logN,
		params.ckkscontext.logQ,
		params.ckkscontext.levels, params.ckkscontext.alpha, params.ckkscontext.beta), func(t *testing.T) {

		values1, _, ciphertext1, err := new_test_vectors(params, -1, 1)
		if err != nil {
			t.Error(err)
		}

		constant := complex(3.1415, -1.4142)

		valuesWant := make([]complex128, params.slots)
		for i := 0; i < len(valuesWant); i++ {
			valuesWant[i] = values1[i] + constant
		}

		if err := params.evaluator.AddConst(ciphertext1, constant, ciphertext1); err != nil {
			t.Error(err)
		}

		if err := verify_test_vectors(params, valuesWant, ciphertext1.Element(), t); err != nil {
			t.Error(err)
		}
	})
}

func test_MulConst(params *CKKSTESTPARAMS, t *testing.T) {

	t.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/a=%d/b=%d/MultCmplx(Ct,complex128)", params.ckkscontext.logN,
		params.ckkscontext.logQ,
		params.ckkscontext.levels, params.ckkscontext.alpha, params.ckkscontext.beta), func(t *testing.T) {

		values1, _, ciphertext1, err := new_test_vectors(params, -5, 5)
		if err != nil {
			t.Error(err)
		}

		constant := complex(1.4142, -3.1415)

		valuesWant := make([]complex128, params.slots)
		for i := 0; i < len(valuesWant); i++ {
			valuesWant[i] = values1[i] * (1 / constant)
		}

		if err = params.evaluator.MultConst(ciphertext1, 1/constant, ciphertext1); err != nil {
			t.Error(err)
		}

		if err := verify_test_vectors(params, valuesWant, ciphertext1.Element(), t); err != nil {
			t.Error(err)
		}
	})

	t.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/a=%d/b=%d/MultByi", params.ckkscontext.logN,
		params.ckkscontext.logQ,
		params.ckkscontext.levels, params.ckkscontext.alpha, params.ckkscontext.beta), func(t *testing.T) {

		values, _, ciphertext1, err := new_test_vectors(params, -1, 1)
		if err != nil {
			t.Error(err)
		}

		valuesWant := make([]complex128, params.slots)

		for i := 0; i < len(valuesWant); i++ {
			valuesWant[i] = values[i] * complex(0, 1)
		}

		if err = params.evaluator.MultByi(ciphertext1, ciphertext1); err != nil {
			t.Error(err)
		}

		if err := verify_test_vectors(params, valuesWant, ciphertext1.Element(), t); err != nil {
			t.Error(err)
		}
	})

	t.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/a=%d/b=%d/DivByi", params.ckkscontext.logN,
		params.ckkscontext.logQ,
		params.ckkscontext.levels, params.ckkscontext.alpha, params.ckkscontext.beta), func(t *testing.T) {

		values, _, ciphertext1, err := new_test_vectors(params, -1, 1)
		if err != nil {
			t.Error(err)
		}

		valuesWant := make([]complex128, params.slots)

		for i := 0; i < len(valuesWant); i++ {
			valuesWant[i] = values[i] / complex(0, 1)
		}

		if err = params.evaluator.DivByi(ciphertext1, ciphertext1); err != nil {
			t.Error(err)
		}

		if err := verify_test_vectors(params, valuesWant, ciphertext1.Element(), t); err != nil {
			t.Error(err)
		}
	})
}

func test_MultByConstAndAdd(params *CKKSTESTPARAMS, t *testing.T) {

	t.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/a=%d/b=%d/MultByCmplxAndAdd(Ct0, complex128, Ct1)", params.ckkscontext.logN,
		params.ckkscontext.logQ,
		params.ckkscontext.levels, params.ckkscontext.alpha, params.ckkscontext.beta), func(t *testing.T) {

		values1, _, ciphertext1, err := new_test_vectors(params, -1, 1)
		if err != nil {
			t.Error(err)
		}

		values2, _, ciphertext2, err := new_test_vectors(params, -1, 1)
		if err != nil {
			t.Error(err)
		}

		constant := complex(3.1415, -1.4142)

		valuesWant := make([]complex128, params.slots)
		for i := 0; i < len(valuesWant); i++ {
			values2[i] += (values1[i] * constant) + (values1[i] * constant)
		}

		if err = params.evaluator.MultByConstAndAdd(ciphertext1, constant, ciphertext2); err != nil {
			t.Error(err)
		}

		params.evaluator.Rescale(ciphertext2, params.ckkscontext.scale, ciphertext2)

		if err = params.evaluator.MultByConstAndAdd(ciphertext1, constant, ciphertext2); err != nil {
			t.Error(err)
		}

		if err := verify_test_vectors(params, values2, ciphertext2.Element(), t); err != nil {
			t.Error(err)
		}
	})
}

func test_ComplexOperations(params *CKKSTESTPARAMS, t *testing.T) {

	t.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/a=%d/b=%d/ExtractImag", params.ckkscontext.logN,
		params.ckkscontext.logQ,
		params.ckkscontext.levels, params.ckkscontext.alpha, params.ckkscontext.beta), func(t *testing.T) {

		values, _, ciphertext, err := new_test_vectors(params, -1, 1)
		if err != nil {
			t.Error(err)
		}

		for i := 0; i < len(values); i++ {
			values[i] = complex(imag(values[i]), 0)
		}

		if err = params.evaluator.ExtractImag(ciphertext, params.rotkey, ciphertext); err != nil {
			t.Error(err)
		}

		if err := verify_test_vectors(params, values, ciphertext.Element(), t); err != nil {
			t.Error(err)
		}
	})

	t.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/a=%d/b=%d/SwapRealImag", params.ckkscontext.logN,
		params.ckkscontext.logQ,
		params.ckkscontext.levels, params.ckkscontext.alpha, params.ckkscontext.beta), func(t *testing.T) {

		values, _, ciphertext, err := new_test_vectors(params, -1, 1)
		if err != nil {
			t.Error(err)
		}

		for i := 0; i < len(values); i++ {
			values[i] = complex(imag(values[i]), real(values[i]))
		}

		if err = params.evaluator.SwapRealImag(ciphertext, params.rotkey, ciphertext); err != nil {
			t.Error(err)
		}

		if err := verify_test_vectors(params, values, ciphertext.Element(), t); err != nil {
			t.Error(err)
		}
	})

	t.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/a=%d/b=%d/RemoveReal", params.ckkscontext.logN,
		params.ckkscontext.logQ,
		params.ckkscontext.levels, params.ckkscontext.alpha, params.ckkscontext.beta), func(t *testing.T) {

		values, _, ciphertext, err := new_test_vectors(params, -1, 1)
		if err != nil {
			t.Error(err)
		}

		for i := 0; i < len(values); i++ {
			values[i] = complex(0, imag(values[i]))
		}

		if err = params.evaluator.RemoveReal(ciphertext, params.rotkey, ciphertext); err != nil {
			t.Error(err)
		}

		if err := verify_test_vectors(params, values, ciphertext.Element(), t); err != nil {
			t.Error(err)
		}
	})

	t.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/a=%d/b=%d/RemoveImag", params.ckkscontext.logN,
		params.ckkscontext.logQ,
		params.ckkscontext.levels, params.ckkscontext.alpha, params.ckkscontext.beta), func(t *testing.T) {

		values, _, ciphertext, err := new_test_vectors(params, -1, 1)
		if err != nil {
			t.Error(err)
		}

		for i := 0; i < len(values); i++ {
			values[i] = complex(real(values[i]), 0)
		}

		if err = params.evaluator.RemoveImag(ciphertext, params.rotkey, ciphertext); err != nil {
			t.Error(err)
		}

		if err := verify_test_vectors(params, values, ciphertext.Element(), t); err != nil {
			t.Error(err)
		}
	})
}

func test_Rescaling(params *CKKSTESTPARAMS, t *testing.T) {

	t.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/a=%d/b=%d/Rescaling", params.ckkscontext.logN,
		params.ckkscontext.logQ,
		params.ckkscontext.levels, params.ckkscontext.alpha, params.ckkscontext.beta), func(t *testing.T) {

		coeffs := make([]*ring.Int, params.ckkscontext.n)
		for i := uint64(0); i < params.ckkscontext.n; i++ {
			coeffs[i] = ring.RandInt(params.ckkscontext.contextQ.ModulusBigint)
			coeffs[i].Div(coeffs[i], ring.NewUint(10))
		}

		coeffsWant := make([]*ring.Int, params.ckkscontext.contextQ.N)
		for i := range coeffs {
			coeffsWant[i] = coeffs[i].Copy()
			coeffsWant[i].DivRound(coeffsWant[i], ring.NewUint(params.ckkscontext.moduli[len(params.ckkscontext.moduli)-1]))
		}

		polTest := params.ckkscontext.contextQ.NewPoly()
		polWant := params.ckkscontext.contextQ.NewPoly()

		params.ckkscontext.contextQ.SetCoefficientsBigint(coeffs, polTest)
		params.ckkscontext.contextQ.SetCoefficientsBigint(coeffsWant, polWant)

		params.ckkscontext.contextQ.NTT(polTest, polTest)
		params.ckkscontext.contextQ.NTT(polWant, polWant)

		rescale(params.evaluator, polTest, polTest)
		state := true
		for i := uint64(0); i < params.ckkscontext.n && state; i++ {
			for j := 0; j < len(params.ckkscontext.moduli)-1 && state; j++ {
				if polWant.Coeffs[j][i] != polTest.Coeffs[j][i] {
					t.Errorf("error : coeff %v Qi%v = %s, want %v have %v", i, j, coeffs[i].String(), polWant.Coeffs[j][i], polTest.Coeffs[j][i])
					state = false
				}
			}
		}
	})

	t.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/a=%d/b=%d/RescalingMultiple", params.ckkscontext.logN,
		params.ckkscontext.logQ,
		params.ckkscontext.levels, params.ckkscontext.alpha, params.ckkscontext.beta), func(t *testing.T) {

		coeffs := make([]*ring.Int, params.ckkscontext.n)
		for i := uint64(0); i < params.ckkscontext.n; i++ {
			coeffs[i] = ring.RandInt(params.ckkscontext.contextQ.ModulusBigint)
			coeffs[i].Div(coeffs[i], ring.NewUint(10))
		}

		nbRescals := 4

		coeffsWant := make([]*ring.Int, params.ckkscontext.contextQ.N)
		for i := range coeffs {
			coeffsWant[i] = coeffs[i].Copy()
			for j := 0; j < nbRescals; j++ {
				coeffsWant[i].DivRound(coeffsWant[i], ring.NewUint(params.ckkscontext.moduli[len(params.ckkscontext.moduli)-1-j]))
			}
		}

		polTest := params.ckkscontext.contextQ.NewPoly()
		polWant := params.ckkscontext.contextQ.NewPoly()

		params.ckkscontext.contextQ.SetCoefficientsBigint(coeffs, polTest)
		params.ckkscontext.contextQ.SetCoefficientsBigint(coeffsWant, polWant)

		params.ckkscontext.contextQ.NTT(polTest, polTest)
		params.ckkscontext.contextQ.NTT(polWant, polWant)

		rescaleMany(params.evaluator, uint64(nbRescals), polTest, polTest)
		state := true
		for i := uint64(0); i < params.ckkscontext.n && state; i++ {
			for j := 0; j < len(params.ckkscontext.moduli)-nbRescals && state; j++ {
				if polWant.Coeffs[j][i] != polTest.Coeffs[j][i] {
					t.Errorf("error : coeff %v Qi%v = %s, want %v have %v", i, j, coeffs[i].String(), polWant.Coeffs[j][i], polTest.Coeffs[j][i])
					state = false
				}
			}
		}
	})
}

func test_Mul(params *CKKSTESTPARAMS, t *testing.T) {

	t.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/a=%d/b=%d/Mul(Ct,Ct)", params.ckkscontext.logN,
		params.ckkscontext.logQ,
		params.ckkscontext.levels, params.ckkscontext.alpha, params.ckkscontext.beta), func(t *testing.T) {

		values1, _, ciphertext1, err := new_test_vectors(params, -1, 1)
		if err != nil {
			t.Error(err)
		}

		values2, _, ciphertext2, err := new_test_vectors(params, -1, 1)
		if err != nil {
			t.Error(err)
		}

		// up to a level equal to 2 modulus
		valuesWant := make([]complex128, params.slots)

		for i := 0; i < len(valuesWant); i++ {
			valuesWant[i] = values1[i] * values2[i]
		}

		if err = params.evaluator.MulRelin(ciphertext1, ciphertext2, nil, ciphertext1); err != nil {
			t.Error(err)
		}

		if err := verify_test_vectors(params, valuesWant, ciphertext1.Element(), t); err != nil {
			t.Error(err)
		}
	})

	t.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/a=%d/b=%d/Relinearize", params.ckkscontext.logN,
		params.ckkscontext.logQ,
		params.ckkscontext.levels, params.ckkscontext.alpha, params.ckkscontext.beta), func(t *testing.T) {

		values1, _, ciphertext1, err := new_test_vectors(params, -1, 1)
		if err != nil {
			t.Error(err)
		}

		values2, _, ciphertext2, err := new_test_vectors(params, -1, 1)
		if err != nil {
			t.Error(err)
		}

		// up to a level equal to 2 modulus
		valuesWant := make([]complex128, params.slots)

		for i := 0; i < len(valuesWant); i++ {
			valuesWant[i] = values1[i] * values2[i]
		}

		if err = params.evaluator.MulRelin(ciphertext1, ciphertext2, nil, ciphertext1); err != nil {
			t.Error(err)
		}

		if err = params.evaluator.Relinearize(ciphertext1, params.rlk, ciphertext1); err != nil {
			t.Error(err)
		}

		if err := verify_test_vectors(params, valuesWant, ciphertext1.Element(), t); err != nil {
			t.Error(err)
		}
	})

	t.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/a=%d/b=%d/MulRelin(Ct,Ct)->Rescale", params.ckkscontext.logN,
		params.ckkscontext.logQ,
		params.ckkscontext.levels, params.ckkscontext.alpha, params.ckkscontext.beta), func(t *testing.T) {

		values1, _, ciphertext1, err := new_test_vectors(params, -1, 1)
		if err != nil {
			t.Error(err)
		}

		values2, _, ciphertext2, err := new_test_vectors(params, -1, 1)
		if err != nil {
			t.Error(err)
		}

		// up to a level equal to 2 modulus
		valuesWant := make([]complex128, params.slots)

		for i := 0; i < len(valuesWant); i++ {
			valuesWant[i] = values1[i]
		}

		for i := uint64(0); i < params.ckkscontext.Levels()-1; i++ {

			for i := 0; i < len(valuesWant); i++ {
				valuesWant[i] *= values2[i]
			}

			if err = params.evaluator.MulRelin(ciphertext1, ciphertext2, params.rlk, ciphertext1); err != nil {
				t.Error(err)
			}

			if err = params.evaluator.Rescale(ciphertext1, params.ckkscontext.scale, ciphertext1); err != nil {
				t.Error(err)
			}
		}

		if err := verify_test_vectors(params, valuesWant, ciphertext1.Element(), t); err != nil {
			t.Error(err)
		}
	})

	t.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/a=%d/b=%d/MulRelin(Ct,Ct)->RescaleMany", params.ckkscontext.logN,
		params.ckkscontext.logQ,
		params.ckkscontext.levels, params.ckkscontext.alpha, params.ckkscontext.beta), func(t *testing.T) {

		values1, _, ciphertext1, err := new_test_vectors(params, -1, 1)
		if err != nil {
			t.Error(err)
		}

		values2, _, ciphertext2, err := new_test_vectors(params, -1, 1)
		if err != nil {
			t.Error(err)
		}

		// up to a level equal to 2 modulus
		valuesWant := make([]complex128, params.slots)

		for i := 0; i < len(valuesWant); i++ {
			valuesWant[i] = values1[i]
		}

		for i := uint64(0); i < params.ckkscontext.Levels()-1; i++ {

			for i := 0; i < len(valuesWant); i++ {
				valuesWant[i] *= values2[i]
			}

			if err = params.evaluator.MulRelin(ciphertext1, ciphertext2, params.rlk, ciphertext1); err != nil {
				t.Error(err)
			}
		}

		if err = params.evaluator.RescaleMany(ciphertext1, params.ckkscontext.Levels()-2, ciphertext1); err != nil {
			t.Error(err)
		}

		if err := verify_test_vectors(params, valuesWant, ciphertext1.Element(), t); err != nil {
			t.Error(err)
		}
	})

	t.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/a=%d/b=%d/MulRelin(Ct,Plain)->Rescale", params.ckkscontext.logN,
		params.ckkscontext.logQ,
		params.ckkscontext.levels, params.ckkscontext.alpha, params.ckkscontext.beta), func(t *testing.T) {

		values1, _, ciphertext1, err := new_test_vectors(params, -1, 1)
		if err != nil {
			t.Error(err)
		}

		values2, plaintext2, _, err := new_test_vectors(params, -1, 1)
		if err != nil {
			t.Error(err)
		}

		// up to a level equal to 2 modulus
		valuesWant := make([]complex128, params.slots)

		for i := 0; i < len(valuesWant); i++ {
			valuesWant[i] = values1[i]
		}

		for i := uint64(0); i < params.ckkscontext.Levels()-1; i++ {

			for i := 0; i < len(valuesWant); i++ {
				valuesWant[i] *= values2[i]
			}

			if err = params.evaluator.MulRelin(ciphertext1, plaintext2, params.rlk, ciphertext1); err != nil {
				t.Error(err)
			}

			if err = params.evaluator.Rescale(ciphertext1, params.ckkscontext.scale, ciphertext1); err != nil {
				t.Error(err)
			}
		}

		if err := verify_test_vectors(params, valuesWant, ciphertext1.Element(), t); err != nil {
			t.Error(err)
		}
	})

	t.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/a=%d/b=%d/MulRelin(Plain,Ct)->Rescale", params.ckkscontext.logN,
		params.ckkscontext.logQ,
		params.ckkscontext.levels, params.ckkscontext.alpha, params.ckkscontext.beta), func(t *testing.T) {

		values1, plaintext1, _, err := new_test_vectors(params, -1, 1)
		if err != nil {
			t.Error(err)
		}

		values2, _, ciphertext2, err := new_test_vectors(params, -1, 1)
		if err != nil {
			t.Error(err)
		}

		// up to a level equal to 2 modulus
		valuesWant := make([]complex128, params.slots)

		for i := 0; i < len(valuesWant); i++ {
			valuesWant[i] = values2[i]
		}

		for i := uint64(0); i < params.ckkscontext.Levels()-1; i++ {

			for i := 0; i < len(valuesWant); i++ {
				valuesWant[i] *= values1[i]
			}

			if err = params.evaluator.MulRelin(plaintext1, ciphertext2, params.rlk, ciphertext2); err != nil {
				t.Error(err)
			}

			if err = params.evaluator.Rescale(ciphertext2, params.ckkscontext.scale, ciphertext2); err != nil {
				t.Error(err)
			}
		}

		if err := verify_test_vectors(params, valuesWant, ciphertext2.Element(), t); err != nil {
			t.Error(err)
		}
	})

}

func test_Functions(params *CKKSTESTPARAMS, t *testing.T) {

	t.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/a=%d/b=%d/PowerOf2", params.ckkscontext.logN,
		params.ckkscontext.logQ,
		params.ckkscontext.levels, params.ckkscontext.alpha, params.ckkscontext.beta), func(t *testing.T) {

		values1, _, ciphertext1, err := new_test_vectors(params, -1, 1)
		if err != nil {
			t.Error(err)
		}

		// up to a level equal to 2 modulus
		valuesWant := make([]complex128, params.slots)

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
			t.Error(err)
		}

		if ciphertext1.Scale() >= 100 {
			params.evaluator.Rescale(ciphertext1, params.ckkscontext.scale, ciphertext1)
		}

		if err := verify_test_vectors(params, valuesWant, ciphertext1.Element(), t); err != nil {
			t.Error(err)
		}
	})

	t.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/a=%d/b=%d/Power", params.ckkscontext.logN,
		params.ckkscontext.logQ,
		params.ckkscontext.levels, params.ckkscontext.alpha, params.ckkscontext.beta), func(t *testing.T) {

		values1, _, ciphertext1, err := new_test_vectors(params, -1, 1)
		if err != nil {
			t.Error(err)
		}

		valuesWant := make([]complex128, params.slots)
		tmp := make([]complex128, params.slots)

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
			t.Error(err)
		}

		if ciphertext1.Scale() >= 100 {
			params.evaluator.Rescale(ciphertext1, params.ckkscontext.scale, ciphertext1)
		}

		if err := verify_test_vectors(params, valuesWant, ciphertext1.Element(), t); err != nil {
			t.Error(err)
		}
	})

	t.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/a=%d/b=%d/Inverse", params.ckkscontext.logN,
		params.ckkscontext.logQ,
		params.ckkscontext.levels, params.ckkscontext.alpha, params.ckkscontext.beta), func(t *testing.T) {

		values, _, ciphertext1, err := new_test_vectors_reals(params, 0.1, 1)
		if err != nil {
			t.Error(err)
		}

		valuesWant := make([]complex128, params.slots)

		for i := 0; i < len(valuesWant); i++ {
			valuesWant[i] = complex(1, 0) / values[i]
		}

		if ciphertext1, err = params.evaluator.InverseNew(ciphertext1, 7, params.rlk); err != nil {
			t.Error(err)
		}

		if err := verify_test_vectors(params, valuesWant, ciphertext1.Element(), t); err != nil {
			t.Error(err)
		}
	})

	t.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/a=%d/b=%d/sin(x) [-1-1i, 1+1i] deg16", params.ckkscontext.logN,
		params.ckkscontext.logQ,
		params.ckkscontext.levels, params.ckkscontext.alpha, params.ckkscontext.beta), func(t *testing.T) {

		values, _, ciphertext1, err := new_test_vectors(params, -1, 1)
		if err != nil {
			t.Error(err)
		}

		valuesWant := make([]complex128, params.slots)

		for i := 0; i < len(valuesWant); i++ {
			valuesWant[i] = cmplx.Sin(values[i])
		}

		cheby := Approximate(cmplx.Sin, complex(-1, -1), complex(1, 1), 16)

		if ciphertext1, err = params.evaluator.EvaluateCheby(ciphertext1, cheby, params.rlk); err != nil {
			t.Error(err)
		}

		if err := verify_test_vectors(params, valuesWant, ciphertext1.Element(), t); err != nil {
			t.Error(err)
		}
	})
}

func test_SwitchKeys(params *CKKSTESTPARAMS, t *testing.T) {

	t.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/a=%d/b=%d/SwitchKeys", params.ckkscontext.logN,
		params.ckkscontext.logQ,
		params.ckkscontext.levels, params.ckkscontext.alpha, params.ckkscontext.beta), func(t *testing.T) {

		valuesWant := make([]complex128, params.slots)

		for i := uint64(0); i < params.slots; i++ {
			valuesWant[i] = randomComplex(0, 5)
		}

		plaintext := params.ckkscontext.NewPlaintext(params.ckkscontext.Levels()-1, params.ckkscontext.Scale())

		if err := params.encoder.Encode(plaintext, valuesWant, params.slots); err != nil {
			t.Error(err)
		}

		ciphertext, err := params.encryptorPk.EncryptNew(plaintext)
		if err != nil {
			t.Error(err)
		}

		sk2 := params.kgen.NewSecretKey()

		switchingkeys, err := params.kgen.NewSwitchingKey(params.sk, sk2)
		if err != nil {
			t.Error(err)
		}

		if err = params.evaluator.SwitchKeys(ciphertext, switchingkeys, ciphertext); err != nil {
			t.Error(err)
		}

		decryptorSk2, err := params.ckkscontext.NewDecryptor(sk2)
		if err != nil {
			t.Error(err)
		}

		plaintextTest := decryptorSk2.DecryptNew(ciphertext)

		if err := verify_test_vectors(params, valuesWant, plaintextTest.Element(), t); err != nil {
			t.Error(err)
		}
	})
}

func test_Conjugate(params *CKKSTESTPARAMS, t *testing.T) {

	values, _, ciphertext, err := new_test_vectors_reals(params, -15, 15)
	if err != nil {
		t.Error(err)
	}

	valuesWant := make([]complex128, len(values))
	for i := range values {
		valuesWant[i] = complex(real(values[i]), -imag(values[i]))
	}

	t.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/a=%d/b=%d/Conjugate(Ct)", params.ckkscontext.logN,
		params.ckkscontext.logQ,
		params.ckkscontext.levels, params.ckkscontext.alpha, params.ckkscontext.beta), func(t *testing.T) {

		if err := params.evaluator.Conjugate(ciphertext, params.rotkey, ciphertext); err != nil {
			t.Error(err)
		}

		if err := verify_test_vectors(params, valuesWant, ciphertext.Element(), t); err != nil {
			t.Error(err)
		}
	})
}

func test_RotColumns(params *CKKSTESTPARAMS, t *testing.T) {

	mask := params.slots - 1

	values, _, ciphertext, err := new_test_vectors_reals(params, 0.1, 1)
	if err != nil {
		t.Error(err)
	}

	valuesWant := make([]complex128, params.slots)

	ciphertextTest := params.ckkscontext.NewCiphertext(1, ciphertext.Level(), ciphertext.Scale())
	ciphertextTest.SetScale(ciphertext.Scale())

	t.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/a=%d/b=%d/RotColumnsPow2(Ct)", params.ckkscontext.logN,
		params.ckkscontext.logQ,
		params.ckkscontext.levels, params.ckkscontext.alpha, params.ckkscontext.beta), func(t *testing.T) {

		// Applies the column rotation to the values
		for i := uint64(0); i < params.slots; i++ {
			valuesWant[i] = values[i]
		}

		if err := params.evaluator.RotateColumns(ciphertext, 0, params.rotkey, ciphertextTest); err != nil {
			t.Error(err)
		}

		if err := verify_test_vectors(params, valuesWant, ciphertextTest.Element(), t); err != nil {
			t.Error(err)
		}

		for n := uint64(1); n < params.slots; n <<= 1 {

			// Applies the column rotation to the values
			for i := uint64(0); i < params.slots; i++ {
				valuesWant[i] = values[(i+n)&mask]
			}

			if err := params.evaluator.RotateColumns(ciphertext, n, params.rotkey, ciphertextTest); err != nil {
				t.Error(err)
			}

			if err := verify_test_vectors(params, valuesWant, ciphertextTest.Element(), t); err != nil {
				t.Error(err)
			}
		}
	})

	t.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/a=%d/b=%d/RotColumnsRandom(Ct)", params.ckkscontext.logN,
		params.ckkscontext.logQ,
		params.ckkscontext.levels, params.ckkscontext.alpha, params.ckkscontext.beta), func(t *testing.T) {

		for n := uint64(0); n < 4; n += 1 {

			rand := ring.RandUniform(params.slots, mask)

			// Applies the column rotation to the values
			for i := uint64(0); i < params.slots; i++ {
				valuesWant[i] = values[(i+rand)&mask]
			}

			if err := params.evaluator.RotateColumns(ciphertext, rand, params.rotkey, ciphertextTest); err != nil {
				t.Error(err)
			}

			if err := verify_test_vectors(params, valuesWant, ciphertextTest.Element(), t); err != nil {
				t.Error(err)
			}
		}
	})
}
