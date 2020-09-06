package ckks

import (
	"flag"
	"fmt"
	"math"
	"math/cmplx"
	"math/rand"
	"testing"
	"time"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"

	"github.com/ldsec/lattigo/utils"
)

var err error
var params = new(testParams)
var defaultParams = DefaultParams[PN12QP109 : PN12QP109+3]
var printPrecisionStats = flag.Bool("print-precision", false, "print precision stats")

func testString(opname string) string {
	return fmt.Sprintf("%slogN=%d/LogSlots=%d/logQP=%d/levels=%d/a=%d/b=%d",
		opname,
		params.params.LogN(),
		params.params.LogSlots(),
		params.params.LogQP(),
		params.params.MaxLevel()+1,
		params.params.Alpha(),
		params.params.Beta())
}

type testParams struct {
	params      *Parameters
	medianprec  float64
	ckkscontext *Context
	prng        utils.PRNG
	encoder     Encoder
	kgen        KeyGenerator
	sk          *SecretKey
	pk          *PublicKey
	rlk         *EvaluationKey
	encryptorPk Encryptor
	encryptorSk Encryptor
	decryptor   Decryptor
	evaluator   Evaluator
}

func TestCKKS(t *testing.T) {

	rand.Seed(time.Now().UnixNano())

	params.medianprec = 15

	for _, defaultParam := range defaultParams {

		if err = genTestParams(defaultParam); err != nil {
			panic(err)
		}

		t.Run("Parameters", testParameters)
		t.Run("Encoder", testEncoder)
		t.Run("Encryptor", testEncryptor)
		t.Run("Evaluator/Add", testEvaluatorAdd)
		t.Run("Evaluator/Sub", testEvaluatorSub)
		t.Run("Evaluator/Rescale", testEvaluatorRescale)
		t.Run("Evaluator/AddConst", testEvaluatorAddConst)
		t.Run("Evaluator/MultByConst", testEvaluatorMultByConst)
		t.Run("Evaluator/MultByConstAndAdd", testEvaluatorMultByConstAndAdd)
		t.Run("Evaluator/Mul", testEvaluatorMul)
		t.Run("Evaluator/Functions", testFunctions)
		t.Run("Evaluator/EvaluatePoly", testEvaluatePoly)
		t.Run("Evaluator/ChebyshevInterpolator", testChebyshevInterpolator)
		t.Run("Evaluator/SwitchKeys", testSwitchKeys)
		t.Run("Evaluator/Conjugate", testConjugate)
		t.Run("Evaluator/RotateColumns", testRotateColumns)
		t.Run("Marshalling", testMarshaller)
	}
}

func genTestParams(defaultParams interface{}) (err error) {

	params = new(testParams)

	switch defaultParams.(type) {
	case *DefaultParam:

		p := defaultParams.(*DefaultParam)

		if params.params, err = NewParametersFromLogModuli(p.LogN, p.LogModuli); err != nil {
			return err
		}

		params.params.SetLogSlots(params.params.LogMaxSlots())
		params.params.SetScale(p.Scale)

		params.kgen = NewKeyGenerator(params.params)
		params.sk, params.pk = params.kgen.GenKeyPair()

	case *BootParams:

		p := defaultParams.(*BootParams)

		if params.params, err = NewParametersFromModuli(p.logN, p.Moduli); err != nil {
			return err
		}

		params.params.SetLogSlots(p.logSlots)
		params.params.SetScale(p.scale)

		params.kgen = NewKeyGenerator(params.params)
		params.sk, params.pk = params.kgen.GenKeyPairSparse(p.H)

	default:
		return fmt.Errorf("Invalid input type for 'genTestParams'")

	}

	params.ckkscontext = newContext(params.params)

	if params.prng, err = utils.NewPRNG(); err != nil {
		return err
	}

	params.encoder = NewEncoder(params.params)

	params.rlk = params.kgen.GenRelinKey(params.sk)

	params.encryptorPk = NewEncryptorFromPk(params.params, params.pk)
	params.encryptorSk = NewEncryptorFromSk(params.params, params.sk)
	params.decryptor = NewDecryptor(params.params, params.sk)

	params.evaluator = NewEvaluator(params.params)

	return nil

}

func newTestVectors(encryptor Encryptor, a float64, t *testing.T) (values []complex128, plaintext *Plaintext, ciphertext *Ciphertext) {

	slots := params.params.Slots()

	values = make([]complex128, slots)

	for i := uint64(0); i < slots; i++ {
		values[i] = randomComplex(-a, a)
	}

	values[0] = complex(0.607538, 0.555668)

	plaintext = NewPlaintext(params.params, params.params.MaxLevel(), params.params.Scale())

	params.encoder.EncodeNTT(plaintext, values, slots)

	if encryptor != nil {
		ciphertext = encryptor.EncryptNew(plaintext)
	}

	return values, plaintext, ciphertext
}

func newTestVectorsReals(encryptor Encryptor, a, b float64, t *testing.T) (values []complex128, plaintext *Plaintext, ciphertext *Ciphertext) {

	slots := params.params.Slots()

	values = make([]complex128, slots)

	for i := uint64(0); i < slots; i++ {
		values[i] = complex(randomFloat(a, b), 0)
	}

	values[0] = complex(0.607538, 0)

	plaintext = NewPlaintext(params.params, params.params.MaxLevel(), params.params.Scale())

	params.encoder.EncodeNTT(plaintext, values, slots)

	if encryptor != nil {
		ciphertext = encryptor.EncryptNew(plaintext)
	}

	return values, plaintext, ciphertext
}

func newTestVectorsSineBoot(encryptor Encryptor, a, b float64, t *testing.T) (values []complex128, plaintext *Plaintext, ciphertext *Ciphertext) {

	slots := params.params.Slots()

	values = make([]complex128, slots)

	for i := uint64(0); i < slots; i++ {
		values[i] = complex(math.Round(randomFloat(a, b))+randomFloat(-1, 1)/1000, 0)
	}

	plaintext = NewPlaintext(params.params, params.params.MaxLevel(), params.params.Scale())

	params.encoder.EncodeNTT(plaintext, values, slots)

	if encryptor != nil {
		ciphertext = encryptor.EncryptNew(plaintext)
	}

	return values, plaintext, ciphertext
}

func verifyTestVectors(decryptor Decryptor, valuesWant []complex128, element interface{}, t *testing.T) {

	var plaintextTest *Plaintext
	var valuesTest []complex128

	switch element.(type) {
	case *Ciphertext:
		plaintextTest = decryptor.DecryptNew(element.(*Ciphertext))
	case *Plaintext:
		plaintextTest = element.(*Plaintext)
	}

	slots := params.params.Slots()

	valuesTest = params.encoder.Decode(plaintextTest, slots)

	var deltaReal, deltaImag float64

	var delta, minprec, maxprec, meanprec, medianprec complex128

	diff := make([]complex128, slots)

	minprec = complex(0, 0)
	maxprec = complex(1, 1)

	meanprec = complex(0, 0)

	distribReal := make(map[uint64]uint64)
	distribImag := make(map[uint64]uint64)

	distribPrec := float64(25)

	for i := range valuesWant {

		delta = valuesTest[i] - valuesWant[i]
		deltaReal = math.Abs(real(delta))
		deltaImag = math.Abs(imag(delta))

		diff[i] += complex(deltaReal, deltaImag)

		meanprec += diff[i]

		if deltaReal > real(minprec) {
			minprec = complex(deltaReal, imag(minprec))
		}

		if deltaImag > imag(minprec) {
			minprec = complex(real(minprec), deltaImag)
		}

		if deltaReal < real(maxprec) {
			maxprec = complex(deltaReal, imag(maxprec))
		}

		if deltaImag < imag(maxprec) {
			maxprec = complex(real(maxprec), deltaImag)
		}

		distribReal[uint64(math.Floor(distribPrec*math.Log2(1/deltaReal)))]++
		distribImag[uint64(math.Floor(distribPrec*math.Log2(1/deltaImag)))]++
	}

	meanprec /= complex(float64(slots), 0)
	medianprec = calcmedian(diff)

	if *printPrecisionStats {
		t.Logf("Minimum precision : (%.2f, %.2f) bits \n", math.Log2(1/real(minprec)), math.Log2(1/imag(minprec)))
		t.Logf("Maximum precision : (%.2f, %.2f) bits \n", math.Log2(1/real(maxprec)), math.Log2(1/imag(maxprec)))
		t.Logf("Mean    precision : (%.2f, %.2f) bits \n", math.Log2(1/real(meanprec)), math.Log2(1/imag(meanprec)))
		t.Logf("Median  precision : (%.2f, %.2f) bits \n", math.Log2(1/real(medianprec)), math.Log2(1/imag(medianprec)))
		t.Log()
	}

	require.GreaterOrEqual(t, math.Log2(1/real(medianprec)), params.medianprec)
	require.GreaterOrEqual(t, math.Log2(1/imag(medianprec)), params.medianprec)
}

func testParameters(t *testing.T) {
	t.Run("NewParametersFromModuli", func(t *testing.T) {
		p, err := NewParametersFromModuli(params.params.LogN(), params.params.Moduli)
		p.SetLogSlots(params.params.LogSlots())
		p.SetScale(params.params.Scale())
		assert.NoError(t, err)
		assert.True(t, p.Equals(params.params))
	})

	t.Run("NewParametersFromLogModuli", func(t *testing.T) {
		p, err := NewParametersFromLogModuli(params.params.LogN(), params.params.LogModuli())
		p.SetLogSlots(params.params.LogSlots())
		p.SetScale(params.params.Scale())
		assert.NoError(t, err)
		assert.True(t, p.Equals(params.params))
	})
}

func testEncoder(t *testing.T) {

	t.Run(testString("Encode/"), func(t *testing.T) {

		values, plaintext, _ := newTestVectors(nil, 1, t)

		verifyTestVectors(params.decryptor, values, plaintext, t)
	})

	t.Run(testString("EncodeCoeffs/"), func(t *testing.T) {

		slots := params.params.N()

		valuesWant := make([]float64, slots)

		for i := uint64(0); i < slots; i++ {
			valuesWant[i] = randomFloat(-1, 1)
		}

		valuesWant[0] = 0.607538

		plaintext := NewPlaintext(params.params, params.params.MaxLevel(), params.params.Scale())

		params.encoder.EncodeCoeffs(valuesWant, plaintext)

		valuesTest := params.encoder.DecodeCoeffs(plaintext)

		var meanprec float64

		for i := range valuesWant {
			meanprec += math.Abs(valuesTest[i] - valuesWant[i])
		}

		meanprec /= float64(slots)

		require.GreaterOrEqual(t, math.Log2(1/meanprec), params.medianprec)
	})

}

func testEncryptor(t *testing.T) {

	t.Run(testString("EncryptFromPk/"), func(t *testing.T) {

		values, _, ciphertext := newTestVectors(params.encryptorPk, 1, t)

		verifyTestVectors(params.decryptor, values, ciphertext, t)
	})

	t.Run(testString("EncryptFromPkFast/"), func(t *testing.T) {

		slots := params.params.Slots()

		values := make([]complex128, slots)

		for i := uint64(0); i < slots; i++ {
			values[i] = randomComplex(-1, 1)
		}

		values[0] = complex(0.607538, 0.555668)

		plaintext := NewPlaintext(params.params, params.params.MaxLevel(), params.params.Scale())

		params.encoder.Encode(plaintext, values, slots)

		verifyTestVectors(params.decryptor, values, params.encryptorPk.EncryptFastNew(plaintext), t)
	})

	t.Run(testString("EncryptFromSk/"), func(t *testing.T) {

		slots := params.params.Slots()

		values := make([]complex128, slots)

		for i := uint64(0); i < slots; i++ {
			values[i] = randomComplex(-1, 1)
		}

		values[0] = complex(0.607538, 0.555668)

		plaintext := NewPlaintext(params.params, params.params.MaxLevel(), params.params.Scale())

		params.encoder.Encode(plaintext, values, slots)

		verifyTestVectors(params.decryptor, values, params.encryptorSk.EncryptNew(plaintext), t)
	})

}

func testEvaluatorAdd(t *testing.T) {

	t.Run(testString("CtCtInPlace/"), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectors(params.encryptorSk, 1, t)
		values2, _, ciphertext2 := newTestVectors(params.encryptorSk, 1, t)

		for i := range values1 {
			values1[i] += values2[i]
		}

		params.evaluator.Add(ciphertext1, ciphertext2, ciphertext1)

		verifyTestVectors(params.decryptor, values1, ciphertext1, t)
	})

	t.Run(testString("CtCtNew/"), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectors(params.encryptorSk, 1, t)
		values2, _, ciphertext2 := newTestVectors(params.encryptorSk, 1, t)

		for i := range values1 {
			values1[i] += values2[i]
		}

		ciphertext3 := params.evaluator.AddNew(ciphertext1, ciphertext2)

		verifyTestVectors(params.decryptor, values1, ciphertext3, t)
	})

	t.Run(testString("CtPlainInPlace/"), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectors(params.encryptorSk, 1, t)
		values2, plaintext2, _ := newTestVectors(params.encryptorSk, 1, t)

		for i := range values1 {
			values1[i] += values2[i]
		}

		params.evaluator.Add(ciphertext1, plaintext2, ciphertext1)

		verifyTestVectors(params.decryptor, values1, ciphertext1, t)

		for i := range values1 {
			values1[i] += values2[i]
		}

		params.evaluator.Add(plaintext2, ciphertext1, ciphertext1)

		verifyTestVectors(params.decryptor, values1, ciphertext1, t)
	})

	t.Run(testString("CtPlainInPlaceNew/"), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectors(params.encryptorSk, 1, t)
		values2, plaintext2, _ := newTestVectors(params.encryptorSk, 1, t)

		for i := range values1 {
			values1[i] += values2[i]
		}

		ciphertext3 := params.evaluator.AddNew(ciphertext1, plaintext2)

		verifyTestVectors(params.decryptor, values1, ciphertext3, t)

		ciphertext3 = params.evaluator.AddNew(plaintext2, ciphertext1)

		verifyTestVectors(params.decryptor, values1, ciphertext3, t)
	})

}

func testEvaluatorSub(t *testing.T) {

	t.Run(testString("CtCtInPlace/"), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectors(params.encryptorSk, 1, t)
		values2, _, ciphertext2 := newTestVectors(params.encryptorSk, 1, t)

		for i := range values1 {
			values1[i] -= values2[i]
		}

		params.evaluator.Sub(ciphertext1, ciphertext2, ciphertext1)

		verifyTestVectors(params.decryptor, values1, ciphertext1, t)
	})

	t.Run(testString("CtCtNew/"), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectors(params.encryptorSk, 1, t)
		values2, _, ciphertext2 := newTestVectors(params.encryptorSk, 1, t)

		for i := range values1 {
			values1[i] -= values2[i]
		}

		ciphertext3 := params.evaluator.SubNew(ciphertext1, ciphertext2)

		verifyTestVectors(params.decryptor, values1, ciphertext3, t)
	})

	t.Run(testString("CtPlainInPlace/"), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectors(params.encryptorSk, 1, t)
		values2, plaintext2, ciphertext2 := newTestVectors(params.encryptorSk, 1, t)

		valuesTest := make([]complex128, len(values1))
		for i := range values1 {
			valuesTest[i] = values1[i] - values2[i]
		}

		params.evaluator.Sub(ciphertext1, plaintext2, ciphertext2)

		verifyTestVectors(params.decryptor, valuesTest, ciphertext2, t)

		for i := range values1 {
			valuesTest[i] = values2[i] - values1[i]
		}

		params.evaluator.Sub(plaintext2, ciphertext1, ciphertext2)

		verifyTestVectors(params.decryptor, valuesTest, ciphertext2, t)
	})

	t.Run(testString("CtPlainNew/"), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectors(params.encryptorSk, 1, t)
		values2, plaintext2, _ := newTestVectors(params.encryptorSk, 1, t)

		valuesTest := make([]complex128, len(values1))
		for i := range values1 {
			valuesTest[i] = values1[i] - values2[i]
		}

		ciphertext3 := params.evaluator.SubNew(ciphertext1, plaintext2)

		verifyTestVectors(params.decryptor, valuesTest, ciphertext3, t)

		for i := range values1 {
			valuesTest[i] = values2[i] - values1[i]
		}

		ciphertext3 = params.evaluator.SubNew(plaintext2, ciphertext1)

		verifyTestVectors(params.decryptor, valuesTest, ciphertext3, t)
	})

}

func testEvaluatorRescale(t *testing.T) {

	t.Run(testString("Single/"), func(t *testing.T) {

		values, _, ciphertext := newTestVectors(params.encryptorSk, 1, t)

		constant := params.ckkscontext.contextQ.Modulus[ciphertext.Level()]

		params.evaluator.MultByConst(ciphertext, constant, ciphertext)

		ciphertext.MulScale(float64(constant))

		params.evaluator.Rescale(ciphertext, params.params.Scale(), ciphertext)

		verifyTestVectors(params.decryptor, values, ciphertext, t)
	})

	t.Run(testString("Many/"), func(t *testing.T) {

		values, _, ciphertext := newTestVectors(params.encryptorSk, 1, t)

		nbRescales := params.params.MaxLevel()

		for i := uint64(0); i < nbRescales; i++ {
			constant := params.ckkscontext.contextQ.Modulus[ciphertext.Level()]
			params.evaluator.MultByConst(ciphertext, constant, ciphertext)
			ciphertext.MulScale(float64(constant))
		}

		params.evaluator.RescaleMany(ciphertext, nbRescales, ciphertext)

		verifyTestVectors(params.decryptor, values, ciphertext, t)
	})
}

func testEvaluatorAddConst(t *testing.T) {

	t.Run(testString(""), func(t *testing.T) {

		values, _, ciphertext := newTestVectors(params.encryptorSk, 1, t)

		constant := complex(3.1415, -1.4142)

		for i := range values {
			values[i] += constant
		}

		params.evaluator.AddConst(ciphertext, constant, ciphertext)

		verifyTestVectors(params.decryptor, values, ciphertext, t)
	})

}

func testEvaluatorMultByConst(t *testing.T) {

	t.Run(testString(""), func(t *testing.T) {

		values, _, ciphertext := newTestVectors(params.encryptorSk, 1, t)

		constant := 1.0 / complex(3.1415, -1.4142)

		for i := range values {
			values[i] *= constant
		}

		params.evaluator.MultByConst(ciphertext, constant, ciphertext)

		verifyTestVectors(params.decryptor, values, ciphertext, t)
	})

}

func testEvaluatorMultByConstAndAdd(t *testing.T) {

	t.Run(testString(""), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectors(params.encryptorSk, 1, t)
		values2, _, ciphertext2 := newTestVectors(params.encryptorSk, 1, t)

		constant := 1.0 / complex(3.1415, -1.4142)

		for i := range values1 {
			values2[i] += (constant * values1[i])
		}

		params.evaluator.MultByConstAndAdd(ciphertext1, constant, ciphertext2)

		verifyTestVectors(params.decryptor, values2, ciphertext2, t)
	})

}

func testEvaluatorMul(t *testing.T) {

	t.Run(testString("CtCtInPlace/"), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectors(params.encryptorSk, 1, t)
		values2, _, ciphertext2 := newTestVectors(params.encryptorSk, 1, t)

		for i := range values1 {
			values2[i] *= values1[i]
		}

		params.evaluator.MulRelin(ciphertext1, ciphertext2, nil, ciphertext2)

		verifyTestVectors(params.decryptor, values2, ciphertext2, t)
	})

	t.Run(testString("CtCtNew/"), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectors(params.encryptorSk, 1, t)
		values2, _, ciphertext2 := newTestVectors(params.encryptorSk, 1, t)

		for i := range values1 {
			values2[i] *= values1[i]
		}

		var ciphertext3 *Ciphertext
		ciphertext3 = params.evaluator.MulRelinNew(ciphertext1, ciphertext2, nil)

		verifyTestVectors(params.decryptor, values2, ciphertext3, t)
	})

	t.Run(testString("CtPlain/"), func(t *testing.T) {

		values1, plaintext1, ciphertext1 := newTestVectors(params.encryptorSk, 1, t)
		values2, plaintext2, ciphertext2 := newTestVectors(params.encryptorSk, 1, t)

		for i := range values1 {
			values1[i] *= values1[i]
			values2[i] *= values2[i]
		}

		params.evaluator.MulRelin(ciphertext1, plaintext1, nil, ciphertext1)

		verifyTestVectors(params.decryptor, values1, ciphertext1, t)

		params.evaluator.MulRelin(plaintext2, ciphertext2, nil, ciphertext2)

		verifyTestVectors(params.decryptor, values2, ciphertext2, t)
	})

	t.Run(testString("Relinearize/"), func(t *testing.T) {

		rlk := params.kgen.GenRelinKey(params.sk)

		values1, _, ciphertext1 := newTestVectors(params.encryptorSk, 1, t)
		values2, _, ciphertext2 := newTestVectors(params.encryptorSk, 1, t)

		for i := range values1 {
			values2[i] *= values1[i]
		}

		params.evaluator.MulRelin(ciphertext1, ciphertext2, nil, ciphertext2)

		params.evaluator.Relinearize(ciphertext2, rlk, ciphertext2)

		require.Equal(t, ciphertext2.Degree(), uint64(1))

		verifyTestVectors(params.decryptor, values2, ciphertext2, t)
	})

}

func testFunctions(t *testing.T) {

	if params.params.MaxLevel() > 2 {

		t.Run(testString("PowerOf2/"), func(t *testing.T) {

			values, _, ciphertext := newTestVectors(params.encryptorSk, 1, t)

			n := uint64(2)

			valuesWant := make([]complex128, len(values))
			for i := 0; i < len(valuesWant); i++ {
				valuesWant[i] = values[i]
			}

			for i := uint64(0); i < n; i++ {
				for j := 0; j < len(valuesWant); j++ {
					valuesWant[j] *= valuesWant[j]
				}
			}

			params.evaluator.PowerOf2(ciphertext, n, params.rlk, ciphertext)

			verifyTestVectors(params.decryptor, valuesWant, ciphertext, t)
		})
	}

	if params.params.MaxLevel() > 3 {

		t.Run(testString("Power/"), func(t *testing.T) {

			values, _, ciphertext := newTestVectors(params.encryptorSk, 1, t)

			n := uint64(3)

			for i := range values {
				values[i] = cmplx.Pow(values[i], complex(float64(n), 0))
			}

			params.evaluator.Power(ciphertext, n, params.rlk, ciphertext)

			verifyTestVectors(params.decryptor, values, ciphertext, t)
		})
	}

	if params.params.MaxLevel() > 6 {
		t.Run(testString("Inverse/"), func(t *testing.T) {

			values, _, ciphertext := newTestVectorsReals(params.encryptorSk, 0.1, 1, t)

			n := uint64(7)

			for i := range values {
				values[i] = 1.0 / values[i]
			}

			ciphertext = params.evaluator.InverseNew(ciphertext, n, params.rlk)

			verifyTestVectors(params.decryptor, values, ciphertext, t)
		})
	}
}

func testEvaluatePoly(t *testing.T) {

	if params.params.MaxLevel() > 3 {

		t.Run(testString("Exp/"), func(t *testing.T) {

			values, _, ciphertext := newTestVectorsReals(params.encryptorSk, -1, 1, t)

			coeffs := []complex128{
				complex(1.0, 0),
				complex(1.0, 0),
				complex(1.0/2, 0),
				complex(1.0/6, 0),
				complex(1.0/24, 0),
				complex(1.0/120, 0),
				complex(1.0/720, 0),
				complex(1.0/5040, 0),
			}

			poly := NewPoly(coeffs)

			for i := range values {
				values[i] = cmplx.Exp(values[i])
			}

			ciphertext = params.evaluator.EvaluatePoly(ciphertext, poly, params.rlk)

			verifyTestVectors(params.decryptor, values, ciphertext, t)
		})
	}
}

func testChebyshevInterpolator(t *testing.T) {

	if params.params.MaxLevel() > 5 {

		rlk := params.kgen.GenRelinKey(params.sk)

		t.Run(testString("Sin/"), func(t *testing.T) {

			values, _, ciphertext := newTestVectorsReals(params.encryptorSk, -1, 1, t)

			cheby := Approximate(cmplx.Sin, complex(-3, 0), complex(3, 0), 15)

			for i := range values {
				values[i] = cmplx.Sin(values[i])
			}

			ciphertext = params.evaluator.EvaluateCheby(ciphertext, cheby, rlk)

			verifyTestVectors(params.decryptor, values, ciphertext, t)
		})
	}
}

func testSwitchKeys(t *testing.T) {

	sk2 := params.kgen.GenSecretKey()
	decryptorSk2 := NewDecryptor(params.params, sk2)
	switchingKey := params.kgen.GenSwitchingKey(params.sk, sk2)

	t.Run(testString("InPlace/"), func(t *testing.T) {

		values, _, ciphertext := newTestVectorsReals(params.encryptorSk, -1, 1, t)

		params.evaluator.SwitchKeys(ciphertext, switchingKey, ciphertext)

		verifyTestVectors(decryptorSk2, values, ciphertext, t)
	})

	t.Run(testString("New/"), func(t *testing.T) {

		values, _, ciphertext := newTestVectorsReals(params.encryptorSk, -1, 1, t)

		ciphertext = params.evaluator.SwitchKeysNew(ciphertext, switchingKey)

		verifyTestVectors(decryptorSk2, values, ciphertext, t)
	})

}

func testConjugate(t *testing.T) {

	rotKey := NewRotationKeys()
	params.kgen.GenRot(Conjugate, params.sk, 0, rotKey)

	t.Run(testString("InPlace/"), func(t *testing.T) {

		values, _, ciphertext := newTestVectorsReals(params.encryptorSk, -1, 1, t)

		for i := range values {
			values[i] = complex(real(values[i]), -imag(values[i]))
		}

		params.evaluator.Conjugate(ciphertext, rotKey, ciphertext)

		verifyTestVectors(params.decryptor, values, ciphertext, t)
	})

	t.Run(testString("New/"), func(t *testing.T) {

		values, _, ciphertext := newTestVectorsReals(params.encryptorSk, -1, 1, t)

		for i := range values {
			values[i] = complex(real(values[i]), -imag(values[i]))
		}

		ciphertext = params.evaluator.ConjugateNew(ciphertext, rotKey)

		verifyTestVectors(params.decryptor, values, ciphertext, t)
	})

}

func testRotateColumns(t *testing.T) {

	rotKey := params.kgen.GenRotationKeysPow2(params.sk)

	t.Run(testString("InPlace/"), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectorsReals(params.encryptorSk, -1, 1, t)

		values2 := make([]complex128, len(values1))
		ciphertext2 := NewCiphertext(params.params, ciphertext1.Degree(), ciphertext1.Level(), ciphertext1.Scale())

		for n := 1; n < len(values1); n <<= 1 {

			// Applies the column rotation to the values
			for i := range values1 {
				values2[i] = values1[(i+n)%len(values1)]
			}

			params.evaluator.RotateColumns(ciphertext1, uint64(n), rotKey, ciphertext2)

			verifyTestVectors(params.decryptor, values2, ciphertext2, t)
		}

	})

	t.Run(testString("New/"), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectorsReals(params.encryptorSk, -1, 1, t)

		values2 := make([]complex128, len(values1))
		ciphertext2 := NewCiphertext(params.params, ciphertext1.Degree(), ciphertext1.Level(), ciphertext1.Scale())

		for n := 1; n < len(values1); n <<= 1 {

			// Applies the column rotation to the values
			for i := range values1 {
				values2[i] = values1[(i+n)%len(values1)]
			}

			ciphertext2 = params.evaluator.RotateColumnsNew(ciphertext1, uint64(n), rotKey)

			verifyTestVectors(params.decryptor, values2, ciphertext2, t)
		}

	})

	t.Run(testString("Random/"), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectorsReals(params.encryptorSk, -1, 1, t)

		values2 := make([]complex128, len(values1))
		ciphertext2 := NewCiphertext(params.params, ciphertext1.Degree(), ciphertext1.Level(), ciphertext1.Scale())

		for n := 1; n < 4; n++ {

			rand := rand.Uint64() % uint64(len(values1))

			// Applies the column rotation to the values
			for i := range values1 {
				values2[i] = values1[(i+int(rand))%len(values1)]
			}

			params.evaluator.RotateColumns(ciphertext1, rand, rotKey, ciphertext2)

			verifyTestVectors(params.decryptor, values2, ciphertext2, t)
		}

	})

	t.Run(testString("Hoisted/"), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectorsReals(params.encryptorSk, -1, 1, t)

		values2 := make([]complex128, len(values1))
		rotations := []uint64{0, 1, 2, 3, 4, 5}
		for _, n := range rotations {
			params.kgen.GenRot(RotationLeft, params.sk, n, rotKey)
		}

		ciphertexts := params.evaluator.RotateHoisted(ciphertext1, rotations, rotKey)

		for _, n := range rotations {

			for i := range values1 {
				values2[i] = values1[(i+int(n))%len(values1)]
			}

			verifyTestVectors(params.decryptor, values2, ciphertexts[n], t)
		}

	})
}

func testMarshaller(t *testing.T) {

	contextQP := params.ckkscontext.contextQP

	t.Run(testString("Ciphertext/"), func(t *testing.T) {
		t.Run(testString("Ciphertext/EndToEnd"), func(t *testing.T) {
			t.Parallel()

			ciphertextWant := NewCiphertextRandom(params.prng, params.params, 2, params.params.MaxLevel(), params.params.Scale())

			marshalledCiphertext, err := ciphertextWant.MarshalBinary()
			require.NoError(t, err)

			ciphertextTest := new(Ciphertext)
			require.NoError(t, ciphertextTest.UnmarshalBinary(marshalledCiphertext))

			require.Equal(t, ciphertextWant.Degree(), ciphertextTest.Degree())
			require.Equal(t, ciphertextWant.Level(), ciphertextTest.Level())
			require.Equal(t, ciphertextWant.Scale(), ciphertextTest.Scale())

			for i := range ciphertextWant.value {
				require.True(t, params.ckkscontext.contextQ.EqualLvl(ciphertextWant.Level(), ciphertextWant.Value()[i], ciphertextTest.Value()[i]))
			}
		})

		t.Run(testString("Ciphertext/Minimal"), func(t *testing.T) {
			t.Parallel()

			ciphertext := NewCiphertextRandom(params.prng, params.params, 0, params.params.MaxLevel(), params.params.Scale())

			marshalledCiphertext, err := ciphertext.MarshalBinary()
			require.NoError(t, err)

			ciphertextTest := new(Ciphertext)
			require.Error(t, ciphertextTest.UnmarshalBinary(nil))
			require.NoError(t, ciphertextTest.UnmarshalBinary(marshalledCiphertext))

			require.Equal(t, ciphertext.Degree(), uint64(0))
			require.Equal(t, ciphertext.Level(), params.params.MaxLevel())
			require.Equal(t, ciphertext.Scale(), params.params.Scale())
			require.Equal(t, len(ciphertext.Value()), 1)
		})
	})

	t.Run(testString("Sk"), func(t *testing.T) {

		marshalledSk, err := params.sk.MarshalBinary()
		require.NoError(t, err)

		sk := new(SecretKey)
		err = sk.UnmarshalBinary(marshalledSk)
		require.NoError(t, err)

		require.True(t, contextQP.Equal(sk.sk, params.sk.sk))

	})

	t.Run(testString("Pk"), func(t *testing.T) {

		marshalledPk, err := params.pk.MarshalBinary()
		require.NoError(t, err)

		pk := new(PublicKey)
		err = pk.UnmarshalBinary(marshalledPk)
		require.NoError(t, err)

		for k := range params.pk.pk {
			require.Truef(t, contextQP.Equal(pk.pk[k], params.pk.pk[k]), "Marshal PublicKey element [%d]", k)
		}
	})

	t.Run(testString("EvaluationKey"), func(t *testing.T) {

		evalKey := params.kgen.GenRelinKey(params.sk)
		data, err := evalKey.MarshalBinary()
		require.NoError(t, err)

		resEvalKey := new(EvaluationKey)
		err = resEvalKey.UnmarshalBinary(data)
		require.NoError(t, err)

		evakeyWant := evalKey.evakey.evakey
		evakeyTest := resEvalKey.evakey.evakey

		for j := range evakeyWant {
			for k := range evakeyWant[j] {
				require.Truef(t, contextQP.Equal(evakeyWant[j][k], evakeyTest[j][k]), "Marshal EvaluationKey element [%d][%d]", j, k)
			}
		}
	})

	t.Run(testString("SwitchingKey"), func(t *testing.T) {

		skOut := params.kgen.GenSecretKey()

		switchingKey := params.kgen.GenSwitchingKey(params.sk, skOut)
		data, err := switchingKey.MarshalBinary()
		require.NoError(t, err)

		resSwitchingKey := new(SwitchingKey)
		err = resSwitchingKey.UnmarshalBinary(data)
		require.NoError(t, err)

		evakeyWant := switchingKey.evakey
		evakeyTest := resSwitchingKey.evakey

		for j := range evakeyWant {
			for k := range evakeyWant[j] {
				require.True(t, contextQP.Equal(evakeyWant[j][k], evakeyTest[j][k]))
			}
		}
	})

	t.Run(testString("RotationKey"), func(t *testing.T) {

		rotationKey := NewRotationKeys()

		params.kgen.GenRot(Conjugate, params.sk, 0, rotationKey)
		params.kgen.GenRot(RotationLeft, params.sk, 1, rotationKey)
		params.kgen.GenRot(RotationLeft, params.sk, 2, rotationKey)
		params.kgen.GenRot(RotationRight, params.sk, 3, rotationKey)
		params.kgen.GenRot(RotationRight, params.sk, 5, rotationKey)

		data, err := rotationKey.MarshalBinary()
		require.NoError(t, err)

		resRotationKey := new(RotationKeys)
		err = resRotationKey.UnmarshalBinary(data)
		require.NoError(t, err)

		for i := uint64(1); i < params.ckkscontext.n>>1; i++ {

			if rotationKey.evakeyRotColLeft[i] != nil {

				evakeyWant := rotationKey.evakeyRotColLeft[i].evakey
				evakeyTest := resRotationKey.evakeyRotColLeft[i].evakey

				evakeyNTTIndexWant := rotationKey.permuteNTTLeftIndex[i]
				evakeyNTTIndexTest := resRotationKey.permuteNTTLeftIndex[i]

				require.True(t, utils.EqualSliceUint64(evakeyNTTIndexWant, evakeyNTTIndexTest))

				for j := range evakeyWant {
					for k := range evakeyWant[j] {
						require.Truef(t, contextQP.Equal(evakeyWant[j][k], evakeyTest[j][k]), "Marshal RotationKey RotateLeft %d element [%d][%d]", i, j, k)
					}
				}
			}

			if rotationKey.evakeyRotColRight[i] != nil {

				evakeyWant := rotationKey.evakeyRotColRight[i].evakey
				evakeyTest := resRotationKey.evakeyRotColRight[i].evakey

				evakeyNTTIndexWant := rotationKey.permuteNTTRightIndex[i]
				evakeyNTTIndexTest := resRotationKey.permuteNTTRightIndex[i]

				require.True(t, utils.EqualSliceUint64(evakeyNTTIndexWant, evakeyNTTIndexTest))

				for j := range evakeyWant {
					for k := range evakeyWant[j] {
						require.Truef(t, contextQP.Equal(evakeyWant[j][k], evakeyTest[j][k]), "Marshal RotationKey RotateRight %d element [%d][%d]", i, j, k)
					}
				}
			}
		}

		if rotationKey.evakeyConjugate != nil {

			evakeyWant := rotationKey.evakeyConjugate.evakey
			evakeyTest := resRotationKey.evakeyConjugate.evakey

			evakeyNTTIndexWant := rotationKey.permuteNTTConjugateIndex
			evakeyNTTIndexTest := resRotationKey.permuteNTTConjugateIndex

			require.True(t, utils.EqualSliceUint64(evakeyNTTIndexWant, evakeyNTTIndexTest))

			for j := range evakeyWant {
				for k := range evakeyWant[j] {
					require.Truef(t, contextQP.Equal(evakeyWant[j][k], evakeyTest[j][k]), "Marshal RotationKey RotateRow element [%d][%d]", j, k)
				}
			}
		}
	})

}
