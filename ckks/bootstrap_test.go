package ckks

import (
	"math"
	"math/cmplx"
	"runtime"
	"testing"

	"github.com/ldsec/lattigo/v2/ckks/bettersine"
	"github.com/ldsec/lattigo/v2/utils"
)

func TestBootstrap(t *testing.T) {

	if !*testBootstrapping {
		t.Skip("skipping bootstrapping test")
	}

	if runtime.GOARCH == "wasm" {
		t.Skip("skipping bootstrapping tests for GOARCH=wasm")
	}

	var testContext = new(testParams)

	paramSet := uint64(0)

	bootstrapParams := DefaultBootstrapParams[paramSet : paramSet+1]

	for paramSet := range bootstrapParams {

		btpParams := bootstrapParams[paramSet]

		// Insecure params for fast testing only
		if !*flagLongTest {
			btpParams.LogN = 14
			btpParams.LogSlots = 10
		}

		params, err := btpParams.Params()
		if err != nil {
			panic(err)
		}

		if testContext, err = genTestParams(params, btpParams.H); err != nil {
			panic(err)
		}

		for _, testSet := range []func(testContext *testParams, btpParams *BootstrappingParameters, t *testing.T){
			testEvalSine,
			testCoeffsToSlots,
			testSlotsToCoeffs,
			testbootstrap,
		} {
			testSet(testContext, btpParams, t)
			runtime.GC()
		}
	}
}

func testEvalSine(testContext *testParams, btpParams *BootstrappingParameters, t *testing.T) {

	t.Run(testString(testContext, "Sin/"), func(t *testing.T) {

		var err error

		eval := testContext.evaluator

		DefaultScale := testContext.params.scale

		SineScale := btpParams.SineEvalModuli.ScalingFactor

		testContext.params.scale = SineScale
		eval.(*evaluator).scale = SineScale

		deg := 127
		K := float64(15)

		values, _, ciphertext := newTestVectorsSineBootstrapp(testContext, btpParams, testContext.encryptorSk, -K+1, K-1, t)
		eval.DropLevel(ciphertext, btpParams.CtSDepth(true)-1)

		cheby := Approximate(sin2pi2pi, -complex(K, 0), complex(K, 0), deg)

		for i := range values {
			values[i] = sin2pi2pi(values[i])
		}

		eval.MultByConst(ciphertext, 2/(cheby.b-cheby.a), ciphertext)
		eval.AddConst(ciphertext, (-cheby.a-cheby.b)/(cheby.b-cheby.a), ciphertext)
		eval.Rescale(ciphertext, eval.(*evaluator).scale, ciphertext)

		if ciphertext, err = eval.EvaluateCheby(ciphertext, cheby, ciphertext.Scale()); err != nil {
			t.Error(err)
		}

		verifyTestVectors(testContext, testContext.decryptor, values, ciphertext, testContext.params.LogSlots(), 0, t)

		testContext.params.scale = DefaultScale
		eval.(*evaluator).scale = DefaultScale
	})

	t.Run(testString(testContext, "Cos1/"), func(t *testing.T) {

		var err error

		eval := testContext.evaluator

		DefaultScale := testContext.params.scale

		SineScale := btpParams.SineEvalModuli.ScalingFactor

		testContext.params.scale = SineScale
		eval.(*evaluator).scale = SineScale

		K := 25
		deg := 63
		dev := btpParams.MessageRatio
		scNum := 2

		scFac := complex(float64(int(1<<scNum)), 0)

		values, _, ciphertext := newTestVectorsSineBootstrapp(testContext, btpParams, testContext.encryptorSk, float64(-K+1), float64(K-1), t)
		eval.DropLevel(ciphertext, btpParams.CtSDepth(true)-1)

		cheby := new(ChebyshevInterpolation)
		cheby.coeffs = bettersine.Approximate(K, deg, dev, scNum)
		cheby.maxDeg = cheby.Degree()
		cheby.a = complex(float64(-K), 0) / scFac
		cheby.b = complex(float64(K), 0) / scFac
		cheby.lead = true

		var sqrt2pi float64
		if btpParams.ArcSineDeg > 0 {
			sqrt2pi = math.Pow(1, 1.0/real(scFac))
		} else {
			sqrt2pi = math.Pow(0.15915494309189535, 1.0/real(scFac))
		}

		for i := range cheby.coeffs {
			cheby.coeffs[i] *= complex(sqrt2pi, 0)
		}

		verifyTestVectors(testContext, testContext.decryptor, values, ciphertext, testContext.params.LogSlots(), 0, t)

		for i := range values {

			values[i] = cmplx.Cos(6.283185307179586 * (1 / scFac) * (values[i] - 0.25))

			for j := 0; j < scNum; j++ {
				values[i] = 2*values[i]*values[i] - 1
			}

			if btpParams.ArcSineDeg == 0 {
				values[i] /= 6.283185307179586
			}
		}

		eval.AddConst(ciphertext, -0.25, ciphertext)

		eval.MultByConst(ciphertext, 2/((cheby.b-cheby.a)*scFac), ciphertext)
		eval.AddConst(ciphertext, (-cheby.a-cheby.b)/(cheby.b-cheby.a), ciphertext)
		eval.Rescale(ciphertext, eval.(*evaluator).scale, ciphertext)

		if ciphertext, err = eval.EvaluateCheby(ciphertext, cheby, ciphertext.Scale()); err != nil {
			t.Error(err)
		}

		for i := 0; i < scNum; i++ {
			sqrt2pi *= sqrt2pi
			eval.MulRelin(ciphertext, ciphertext, ciphertext)
			eval.Add(ciphertext, ciphertext, ciphertext)
			eval.AddConst(ciphertext, -sqrt2pi, ciphertext)
			eval.Rescale(ciphertext, eval.(*evaluator).scale, ciphertext)
		}

		verifyTestVectors(testContext, testContext.decryptor, values, ciphertext, testContext.params.LogSlots(), 0, t)

		testContext.params.scale = DefaultScale
		eval.(*evaluator).scale = DefaultScale

	})

	t.Run(testString(testContext, "Cos2/"), func(t *testing.T) {

		if len(btpParams.SineEvalModuli.Qi) < 12 {
			t.Skip()
		}

		var err error

		eval := testContext.evaluator

		DefaultScale := testContext.params.scale

		SineScale := btpParams.SineEvalModuli.ScalingFactor

		testContext.params.scale = SineScale
		eval.(*evaluator).scale = SineScale

		K := 325
		deg := 255
		scNum := 4

		scFac := complex(float64(int(1<<scNum)), 0)

		values, _, ciphertext := newTestVectorsSineBootstrapp(testContext, btpParams, testContext.encryptorSk, float64(-K+1), float64(K-1), t)
		eval.DropLevel(ciphertext, btpParams.CtSDepth(true)-1)

		cheby := Approximate(cos2pi, -complex(float64(K), 0)/scFac, complex(float64(K), 0)/scFac, deg)

		sqrt2pi := math.Pow(0.15915494309189535, 1.0/real(scFac))

		for i := range cheby.coeffs {
			cheby.coeffs[i] *= complex(sqrt2pi, 0)
		}

		for i := range values {

			values[i] = cmplx.Cos(6.283185307179586 * (1 / scFac) * (values[i] - 0.25))

			for j := 0; j < scNum; j++ {
				values[i] = 2*values[i]*values[i] - 1
			}

			values[i] /= 6.283185307179586
		}

		testContext.evaluator.AddConst(ciphertext, -0.25, ciphertext)

		eval.MultByConst(ciphertext, 2/((cheby.b-cheby.a)*scFac), ciphertext)
		eval.AddConst(ciphertext, (-cheby.a-cheby.b)/(cheby.b-cheby.a), ciphertext)
		eval.Rescale(ciphertext, eval.(*evaluator).scale, ciphertext)

		if ciphertext, err = eval.EvaluateCheby(ciphertext, cheby, ciphertext.Scale()); err != nil {
			t.Error(err)
		}

		for i := 0; i < scNum; i++ {
			sqrt2pi *= sqrt2pi
			eval.MulRelin(ciphertext, ciphertext, ciphertext)
			eval.Add(ciphertext, ciphertext, ciphertext)
			eval.AddConst(ciphertext, -sqrt2pi, ciphertext)
			eval.Rescale(ciphertext, eval.(*evaluator).scale, ciphertext)
		}

		verifyTestVectors(testContext, testContext.decryptor, values, ciphertext, testContext.params.LogSlots(), 0, t)

		testContext.params.scale = DefaultScale
		eval.(*evaluator).scale = DefaultScale

	})
}

func testCoeffsToSlots(testContext *testParams, btpParams *BootstrappingParameters, t *testing.T) {
	t.Run(testString(testContext, "CoeffsToSlots/"), func(t *testing.T) {

		params := testContext.params

		// Generates the encoding matrices
		CoeffsToSlotMatrices := btpParams.GenCoeffsToSlotsMatrix(1.0, testContext.encoder)

		rotations := []int{}

		// Compute what rotations are needed for each matrix
		for i := range CoeffsToSlotMatrices {
			rotations = AddMatrixRotToList(CoeffsToSlotMatrices[i], rotations, params.Slots(), false)
		}

		// rotation for repacking sparse plaintexts
		if params.LogSlots() < params.LogN()-1 {
			rotations = append(rotations, params.Slots())
		}

		// Generates the rotation keys
		rotKey := testContext.kgen.GenRotationKeysForRotations(rotations, true, testContext.sk)

		// Generates a random test vectors
		values := make([]complex128, params.Slots())
		for i := range values {
			values[i] = complex(float64(i+1), 0)
		}

		// Encodes and encrypts the test vector
		plaintext := NewPlaintext(params, params.MaxLevel(), params.Scale())
		testContext.encoder.Encode(plaintext, values, params.logSlots)
		ciphertext := testContext.encryptorPk.EncryptNew(plaintext)

		// Creates an evaluator with the rotation keys
		eval := testContext.evaluator.WithKey(EvaluationKey{testContext.rlk, rotKey})

		// Applies the homomorphic DFT
		ct0, ct1 := CoeffsToSlots(ciphertext, CoeffsToSlotMatrices, eval)

		// Applies the same on the plaintext
		encoder := testContext.encoder

		// Data is not bit-reversed
		//sliceBitReverseInPlaceComplex128(values, params.Slots())
		invfft(values, params.Slots(), encoder.(*encoderComplex128).m, encoder.(*encoderComplex128).rotGroup, encoder.(*encoderComplex128).roots)
		sliceBitReverseInPlaceComplex128(values, params.Slots())

		// Verify the output values, and switch depending on if the original plaintext was sparse or not
		if params.LogSlots() < params.LogN()-1 {
			logSlots := params.LogSlots() + 1
			// Split the real and imaginary parts, puts the real part in the first 2*slots/2 slots and the imaginary part in the last 2*slots/2 slots.
			valuesFloat := make([]complex128, 1<<logSlots)
			for i, idx, jdx := 0, 0, 1<<(logSlots-1); i < 1<<(logSlots-1); i, jdx, idx = i+1, jdx+1, idx+1 {
				valuesFloat[idx] = complex(real(values[i]), 0)
				valuesFloat[jdx] = complex(imag(values[i]), 0)
			}

			ct0.MulScale(float64(int(1 << logSlots)))

			valuesTest := testContext.encoder.DecodePublic(testContext.decryptor.DecryptNew(ct0), logSlots, 0)

			/*
				for i := range valuesFloat{
					fmt.Println(i, valuesFloat[i], valuesTest[i])
				}
				fmt.Println()
			*/

			verifyTestVectors(testContext, testContext.decryptor, valuesFloat, valuesTest, logSlots, 0, t)
		} else {
			logSlots := params.LogSlots()

			/*
				for i := 0; i < 1<<(logSlots); i++{
					fmt.Println(i, values[i])
				}
				fmt.Println()
			*/

			// Splits the real and imaginary parts into two different slices.
			valuesFloat0 := make([]complex128, 1<<logSlots)
			valuesFloat1 := make([]complex128, 1<<logSlots)
			for i, idx, jdx := 0, 0, 1<<logSlots; i < 1<<logSlots; i, jdx, idx = i+1, jdx+1, idx+1 {
				valuesFloat0[idx] = complex(real(values[i]), 0)
				valuesFloat1[idx] = complex(imag(values[i]), 0)
			}

			ct0.MulScale(2 * float64(int(1<<logSlots)))
			ct1.MulScale(2 * float64(int(1<<logSlots)))

			valuesTest0 := testContext.encoder.DecodePublic(testContext.decryptor.DecryptNew(ct0), logSlots, 0)
			valuesTest1 := testContext.encoder.DecodePublic(testContext.decryptor.DecryptNew(ct1), logSlots, 0)

			/*
				for i := 0 ; i < 1<<logSlots; i++ {
					fmt.Println(i, valuesTest0[i], valuesTest1[i])
				}
			*/

			verifyTestVectors(testContext, testContext.decryptor, valuesFloat0, valuesTest0, logSlots, 0, t)
			verifyTestVectors(testContext, testContext.decryptor, valuesFloat1, valuesTest1, logSlots, 0, t)
		}
	})
}

func testSlotsToCoeffs(testContext *testParams, btpParams *BootstrappingParameters, t *testing.T) {
	t.Run(testString(testContext, "SlotsToCoeffs/"), func(t *testing.T) {

		params := testContext.params

		// Generates the encoding matrices
		SlotsToCoeffsMatrix := btpParams.GenSlotsToCoeffsMatrix(1.0, testContext.encoder)

		rotations := []int{}

		// Compute what rotations are needed for each matrix
		for i := range SlotsToCoeffsMatrix {
			rotations = AddMatrixRotToList(SlotsToCoeffsMatrix[i], rotations, params.Slots(), (i == 0) && (params.LogSlots() < params.LogN()-1))
		}

		// rotation for repacking sparse plaintexts
		if params.LogSlots() < params.LogN()-1 {
			rotations = append(rotations, params.Slots())
		}

		// Generates the rotation keys
		rotKey := testContext.kgen.GenRotationKeysForRotations(rotations, true, testContext.sk)

		// Creates an evaluator with the rotation keys
		eval := testContext.evaluator.WithKey(EvaluationKey{testContext.rlk, rotKey})

		// Generates a random test vectors that simulates the encoding of a real vector
		values0 := make([]complex128, params.Slots())
		values1 := make([]complex128, params.Slots())
		for i := range values0 {
			values0[i] = complex(utils.RandFloat64(-1, 1), 0)
		}

		for i := range values1[1:] {
			values1[i+1] = -values0[len(values0)-i-1]
		}

		// If sparse, puts the second vector in the imaginary part of the first one
		if params.LogSlots() < params.LogN()-1 {
			for i := range values0 {
				values0[i] += complex(0, real(values1[i]))
			}
		}

		// Ouputs of the homomorphic FFT^-1 is bit-reversed
		sliceBitReverseInPlaceComplex128(values0, params.Slots())
		sliceBitReverseInPlaceComplex128(values1, params.Slots())

		// Encodes and encrypts the test vectors
		logSlots := params.LogSlots()
		if params.LogSlots() < params.LogN()-1 {
			logSlots++
		}
		encoder := testContext.encoder
		plaintext := NewPlaintext(params, params.MaxLevel(), params.Scale())
		encoder.Encode(plaintext, values0, logSlots)
		ct0 := testContext.encryptorPk.EncryptNew(plaintext)
		var ct1 *Ciphertext
		if params.LogSlots() == params.LogN()-1 {
			testContext.encoder.Encode(plaintext, values1, logSlots)
			ct1 = testContext.encryptorPk.EncryptNew(plaintext)
		}

		// Applies the homomorphic DFT
		res := SlotsToCoeffs(ct0, ct1, SlotsToCoeffsMatrix, eval)

		eval.Rotate(res, 1, res)

		// Applies the DFT on the plaintext
		// If not sparse, puts the second vector in the imaginary part of the first one
		if params.LogSlots() == params.LogN()-1 {
			for i := range values0 {
				values0[i] += complex(0, real(values1[i]))
			}
		}

		/*
			for i := range values0{
				fmt.Println(i, values0[i])
			}
			fmt.Println()
		*/

		sliceBitReverseInPlaceComplex128(values0, params.Slots())
		fft(values0, params.Slots(), encoder.(*encoderComplex128).m, encoder.(*encoderComplex128).rotGroup, encoder.(*encoderComplex128).roots)
		//sliceBitReverseInPlaceComplex128(values0, params.Slots())

		values0 = utils.RotateComplex128Slice(values0, 1)

		/*
			for i := range values0{
				fmt.Println(i, values0[i])
			}
			fmt.Println()
		*/

		valuesTest := testContext.encoder.DecodePublic(testContext.decryptor.DecryptNew(res), params.LogSlots(), 0)

		/*
			for i := range values0{
				fmt.Println(i, valuesTest[i], values0[i])
			}
			fmt.Println()
		*/

		verifyTestVectors(testContext, testContext.decryptor, values0, valuesTest, params.LogSlots(), 0, t)

	})
}

func testbootstrap(testContext *testParams, btpParams *BootstrappingParameters, t *testing.T) {

	t.Run(testString(testContext, "Bootstrapping/FullCircuit/"), func(t *testing.T) {

		params := testContext.params

		rotations := testContext.kgen.GenRotationIndexesForBootstrapping(testContext.params.logSlots, btpParams)
		rotkeys := testContext.kgen.GenRotationKeysForRotations(rotations, true, testContext.sk)
		btpKey := BootstrappingKey{testContext.rlk, rotkeys}

		btp, err := NewBootstrapper(testContext.params, btpParams, btpKey)
		if err != nil {
			panic(err)
		}

		values := make([]complex128, 1<<params.logSlots)
		for i := range values {
			values[i] = utils.RandComplex128(-1, 1)
		}

		values[0] = complex(0.9238795325112867, 0.3826834323650898)
		values[1] = complex(0.9238795325112867, 0.3826834323650898)
		if 1<<params.logSlots > 2 {
			values[2] = complex(0.9238795325112867, 0.3826834323650898)
			values[3] = complex(0.9238795325112867, 0.3826834323650898)
		}

		plaintext := NewPlaintext(params, params.MaxLevel(), params.Scale())
		testContext.encoder.Encode(plaintext, values, params.logSlots)

		ciphertext := testContext.encryptorPk.EncryptNew(plaintext)

		eval := testContext.evaluator
		for ciphertext.Level() != 0 {
			eval.DropLevel(ciphertext, 1)
		}

		for i := 0; i < 1; i++ {

			ciphertext = btp.Bootstrapp(ciphertext)
			//testContext.evaluator.SetScale(ciphertext, testContext.params.scale)
			verifyTestVectors(testContext, testContext.decryptor, values, ciphertext, testContext.params.LogSlots(), 0, t)
		}

	})
}

func newTestVectorsSineBootstrapp(testContext *testParams, btpParams *BootstrappingParameters, encryptor Encryptor, a, b float64, t *testing.T) (values []complex128, plaintext *Plaintext, ciphertext *Ciphertext) {

	logSlots := testContext.params.LogSlots()

	values = make([]complex128, 1<<logSlots)

	ratio := btpParams.MessageRatio

	for i := uint64(0); i < 1<<logSlots; i++ {
		values[i] = complex(math.Round(utils.RandFloat64(a, b))+utils.RandFloat64(-1, 1)/ratio, 0)
	}

	plaintext = NewPlaintext(testContext.params, testContext.params.MaxLevel(), testContext.params.Scale())

	testContext.encoder.EncodeNTT(plaintext, values, logSlots)

	if encryptor != nil {
		ciphertext = encryptor.EncryptNew(plaintext)
	}

	return values, plaintext, ciphertext
}
