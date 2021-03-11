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

		params, err := btpParams.Params()
		if err != nil {
			panic(err)
		}

		// Insecure params for fast testing only
		if !*flagLongTest {
			params.logN = 14
			params.logSlots = 13
		}

		if testContext, err = genTestParams(params, btpParams.H); err != nil {
			panic(err)
		}

		for _, testSet := range []func(testContext *testParams, btpParams *BootstrappingParameters, t *testing.T){
			testSin,
			testCos1,
			testCos2,
			testbootstrap,
		} {
			testSet(testContext, btpParams, t)
			runtime.GC()
		}
	}
}

func testSin(testContext *testParams, btpParams *BootstrappingParameters, t *testing.T) {
	t.Run(testString(testContext, "Sin/"), func(t *testing.T) {

		var err error

		eval := testContext.evaluator

		DefaultScale := testContext.params.scale

		SineScale := btpParams.SineEvalModuli.ScalingFactor

		testContext.params.scale = SineScale
		eval.(*evaluator).scale = SineScale

		deg := 127
		K := float64(15)

		values, _, ciphertext := newTestVectorsSineBootstrapp(testContext, testContext.encryptorSk, -K+1, K-1, t)
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

		verifyTestVectors(testContext, testContext.decryptor, values, ciphertext, t, 0)

		testContext.params.scale = DefaultScale
		eval.(*evaluator).scale = DefaultScale
	})
}

func testCos1(testContext *testParams, btpParams *BootstrappingParameters, t *testing.T) {
	t.Run(testString(testContext, "Cos1/"), func(t *testing.T) {

		var err error

		eval := testContext.evaluator

		DefaultScale := testContext.params.scale

		SineScale := btpParams.SineEvalModuli.ScalingFactor

		testContext.params.scale = SineScale
		eval.(*evaluator).scale = SineScale

		K := 25
		deg := 62
		dev := float64(testContext.params.qi[0]) / DefaultScale
		scNum := 2

		scFac := complex(float64(int(1<<scNum)), 0)

		values, _, ciphertext := newTestVectorsSineBootstrapp(testContext, testContext.encryptorSk, float64(-K+1), float64(K-1), t)
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

		
		verifyTestVectors(testContext, testContext.decryptor, values, ciphertext, t, 0)
		
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

		verifyTestVectors(testContext, testContext.decryptor, values, ciphertext, t, 0)

	})
}

func testCos2(testContext *testParams, btpParams *BootstrappingParameters, t *testing.T) {
	t.Run(testString(testContext, "Cos2/"), func(t *testing.T) {

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

		values, _, ciphertext := newTestVectorsSineBootstrapp(testContext, testContext.encryptorSk, float64(-K+1), float64(K-1), t)
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

		verifyTestVectors(testContext, testContext.decryptor, values, ciphertext, t, 0)

		testContext.params.scale = DefaultScale
		eval.(*evaluator).scale = DefaultScale

	})
}

func testbootstrap(testContext *testParams, btpParams *BootstrappingParameters, t *testing.T) {
	t.Run(testString(testContext, "Bootstrap/"), func(t *testing.T) {

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
			verifyTestVectors(testContext, testContext.decryptor, values, ciphertext, t, 0)
		}

	})
}

func newTestVectorsSineBootstrapp(testContext *testParams, encryptor Encryptor, a, b float64, t *testing.T) (values []complex128, plaintext *Plaintext, ciphertext *Ciphertext) {

	logSlots := testContext.params.LogSlots()

	values = make([]complex128, 1<<logSlots)

	ratio := float64(testContext.params.Qi()[0]) / float64(35184372088832)

	for i := uint64(0); i < 1<<logSlots; i++ {
		values[i] = complex(math.Round(utils.RandFloat64(a, b))+utils.RandFloat64(-1, 1)/ratio, 0)
	}

	values[0] = complex(3+27/ratio, 0)

	plaintext = NewPlaintext(testContext.params, testContext.params.MaxLevel(), testContext.params.Scale())

	testContext.encoder.EncodeNTT(plaintext, values, logSlots)

	if encryptor != nil {
		ciphertext = encryptor.EncryptNew(plaintext)
	}

	return values, plaintext, ciphertext
}
