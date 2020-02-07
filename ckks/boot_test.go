package ckks

import (
	"log"
	"math/cmplx"
	"math/rand"
	"testing"
	"time"
	"fmt"

	"github.com/ldsec/lattigo/ckks/bettersine"
)

func TestBootstrapp(t *testing.T) {

	rand.Seed(time.Now().UnixNano())

	var DefaultScale, LTScale, SineScale float64

	DefaultScale = 1 << 40
	LTScale = 1 << 45
	_ = LTScale
	SineScale = 1 << 55

	logSlots := uint64(13)
	ctsDepth := uint64(3)
	stcDepth := uint64(3)

	bootParams := new(Parameters)
	bootParams.LogN = 16
	bootParams.LogSlots = logSlots
	bootParams.Scale = DefaultScale
	bootParams.LogQi = []uint64{55, 40, 40, 40, 40, 40, 40, 40, 40, 40, 45, 45, 45, 55, 55, 55, 55, 55, 55, 55, 55, 55, 45, 45, 45}
	bootParams.LogPi = []uint64{55, 55, 55, 55}
	bootParams.Sigma = 3.2

	bootParams.GenFromLogModuli()

	params := genCkksParams(bootParams)

	slots := uint64(1 << logSlots)

	rlk := params.kgen.GenRelinKey(params.sk)

	t.Run(testString("SineOriginal/", bootParams), func(t *testing.T) {

		params.params.Scale = SineScale

		evaluator := NewEvaluator(bootParams)

		deg := 127
		K := float64(15)

		values, _, ciphertext := newTestVectorsSineBoot(params, params.encryptorSk, -K+1, K-1, t)
		evaluator.DropLevel(ciphertext, ctsDepth)

		cheby := Approximate(sin2pi2pi, -complex(K, 0), complex(K, 0), deg)

		for i := range values {
			values[i] = sin2pi2pi(values[i])
		}

		ciphertext = params.evaluator.EvaluateChebyFast(ciphertext, cheby, rlk)

		verifyTestVectors(params, params.decryptor, values, ciphertext, t)

		params.params.Scale = DefaultScale
	})

	t.Run(testString("SineFaster/", bootParams), func(t *testing.T) {

		params.params.Scale = SineScale

		evaluator := NewEvaluator(bootParams)

		K := 12
		deg := 137
		dev := 10
		sc_num := 0

		sc_fac := complex(float64(int(1<<sc_num)), 0)

		values, _, ciphertext := newTestVectorsSineBoot(params, params.encryptorSk, float64(-K+1), float64(K-1), t)
		evaluator.DropLevel(ciphertext, ctsDepth)

		cheby := new(ChebyshevInterpolation)
		cheby.coeffs = bettersine.Approximate(K, deg, dev, sc_num)
		cheby.maxDeg = uint64(len(cheby.coeffs)-1)
		cheby.a = complex(float64(-K), 0) / sc_fac
		cheby.b = complex(float64(K), 0) / sc_fac

		fmt.Println(len(cheby.coeffs)-1, cheby.coeffs[len(cheby.coeffs)-1], max(cheby.coeffs), min(cheby.coeffs))

		if sc_num == 0 {
			for i := range cheby.coeffs {
				cheby.coeffs[i] *= 0.15915494309189535
			}
		}
		
		if sc_num == 1 {
			for i := range cheby.coeffs {
				cheby.coeffs[i] *= 0.5641895835477563
			}
		}

		if sc_num == 2{
			for i := range cheby.coeffs {
				cheby.coeffs[i] *= 0.7511255444649425
			}
		}

		for i := range values {

			values[i] = cmplx.Cos(6.283185307179586 * (1 / sc_fac) * (values[i] - 0.25))

			for j := 0; j < sc_num; j++ {
				values[i] = 2*values[i]*values[i] - 1
			}

			values[i] /= 6.283185307179586
		}

		params.evaluator.AddConst(ciphertext, -0.25, ciphertext)

		fmt.Println(ciphertext.Level())
		ciphertext = params.evaluator.EvaluateChebyFastSpecial(ciphertext, sc_fac, cheby, rlk)
		fmt.Println(ciphertext.Level())
		
		

		
		if sc_num == 1 {
			params.evaluator.MulRelin(ciphertext, ciphertext, rlk, ciphertext)
			params.evaluator.AddConst(ciphertext, -1.0/6.283185307179586, ciphertext)
			params.evaluator.Rescale(ciphertext, params.params.Scale, ciphertext)
		}

		if sc_num == 2 {

			params.evaluator.MulRelin(ciphertext, ciphertext, rlk, ciphertext)
			params.evaluator.Rescale(ciphertext, params.params.Scale, ciphertext)
			y := params.evaluator.AddConstNew(ciphertext, -0.5641895835477563)

			params.evaluator.MulRelin(ciphertext, y, rlk, ciphertext)
			params.evaluator.MultByConst(ciphertext, 4, ciphertext)
			params.evaluator.AddConst(ciphertext, 1.0/6.283185307179586, ciphertext)

			params.evaluator.Rescale(ciphertext, params.params.Scale, ciphertext)
		}


		if sc_num == 3 {
			for i:= 0 ; i < sc_num; i++ {
				params.evaluator.MulRelin(ciphertext, ciphertext, rlk, ciphertext)
				params.evaluator.MultByConst(ciphertext, 2, ciphertext)
				params.evaluator.AddConst(ciphertext, -1, ciphertext)
				params.evaluator.Rescale(ciphertext, params.params.Scale, ciphertext)
			}

			params.evaluator.MultByConst(ciphertext, 1.0 / 6.283185307179586, ciphertext)
		}


		

		verifyTestVectors(params, params.decryptor, values, ciphertext, t)

		params.params.Scale = DefaultScale

	})

	t.Run(testString("Bootstrapp/", bootParams), func(t *testing.T) {

		var bootcontext *BootContext
		var err error

		if bootcontext, err = NewBootContext(bootParams, params.sk, ctsDepth, stcDepth); err != nil {
			log.Fatal(err)
		}

		values := make([]complex128, slots)
		for i := range values {
			values[i] = complex(randomFloat(-1, 1), 0)
		}

		values[0] = complex(0.516015, 0)
		values[1] = complex(0.772621, 0)
		if slots > 2 {
			values[2] = complex(0.939175, 0)
			values[3] = complex(0.345987, 0)
		}

		plaintext := NewPlaintext(bootParams, bootParams.MaxLevel(), bootParams.Scale)
		params.encoder.Encode(plaintext, values, slots)

		ciphertext := params.encryptorPk.EncryptNew(plaintext)

		for i := 0; i < 1; i++ {

			ciphertext = params.evaluator.Bootstrapp(ciphertext, bootcontext)

			//if err = evaluator.SetScale(ciphertext, params.Scale); err != nil {
			//	log.Fatal(err)
			//}

			verifyTestVectors(params, params.decryptor, values, ciphertext, t)
		}

	})

	t.Run(testString("BootstrappBetterSine/", bootParams), func(t *testing.T) {

		var bootcontext *BootContext
		var err error

		if bootcontext, err = NewBootContext(bootParams, params.sk, ctsDepth, stcDepth); err != nil {
			log.Fatal(err)
		}

		values := make([]complex128, slots)
		for i := range values {
			values[i] = complex(randomFloat(-1, 1), 0)
		}

		values[0] = complex(0.516015, 0)
		values[1] = complex(0.772621, 0)
		if slots > 2 {
			values[2] = complex(0.939175, 0)
			values[3] = complex(0.345987, 0)
		}

		plaintext := NewPlaintext(bootParams, bootParams.MaxLevel(), bootParams.Scale)
		params.encoder.Encode(plaintext, values, slots)

		ciphertext := params.encryptorPk.EncryptNew(plaintext)

		for i := 0; i < 1; i++ {

			ciphertext = params.evaluator.BootstrappBetterSine(ciphertext, bootcontext)

			//if err = evaluator.SetScale(ciphertext, params.Scale); err != nil {
			//	log.Fatal(err)
			//}

			verifyTestVectors(params, params.decryptor, values, ciphertext, t)
		}

	})
}
