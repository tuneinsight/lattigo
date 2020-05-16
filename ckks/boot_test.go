package ckks

import (
	"fmt"
	"math/cmplx"
	"math/rand"
	"testing"
	"time"

	"github.com/ldsec/lattigo/ckks/bettersine"
)

func TestBootstrapp(t *testing.T) {

	rand.Seed(time.Now().UnixNano())

	var SineScale float64

	SineScale = 1 << 55

	bootparams := BootstrappParams[3]

	parameters := &bootparams.Parameters

	bootparams.Gen()

	//for i, qi := range bootparams.Qi{
	//	fmt.Println(i, qi)
	//}

	params := genCkksParams(parameters)

	slots := uint64(1 << bootparams.Parameters.LogSlots)

	rlk := params.kgen.GenRelinKey(params.sk)

	t.Run(testString("ChebySin/", parameters), func(t *testing.T) {

		DefaultScale := parameters.Scale

		params.params.Scale = SineScale
		parameters.Scale = SineScale

		evaluator := NewEvaluator(parameters)

		deg := 123
		K := float64(15)

		values, _, ciphertext := newTestVectorsSineBoot(params, params.encryptorSk, -K+1, K-1, t)
		evaluator.DropLevel(ciphertext, uint64(len(bootparams.CtSLevel)))

		cheby := Approximate(sin2pi2pi, -complex(K, 0), complex(K, 0), deg)

		for i := range values {
			values[i] = sin2pi2pi(values[i])
		}

		//fmt.Println(ciphertext.Level())
		start := time.Now()
		ciphertext = params.evaluator.EvaluateChebyFast(ciphertext, cheby, rlk)
		fmt.Printf("Elapsed : %s \n", time.Since(start))
		//fmt.Println(ciphertext.Level())

		verifyTestVectors(params, params.decryptor, values, ciphertext, t)
		t.Fail()

		params.params.Scale = DefaultScale
		parameters.Scale = DefaultScale
	})

	t.Run(testString("ChebyCos/", parameters), func(t *testing.T) {

		DefaultScale := parameters.Scale

		params.params.Scale = SineScale
		parameters.Scale = SineScale

		evaluator := NewEvaluator(parameters)

		K := 16
		deg := 44
		dev := 10
		sc_num := 2

		sc_fac := complex(float64(int(1<<sc_num)), 0)

		values, _, ciphertext := newTestVectorsSineBoot(params, params.encryptorSk, float64(-K+1), float64(K-1), t)
		evaluator.DropLevel(ciphertext, uint64(len(bootparams.CtSLevel)))

		cheby := new(ChebyshevInterpolation)
		cheby.coeffs = bettersine.Approximate(K, deg, dev, sc_num)
		cheby.maxDeg = uint64(len(cheby.coeffs) - 1)
		cheby.a = complex(float64(-K), 0) / sc_fac
		cheby.b = complex(float64(K), 0) / sc_fac

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

		if sc_num == 2 {
			for i := range cheby.coeffs {
				cheby.coeffs[i] *= 0.7511255444649425
			}
		}

		if sc_num == 3 {
			for i := range cheby.coeffs {
				cheby.coeffs[i] *= 0.8666749935615672
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

		fmt.Println(ciphertext.Level(), ciphertext.Scale())
		start := time.Now()
		ciphertext = params.evaluator.EvaluateChebyFastSpecial(ciphertext, sc_fac, cheby, rlk)
		fmt.Println(ciphertext.Level(), ciphertext.Scale())

		if sc_num == 1 {
			params.evaluator.MulRelin(ciphertext, ciphertext, rlk, ciphertext)
			params.evaluator.AddConst(ciphertext, -1.0/6.283185307179586, ciphertext)
			params.evaluator.Rescale(ciphertext, parameters.Scale, ciphertext)
		}

		if sc_num == 2 {

			params.evaluator.MulRelin(ciphertext, ciphertext, rlk, ciphertext)
			params.evaluator.Rescale(ciphertext, parameters.Scale, ciphertext)
			y := params.evaluator.AddConstNew(ciphertext, -0.5641895835477563)

			params.evaluator.MulRelin(ciphertext, y, rlk, ciphertext)
			params.evaluator.MultByConst(ciphertext, 4, ciphertext)
			params.evaluator.AddConst(ciphertext, 1.0/6.283185307179586, ciphertext)

			params.evaluator.Rescale(ciphertext, parameters.Scale, ciphertext)
		}

		if sc_num == 3 {

			// r = 16*(y4 * (a*y4 - b*y2 + c) - d*y2) + 1/(2*pi)

			a := 4.0
			b := -6.00900435571954
			c := 2.8209479177387813
			d := -0.42377720812375763

			y2 := evaluator.MulRelinNew(ciphertext, ciphertext, rlk)
			evaluator.Rescale(y2, parameters.Scale, y2)

			y4 := evaluator.MulRelinNew(y2, y2, rlk)

			ciphertext = y4.CopyNew().Ciphertext()

			evaluator.MultByConst(ciphertext, a, ciphertext)
			evaluator.MultByConstAndAdd(y2, b, ciphertext)
			evaluator.AddConst(ciphertext, c, ciphertext)
			evaluator.Rescale(ciphertext, parameters.Scale, ciphertext)
			evaluator.Rescale(y4, parameters.Scale, y4)

			evaluator.MulRelin(ciphertext, y4, rlk, ciphertext)
			evaluator.MultByConstAndAdd(y2, d, ciphertext)

			evaluator.MultByConst(ciphertext, 16, ciphertext)

			evaluator.AddConst(ciphertext, 1.0/6.283185307179586, ciphertext)

			evaluator.Rescale(ciphertext, parameters.Scale, ciphertext)

		}
		fmt.Printf("Elapsed : %s \n", time.Since(start))
		fmt.Println(ciphertext.Level(), ciphertext.Scale())
		verifyTestVectors(params, params.decryptor, values, ciphertext, t)

		t.Fail()

		params.params.Scale = DefaultScale
		parameters.Scale = DefaultScale

	})

	t.Run(testString("Bootstrapp/", parameters), func(t *testing.T) {

		bootcontext := NewBootContext(bootparams)
		bootcontext.GenBootKeys(params.sk)

		rlk, rotkey := bootcontext.ExportKeys()

		bootcontext.ImportKeys(rlk, rotkey)
		if err := bootcontext.CheckKeys(); err != nil {
			panic(err)
		}

		values := make([]complex128, slots)
		for i := range values {
			values[i] = complex(randomFloat(-1, 1), randomFloat(-1, 1))
		}

		values[0] = complex(0.9238795325112867, 0.3826834323650898)
		values[1] = complex(0.9238795325112867, 0.3826834323650898)
		if slots > 2 {
			values[2] = complex(0.9238795325112867, 0.3826834323650898)
			values[3] = complex(0.9238795325112867, 0.3826834323650898)
		}

		plaintext := NewPlaintext(parameters, parameters.MaxLevel, parameters.Scale)
		params.encoder.Encode(plaintext, values, slots)

		ciphertext := params.encryptorPk.EncryptNew(plaintext)

		for i := 0; i < 1; i++ {

			ciphertext = bootcontext.Bootstrapp(ciphertext)

			//params.evaluator.SetScale(ciphertext, parameters.Scale)

			verifyTestVectors(params, params.decryptor, values, ciphertext, t)
			//fmt.Println()
		}
		t.Fail()

	})
}
