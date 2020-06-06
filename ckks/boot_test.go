package ckks

import (
	"fmt"
	"math"
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

	bootparams := BootstrappParams[1]

	parameters := &bootparams.Parameters

	bootparams.Gen()

	/*
		for i, qi := range bootparams.Qi{
			fmt.Printf("%d : 0x%x\n", i, qi)
		}

		for i, pj := range bootparams.Pi{
			fmt.Printf("%d : 0x%x\n", i, pj)
		}
	*/

	params := genCkksParams(parameters)

	slots := uint64(1 << bootparams.Parameters.LogSlots)

	rlk := params.kgen.GenRelinKey(params.sk)

	t.Run(testString("ChebySin/", parameters), func(t *testing.T) {

		eval := params.evaluator

		DefaultScale := parameters.Scale

		params.params.Scale = SineScale
		parameters.Scale = SineScale
		eval.(*evaluator).ckksContext.scale = SineScale

		deg := 119
		K := float64(15)

		values, _, ciphertext := newTestVectorsSineBoot(params, params.encryptorSk, -K+1, K-1, t)
		eval.DropLevel(ciphertext, uint64(len(bootparams.CtSLevel))-1)

		cheby := Approximate(sin2pi2pi, -complex(K, 0), complex(K, 0), deg)

		for i := range values {
			values[i] = sin2pi2pi(values[i])
		}

		fmt.Println(ciphertext.Level() - 1)
		start := time.Now()
		ciphertext = params.evaluator.EvaluateCheby(ciphertext, cheby, rlk)
		fmt.Printf("Elapsed : %s \n", time.Since(start))
		fmt.Println(ciphertext.Level())

		verifyTestVectors(parameters, params.encoder, params.decryptor, values, ciphertext, t)

		params.params.Scale = DefaultScale
		parameters.Scale = DefaultScale
		eval.(*evaluator).params.Scale = DefaultScale
	})

	t.Run(testString("ChebyCos/", parameters), func(t *testing.T) {

		eval := params.evaluator

		DefaultScale := parameters.Scale

		params.params.Scale = SineScale
		parameters.Scale = SineScale
		eval.(*evaluator).ckksContext.scale = SineScale

		K := 16
		deg := 40
		dev := 10
		sc_num := 2

		sc_fac := complex(float64(int(1<<sc_num)), 0)

		values, _, ciphertext := newTestVectorsSineBoot(params, params.encryptorSk, float64(-K+1), float64(K-1), t)
		eval.DropLevel(ciphertext, uint64(len(bootparams.CtSLevel))-1)

		cheby := new(ChebyshevInterpolation)
		cheby.coeffs = bettersine.Approximate(K, deg, dev, sc_num)
		cheby.maxDeg = cheby.degree()
		cheby.a = complex(float64(-K), 0) / sc_fac
		cheby.b = complex(float64(K), 0) / sc_fac

		sqrt2pi := math.Pow(0.15915494309189535, 1.0/real(sc_fac))

		for i := range cheby.coeffs {
			cheby.coeffs[i] *= complex(sqrt2pi, 0)
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
		ciphertext = params.evaluator.EvaluateChebySpecial(ciphertext, sc_fac, cheby, rlk)
		fmt.Println(ciphertext.Level(), ciphertext.Scale())

		for i := 0; i < sc_num; i++ {
			sqrt2pi *= sqrt2pi
			params.evaluator.MulRelin(ciphertext, ciphertext, rlk, ciphertext)
			params.evaluator.Add(ciphertext, ciphertext, ciphertext)
			params.evaluator.AddConst(ciphertext, -sqrt2pi, ciphertext)
			params.evaluator.Rescale(ciphertext, parameters.Scale, ciphertext)
		}

		fmt.Printf("Elapsed : %s \n", time.Since(start))
		fmt.Println(ciphertext.Level(), ciphertext.Scale())
		verifyTestVectors(parameters, params.encoder, params.decryptor, values, ciphertext, t)

		params.params.Scale = DefaultScale
		parameters.Scale = DefaultScale
		eval.(*evaluator).params.Scale = DefaultScale

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
			values[i] = complex(RandomFloat(-1, 1), RandomFloat(-1, 1))
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

			verifyTestVectors(parameters, params.encoder, params.decryptor, values, ciphertext, t)
			//fmt.Println()
		}

	})
}
