package ckks

import (
	//"fmt"
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

	var shemeParams []*Parameters
	var bootstrappParams []*BootstrappParams

	if testing.Short() {
		shemeParams = DefaultBootstrappSchemeParamsShort
		bootstrappParams = DefaultBootstrappParamsShort
	} else {
		shemeParams = DefaultBootstrappSchemeParams
		bootstrappParams = DefaultBootstrappParams
	}

	for paramSet := range shemeParams {

		btpParams := bootstrappParams[paramSet]

		if err := genTestParams(shemeParams[paramSet], btpParams.H); err != nil {
			panic(err)
		}

		slots := params.params.Slots()

		t.Run(testString("ChebySin/"), func(t *testing.T) {

			eval := params.evaluator

			DefaultScale := params.params.scale

			params.params.scale = SineScale
			eval.(*evaluator).scale = SineScale

			deg := 127
			K := float64(15)

			values, _, ciphertext := newTestVectorsSineBoot(params.encryptorSk, -K+1, K-1, t)
			eval.DropLevel(ciphertext, uint64(len(btpParams.CtSLevel))-1)

			cheby := Approximate(sin2pi2pi, -complex(K, 0), complex(K, 0), deg)

			for i := range values {
				values[i] = sin2pi2pi(values[i])
			}

			//fmt.Println(ciphertext.Level() - 1)
			//start := time.Now()
			ciphertext = params.evaluator.EvaluateCheby(ciphertext, cheby, params.rlk)
			//fmt.Printf("Elapsed : %s \n", time.Since(start))
			//fmt.Println(ciphertext.Level())

			verifyTestVectors(params.decryptor, values, ciphertext, t)

			params.params.scale = DefaultScale
			eval.(*evaluator).scale = DefaultScale
		})

		t.Run(testString("ChebyCos/"), func(t *testing.T) {

			eval := params.evaluator

			DefaultScale := params.params.scale

			params.params.scale = SineScale
			eval.(*evaluator).scale = SineScale

			K := 26
			deg := 63
			dev := float64(params.params.qi[0]) / DefaultScale
			sc_num := 2

			sc_fac := complex(float64(int(1<<sc_num)), 0)

			values, _, ciphertext := newTestVectorsSineBoot(params.encryptorSk, float64(-K+1), float64(K-1), t)
			eval.DropLevel(ciphertext, uint64(len(btpParams.CtSLevel))-1)

			cheby := new(ChebyshevInterpolation)
			cheby.coeffs = bettersine.Approximate(K, deg, dev, sc_num)
			cheby.maxDeg = cheby.Degree()
			cheby.a = complex(float64(-K), 0) / sc_fac
			cheby.b = complex(float64(K), 0) / sc_fac
			cheby.lead = true

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

			//fmt.Println(ciphertext.Level(), ciphertext.Scale())
			//start := time.Now()
			ciphertext = params.evaluator.EvaluateChebySpecial(ciphertext, sc_fac, cheby, params.rlk)
			//fmt.Println(ciphertext.Level(), ciphertext.Scale())

			for i := 0; i < sc_num; i++ {
				sqrt2pi *= sqrt2pi
				params.evaluator.MulRelin(ciphertext, ciphertext, params.rlk, ciphertext)
				params.evaluator.Add(ciphertext, ciphertext, ciphertext)
				params.evaluator.AddConst(ciphertext, -sqrt2pi, ciphertext)
				params.evaluator.Rescale(ciphertext, eval.(*evaluator).scale, ciphertext)
			}

			//fmt.Printf("Elapsed : %s \n", time.Since(start))
			//fmt.Println(ciphertext.Level(), ciphertext.Scale())
			verifyTestVectors(params.decryptor, values, ciphertext, t)

			params.params.scale = DefaultScale
			eval.(*evaluator).scale = DefaultScale

		})

		t.Run(testString("Bootstrapp/"), func(t *testing.T) {

			btp, err := NewBootstrapper(params.params, btpParams)
			if err != nil {
				panic(err)
			}
			btp.GenKeys(params.sk)

			rlk, rotkey := btp.ExportKeys()

			btp.ImportKeys(rlk, rotkey)
			if err := btp.CheckKeys(); err != nil {
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

			plaintext := NewPlaintext(params.params, params.params.MaxLevel(), params.params.scale)
			params.encoder.Encode(plaintext, values, slots)

			ciphertext := params.encryptorPk.EncryptNew(plaintext)

			for i := 0; i < 1; i++ {

				ciphertext = btp.Bootstrapp(ciphertext)

				//params.evaluator.SetScale(ciphertext, params.params.scale)

				verifyTestVectors(params.decryptor, values, ciphertext, t)
			}

		})
	}
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
