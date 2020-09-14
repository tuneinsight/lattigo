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

	var err error
	var testContext = new(testParams)
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

		if testContext, err = genTestParams(shemeParams[paramSet], btpParams.H); err != nil {
			panic(err)
		}

		slots := testContext.params.Slots()

		t.Run(testString(testContext, "ChebySin/"), func(t *testing.T) {

			eval := testContext.evaluator

			DefaultScale := testContext.params.scale

			testContext.params.scale = SineScale
			eval.(*evaluator).scale = SineScale

			deg := 127
			K := float64(15)

			values, _, ciphertext := newTestVectorsSineBootstrapp(testContext, testContext.encryptorSk, -K+1, K-1, t)
			eval.DropLevel(ciphertext, uint64(len(btpParams.CtSLevel))-1)

			cheby := Approximate(sin2pi2pi, -complex(K, 0), complex(K, 0), deg)

			for i := range values {
				values[i] = sin2pi2pi(values[i])
			}

			//fmt.Println(ciphertext.Level() - 1)
			//start := time.Now()
			ciphertext = testContext.evaluator.EvaluateCheby(ciphertext, cheby, testContext.rlk)
			//fmt.Printf("Elapsed : %s \n", time.Since(start))
			//fmt.Println(ciphertext.Level())

			verifyTestVectors(testContext, testContext.decryptor, values, ciphertext, t)

			testContext.params.scale = DefaultScale
			eval.(*evaluator).scale = DefaultScale
		})

		t.Run(testString(testContext, "ChebyCos/"), func(t *testing.T) {

			eval := testContext.evaluator

			DefaultScale := testContext.params.scale

			testContext.params.scale = SineScale
			eval.(*evaluator).scale = SineScale

			K := 26
			deg := 63
			dev := float64(testContext.params.qi[0]) / DefaultScale
			scNum := 2

			scFac := complex(float64(int(1<<scNum)), 0)

			values, _, ciphertext := newTestVectorsSineBootstrapp(testContext, testContext.encryptorSk, float64(-K+1), float64(K-1), t)
			eval.DropLevel(ciphertext, uint64(len(btpParams.CtSLevel))-1)

			cheby := new(ChebyshevInterpolation)
			cheby.coeffs = bettersine.Approximate(K, deg, dev, scNum)
			cheby.maxDeg = cheby.Degree()
			cheby.a = complex(float64(-K), 0) / scFac
			cheby.b = complex(float64(K), 0) / scFac
			cheby.lead = true

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

			//fmt.Println(ciphertext.Level(), ciphertext.Scale())
			//start := time.Now()
			ciphertext = testContext.evaluator.EvaluateChebySpecial(ciphertext, scFac, cheby, testContext.rlk)
			//fmt.Println(ciphertext.Level(), ciphertext.Scale())

			for i := 0; i < scNum; i++ {
				sqrt2pi *= sqrt2pi
				testContext.evaluator.MulRelin(ciphertext, ciphertext, testContext.rlk, ciphertext)
				testContext.evaluator.Add(ciphertext, ciphertext, ciphertext)
				testContext.evaluator.AddConst(ciphertext, -sqrt2pi, ciphertext)
				testContext.evaluator.Rescale(ciphertext, eval.(*evaluator).scale, ciphertext)
			}

			//fmt.Printf("Elapsed : %s \n", time.Since(start))
			//fmt.Println(ciphertext.Level(), ciphertext.Scale())
			verifyTestVectors(testContext, testContext.decryptor, values, ciphertext, t)

			testContext.params.scale = DefaultScale
			eval.(*evaluator).scale = DefaultScale

		})

		t.Run(testString(testContext, "Bootstrapp/"), func(t *testing.T) {

			btpKey := testContext.kgen.GenBootstrappingKey(testContext.params.logSlots, btpParams, testContext.sk)
			btp, err := NewBootstrapper(testContext.params, btpParams, btpKey)
			if err != nil {
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

			plaintext := NewPlaintext(testContext.params, testContext.params.MaxLevel(), testContext.params.scale)
			testContext.encoder.Encode(plaintext, values, slots)

			ciphertext := testContext.encryptorPk.EncryptNew(plaintext)

			for i := 0; i < 1; i++ {

				ciphertext = btp.Bootstrapp(ciphertext)

				//testContext.evaluator.SetScale(ciphertext, testContext.params.scale)

				verifyTestVectors(testContext, testContext.decryptor, values, ciphertext, t)
			}

		})
	}
}

func newTestVectorsSineBootstrapp(testContext *testParams, encryptor Encryptor, a, b float64, t *testing.T) (values []complex128, plaintext *Plaintext, ciphertext *Ciphertext) {

	slots := testContext.params.Slots()

	values = make([]complex128, slots)

	for i := uint64(0); i < slots; i++ {
		values[i] = complex(math.Round(randomFloat(a, b))+randomFloat(-1, 1)/1000, 0)
	}

	plaintext = NewPlaintext(testContext.params, testContext.params.MaxLevel(), testContext.params.Scale())

	testContext.encoder.EncodeNTT(plaintext, values, slots)

	if encryptor != nil {
		ciphertext = encryptor.EncryptNew(plaintext)
	}

	return values, plaintext, ciphertext
}

// DefaultBootstrappSchemeParamsShort are insecure params
// for quick correctness testing of the bootstrapping
var DefaultBootstrappSchemeParamsShort = []*Parameters{

	{
		logN:     14,
		logSlots: 13,
		Moduli: Moduli{
			qi: []uint64{
				0x80000000080001,  // 55 Q0
				0x2000000a0001,    // 45
				0x2000000e0001,    // 45
				0x1fffffc20001,    // 45
				0x200000440001,    // 45
				0x200000500001,    // 45
				0x200000620001,    // 45
				0x1fffff980001,    // 45
				0x2000006a0001,    // 45
				0x1fffff7e0001,    // 45
				0x200000860001,    // 45
				0x100000000060001, // 56 StC (28 + 28)
				0xffa0001,         // 28 StC
				0x80000000440001,  // 55 Sine (double angle)
				0x7fffffffba0001,  // 55 Sine (double angle)
				0x80000000500001,  // 55 Sine
				0x7fffffffaa0001,  // 55 Sine
				0x800000005e0001,  // 55 Sine
				0x7fffffff7e0001,  // 55 Sine
				0x7fffffff380001,  // 55 Sine
				0x80000000ca0001,  // 55 Sine
				0x200000000e0001,  // 53 CtS
				0x20000000140001,  // 53 CtS
				0x20000000280001,  // 53 CtS
				0x1fffffffd80001,  // 53 CtS
			},
			pi: []uint64{
				0xfffffffff00001,  // 56
				0xffffffffd80001,  // 56
				0x1000000002a0001, // 56
				0xffffffffd20001,  // 56
				0x100000000480001, // 56
			},
		},
		scale: 1 << 45,
		sigma: DefaultSigma,
	},

	{
		logN:     14,
		logSlots: 12,
		Moduli: Moduli{
			qi: []uint64{
				0x80000000080001,  // 55 Q0
				0x2000000a0001,    // 45
				0x2000000e0001,    // 45
				0x1fffffc20001,    // 45
				0x200000440001,    // 45
				0x200000500001,    // 45
				0x200000620001,    // 45
				0x1fffff980001,    // 45
				0x2000006a0001,    // 45
				0x1fffff7e0001,    // 45
				0x200000860001,    // 45
				0x100000000060001, // 56 StC (28 + 28)
				0xffa0001,         // 28 StC
				0x80000000440001,  // 55 Sine (double angle)
				0x7fffffffba0001,  // 55 Sine (double angle)
				0x80000000500001,  // 55 Sine
				0x7fffffffaa0001,  // 55 Sine
				0x800000005e0001,  // 55 Sine
				0x7fffffff7e0001,  // 55 Sine
				0x7fffffff380001,  // 55 Sine
				0x80000000ca0001,  // 55 Sine
				0x200000000e0001,  // 53 CtS
				0x20000000140001,  // 53 CtS
				0x20000000280001,  // 53 CtS
				0x1fffffffd80001,  // 53 CtS
			},
			pi: []uint64{
				0xfffffffff00001,  // 56
				0xffffffffd80001,  // 56
				0x1000000002a0001, // 56
				0xffffffffd20001,  // 56
				0x100000000480001, // 56
			},
		},
		scale: 1 << 45,
		sigma: DefaultSigma,
	},
}

// DefaultBootstrappParamsShort are default bootstrapping params for the
// DefaultBootstrappSchemeParamsShort scheme params
var DefaultBootstrappParamsShort = []*BootstrappParams{

	// SET II
	// 1525 Cos - 550
	{
		H:            196,
		SinType:      Cos,
		SinRange:     21,
		SinDeg:       52,
		SinRescal:    2,
		CtSLevel:     []uint64{24, 23, 22, 21},
		StCLevel:     []uint64{12, 11, 11},
		MaxN1N2Ratio: 16.0,
	},

	// SET II
	// 1525 Cos - 550
	{
		H:            196,
		SinType:      Cos,
		SinRange:     21,
		SinDeg:       52,
		SinRescal:    2,
		CtSLevel:     []uint64{24, 23, 22, 21},
		StCLevel:     []uint64{12, 11, 11},
		MaxN1N2Ratio: 16.0,
	},
}
