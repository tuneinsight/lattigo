package ckks

import (
	"log"
	"math/rand"
	"testing"
	"time"
)

func Test_Bootstrapp(t *testing.T) {

	rand.Seed(time.Now().UnixNano())

	bootParams := new(Parameters)
	bootParams.LogN = 15
	bootParams.LogSlots = 12
	bootParams.Scale = 1 << 40
	bootParams.LogQi = []uint64{55, 45, 45, 45, 55, 55, 55, 55, 55, 55, 55, 55, 55, 45, 45, 45}
	bootParams.LogPi = []uint64{55, 55, 55, 55}
	bootParams.Sigma = 3.2

	bootParams.GenFromLogModuli()

	params := genCkksParams(bootParams)

	//medianprec := float64(20) // target median precision in log2 among all the coeffs, determines the success/failure of a test

	ctsDepth := uint64(3)
	stcDepth := uint64(3)

	slots := uint64(1 << bootParams.LogSlots)

	t.Run(testString("TestBoot/", bootParams), func(t *testing.T) {

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
}
