package ckks

import (
	//"github.com/ldsec/lattigo/ring"
	"fmt"
	"log"
	"math/rand"
	"testing"
	"time"
)

func Test_Bootstrapp(t *testing.T) {

	rand.Seed(time.Now().UnixNano())

	params := genCkksParams(&Parameters{14, []uint8{55, 45, 45, 45, 55, 55, 55, 55, 55, 55, 55, 55, 55, 45, 45, 45}, []uint8{55, 55, 55, 55}, 1 << 40, 3.2})

	//medianprec := float64(20) // target median precision in log2 among all the coeffs, determines the success/failure of a test

	ctsDepth := uint64(3)
	stcDepth := uint64(3)

	slots := uint64(1 << 10)

	t.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/a=%d/b=%d/TestBoot", params.ckkscontext.logN,
		params.ckkscontext.logQ,
		params.ckkscontext.levels,
		params.ckkscontext.alpha,
		params.ckkscontext.beta), func(t *testing.T) {

		var bootcontext *BootContext
		var err error

		if bootcontext, err = params.ckkscontext.NewBootContext(slots, params.sk, ctsDepth, stcDepth); err != nil {
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

		plaintext := params.ckkscontext.NewPlaintext(params.ckkscontext.levels-1, params.ckkscontext.scale)
		params.encoder.Encode(plaintext, values, slots)

		ciphertext := params.encryptorPk.EncryptNew(plaintext)

		for i := 0; i < 1; i++ {

			ciphertext = params.evaluator.Bootstrapp(ciphertext, bootcontext)

			//if err = evaluator.SetScale(ciphertext, params.Scale); err != nil {
			//	log.Fatal(err)
			//}

			verify_test_vectors(params, params.decryptor, values, ciphertext, t)
		}

	})
}
