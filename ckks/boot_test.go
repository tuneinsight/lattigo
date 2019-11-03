package ckks

import (
	//"github.com/ldsec/lattigo/ring"
	"fmt"
	"log"
	"testing"
)

func test_Bootstrapp(params *CKKSTESTPARAMS, t *testing.T) {

	repack := true
	ctsDepth := uint64(2)
	stcDepth := uint64(2)

	t.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/a=%d/b=%d/TestBoot", params.ckkscontext.logN,
		params.ckkscontext.logQ,
		params.ckkscontext.levels,
		params.ckkscontext.alpha,
		params.ckkscontext.beta), func(t *testing.T) {

		var bootcontext *BootContext
		var err error

		if bootcontext, err = params.ckkscontext.NewBootContext(params.slots, params.sk, ctsDepth, stcDepth, repack); err != nil {
			log.Fatal(err)
		}

		values := make([]complex128, params.slots)
		for i := range values {
			values[i] = complex(randomFloat(-1, 1), 0)
		}

		values[0] = complex(1, 0) //complex(0.516015, 0)
		values[1] = complex(2, 0) //complex(0.772621, 0)
		if params.slots > 2 {
			values[2] = complex(3, 0) //complex(0.939175, 0)
			values[3] = complex(4, 0) //complex(0.345987, 0)
		}

		plaintext := params.ckkscontext.NewPlaintext(params.ckkscontext.levels-1, params.scale)
		if err = params.encoder.Encode(plaintext, values, params.slots); err != nil {
			log.Fatal(err)
		}

		ciphertext, err := params.encryptorPk.EncryptNew(plaintext)
		if err != nil {
			log.Fatal(err)
		}

		if ciphertext, err = params.evaluator.Bootstrapp(ciphertext, bootcontext); err != nil {
			log.Fatal(err)
		}

		if err := verify_test_vectors(params, values, ciphertext.Element(), t); err != nil {
			log.Fatal(err)
		}
	})
}
