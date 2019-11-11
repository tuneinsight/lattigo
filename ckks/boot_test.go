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

	var err error

	medianprec := float64(20) // target median precision in log2 among all the coeffs, determines the success/failure of a test

	ctsDepth := uint64(3)
	stcDepth := uint64(3)

	params := Parameters{10, []uint8{55, 45, 45, 45, 55, 55, 55, 55, 55, 55, 55, 55, 55, 45, 45, 45}, []uint8{55, 55, 55, 55}, 1 << 40, 3.2}

	ckksTest := new(CKKSTESTPARAMS)

	ckksTest.medianprec = medianprec

	ckksTest.slots = 1 << 8

	ckksTest.levels = uint64(len(params.Modulichain))
	ckksTest.scale = params.Scale

	if ckksTest.ckkscontext, err = NewCkksContext(&params); err != nil {
		t.Error(err)
	}

	log.Printf("Generated CkksContext for logN=%d/logQ=%d/levels=%d/a=%d/b=%d/sigma=%f",
		ckksTest.ckkscontext.LogN(),
		ckksTest.ckkscontext.LogQ(),
		ckksTest.ckkscontext.Levels(),
		ckksTest.ckkscontext.alpha,
		ckksTest.ckkscontext.beta,
		ckksTest.ckkscontext.Sigma())

	ckksTest.kgen = ckksTest.ckkscontext.NewKeyGenerator()

	ckksTest.sk, ckksTest.pk = ckksTest.kgen.NewKeyPairSparse(128)

	ckksTest.encoder = ckksTest.ckkscontext.NewEncoder()

	if ckksTest.encryptorPk, err = ckksTest.ckkscontext.NewEncryptorFromPk(ckksTest.pk); err != nil {
		t.Error(err)
	}

	if ckksTest.encryptorSk, err = ckksTest.ckkscontext.NewEncryptorFromSk(ckksTest.sk); err != nil {
		t.Error(err)
	}

	if ckksTest.decryptor, err = ckksTest.ckkscontext.NewDecryptor(ckksTest.sk); err != nil {
		t.Error(err)
	}

	ckksTest.evaluator = ckksTest.ckkscontext.NewEvaluator()

	t.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/a=%d/b=%d/TestBoot", ckksTest.ckkscontext.logN,
		ckksTest.ckkscontext.logQ,
		ckksTest.ckkscontext.levels,
		ckksTest.ckkscontext.alpha,
		ckksTest.ckkscontext.beta), func(t *testing.T) {

		var bootcontext *BootContext
		var err error

		if bootcontext, err = ckksTest.ckkscontext.NewBootContext(ckksTest.slots, ckksTest.sk, ctsDepth, stcDepth); err != nil {
			log.Fatal(err)
		}

		values := make([]complex128, ckksTest.slots)
		for i := range values {
			values[i] = complex(randomFloat(-1, 1), 0)
		}

		values[0] = complex(0.516015, 0)
		values[1] = complex(0.772621, 0)
		if ckksTest.slots > 2 {
			values[2] = complex(0.939175, 0)
			values[3] = complex(0.345987, 0)
		}

		plaintext := ckksTest.ckkscontext.NewPlaintext(ckksTest.ckkscontext.levels-1, ckksTest.scale)
		if err = ckksTest.encoder.Encode(plaintext, values, ckksTest.slots); err != nil {
			log.Fatal(err)
		}

		ciphertext, err := ckksTest.encryptorPk.EncryptNew(plaintext)
		if err != nil {
			log.Fatal(err)
		}

		for i := 0; i < 1; i++ {

			if ciphertext, err = ckksTest.evaluator.Bootstrapp(ciphertext, bootcontext); err != nil {
				log.Fatal(err)
			}

			//if err = ckksTest.evaluator.SetScale(ciphertext, params.Scale); err != nil {
			//	log.Fatal(err)
			//}

			if err := verify_test_vectors(ckksTest, values, ciphertext.Element(), t); err != nil {
				log.Fatal(err)
			}
		}

	})
}
