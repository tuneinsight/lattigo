package ckks

import (
	//"github.com/ldsec/lattigo/ring"
	"fmt"
	"log"
	"testing"
)

func test_Bootstrapp(params *CKKSTESTPARAMS, ctsDepth, stcDepth uint64, repack bool, verbose bool, t *testing.T) {

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
	/*
		t.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/TestLT", params.ckkscontext.logN,
			params.ckkscontext.logQ,
			params.ckkscontext.levels), func(t *testing.T) {


			ciphertext, err := params.encryptor.EncryptFromPkNew(plaintext)
			if err != nil {
				log.Fatal(err)
			}

			showcoeffs(params.decryptor, params.encoder, ciphertext, "Before LT : ")

			vec := make([]int64, params.ckkscontext.n)
			for i := range vec {
				vec[i] = int64(ring.RandUniform(6, 7) * params.ckkscontext.moduli[0])
				sign := ring.RandUniform(2, 3)
				if sign == 2 {
					vec[i] = -vec[i]
				}
			}

			ciphertext.InvNTT(params.ckkscontext, ciphertext)
			for i := uint64(0) ; i < params.ckkscontext.n ; i++ {
				for j, qi := range params.ckkscontext.moduli[:ciphertext.Level()+1] {
					ciphertext.Value()[0].Coeffs[j][i] = ring.CRed(ciphertext.Value()[0].Coeffs[j][i] + uint64((vec[i]%int64(qi)) + int64(qi))%qi, qi)
				}
			}
			ciphertext.NTT(params.ckkscontext, ciphertext)


			var vec0, vec1 *Ciphertext

			if vec0, vec1, err = bootcontext.coeffsToSlots(params.evaluator, ciphertext) ; err != nil {
				log.Fatal(err)
			}

			showcoeffs(params.decryptor, params.encoder, vec0, "vec0      : ")
			showcoeffs(params.decryptor, params.encoder, vec1, "vec1      : ")

			if ciphertext, err = bootcontext.slotsToCoeffs(params.evaluator, vec0, vec1) ; err != nil {
				log.Fatal(err)
			}


			ciphertext.InvNTT(params.ckkscontext, ciphertext)
			for i := uint64(0) ; i < params.ckkscontext.n ; i++ {
				for j, qi := range params.ckkscontext.moduli[:ciphertext.Level()+1] {
					ciphertext.Value()[0].Coeffs[j][i] = ring.CRed(ciphertext.Value()[0].Coeffs[j][i] + (qi - (uint64((vec[i]%int64(qi)) + int64(qi))%qi)), qi)
				}
			}
			ciphertext.NTT(params.ckkscontext, ciphertext)

			params.evaluator.RemoveImag(ciphertext, bootcontext.rotkeys, ciphertext)

			showcoeffs(params.decryptor, params.encoder, ciphertext, "After  LT : ")


			if err := verify_test_vectors(params, values, ciphertext, t); err != nil {
				log.Fatal(err)
			}

		})
	*/
	t.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/a=%d/b=%d/TestBoot", params.ckkscontext.logN,
		params.ckkscontext.logQ,
		params.ckkscontext.levels,
		params.ckkscontext.alpha,
		params.ckkscontext.beta), func(t *testing.T) {

		ciphertext, err := params.encryptorPk.EncryptNew(plaintext)
		if err != nil {
			log.Fatal(err)
		}

		/*
			ciphertext.InvNTT(ciphertext)
			for i := uint64(0) ; i < params.ckkscontext.n ; i++ {
				q := ring.RandUniform(12)*params.ckkscontext.modulie[0]
				for j, qi := range params.ckkscontext.modulie {
					ciphertext.Value()[0].Coeffs[j][i] = ring.BRedAdd(ciphertext.Value()[0].Coeffs[j][i] + q, qi, params.ckkscontext.contextLevel[params.ckkscontext.levels-1].GetBredParams()[j])
				}
			}
			ciphertext.NTT(ciphertext)
		*/

		if ciphertext, err = params.evaluator.Bootstrapp(ciphertext, bootcontext, params.decryptor, params.encoder, verbose); err != nil {
			log.Fatal(err)
		}

		if err := verify_test_vectors(params, values, ciphertext.Element(), t); err != nil {
			log.Fatal(err)
		}
	})
}
