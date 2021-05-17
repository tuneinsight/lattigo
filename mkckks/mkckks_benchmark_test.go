package mkckks

import (
	"fmt"
	"testing"

	"github.com/ldsec/lattigo/v2/ckks"
	"github.com/ldsec/lattigo/v2/mkrlwe"
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/utils"
)

func testString(opname string, parties uint64, params *ckks.Parameters) string {
	return fmt.Sprintf("%sparties=%d/logN=%d/logQ=%d/levels=%d/alpha=%d/beta=%d",
		opname,
		parties,
		params.LogN(),
		params.LogQP(),
		params.MaxLevel()+1,
		params.Alpha(),
		params.Beta())
}

func BenchmarkMKCKKS(b *testing.B) {

	for _, paramLit := range ckks.DefaultParams {

		params, err := ckks.NewParametersFromLiteral(paramLit)
		if err != nil {
			panic(err)
		}

		p := &params

		benchKeyGen(b, p)
		benchAddTwoCiphertexts(b, p)
		benchEncrypt(b, p)
		benchDecrypt(b, p)
		benchPartialDecrypt(b, p)
		benchMultTwoCiphertexts(b, p)
		benchMulAndRelin(b, p)
		benchRotate(b, p)
		benchMemoryConsumption(b, p)

		for i := uint64(2); i < 20; i++ {
			benchDecryptionIncreasingParticipants(i, b, p)
			benchRotIncreasingParticipants(i, b, p)
			benchAddIncreasingParticipants(i, b, p)
			benchMultIncreasingParticipants(i, b, p)
		}
	}
}

func benchKeyGen(b *testing.B, params *ckks.Parameters) {

	prng, err := utils.NewKeyedPRNG([]byte{'l', 'a', 't', 't', 'i', 'g', 'o'})

	if err != nil {
		panic(err)
	}

	crs := mkrlwe.GenCommonPublicParam(&params.Parameters, prng)

	b.Run(testString("KeyGen/", 1, params), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			mkrlwe.KeyGen(&params.Parameters, crs)
		}
	})
}

func benchAddTwoCiphertexts(b *testing.B, params *ckks.Parameters) {

	participants := setupPeers(2, params, 6.0)

	value1 := newTestValue(params, complex(-1, -1), complex(1, 1))
	value2 := newTestValue(params, complex(-1, -1), complex(1, 1))

	cipher1 := participants[0].Encrypt(value1)
	cipher2 := participants[1].Encrypt(value2)

	evaluator := NewMKEvaluator(params)

	b.Run(testString("Add/", 2, params), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			evaluator.Add(cipher1, cipher2)
		}
	})
}

func benchRotate(b *testing.B, params *ckks.Parameters) {

	participants := setupPeers(1, params, 6.0)

	value1 := newTestValue(params, complex(-1, -1), complex(1, 1))

	cipher1 := participants[0].Encrypt(value1)

	evaluator := NewMKEvaluator(params)

	rotKey := participants[0].GetRotationKeys(15)

	b.Run(testString("Rotate/", 1, params), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			evaluator.Rotate(cipher1, 15, []*mkrlwe.MKEvalGalKey{rotKey})
		}
	})
}

func benchMultTwoCiphertexts(b *testing.B, params *ckks.Parameters) {

	participants := setupPeers(2, params, 6.0)

	value1 := newTestValue(params, complex(-1, -1), complex(1, 1))
	value2 := newTestValue(params, complex(-1, -1), complex(1, 1))

	cipher1 := participants[0].Encrypt(value1)
	cipher2 := participants[1].Encrypt(value2)

	evaluator := NewMKEvaluator(params)

	b.Run(testString("Mul/", 2, params), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			evaluator.Mul(cipher1, cipher2)
		}
	})
}

func benchMulAndRelin(b *testing.B, params *ckks.Parameters) {

	participants := setupPeers(2, params, 6.0)

	value1 := newTestValue(params, complex(-1, -1), complex(1, 1))
	value2 := newTestValue(params, complex(-1, -1), complex(1, 1))

	cipher1 := participants[0].Encrypt(value1)
	cipher2 := participants[1].Encrypt(value2)

	evaluator := NewMKEvaluator(params)
	evalKeys := []*mkrlwe.MKEvaluationKey{participants[0].GetEvaluationKey(), participants[1].GetEvaluationKey()}
	publicKeys := []*mkrlwe.MKPublicKey{participants[0].GetPublicKey(), participants[1].GetPublicKey()}

	b.Run(testString("Mul and Relin/", 2, params), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			res := evaluator.Mul(cipher1, cipher2)
			evaluator.RelinInPlace(res, evalKeys, publicKeys)
		}
	})
}

func benchEncrypt(b *testing.B, params *ckks.Parameters) {

	participants := setupPeers(1, params, 6.0)

	value1 := newTestValue(params, complex(-1, -1), complex(1, 1))

	b.Run(testString("Encrypt/", 2, params), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			participants[0].Encrypt(value1)
		}
	})
}

func benchDecrypt(b *testing.B, params *ckks.Parameters) {

	participants := setupPeers(1, params, 6.0)

	value1 := newTestValue(params, complex(-1, -1), complex(1, 1))

	cipher1 := participants[0].Encrypt(value1)
	partialDec := participants[0].GetPartialDecryption(cipher1)

	b.Run(testString("Decrypt/", 1, params), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			participants[0].Decrypt(cipher1, []*ring.Poly{partialDec})
		}
	})
}

func benchPartialDecrypt(b *testing.B, params *ckks.Parameters) {

	participants := setupPeers(1, params, 6.0)

	value1 := newTestValue(params, complex(-1, -1), complex(1, 1))
	cipher1 := participants[0].Encrypt(value1)

	b.Run(testString("Partial decryption/", 1, params), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			participants[0].GetPartialDecryption(cipher1)
		}
	})
}

func benchMultIncreasingParticipants(nbrParticipants uint64, b *testing.B, params *ckks.Parameters) {

	participants := setupPeers(2*nbrParticipants, params, 6.0)

	ciphers1 := make([]*MKCiphertext, nbrParticipants)
	ciphers2 := make([]*MKCiphertext, nbrParticipants)

	for i := uint64(0); i < nbrParticipants; i++ {
		ciphers1[i] = participants[2*i].Encrypt(newTestValue(params, complex(-1, -1), complex(1, 1)))
		ciphers2[i] = participants[2*i+1].Encrypt(newTestValue(params, complex(-1, -1), complex(1, 1)))
	}

	evaluator := NewMKEvaluator(params)

	evalKeys := make([]*mkrlwe.MKEvaluationKey, 2*nbrParticipants)
	pubKeys := make([]*mkrlwe.MKPublicKey, 2*nbrParticipants)

	// perform additions until ciphertexts concerns all participants and then Square + Relin
	resCipher1 := ciphers1[0]
	resCipher2 := ciphers2[0]
	evalKeys[0] = participants[0].GetEvaluationKey()
	pubKeys[0] = participants[0].GetPublicKey()
	evalKeys[1] = participants[1].GetEvaluationKey()
	pubKeys[1] = participants[1].GetPublicKey()

	for i := uint64(1); i < nbrParticipants; i++ {
		resCipher1 = evaluator.Add(resCipher1, ciphers1[i])
		resCipher2 = evaluator.Add(resCipher2, ciphers2[i])

		// prepare public material
		evalKeys[2*i] = participants[2*i].GetEvaluationKey()
		evalKeys[2*i+1] = participants[2*i+1].GetEvaluationKey()
		pubKeys[2*i] = participants[2*i].GetPublicKey()
		pubKeys[2*i+1] = participants[2*i+1].GetPublicKey()
	}

	b.Run(testString("Mul + Relin Increasing number of participants/", nbrParticipants, params), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			res := evaluator.Mul(resCipher1, resCipher2)
			evaluator.RelinInPlace(res, evalKeys, pubKeys)
		}
	})

}

func benchAddIncreasingParticipants(nbrParticipants uint64, b *testing.B, params *ckks.Parameters) {

	participants := setupPeers(2*nbrParticipants, params, 6.0)

	ciphers1 := make([]*MKCiphertext, nbrParticipants)
	ciphers2 := make([]*MKCiphertext, nbrParticipants)

	for i := uint64(0); i < nbrParticipants; i++ {
		ciphers1[i] = participants[2*i].Encrypt(newTestValue(params, complex(-1, -1), complex(1, 1)))
		ciphers2[i] = participants[2*i+1].Encrypt(newTestValue(params, complex(-1, -1), complex(1, 1)))
	}

	evaluator := NewMKEvaluator(params)

	// perform additions until ciphertexts concerns all participants and then Add both ciphertexts
	resCipher1 := ciphers1[0]
	resCipher2 := ciphers2[0]

	for i := uint64(1); i < nbrParticipants; i++ {
		resCipher1 = evaluator.Add(resCipher1, ciphers1[i])
		resCipher2 = evaluator.Add(resCipher2, ciphers2[i])
	}

	b.Run(testString("Add Increasing number of participants/", nbrParticipants, params), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			evaluator.Add(resCipher1, resCipher2)
		}
	})

}

func benchRotIncreasingParticipants(nbrParticipants uint64, b *testing.B, params *ckks.Parameters) {
	participants := setupPeers(nbrParticipants, params, 6.0)

	ciphers := make([]*MKCiphertext, nbrParticipants)

	for i := uint64(0); i < nbrParticipants; i++ {
		ciphers[i] = participants[i].Encrypt(newTestValue(params, complex(-1, -1), complex(1, 1)))
	}

	evaluator := NewMKEvaluator(params)

	galKeys := make([]*mkrlwe.MKEvalGalKey, nbrParticipants)

	// perform additions until ciphertexts concerns all participants and then Square + Relin
	resCipher := ciphers[0]
	galKeys[0] = participants[0].GetRotationKeys(15)

	for i := uint64(1); i < nbrParticipants; i++ {
		resCipher = evaluator.Add(resCipher, ciphers[i])

		// prepare public material
		galKeys[i] = participants[i].GetRotationKeys(15)
	}

	b.Run(testString("Rotation Increasing number of participants/", nbrParticipants, params), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			evaluator.Rotate(resCipher, 15, galKeys)

		}
	})

}

func benchDecryptionIncreasingParticipants(nbrParticipants uint64, b *testing.B, params *ckks.Parameters) {
	participants := setupPeers(nbrParticipants, params, 6.0)

	ciphers1 := make([]*MKCiphertext, nbrParticipants)

	for i := uint64(0); i < nbrParticipants; i++ {
		ciphers1[i] = participants[i].Encrypt(newTestValue(params, complex(-1, -1), complex(1, 1)))
	}

	evaluator := NewMKEvaluator(params)

	partialDec := make([]*ring.Poly, nbrParticipants)

	// perform additions until ciphertexts concerns all participants and then Square + Relin
	resCipher := ciphers1[0]

	for i := uint64(1); i < nbrParticipants; i++ {
		resCipher = evaluator.Add(resCipher, ciphers1[i])
	}

	for i := uint64(0); i < nbrParticipants; i++ {
		partialDec[i] = participants[i].GetPartialDecryption(resCipher)

	}

	b.Run(testString("Decryption Increasing number of participants/", nbrParticipants, params), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			participants[0].Decrypt(resCipher, partialDec)
		}
	})

}

func benchMemoryConsumption(b *testing.B, params *ckks.Parameters) {

	participants := setupPeers(1, params, 6.0)

	pk := participants[0].GetPublicKey()
	evalKey := participants[0].GetEvaluationKey()
	evalGalKey := participants[0].GetRotationKeys(12)

	b.Run("Measure Memory Public Key", func(b *testing.B) {

		data, err := pk.MarshalBinary()

		if err != nil {
			b.Error("Couldn't marshal public key")
		}

		b.Logf("Size of public key: %d bytes", len(data[0])+len(data[1])+len(data[2]))
	})

	b.Run("Measure Memory Evaluation Key", func(b *testing.B) {

		data, err := evalKey.MarshalBinary()

		if err != nil {
			b.Error("Couldn't marshal evaluation key")
		}

		b.Logf("Size of evaluation key: %d bytes", len(data[0])+len(data[1])+len(data[2])+len(data[3]))
	})

	b.Run("Measure Memory Galois Evaluation Key", func(b *testing.B) {

		data, err := evalGalKey.MarshalBinary()

		if err != nil {
			b.Error("Couldn't marshal galois evaluation key")
		}

		b.Logf("Size of galois evaluation key: %d bytes", len(data[0])+len(data[1])+len(data[2]))
	})

}
