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

	for i, paramLit := range ckks.DefaultParams {

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
		benchRelin(b, p)
		benchRotate(b, p)

		if testing.Verbose() {
			benchMemoryConsumption(b, p)
		}

		if i != 4 {
			for i := uint64(1); i < 11; i++ {
				benchAddIncreasingParticipants(i, b, p)
				benchMultIncreasingParticipants(i, b, p)
			}

			for i := uint64(1); i < 21; i++ {
				benchDecryptionIncreasingParticipants(i, b, p)
				benchRotIncreasingParticipants(i, b, p)
			}
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
	ciphers := evaluator.ConvertToMKCiphertext([]*ckks.Ciphertext{cipher1, cipher2}, []uint64{1, 2})

	b.Run(testString("Add/", 2, params), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			evaluator.Add(ciphers[0], ciphers[1])
		}
	})
}

func benchRotate(b *testing.B, params *ckks.Parameters) {

	participants := setupPeers(1, params, 6.0)

	value1 := newTestValue(params, complex(-1, -1), complex(1, 1))

	cipher1 := participants[0].Encrypt(value1)

	evaluator := NewMKEvaluator(params)
	ciphers := evaluator.ConvertToMKCiphertext([]*ckks.Ciphertext{cipher1}, []uint64{1})
	rotKey := participants[0].GetRotationKeys(15)
	rotKey.PeerID = 1

	b.Run(testString("Rotate/", 1, params), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			evaluator.Rotate(ciphers[0], 15, []*mkrlwe.MKRotationKey{rotKey})
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
	ciphers := evaluator.ConvertToMKCiphertext([]*ckks.Ciphertext{cipher1, cipher2}, []uint64{1, 2})

	b.Run(testString("Mul/", 2, params), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			evaluator.Mul(ciphers[0], ciphers[1])
		}
	})
}

func benchRelin(b *testing.B, params *ckks.Parameters) {

	participants := setupPeers(2, params, 6.0)

	value1 := newTestValue(params, complex(-1, -1), complex(1, 1))
	value2 := newTestValue(params, complex(-1, -1), complex(1, 1))

	cipher1 := participants[0].Encrypt(value1)
	cipher2 := participants[1].Encrypt(value2)

	evaluator := NewMKEvaluator(params)
	ciphers := evaluator.ConvertToMKCiphertext([]*ckks.Ciphertext{cipher1, cipher2}, []uint64{1, 2})

	evk1 := participants[0].GetEvaluationKey()
	evk1.PeerID = 1
	evk2 := participants[1].GetEvaluationKey()
	evk2.PeerID = 2

	evalKeys := []*mkrlwe.MKRelinearizationKey{participants[0].GetEvaluationKey(), participants[1].GetEvaluationKey()}

	pk1 := participants[0].GetPublicKey()
	pk2 := participants[1].GetPublicKey()
	pk1.PeerID = 1
	pk2.PeerID = 2

	publicKeys := []*mkrlwe.MKPublicKey{pk1, pk2}

	b.Run(testString("Relin/", 2, params), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			b.StopTimer()
			res := evaluator.Mul(ciphers[0], ciphers[1])
			b.StartTimer()
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

	ckksCipher1 := make([]*ckks.Ciphertext, nbrParticipants)
	ckksCipher2 := make([]*ckks.Ciphertext, nbrParticipants)

	ids1 := make([]uint64, nbrParticipants)
	ids2 := make([]uint64, nbrParticipants)

	for i := uint64(0); i < nbrParticipants; i++ {
		ckksCipher1[i] = participants[2*i].Encrypt(newTestValue(params, complex(-1, -1), complex(1, 1)))
		ckksCipher2[i] = participants[2*i+1].Encrypt(newTestValue(params, complex(-1, -1), complex(1, 1)))
		ids1[i] = 2 * i
		ids2[i] = 2*i + 1
	}

	evaluator := NewMKEvaluator(params)

	ciphers1 := evaluator.ConvertToMKCiphertext(ckksCipher1, ids1)
	ciphers2 := evaluator.ConvertToMKCiphertext(ckksCipher2, ids2)

	evalKeys := make([]*mkrlwe.MKRelinearizationKey, 2*nbrParticipants)
	pubKeys := make([]*mkrlwe.MKPublicKey, 2*nbrParticipants)

	// perform additions until ciphertexts concerns all participants and then Square + Relin
	resCipher1 := ciphers1[0]
	resCipher2 := ciphers2[0]
	evalKeys[0] = participants[0].GetEvaluationKey()
	evalKeys[0].PeerID = 0
	pubKeys[0] = participants[0].GetPublicKey()
	pubKeys[0].PeerID = 0
	evalKeys[1] = participants[1].GetEvaluationKey()
	evalKeys[1].PeerID = 1
	pubKeys[1] = participants[1].GetPublicKey()
	pubKeys[1].PeerID = 1

	for i := uint64(1); i < nbrParticipants; i++ {
		resCipher1 = evaluator.Add(resCipher1, ciphers1[i])
		resCipher2 = evaluator.Add(resCipher2, ciphers2[i])

		// prepare public material
		evalKeys[2*i] = participants[2*i].GetEvaluationKey()
		evalKeys[2*i].PeerID = 2 * i
		evalKeys[2*i+1] = participants[2*i+1].GetEvaluationKey()
		evalKeys[2*i+1].PeerID = 2*i + 1
		pubKeys[2*i] = participants[2*i].GetPublicKey()
		pubKeys[2*i].PeerID = 2 * i
		pubKeys[2*i+1] = participants[2*i+1].GetPublicKey()
		pubKeys[2*i+1].PeerID = 2*i + 1
	}

	b.Run(testString("Mul Increasing number of participants/", nbrParticipants, params), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			evaluator.Mul(resCipher1, resCipher2)
		}
	})

	b.Run(testString("Relin Increasing number of participants/", nbrParticipants, params), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			b.StopTimer()
			res := evaluator.Mul(resCipher1, resCipher2)
			b.StartTimer()
			evaluator.RelinInPlace(res, evalKeys, pubKeys)
		}
	})

}

func benchAddIncreasingParticipants(nbrParticipants uint64, b *testing.B, params *ckks.Parameters) {

	participants := setupPeers(2*nbrParticipants, params, 6.0)

	ckksCipher1 := make([]*ckks.Ciphertext, nbrParticipants)
	ckksCipher2 := make([]*ckks.Ciphertext, nbrParticipants)

	ids1 := make([]uint64, nbrParticipants)
	ids2 := make([]uint64, nbrParticipants)

	for i := uint64(0); i < nbrParticipants; i++ {
		ckksCipher1[i] = participants[2*i].Encrypt(newTestValue(params, complex(-1, -1), complex(1, 1)))
		ckksCipher2[i] = participants[2*i+1].Encrypt(newTestValue(params, complex(-1, -1), complex(1, 1)))
		ids1[i] = 2 * i
		ids2[i] = 2*i + 1
	}

	evaluator := NewMKEvaluator(params)

	ciphers1 := evaluator.ConvertToMKCiphertext(ckksCipher1, ids1)
	ciphers2 := evaluator.ConvertToMKCiphertext(ckksCipher2, ids2)

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

	ciphersCKKS := make([]*ckks.Ciphertext, nbrParticipants)

	ids := make([]uint64, nbrParticipants)

	for i := uint64(0); i < nbrParticipants; i++ {
		ciphersCKKS[i] = participants[i].Encrypt(newTestValue(params, complex(-1, -1), complex(1, 1)))
		ids[i] = i
	}

	evaluator := NewMKEvaluator(params)

	ciphers := evaluator.ConvertToMKCiphertext(ciphersCKKS, ids)

	galKeys := make([]*mkrlwe.MKRotationKey, nbrParticipants)

	// perform additions until ciphertexts concerns all participants and then Square + Relin
	resCipher := ciphers[0]
	galKeys[0] = participants[0].GetRotationKeys(15)
	galKeys[0].PeerID = 0

	for i := uint64(1); i < nbrParticipants; i++ {
		resCipher = evaluator.Add(resCipher, ciphers[i])

		// prepare public material
		galKeys[i] = participants[i].GetRotationKeys(15)
		galKeys[i].PeerID = i
	}

	b.Run(testString("Rotation Increasing number of participants/", nbrParticipants, params), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			evaluator.Rotate(resCipher, 15, galKeys)

		}
	})

}

func benchDecryptionIncreasingParticipants(nbrParticipants uint64, b *testing.B, params *ckks.Parameters) {
	participants := setupPeers(nbrParticipants, params, 6.0)

	ciphersCKKS := make([]*ckks.Ciphertext, nbrParticipants)

	ids := make([]uint64, nbrParticipants)

	for i := uint64(0); i < nbrParticipants; i++ {
		ciphersCKKS[i] = participants[i].Encrypt(newTestValue(params, complex(-1, -1), complex(1, 1)))
		ids[i] = i
	}

	evaluator := NewMKEvaluator(params)

	ciphers := evaluator.ConvertToMKCiphertext(ciphersCKKS, ids)

	partialDec := make([]*ring.Poly, nbrParticipants)

	// perform additions until ciphertexts concerns all participants and then Square + Relin
	resCipher := ciphers[0]

	for i := uint64(1); i < nbrParticipants; i++ {
		resCipher = evaluator.Add(resCipher, ciphers[i])
	}

	resCKKS := evaluator.ConvertToCKKSCiphertext(resCipher)

	for i := uint64(0); i < nbrParticipants; i++ {
		partialDec[i] = participants[i].GetPartialDecryption(resCKKS[i])

	}

	b.Run(testString("Decryption Increasing number of participants/", nbrParticipants, params), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			participants[0].Decrypt(resCKKS[0], partialDec)
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

		b.Logf("Size of public key: %d bytes", len(data[8:]))
	})

	b.Run("Measure Memory Evaluation Key", func(b *testing.B) {

		data, err := evalKey.MarshalBinary()

		if err != nil {
			b.Error("Couldn't marshal evaluation key")
		}

		b.Logf("Size of evaluation key: %d bytes", len(data[8:]))
	})

	b.Run("Measure Memory Galois Evaluation Key", func(b *testing.B) {

		data, err := evalGalKey.MarshalBinary()

		if err != nil {
			b.Error("Couldn't marshal galois evaluation key")
		}

		b.Logf("Size of galois evaluation key: %d bytes", len(data))
	})

	b.Run("Measure Memory Ciphertext", func(b *testing.B) {

		value := newTestValue(params, complex(-1, -1), complex(1, 1))

		cipher := participants[0].Encrypt(value)
		mkCipher := &MKCiphertext{Ciphertexts: cipher, PeerID: []uint64{1}}

		data, err := mkCipher.MarshalBinary()
		if err != nil {
			b.Error("ciphertext not correctly marshalled")
		}

		b.Logf("Size of ciphertext: %d bytes", len(data[8:]))
	})

}
