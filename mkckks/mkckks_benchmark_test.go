package mkckks

import (
	"fmt"
	"testing"

	"github.com/ldsec/lattigo/v2/ckks"
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

	for _, p := range ckks.DefaultParams {
		benchKeyGen(b, p)
		benchAddTwoCiphertexts(b, p)
		benchEncrypt(b, p)
		benchDecrypt(b, p)
		benchPartialDecrypt(b, p)
		//benchMultTwoCiphertexts(b, p)
		//benchRelin(b, p)
	}
}

func benchKeyGen(b *testing.B, params *ckks.Parameters) {

	prng, err := utils.NewKeyedPRNG([]byte{'l', 'a', 't', 't', 'i', 'g', 'o'})

	if err != nil {
		panic(err)
	}

	crs := GenCommonPublicParam(params, prng)

	b.Run(testString("KeyGen/", 1, params), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			KeyGen(params, crs)
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

func benchRelin(b *testing.B, params *ckks.Parameters) {

	participants := setupPeers(2, params, 6.0)

	value1 := newTestValue(params, complex(-1, -1), complex(1, 1))
	value2 := newTestValue(params, complex(-1, -1), complex(1, 1))

	cipher1 := participants[0].Encrypt(value1)
	cipher2 := participants[1].Encrypt(value2)

	evaluator := NewMKEvaluator(params)
	evalKeys := []*MKEvaluationKey{participants[0].GetEvaluationKey(), participants[1].GetEvaluationKey()}
	publicKeys := []*MKPublicKey{participants[0].GetPublicKey(), participants[1].GetPublicKey()}

	res := evaluator.Mul(cipher1, cipher2)

	b.Run(testString("Relin/", 2, params), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
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
