package mkbfv

import (
	"fmt"
	"testing"

	"github.com/ldsec/lattigo/v2/bfv"
	"github.com/ldsec/lattigo/v2/mkrlwe"
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/utils"
)

func testString(opname string, parties uint64, params *bfv.Parameters) string {
	return fmt.Sprintf("%sparties=%d/LogN=%d/logQ=%d", opname, parties, params.LogN(), params.LogQP())
}

func Benchmark_MKBFV(b *testing.B) {

	for _, paramLit := range bfv.DefaultParams {

		params, err := bfv.NewParametersFromLiteral(paramLit)
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
	}
}

func benchKeyGen(b *testing.B, params *bfv.Parameters) {

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

func benchAddTwoCiphertexts(b *testing.B, params *bfv.Parameters) {

	participants := setupPeers(2, params, 6.0)

	ringT := getRingT(params)

	value1 := getRandomPlaintextValue(ringT, params)
	value2 := getRandomPlaintextValue(ringT, params)

	cipher1 := participants[0].Encrypt(value1)
	cipher2 := participants[1].Encrypt(value2)

	evaluator := NewMKEvaluator(params)

	b.Run(testString("Add/", 2, params), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			evaluator.Add(cipher1, cipher2)
		}
	})
}

func benchRotate(b *testing.B, params *bfv.Parameters) {

	participants := setupPeers(1, params, 6.0)

	ringT := getRingT(params)

	value1 := getRandomPlaintextValue(ringT, params)

	cipher1 := participants[0].Encrypt(value1)

	evaluator := NewMKEvaluator(params)

	rotKey := participants[0].GetRotationKeys(15)

	b.Run(testString("Rotate/", 1, params), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			evaluator.Rotate(cipher1, 15, []*mkrlwe.MKEvalGalKey{rotKey})
		}
	})
}

func benchMultTwoCiphertexts(b *testing.B, params *bfv.Parameters) {

	participants := setupPeers(2, params, 6.0)

	ringT := getRingT(params)

	value1 := getRandomPlaintextValue(ringT, params)
	value2 := getRandomPlaintextValue(ringT, params)

	cipher1 := participants[0].Encrypt(value1)
	cipher2 := participants[1].Encrypt(value2)

	evaluator := NewMKEvaluator(params)

	b.Run(testString("Mul/", 2, params), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			evaluator.Mul(cipher1, cipher2)
		}
	})
}

func benchMulAndRelin(b *testing.B, params *bfv.Parameters) {

	participants := setupPeers(2, params, 6.0)

	ringT := getRingT(params)

	value1 := getRandomPlaintextValue(ringT, params)
	value2 := getRandomPlaintextValue(ringT, params)

	cipher1 := participants[0].Encrypt(value1)
	cipher2 := participants[1].Encrypt(value2)

	evaluator := NewMKEvaluator(params)

	evalKeys := []*mkrlwe.MKEvaluationKey{participants[0].GetEvaluationKey(), participants[1].GetEvaluationKey()}
	pubKeys := []*mkrlwe.MKPublicKey{participants[0].GetPublicKey(), participants[1].GetPublicKey()}

	b.Run(testString("Mul and Relin/", 2, params), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			res := evaluator.Mul(cipher1, cipher2)
			evaluator.RelinInPlace(res, evalKeys, pubKeys)
		}
	})
}

func benchEncrypt(b *testing.B, params *bfv.Parameters) {

	participants := setupPeers(1, params, 6.0)

	ringT := getRingT(params)

	value1 := getRandomPlaintextValue(ringT, params)

	b.Run(testString("Encrypt/", 2, params), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			participants[0].Encrypt(value1)
		}
	})
}

func benchDecrypt(b *testing.B, params *bfv.Parameters) {

	participants := setupPeers(1, params, 6.0)

	ringT := getRingT(params)

	value1 := getRandomPlaintextValue(ringT, params)

	cipher1 := participants[0].Encrypt(value1)
	partialDec := participants[0].GetPartialDecryption(cipher1)

	b.Run(testString("Decrypt/", 1, params), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			participants[0].Decrypt(cipher1, []*ring.Poly{partialDec})
		}
	})
}

func benchPartialDecrypt(b *testing.B, params *bfv.Parameters) {

	participants := setupPeers(1, params, 6.0)

	ringT := getRingT(params)

	value1 := getRandomPlaintextValue(ringT, params)

	cipher1 := participants[0].Encrypt(value1)

	b.Run(testString("Partial decryption/", 1, params), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			participants[0].GetPartialDecryption(cipher1)
		}
	})
}
