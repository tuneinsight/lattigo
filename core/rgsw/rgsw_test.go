package rgsw

import (
	"math/big"
	"testing"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring/ringqp"
	"github.com/tuneinsight/lattigo/v6/utils/bignum"
	"github.com/tuneinsight/lattigo/v6/utils/buffer"

	"github.com/stretchr/testify/require"
)

func TestRGSW(t *testing.T) {

	// <<<<!Insecure parameters!>>>>
	params, err := rlwe.NewParametersFromLiteral(rlwe.ParametersLiteral{
		LogN:    10,
		LogQ:    []int{35, 20},
		LogP:    []int{61, 61},
		NTTFlag: true,
	})

	require.NoError(t, err)

	kgen := rlwe.NewKeyGenerator(params)
	sk, pk := kgen.GenKeyPairNew()

	bound := 10.0

	// plaintext [-1, 0, 1]
	pt := rlwe.NewPlaintext(params, params.MaxLevel())
	kgen.GenSecretKey(&rlwe.SecretKey{Value: ringqp.Poly{Q: pt.Value, P: params.RingP().NewPoly()}})
	pt.IsMontgomery = true

	t.Run("Encryptor/SK", func(t *testing.T) {

		ct := NewCiphertext(params, params.MaxLevelQ(), params.MaxLevelP(), 0)

		NewEncryptor(params, sk).Encrypt(pt, ct)

		left, right := NoiseRGSWCiphertext(ct, pt.Value, sk, params)

		require.GreaterOrEqual(t, bound, left)
		require.GreaterOrEqual(t, bound, right)
	})

	t.Run("Encryptor/PK", func(t *testing.T) {

		ct := NewCiphertext(params, params.MaxLevelQ(), params.MaxLevelP(), 0)

		NewEncryptor(params, pk).Encrypt(pt, ct)

		left, right := NoiseRGSWCiphertext(ct, pt.Value, sk, params)

		require.GreaterOrEqual(t, bound, left)
		require.GreaterOrEqual(t, bound, right)
	})

	t.Run("Evaluator/ExternalProduct", func(t *testing.T) {

		ptRGSW := rlwe.NewPlaintext(params, params.MaxLevel())
		ptRLWE := rlwe.NewPlaintext(params, params.MaxLevel())

		k0 := 0
		k1 := 1

		setPlaintext(params, ptRGSW, k0) // X^{k0}
		setPlaintext(params, ptRLWE, k1) // X^{k1}

		scale := new(big.Int).SetUint64(params.Q()[0])

		// Scale * X^{k1}
		params.RingQ().MulScalarBigint(ptRLWE.Value, scale, ptRLWE.Value)

		ctRGSW := NewCiphertext(params, params.MaxLevelQ(), params.MaxLevelP(), 0)
		ctRLWE := rlwe.NewCiphertext(params, 1, params.MaxLevelQ())

		NewEncryptor(params, sk).Encrypt(ptRGSW, ctRGSW)
		rlwe.NewEncryptor(params, sk).Encrypt(ptRLWE, ctRLWE)

		// X^{k0} * Scale * X^{k1}
		NewEvaluator(params, nil).ExternalProduct(ctRLWE, ctRGSW, ctRLWE)

		ptHave := rlwe.NewDecryptor(params, sk).DecryptNew(ctRLWE)

		params.RingQ().INTT(ptHave.Value, ptHave.Value)

		coeffs := make([]*big.Int, params.N())

		for i := range coeffs {
			coeffs[i] = new(big.Int)
		}

		params.RingQ().PolyToBigintCentered(ptHave.Value, 1, coeffs)

		// X^{k0} * Scale * X^{k1} / Scale
		for i := range coeffs {
			bignum.DivRound(coeffs[i], scale, coeffs[i])
		}

		have := make([]uint64, params.N())
		want := make([]uint64, params.N())

		for i := range coeffs {
			have[i] = coeffs[i].Uint64()
		}

		want[k0+k1] = 1

		require.Equal(t, have, want)
	})

	t.Run("WriteAndRead", func(t *testing.T) {
		ct := NewCiphertext(params, params.MaxLevelQ(), params.MaxLevelP(), 0)
		NewEncryptor(params, pk).Encrypt(nil, ct)
		buffer.RequireSerializerCorrect(t, ct)
	})
}

func setPlaintext(params rlwe.Parameters, pt *rlwe.Plaintext, k int) {
	r := params.RingQ()

	for i := range r.SubRings {
		pt.Value.Coeffs[i][k] = 1
	}

	r.NTT(pt.Value, pt.Value)
}
