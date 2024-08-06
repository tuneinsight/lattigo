package blindrot

import (
	"testing"

	"github.com/stretchr/testify/require"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/utils"
	"github.com/tuneinsight/lattigo/v6/utils/sampling"
)

func BenchmarkHEBin(b *testing.B) {

	b.Run("BlindRotateCore/LogN=(9, 10)/LogQ=(13.6,26.99)/Gadget=2^7", func(b *testing.B) {

		// RLWE parameters of the BlindRotation
		// N=1024, Q=0x7fff801 -> 131 bit secure
		paramsBR, err := rlwe.NewParametersFromLiteral(rlwe.ParametersLiteral{
			LogN:    10,
			Q:       []uint64{0x7fff801},
			NTTFlag: NTTFlag,
		})

		require.NoError(b, err)

		// RLWE parameters of the samples
		// N=512, Q=0x3001 -> 135 bit secure
		paramsLWE, err := rlwe.NewParametersFromLiteral(rlwe.ParametersLiteral{
			LogN:    9,
			Q:       []uint64{0x3001},
			NTTFlag: NTTFlag,
		})

		require.NoError(b, err)

		evkParams := rlwe.EvaluationKeyParameters{BaseTwoDecomposition: utils.Pointy(7)}

		// RLWE secret for the samples
		skLWE := rlwe.NewKeyGenerator(paramsLWE).GenSecretKeyNew()

		// Secret of the RGSW ciphertexts encrypting the bits of skLWE
		skBR := rlwe.NewKeyGenerator(paramsBR).GenSecretKeyNew()

		// Collection of RGSW ciphertexts encrypting the bits of skLWE under skBR
		BRK := GenEvaluationKeyNew(paramsBR, skBR, paramsLWE, skLWE, evkParams)

		// Random LWE mask mod 2N with odd coefficients
		a := make([]uint64, paramsLWE.N())
		mask := uint64(2*paramsLWE.N() - 1)
		for i := range a {
			ai := sampling.RandUint64() & mask
			if ai&1 == 0 && ai != 0 {
				ai ^= 1
			}
			a[i] = ai
		}

		acc := rlwe.NewCiphertext(paramsBR, 1, paramsBR.MaxLevel())

		// Evaluator for the Blind Rotation evaluation
		eval := NewEvaluator(paramsBR, paramsLWE)

		b.ResetTimer()

		for i := 0; i < b.N; i++ {
			if err := eval.BlindRotateCore(a, acc, BRK); err != nil {
				panic(err)
			}
		}
	})
}
