package polynomial

import (
	"testing"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/utils/bignum"
	"github.com/tuneinsight/lattigo/v6/utils/buffer"
	"github.com/tuneinsight/lattigo/v6/utils/sampling"
)

func TestPowerBasis(t *testing.T) {
	t.Run("WriteAndRead", func(t *testing.T) {
		var err error
		var params rlwe.Parameters
		if params, err = rlwe.NewParametersFromLiteral(rlwe.ParametersLiteral{
			LogN: 10,
			Q:    []uint64{0x200000440001, 0x7fff80001},
			P:    []uint64{0x3ffffffb80001, 0x4000000800001},
		}); err != nil {
			t.Fatal(err)
		}

		levelQ := params.MaxLevelQ()

		prng, _ := sampling.NewPRNG()

		ct := rlwe.NewCiphertextRandom(prng, params, 1, levelQ)

		basis := NewPowerBasis(ct, bignum.Chebyshev)

		basis.Value[2] = rlwe.NewCiphertextRandom(prng, params, 1, levelQ)
		basis.Value[3] = rlwe.NewCiphertextRandom(prng, params, 2, levelQ)
		basis.Value[4] = rlwe.NewCiphertextRandom(prng, params, 1, levelQ)
		basis.Value[8] = rlwe.NewCiphertextRandom(prng, params, 1, levelQ)

		buffer.RequireSerializerCorrect(t, &basis)
	})
}
