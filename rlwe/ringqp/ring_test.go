package ringqp

import (
	"testing"

	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/utils/buffer"
	"github.com/tuneinsight/lattigo/v4/utils/sampling"

	"github.com/stretchr/testify/require"
)

func TestRingQP(t *testing.T) {
	LogN := 10
	ringQ, err := ring.NewRing(1<<LogN, ring.Qi60[:4])
	require.NoError(t, err)

	ringP, err := ring.NewRing(1<<LogN, ring.Pi60[:4])
	require.NoError(t, err)

	ringQP := Ring{ringQ, ringP}

	prng, err := sampling.NewPRNG()
	require.NoError(t, err)

	usampler := NewUniformSampler(prng, ringQP)

	t.Run("Binary/Poly", func(t *testing.T) {
		p := usampler.ReadNew()
		buffer.TestInterfaceWriteAndRead(t, &p)
	})

	t.Run("Binary/PolyVector", func(t *testing.T) {

		polys := make([]Poly, 4)

		for i := range polys {
			polys[i] = usampler.ReadNew()
		}

		pv := new(PolyVector)
		pv.Set(polys)

		buffer.TestInterfaceWriteAndRead(t, pv)
	})

	t.Run("Binary/PolyMatrix", func(t *testing.T) {

		polys := make([][]Poly, 4)

		for i := range polys {
			polys[i] = make([]Poly, 4)

			for j := range polys {
				polys[i][j] = usampler.ReadNew()
			}
		}

		pm := new(PolyMatrix)
		pm.Set(polys)

		buffer.TestInterfaceWriteAndRead(t, pm)
	})

}
