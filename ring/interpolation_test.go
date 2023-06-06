package ring

import (
	"testing"

	"github.com/stretchr/testify/assert"
)

func TestInterpolation(t *testing.T) {
	t.Run("Roots", func(t *testing.T) {
		var T uint64 = 65537
		roots := make([]uint64, 22)
		for i := range roots {
			roots[i] = uint64(i)
		}

		itp, err := NewInterpolator(len(roots), T)
		assert.Nil(t, err)

		coeffs := itp.Interpolate(roots)
		for _, alpha := range roots {
			assert.Equal(t, uint64(0), EvalPolyModP(alpha, coeffs, T))
		}
	})

	t.Run("Lagrange", func(t *testing.T) {
		var T uint64 = 65537
		n := 512
		x := make([]uint64, n+1)
		y := make([]uint64, n+1)

		for i := 0; i < n>>1; i++ {
			x[i] = T - uint64(n>>1-i)
			y[i] = 0
		}

		y[n>>1] = 1

		for i := 1; i < n>>1+1; i++ {
			x[i+n>>1] = uint64(i)
			y[i+n>>1] = 0
		}

		itp, err := NewInterpolator(len(x), T)
		assert.Nil(t, err)

		coeffs, err := itp.Lagrange(x, y)
		assert.Nil(t, err)

		for i := range x {
			assert.Equal(t, y[i], EvalPolyModP(x[i], coeffs, T))
		}
	})
}
