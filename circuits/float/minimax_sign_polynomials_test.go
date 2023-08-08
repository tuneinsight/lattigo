package float

import (
	"math"
	"sort"
	"testing"

	"github.com/stretchr/testify/require"
	"github.com/tuneinsight/lattigo/v4/utils/bignum"
)

func TestMinimaxCompositeSignPolys30bits(t *testing.T) {

	keys := make([]int, len(SignPolys30))

	idx := 0
	for k := range SignPolys30 {
		keys[idx] = k
		idx++
	}

	sort.Ints(keys)

	for _, alpha := range keys[:] {

		polys, err := GetSignPoly30Polynomials(alpha)
		require.NoError(t, err)

		xPos := bignum.NewFloat(math.Exp2(-float64(alpha)), 53)
		xNeg := bignum.NewFloat(-math.Exp2(-float64(alpha)), 53)

		for _, poly := range polys {
			xPos = poly.Evaluate(xPos)[0]
			xNeg = poly.Evaluate(xNeg)[0]
		}

		xPosF64, _ := xPos.Float64()
		xNegF64, _ := xNeg.Float64()

		require.Greater(t, -30.0, math.Log2(1-xPosF64))
		require.Greater(t, -30.0, math.Log2(1+xNegF64))
	}
}
