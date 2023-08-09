package float

import (
	"fmt"
	"math"
	"sort"
	"testing"

	"github.com/stretchr/testify/require"
	"github.com/tuneinsight/lattigo/v4/utils/bignum"
)

/*
func TestMinimaxApprox(t *testing.T) {
	// Precision of the floating point arithmetic
	prec := uint(512)

	// 2^{-logalpha} distinguishing ability
	logalpha := int(10)

	// Degrees of each minimax polynomial
	deg := []int{8, 8, 18, 32}

	GenSignPoly(prec, logalpha, deg)
}
*/

func TestMinimaxCompositeSignPolys30bits(t *testing.T) {

	keys := make([]int, len(SingPoly30String))

	idx := 0
	for k := range SingPoly30String {
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

func TestMinimaxCompositeSignPolys20bits(t *testing.T) {

	keys := make([]int, len(SingPoly20String))

	idx := 0
	for k := range SingPoly20String {
		keys[idx] = k
		idx++
	}

	sort.Ints(keys)

	for _, alpha := range keys[:] {

		polys, err := GetSignPoly20Polynomials(alpha)
		require.NoError(t, err)

		xPos := bignum.NewFloat(math.Exp2(-float64(alpha)), 53)
		xNeg := bignum.NewFloat(-math.Exp2(-float64(alpha)), 53)

		for _, poly := range polys {
			xPos = poly.Evaluate(xPos)[0]
			xNeg = poly.Evaluate(xNeg)[0]
		}

		xPosF64, _ := xPos.Float64()
		xNegF64, _ := xNeg.Float64()

		fmt.Println(alpha, math.Log2(1-xPosF64), math.Log2(1+xNegF64))

		require.Greater(t, -20.0, math.Log2(1-xPosF64))
		require.Greater(t, -20.0, math.Log2(1+xNegF64))
	}
}
