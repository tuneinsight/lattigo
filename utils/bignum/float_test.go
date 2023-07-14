package bignum

import (
	"math"
	"math/big"
	"testing"

	"github.com/stretchr/testify/require"
)

func TestFloat(t *testing.T) {
	testFunc1("Sin", 1.4142135623730951, math.Sin, Sin, 1e-15, t)
	testFunc1("Cos", 1.4142135623730951, math.Cos, Cos, 1e-15, t)
	testFunc1("Log", 1.4142135623730951, math.Log, Log, 1e-15, t)
	testFunc1("Exp", 1.4142135623730951, math.Exp, Exp, 1e-15, t)
	testFunc2("Pow", 2, 1.4142135623730951, math.Pow, Pow, 1e-15, t)
	testFunc1("SinH", 1.4142135623730951, math.Sinh, SinH, 1e-15, t)
	testFunc1("TanH", 1.4142135623730951, math.Tanh, TanH, 1e-15, t)
}

func testFunc1(name string, x float64, f func(x float64) (y float64), g func(x *big.Float) (y *big.Float), delta float64, t *testing.T) {
	t.Run(name, func(t *testing.T) {
		y, _ := g(NewFloat(x, 53)).Float64()
		require.InDelta(t, f(x), y, delta)
	})
}

func testFunc2(name string, x, e float64, f func(x, e float64) (y float64), g func(x, e *big.Float) (y *big.Float), delta float64, t *testing.T) {
	t.Run(name, func(t *testing.T) {
		y, _ := g(NewFloat(x, 53), NewFloat(e, 53)).Float64()
		require.InDelta(t, f(x, e), y, delta)
	})
}
