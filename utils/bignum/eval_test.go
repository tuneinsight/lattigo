package bignum

import (
	"math/big"
	"testing"

	"github.com/stretchr/testify/require"
)

func TestMonomialEval(t *testing.T) {
	newFloat := func(f float64) *big.Float { return new(big.Float).SetPrec(53).SetFloat64(f) }

	// 1 + 2x at x = 3 = 7: the highest-degree coefficient must not be dropped.
	got := MonomialEval(newFloat(3), []*big.Float{newFloat(1), newFloat(2)})
	require.Zerof(t, got.Cmp(newFloat(7)), "1 + 2x at x=3: got %v, want 7", got)

	// 1 + 2x + 3x^2 at x = 2 = 17.
	got = MonomialEval(newFloat(2), []*big.Float{newFloat(1), newFloat(2), newFloat(3)})
	require.Zerof(t, got.Cmp(newFloat(17)), "1 + 2x + 3x^2 at x=2: got %v, want 17", got)

	// A constant polynomial must return the constant without panicking.
	got = MonomialEval(newFloat(9), []*big.Float{newFloat(5)})
	require.Zerof(t, got.Cmp(newFloat(5)), "constant poly: got %v, want 5", got)
}
