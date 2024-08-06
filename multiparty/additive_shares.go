package multiparty

import (
	"math/big"

	"github.com/tuneinsight/lattigo/v6/ring"
)

// AdditiveShare is a type for storing additively shared values in Z_Q[X] (RNS domain).
type AdditiveShare struct {
	Value ring.Poly
}

// AdditiveShareBigint is a type for storing additively shared values
// in Z (positional domain).
type AdditiveShareBigint struct {
	Value []*big.Int
}

// NewAdditiveShare instantiates a new additive share struct for the ring defined
// by the given parameters at maximum level.
func NewAdditiveShare(r *ring.Ring) AdditiveShare {
	return AdditiveShare{Value: r.NewPoly()}
}

// NewAdditiveShareBigint instantiates a new additive share struct composed of n [big.Int] elements.
func NewAdditiveShareBigint(n int) AdditiveShareBigint {
	v := make([]*big.Int, n)
	for i := range v {
		v[i] = new(big.Int)
	}
	return AdditiveShareBigint{Value: v}
}
