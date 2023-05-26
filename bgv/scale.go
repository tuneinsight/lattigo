package bgv

import (
	"math/big"

	"github.com/tuneinsight/lattigo/v4/rlwe"
)

func MulScale(params Parameters, a, b rlwe.Scale, level int, invariant bool) (c rlwe.Scale) {
	c = a.Mul(b)
	if invariant {
		qModTNeg := new(big.Int).Mod(params.RingQ().ModulusAtLevel[level], new(big.Int).SetUint64(params.T())).Uint64()
		qModTNeg = params.T() - qModTNeg
		c = c.Div(params.NewScale(qModTNeg))
	}

	return
}
