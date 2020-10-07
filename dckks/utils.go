package dckks

import (
	"github.com/ldsec/lattigo/ring"
	"github.com/ldsec/lattigo/utils"
)

func randomFloat(prng utils.PRNG, max float64) float64 {
	return float64(ring.RandUniform(prng, 0x20000000000000, 0x3fffffffffffff))*(2*max/float64(0x20000000000000)) - max
}

func randomComplex(prng utils.PRNG, max float64) complex128 {
	return complex(randomFloat(prng, max), randomFloat(prng, max))
}
