package dckks

import (
	"github.com/tuneinsight/lattigo/ring"
)

func randomFloat(max float64) float64 {
	return float64(ring.RandUniform(0x20000000000000, 0x3fffffffffffff))*(2*max/float64(0x20000000000000)) - max
}

func randomComplex(max float64) complex128 {
	return complex(randomFloat(max), randomFloat(max))
}
