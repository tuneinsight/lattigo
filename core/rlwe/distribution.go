package rlwe

import (
	"math"

	"github.com/tuneinsight/lattigo/v6/ring"
)

type Distribution struct {
	ring.DistributionParameters
	Std      float64
	AbsBound float64
}

func NewDistribution(params ring.DistributionParameters, logN int) (d Distribution) {
	d.DistributionParameters = params
	switch params := params.(type) {
	case ring.DiscreteGaussian:
		d.Std = params.Sigma
		d.AbsBound = params.Bound
	case ring.Ternary:
		if params.P != 0 {
			d.Std = math.Sqrt(1 - params.P)
		} else {
			d.Std = math.Sqrt(float64(params.H) / (math.Exp2(float64(logN)) - 1))
		}
		d.AbsBound = 1
	default:
		// Sanity check
		panic("invalid dist")
	}
	return
}
