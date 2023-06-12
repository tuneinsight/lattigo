package rlwe

import (
	"math"

	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/utils"
)

type distribution struct {
	params   ring.DistributionParameters
	std      float64
	bounds   [2]float64
	absBound float64
	density  float64
}

func newDistribution(params ring.DistributionParameters, logN int, logQP float64) (d distribution) {
	d.params = params
	switch params := params.(type) {
	case ring.DiscreteGaussian:
		d.std = params.Sigma
		d.bounds = [2]float64{-params.Bound, params.Bound}
		d.absBound = params.Bound
		d.density = 1 - utils.Min(1/math.Sqrt(2*math.Pi)*params.Sigma, 1)
	case ring.Ternary:
		N := math.Exp2(float64(logN))
		if params.P != 0 {
			d.std = math.Sqrt(1 - params.P)
			d.density = params.P
		} else {
			d.std = math.Sqrt(float64(params.H) / (math.Exp2(float64(logN)) - 1))
			d.density = float64(params.H) / N
		}
		d.bounds = [2]float64{-1, 1}
		d.absBound = 1
	case ring.Uniform:
		d.std = math.Exp2(logQP) / math.Sqrt(12.0)
		d.bounds = [2]float64{-math.Exp2(logQP - 1), math.Exp2(logQP - 1)}
		d.density = 1 - (1 / (math.Exp2(logQP) + 1))
	default:
		panic("invalid dist")
	}
	return
}
