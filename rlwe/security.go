package rlwe

import (
	"github.com/tuneinsight/lattigo/v4/ring/distribution"
)

const (
	// XsUniformTernary is the standard deviation of a ternary key with uniform distribution
	XsUniformTernary = 0.816496580927726 //Sqrt(2/3)

	// DefaultNoise is the default standard deviation of the error
	DefaultNoise = 3.2

	// DefaultNoiseBound is the default bound (in number of standar deviation) of the noise bound
	DefaultNoiseBound = 19.2 // 6*3.2
)

// DefaultXe is the default discret Gaussian distribution.
var DefaultXe = distribution.DiscreteGaussian{Sigma: DefaultNoise, Bound: DefaultNoiseBound}

var DefaultXs = distribution.Ternary{P: 1 / 3.0}
