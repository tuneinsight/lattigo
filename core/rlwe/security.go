package rlwe

import "github.com/tuneinsight/lattigo/v6/ring"

const (
	// XsUniformTernary is the standard deviation of a ternary key with uniform distribution
	XsUniformTernary = 0.816496580927726 //Sqrt(2/3)

	// DefaultNoise is the default standard deviation of the error
	DefaultNoise = 3.2

	// DefaultNoiseBound is the default bound (in number of standard deviation) of the noise bound
	DefaultNoiseBound = 19.2 // 6*3.2
)

// DefaultXe is the default discrete Gaussian distribution.
var DefaultXe = ring.DiscreteGaussian{Sigma: DefaultNoise, Bound: DefaultNoiseBound}

var DefaultXs = ring.Ternary{P: 2 / 3.0}
