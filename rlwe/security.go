package rlwe

import (
	"fmt"

	"github.com/tuneinsight/lattigo/v4/ring/distribution"
)

const (
	// XsUniformTernary is the standard deviation of a ternary key with uniform distribution
	XsUniformTernary = distribution.StandardDeviation(0.816496580927726)

	// DefaultNoise is the default standard deviation of the error
	DefaultNoise = distribution.StandardDeviation(3.2)

	// DefaultNoiseBound is the default bound (in number of standar deviation) of the noise bound
	DefaultNoiseBound = 19
)

// DefaultXe is the default discret Gaussian distribution.
var DefaultXe = distribution.DiscreteGaussian{Sigma: DefaultNoise, Bound: DefaultNoiseBound}

var DefaultXs = distribution.Ternary{P: 1 / 3.0}

// LWEParameters is a struct
type LWEParameters struct {
	LogN  int
	LogQP float64
	Xe    distribution.StandardDeviation
	Xs    distribution.StandardDeviation
}

func (p *LWEParameters) String() string {
	return fmt.Sprintf("empty\n, %d", 0)
}
