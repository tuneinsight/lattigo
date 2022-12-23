package rlwe

import (
	"fmt"

	"github.com/tuneinsight/lattigo/v4/ring"
)

const (
	// XsUniformTernary is the standard deviation of a ternary key with uniform distribution
	XsUniformTernary = ring.StandardDeviation(0.816496580927726)

	// DefaultNoise is the default standard deviation of the error
	DefaultNoise = ring.StandardDeviation(3.2)

	// DefaultNoiseBound is the default bound (in number of standar deviation) of the noise bound
	DefaultNoiseBound = 19
)

// DefaultXe is the default discret Gaussian distribution.
var DefaultXe = ring.DiscreteGaussianDistribution{Sigma: DefaultNoise, Bound: DefaultNoiseBound}

var DefaultXs = ring.TernaryDistribution{P: 1 / 3.0}

// LWEParameters is a struct
type LWEParameters struct {
	LogN  int
	LogQP float64
	Xe    ring.StandardDeviation
	Xs    ring.StandardDeviation
}

// HomomorphicStandardUSVP128 stores 128-bit secures parameters according to the
// homomorphic encryption standard
var HomomorphicStandardUSVP128 = map[int]LWEParameters{
	10: LWEParameters{LogN: 10, LogQP: 27, Xs: XsUniformTernary, Xe: DefaultNoise},
	11: LWEParameters{LogN: 11, LogQP: 54, Xs: XsUniformTernary, Xe: DefaultNoise},
	12: LWEParameters{LogN: 12, LogQP: 109, Xs: XsUniformTernary, Xe: DefaultNoise},
	13: LWEParameters{LogN: 13, LogQP: 218, Xs: XsUniformTernary, Xe: DefaultNoise},
	14: LWEParameters{LogN: 14, LogQP: 438, Xs: XsUniformTernary, Xe: DefaultNoise},
	15: LWEParameters{LogN: 15, LogQP: 881, Xs: XsUniformTernary, Xe: DefaultNoise},
}

// CheckSecurityForHomomorphicStandardUSVP128 checks if parameters are compliant with the
// HomomorphicStandardUSVP128 security parameters.
func CheckSecurityForHomomorphicStandardUSVP128(params LWEParameters) (err error) {

	if refParams, ok := HomomorphicStandardUSVP128[params.LogN]; ok {
		// We allow a small slack
		if params.LogQP > refParams.LogQP+0.1 {
			return fmt.Errorf("warning: parameters do not comply with the HE Standard 128-bit security: LogQP %f > %f for LogN %d", params.LogQP, refParams.LogQP, params.LogN)
		}

		if params.Xs < refParams.Xs {
			return fmt.Errorf("warning: parameters do not comply with the HE Standard 128-bit security: Xs %f < %f", params.Xs, refParams.Xs)
		}

		if params.Xe < refParams.Xe {
			return fmt.Errorf("warning: parameters do not comply with the HE Standard 128-bit security: Xe %f < %f", params.Xs, refParams.Xe)
		}

	} else {
		return fmt.Errorf("warning: parameters do not comply with the HE Standard 128-bit security: LogN=%d is not supported", params.LogN)
	}

	return nil
}
