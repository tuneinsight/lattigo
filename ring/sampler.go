package ring

import (
	"encoding/json"
	"fmt"

	"github.com/tuneinsight/lattigo/v6/utils/sampling"
)

const (
	discreteGaussianName = "DiscreteGaussian"
	ternaryDistName      = "Ternary"
	uniformDistName      = "Uniform"
)

// Sampler is an interface for random polynomial samplers.
// It has a single Read method which takes as argument the polynomial to be
// populated according to the Sampler's distribution.
type Sampler interface {
	Read(pol Poly)
	ReadNew() (pol Poly)
	ReadAndAdd(pol Poly)
	AtLevel(level int) Sampler
}

// DistributionParameters is an interface for distribution
// parameters in the ring.
// There are three implementation of this interface:
//   - DiscreteGaussian for sampling polynomials with discretized
//     gaussian coefficient of given standard deviation and bound.
//   - Ternary for sampling polynomials with coefficients in [-1, 1].
//   - Uniform for sampling polynomial with uniformly random
//     coefficients in the ring.
type DistributionParameters interface {
	// Type returns a string representation of the distribution name.
	Type() string
	mustBeDist()
}

// DiscreteGaussian represents the parameters of a
// discrete Gaussian distribution with standard
// deviation Sigma and bounds [-Bound, Bound].
type DiscreteGaussian struct {
	Sigma float64
	Bound float64
}

// Ternary represent the parameters of a distribution with coefficients
// in [-1, 0, 1]. Only one of its field must be set to a non-zero value:
//
//   - If P is set, each coefficient in the polynomial is sampled in [-1, 0, 1]
//     with probabilities [0.5*P, 1-P, 0.5*P].
//   - if H is set, the coefficients are sampled uniformly in the set of ternary
//     polynomials with H non-zero coefficients (i.e., of hamming weight H).
type Ternary struct {
	P float64
	H int
}

// Uniform represents the parameters of a uniform distribution
// i.e., with coefficients uniformly distributed in the given ring.
type Uniform struct{}

func NewSampler(prng sampling.PRNG, baseRing *Ring, X DistributionParameters, montgomery bool) (Sampler, error) {
	switch X := X.(type) {
	case DiscreteGaussian:
		return NewGaussianSampler(prng, baseRing, X, montgomery), nil
	case Ternary:
		return NewTernarySampler(prng, baseRing, X, montgomery)
	case Uniform:
		return NewUniformSampler(prng, baseRing), nil
	default:
		return nil, fmt.Errorf("invalid distribution: want ring.DiscreteGaussianDistribution, ring.TernaryDistribution or ring.UniformDistribution but have %T", X)
	}
}

type baseSampler struct {
	prng     sampling.PRNG
	baseRing *Ring
}

type randomBuffer struct {
	randomBufferN []byte
	ptr           int
}

func newRandomBuffer() *randomBuffer {
	return &randomBuffer{
		randomBufferN: make([]byte, 1024),
	}
}

// AtLevel returns an instance of the target base sampler that operates at the target level.
// This instance is not thread safe and cannot be used concurrently to the base instance.
func (b baseSampler) AtLevel(level int) *baseSampler {
	return &baseSampler{
		prng:     b.prng,
		baseRing: b.baseRing.AtLevel(level),
	}
}

func (d DiscreteGaussian) Type() string {
	return discreteGaussianName
}

func (d DiscreteGaussian) MarshalJSON() ([]byte, error) {
	return json.Marshal(struct {
		Type         string
		Sigma, Bound float64 `json:",omitempty"`
	}{d.Type(), d.Sigma, d.Bound})
}

func (d DiscreteGaussian) mustBeDist() {}

func (d Ternary) Type() string {
	return ternaryDistName
}

func (d Ternary) MarshalJSON() ([]byte, error) {
	return json.Marshal(struct {
		Type string
		P    float64 `json:",omitempty"`
		H    int     `json:",omitempty"`
	}{Type: d.Type(), P: d.P, H: d.H})
}

func (d Ternary) mustBeDist() {}

func (d Uniform) Type() string {
	return uniformDistName
}

func (d Uniform) MarshalJSON() ([]byte, error) {
	return json.Marshal(struct {
		Type string
	}{Type: d.Type()})
}

func (d Uniform) mustBeDist() {}

func getFloatFromMap(distDef map[string]interface{}, key string) (float64, error) {
	val, hasVal := distDef[key]
	if !hasVal {
		return 0, fmt.Errorf("map specifies no value for %s", key)
	}
	f, isFloat := val.(float64)
	if !isFloat {
		return 0, fmt.Errorf("value for key %s in map should be of type float", key)
	}
	return f, nil
}

func getIntFromMap(distDef map[string]interface{}, key string) (int, error) {
	val, hasVal := distDef[key]
	if !hasVal {
		return 0, fmt.Errorf("map specifies no value for %s", key)
	}
	f, isNumeric := val.(float64)
	if !isNumeric && f == float64(int(f)) {
		return 0, fmt.Errorf("value for key %s in map should be an integer", key)
	}
	return int(f), nil
}

func ParametersFromMap(distDef map[string]interface{}) (DistributionParameters, error) {
	distTypeVal, specified := distDef["Type"]
	if !specified {
		return nil, fmt.Errorf("map specifies no distribution type")
	}
	distTypeStr, isString := distTypeVal.(string)
	if !isString {
		return nil, fmt.Errorf("value for key Type of map should be of type string")
	}
	switch distTypeStr {
	case uniformDistName:
		return Uniform{}, nil
	case ternaryDistName:
		_, hasP := distDef["P"]
		_, hasH := distDef["H"]

		var (
			p   float64
			h   int
			err error
		)

		// a zero value for both P and H is interpreted as an unset value
		if hasP {
			if p, err = getFloatFromMap(distDef, "P"); err != nil {
				return nil, fmt.Errorf("unable to parse ternary parameters P: %w", err)
			}
			hasP = (p != 0)
		}
		if hasH {
			if h, err = getIntFromMap(distDef, "H"); err != nil {
				return nil, fmt.Errorf("unable to parse ternary parameters H: %w", err)
			}
			hasH = (h != 0)
		}
		if (hasP && hasH) || (!hasP && !hasH) {
			return nil, fmt.Errorf("exactly one of the fields P or H need to be set")
		}

		return Ternary{P: p, H: h}, nil
	case discreteGaussianName:
		sigma, errSigma := getFloatFromMap(distDef, "Sigma")
		if errSigma != nil {
			return nil, errSigma
		}
		bound, errBound := getFloatFromMap(distDef, "Bound")
		if errBound != nil {
			return nil, errBound
		}
		return DiscreteGaussian{Sigma: sigma, Bound: bound}, nil
	default:
		return nil, fmt.Errorf("distribution type %s does not exist", distTypeStr)
	}
}
