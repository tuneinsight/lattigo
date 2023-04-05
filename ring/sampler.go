package ring

import (
	"fmt"

	"github.com/tuneinsight/lattigo/v4/ring/distribution"
	"github.com/tuneinsight/lattigo/v4/utils/sampling"
)

const precision = uint64(56)

type baseSampler struct {
	prng     sampling.PRNG
	baseRing *Ring
}

// AtLevel returns an instance of the target base sampler that operates at the target level.
// This instance is not thread safe and cannot be used concurrently to the base instance.
func (b *baseSampler) AtLevel(level int) baseSampler {
	return baseSampler{
		prng:     b.prng,
		baseRing: b.baseRing.AtLevel(level),
	}
}

// Sampler is an interface for random polynomial samplers.
// It has a single Read method which takes as argument the polynomial to be
// populated according to the Sampler's distribution.
type Sampler interface {
	Read(pol *Poly)
	ReadNew() (pol *Poly)
	ReadAndAdd(pol *Poly)
	AtLevel(level int) Sampler
}

func NewSampler(prng sampling.PRNG, baseRing *Ring, X distribution.Distribution, montgomery bool) Sampler {
	switch X := X.(type) {
	case *distribution.DiscreteGaussian:
		return NewGaussianSampler(prng, baseRing, *X, montgomery)
	case *distribution.Ternary:
		return NewTernarySampler(prng, baseRing, *X, montgomery)
	case *distribution.Uniform:
		return NewUniformSampler(prng, baseRing)
	default:
		panic(fmt.Sprintf("Invalid distribution: want *ring.DiscreteGaussianDistribution, *ring.TernaryDistribution or *ring.UniformDistribution but have %T", X))
	}
}
