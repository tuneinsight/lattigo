package ring

import (
	"fmt"
	"reflect"

	"github.com/tuneinsight/lattigo/v4/utils"
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
	Read(pOut *Poly)
	AtLevel(level int) Sampler
	Read(pol *Poly)
	ReadLvl(level int, pol *Poly)
	ReadNew() (pol *Poly)
	ReadLvlNew(level int) (pol *Poly)
	ReadAndAddLvl(level int, pol *Poly)
}

func NewSampler(prng utils.PRNG, baseRing *Ring, X Distribution, montgomery bool) Sampler {
	switch X := X.(type) {
	case *DiscreteGaussianDistribution:
		return NewGaussianSampler(prng, baseRing, X, montgomery)
	case *TernaryDistribution:
		return NewTernarySampler(prng, baseRing, X, montgomery)
	case *UniformDistribution:
		return NewUniformSampler(prng, baseRing)
	default:
		panic(fmt.Sprintf("Invalid distribution: want *ring.DiscretGaussian, *ring.UniformTernary, *ring.SparseTernary or *ring.Uniform but have %s", reflect.TypeOf(X)))
	}
}
