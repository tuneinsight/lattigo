package ring

import (
	"github.com/tuneinsight/lattigo/v4/utils"
)

const precision = uint64(56)

type baseSampler struct {
	prng     utils.PRNG
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
}
