package dckks

import (
	"github.com/ldsec/lattigo/v2/ckks"
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/utils"
)

type dckksContext struct {
	params *ckks.Parameters

	n uint64

	ringQ  *ring.Ring
	ringP  *ring.Ring
	ringQP *ring.Ring

	alpha uint64
	beta  uint64
}

func newDckksContext(params *ckks.Parameters) (context *dckksContext) {

	context = new(dckksContext)

	context.params = params.Copy()

	context.n = params.N()

	context.alpha = params.Alpha()
	context.beta = params.Beta()

	var err error
	if context.ringQ, err = ring.NewRing(params.N(), params.Qi()); err != nil {
		panic(err)
	}

	if context.ringP, err = ring.NewRing(params.N(), params.Pi()); err != nil {
		panic(err)
	}

	if context.ringQP, err = ring.NewRing(params.N(), append(params.Qi(), params.Pi()...)); err != nil {
		panic(err)
	}

	return
}

// NewCRPGenerator creates a new deterministic random polynomial generator.
func NewCRPGenerator(params *ckks.Parameters, key []byte) *ring.UniformSampler {
	ctx := newDckksContext(params)
	prng, err := utils.NewKeyedPRNG(key)
	if err != nil {
		panic(err)
	}
	return ring.NewUniformSampler(prng, ctx.ringQP)
}
