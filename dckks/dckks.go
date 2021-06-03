package dckks

import (
	"github.com/ldsec/lattigo/v2/ckks"
	"github.com/ldsec/lattigo/v2/ring"
)

type dckksContext struct {
	params *ckks.Parameters

	n int

	ringQ  *ring.Ring
	ringP  *ring.Ring
	ringQP *ring.Ring

	alpha int
	beta  int
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
