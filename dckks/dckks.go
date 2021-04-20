package dckks

import (
	"github.com/ldsec/lattigo/v2/ckks"
	"github.com/ldsec/lattigo/v2/ring"
)

type dckksContext struct {
	params *ckks.ParametersStruct

	n int

	ringQ  *ring.Ring
	ringP  *ring.Ring
	ringQP *ring.Ring

	alpha int
	beta  int
}

func newDckksContext(params *ckks.ParametersStruct) (context *dckksContext) {

	context = new(dckksContext)

	context.params = params.Copy()

	context.n = params.N()

	context.alpha = params.Alpha()
	context.beta = params.Beta()

	var err error
	if context.ringQ, err = ring.NewRing(params.N(), params.Q()); err != nil {
		panic(err)
	}

	if context.ringP, err = ring.NewRing(params.N(), params.P()); err != nil {
		panic(err)
	}

	if context.ringQP, err = ring.NewRing(params.N(), append(params.Q(), params.P()...)); err != nil {
		panic(err)
	}

	return
}
