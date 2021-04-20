package dckks

import (
	"github.com/ldsec/lattigo/v2/ckks"
	"github.com/ldsec/lattigo/v2/ring"
)

type dckksContext struct {
	params ckks.Parameters

	n int

	ringQ  *ring.Ring
	ringP  *ring.Ring
	ringQP *ring.Ring

	alpha int
	beta  int
}

func newDckksContext(params ckks.Parameters) (context *dckksContext) {

	context = new(dckksContext)

	context.params = params

	context.n = params.N()

	context.alpha = params.PCount()
	context.beta = params.Beta()

	context.ringQ = params.RingQ()
	context.ringP = params.RingP()
	context.ringQP = params.RingQP()

	return
}
