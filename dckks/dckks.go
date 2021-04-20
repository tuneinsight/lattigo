package dckks

import (
	"github.com/ldsec/lattigo/v2/ckks"
	"github.com/ldsec/lattigo/v2/ring"
)

type dckksContext struct {
	params ckks.Parameters

	n uint64

	ringQ  *ring.Ring
	ringP  *ring.Ring
	ringQP *ring.Ring

	alpha uint64
	beta  uint64
}

func newDckksContext(params ckks.Parameters) (context *dckksContext) {

	context = new(dckksContext)

	context.params = params

	context.n = params.N()

	context.alpha = params.Alpha()
	context.beta = params.Beta()

	context.ringQ = params.RingQ()
	context.ringP = params.RingP()
	context.ringQP = params.RingQP()

	return
}
