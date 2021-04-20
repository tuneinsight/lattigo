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

	context.params = params

	context.n = params.N()

	context.alpha = params.Alpha()
	context.beta = params.Beta()

	context.ringQ = params.RingQ()
	context.ringP = params.RingP()
	context.ringQP = params.RingQP()

	return
}
