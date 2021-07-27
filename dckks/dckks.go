package dckks

import (
	"github.com/ldsec/lattigo/v2/ckks"
	"github.com/ldsec/lattigo/v2/ring"
)

type Context struct {
	params *ckks.Parameters

	n int

	RingQ  *ring.Ring
	RingP  *ring.Ring
	RingQP *ring.Ring

	Alpha int
	Beta  int
}

func NewContext(params *ckks.Parameters) (context *Context) {

	context = new(Context)

	context.params = params.Copy()

	context.n = params.N()

	context.Alpha = params.Alpha()
	context.Beta = params.Beta()

	var err error
	if context.RingQ, err = ring.NewRing(params.N(), params.Qi()); err != nil {
		panic(err)
	}

	if context.RingP, err = ring.NewRing(params.N(), params.Pi()); err != nil {
		panic(err)
	}

	if context.RingQP, err = ring.NewRing(params.N(), append(params.Qi(), params.Pi()...)); err != nil {
		panic(err)
	}

	return
}
