package dckks

import (
	"github.com/ldsec/lattigo/ckks"
	"github.com/ldsec/lattigo/ring"
	"math"
)

type dckksContext struct {
	params *ckks.Parameters

	n uint64

	gaussianSampler *ring.KYSampler

	contextQ  *ring.Context
	contextP  *ring.Context
	contextQP *ring.Context

	alpha uint64
	beta  uint64
}

func newDckksContext(params *ckks.Parameters) (context *dckksContext) {

	context = new(dckksContext)

	context.params = params.Copy()

	n := uint64(1 << params.LogN)

	moduliQ, moduliP := ckks.GenModuli(params)

	context.n = n

	context.alpha = uint64(len(moduliP))
	context.beta = uint64(math.Ceil(float64(len(moduliQ)) / float64(context.alpha)))

	if context.contextQ, err = ring.NewContextWithParams(n, moduliQ); err != nil {
		panic(err)
	}

	if context.contextP, err = ring.NewContextWithParams(n, moduliP); err != nil {
		panic(err)
	}

	if context.contextQP, err = ring.NewContextWithParams(n, append(moduliQ, moduliP...)); err != nil {
		panic(err)
	}

	context.gaussianSampler = context.contextQP.NewKYSampler(params.Sigma, int(params.Sigma*6))

	return
}

func NewCRPGenerator(params *ckks.Parameters, key []byte) *ring.CRPGenerator {
	ctx := newDckksContext(params)
	return ring.NewCRPGenerator(key, ctx.contextQP)
}
