package dckks

import (
	"math"

	"github.com/tuneinsight/lattigo/ckks"
	"github.com/tuneinsight/lattigo/ring"
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

	if !params.IsValid() {
		panic("cannot newDckksContext : params not valid (check if they where generated properly)")
	}

	context = new(dckksContext)

	context.params = params.Copy()

	n := uint64(1 << params.LogN)

	context.n = n

	context.alpha = uint64(len(params.Pi))
	context.beta = uint64(math.Ceil(float64(len(params.Qi)) / float64(context.alpha)))

	if context.contextQ, err = ring.NewContextWithParams(n, params.Qi); err != nil {
		panic(err)
	}

	if context.contextP, err = ring.NewContextWithParams(n, params.Pi); err != nil {
		panic(err)
	}

	if context.contextQP, err = ring.NewContextWithParams(n, append(params.Qi, params.Pi...)); err != nil {
		panic(err)
	}

	context.gaussianSampler = context.contextQP.NewKYSampler(params.Sigma, int(params.Sigma*6))

	return
}

func NewCRPGenerator(params *ckks.Parameters, key []byte) *ring.CRPGenerator {
	ctx := newDckksContext(params)
	return ring.NewCRPGenerator(key, ctx.contextQP)
}
