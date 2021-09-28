package dbfv

import (
	"github.com/ldsec/lattigo/bfv"
	"github.com/ldsec/lattigo/ring"
)

type dbfvContext struct {
	params *bfv.Parameters

	// Polynomial degree
	n uint64

	// floor(Q/T) mod each Qi in Montgomery form
	deltaMont []uint64

	// Ternary and Gaussian samplers
	gaussianSampler *ring.KYSampler

	// Polynomial contexts
	contextT  *ring.Context
	contextQ  *ring.Context
	contextP  *ring.Context
	contextQP *ring.Context
}

func newDbfvContext(params *bfv.Parameters) *dbfvContext {

	if !params.IsValid() {
		panic("cannot newDbfvContext : params not valid (check if they where generated properly)")
	}

	LogN := params.LogN
	n := uint64(1 << LogN)

	contextT, err := ring.NewContextWithParams(n, []uint64{params.T})
	if err != nil {
		panic(err)
	}

	contextQ, err := ring.NewContextWithParams(n, params.Qi)
	if err != nil {
		panic(err)
	}

	contextP, err := ring.NewContextWithParams(n, params.Pi)
	if err != nil {
		panic(err)
	}

	contextQP, err := ring.NewContextWithParams(n, append(params.Qi, params.Pi...))
	if err != nil {
		panic(err)
	}

	deltaMont := bfv.GenLiftParams(contextQ, params.T)

	gaussianSampler := contextQP.NewKYSampler(params.Sigma, int(6*params.Sigma))

	return &dbfvContext{
		params:          params.Copy(),
		n:               n,
		deltaMont:       deltaMont,
		gaussianSampler: gaussianSampler,
		contextT:        contextT,
		contextQ:        contextQ,
		contextP:        contextP,
		contextQP:       contextQP,
	}
}

func NewCRPGenerator(params *bfv.Parameters, key []byte) *ring.CRPGenerator {
	ctx := newDbfvContext(params)
	return ring.NewCRPGenerator(key, ctx.contextQP)
}
