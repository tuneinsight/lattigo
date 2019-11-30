package dbfv

import (
	"math"

	"github.com/ldsec/lattigo/bfv"
	"github.com/ldsec/lattigo/ring"
)

type dbfvContext struct {
	params *bfv.Parameters

	// Polynomial degree
	n uint64

	// Plaintext Modulus
	t uint64

	// floor(Q/T) mod each Qi in Montgomery form
	deltaMont []uint64

	// Ternary and Gaussian samplers
	gaussianSampler *ring.KYSampler

	// Polynomial contexts
	contextT   *ring.Context
	contextQ1  *ring.Context
	contextQ2  *ring.Context
	contextP   *ring.Context
	contextQ1P *ring.Context
	alpha      uint64
	beta       uint64
}

func newDbfvContext(params *bfv.Parameters) *dbfvContext {
	LogN := params.LogN
	n := uint64(1 << LogN)
	t := params.T

	moduliQ1, moduliP, moduliQ2 := bfv.GenModuli(params)

	contextT, err := ring.NewContextWithParams(n, []uint64{t})
	if err != nil {
		panic(err)
	}

	contextQ1, err := ring.NewContextWithParams(n, moduliQ1)
	if err != nil {
		panic(err)
	}

	contextQ2, err := ring.NewContextWithParams(n, moduliQ2)
	if err != nil {
		panic(err)
	}

	contextP, err := ring.NewContextWithParams(n, moduliP)
	if err != nil {
		panic(err)
	}

	contextQ1P, err := ring.NewContextWithParams(n, append(moduliQ1, moduliP...))
	if err != nil {
		panic(err)
	}

	deltaMont := bfv.GenLiftParams(contextQ1, t)

	alpha := uint64(len(moduliP))
	beta := uint64(math.Ceil(float64(len(moduliQ1)) / float64(alpha)))

	gaussianSampler := contextQ1P.NewKYSampler(params.Sigma, int(6*params.Sigma))

	return &dbfvContext{
		params:          params,
		n:               n,
		t:               t,
		deltaMont:       deltaMont,
		gaussianSampler: gaussianSampler,
		contextT:        contextT,
		contextQ1:       contextQ1,
		contextQ2:       contextQ2,
		contextP:        contextP,
		contextQ1P:      contextQ1P,
		alpha:           alpha,
		beta:            beta,
	}
}

func NewCRPGenerator(params *bfv.Parameters, key []byte) *ring.CRPGenerator {
	ctx := newDbfvContext(params)
	return ring.NewCRPGenerator(key, ctx.contextQ1P)
}