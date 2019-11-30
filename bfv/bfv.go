// Package bfv implements a RNS-accelerated Fan-Vercauteren version of Brakerski's scale invariant homomorphic encryption scheme. It provides modular arithmetic over the integers.
package bfv

import (
	"github.com/ldsec/lattigo/ring"
	"math"
)

// GaloisGen is ... [FIXME]
const GaloisGen uint64 = 5

// bfvContext is a struct which contains all the elements required to instantiate the BFV Scheme. This includes the parameters (N, plaintext modulus, ciphertext modulus,
// sampling, polynomial contexts and other parameters required for the homomorphic operations).
type bfvContext struct {
	params *Parameters

	// Polynomial degree
	n uint64

	// Plaintext Modulus
	t uint64

	logQ uint64

	// Ternary and Gaussian samplers
	sigma           float64
	gaussianSampler *ring.KYSampler

	// Polynomial contexts
	contextT   *ring.Context
	contextQ1  *ring.Context
	contextQ2  *ring.Context
	contextP   *ring.Context
	contextQ1P *ring.Context

	alpha uint64
	beta  uint64

	galElRotRow      uint64
	galElRotColLeft  []uint64
	galElRotColRight []uint64
}

func newBFVContext(params *Parameters) (context *bfvContext) {
	context = new(bfvContext)
	var err error

	LogN := params.LogN
	N := uint64(1 << LogN)
	t := params.T

	ModuliQ1, ModuliP, ModuliQ2 := GenModuli(params)
	sigma := params.Sigma

	context.n = N
	context.t = t

	if context.contextT, err = ring.NewContextWithParams(N, []uint64{t}); err != nil {
		panic(err)
	}

	if context.contextQ1, err = ring.NewContextWithParams(N, ModuliQ1); err != nil {
		panic(err)
	}

	if context.contextQ2, err = ring.NewContextWithParams(N, ModuliQ2); err != nil {
		panic(err)
	}

	if context.contextP, err = ring.NewContextWithParams(N, ModuliP); err != nil {
		panic(err)
	}

	if context.contextQ1P, err = ring.NewContextWithParams(N, append(ModuliQ1, ModuliP...)); err != nil {
		panic(err)
	}

	context.alpha = uint64(len(ModuliP))
	context.beta = uint64(math.Ceil(float64(len(ModuliQ1)) / float64(context.alpha)))

	context.logQ = uint64(context.contextQ1P.ModulusBigint.BitLen())

	context.sigma = sigma

	context.gaussianSampler = context.contextQ1P.NewKYSampler(sigma, int(6*sigma))

	context.galElRotColLeft = ring.GenGaloisParams(context.n, GaloisGen)
	context.galElRotColRight = ring.GenGaloisParams(context.n, ring.ModExp(GaloisGen, 2*context.n-1, 2*context.n))
	context.galElRotRow = 2*context.n - 1
	return
}