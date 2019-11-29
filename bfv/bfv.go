// Package bfv implements a RNS-accelerated Fan-Vercauteren version of Brakerski's scale invariant homomorphic encryption scheme. It provides modular arithmetic over the integers.
package bfv

import (
	"github.com/ldsec/lattigo/ring"
	"math"
)

const GaloisGen uint64 = 5

// Context is a struct which contains all the elements required to instantiate the BFV Scheme. This includes the parameters (N, plaintext modulus, ciphertext modulus,
// sampling, polynomial contexts and other parameters required for the homomorphic operations).
type Context struct {
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

// NewContextWithParam creates a new Context with the given parameters. Returns an error if one of the parameters would not ensure the
// correctness of the scheme (however it doesn't check for security).
//
// Parameters :
//
// - N        : the ring degree (must be a power of 2).
//
// - t        : the plaintext modulus (must be a prime congruent to 1 mod 2N to enable batching).
//
// - ModuliQ1 : the ciphertext modulus composed primes congruent to 1 mod 2N.
//
// - ModuliP : the secondary ciphertext modulus used during the multiplication, composed of primes congruent to 1 mod 2N. Must be bigger than ModuliQ1 by a margin of ~20 bits.
//
// - sigma    : the variance of the gaussian sampling.
func NewContext(params *Parameters) (newContext *Context) {
	newContext = new(Context)
	newContext.SetParameters(params)
	return
}

// SetParameters populates a new Context with the given parameters. Returns an error if one of the parameters would not ensure the
// correctness of the scheme (however it doesn't check for security).
//
// Parameters :
//
// - N        : the ring degree (must be a power of 2).
//
// - t        : the plaintext modulus (must be a prime congruent to 1 mod 2N to enable batching).
//
// - ModuliQ1 : the ciphertext modulus composed primes congruent to 1 mod 2N.
//
// - ModuliP : the secondary ciphertext modulus used during the multiplication, composed of primes congruent to 1 mod 2N. Must be bigger than ModuliQ1 by a margin of ~20 bits.
//
// - sigma    : the variance of the gaussian sampling.
func (context *Context) SetParameters(params *Parameters) {

	var err error

	LogN := params.LogN
	N := uint64(1 << LogN)
	t := params.T

	ModuliQ1, ModuliP, ModuliQ2 := GenModuli(params)
	sigma := params.Sigma

	context.n = N
	context.t = t

	var err error
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
}

// N returns N which is the degree of the ring, of the target context.
func (context *Context) N() uint64 {
	return context.n
}

// LogQ returns logQ which is the total bitzise of the ciphertext modulus.
func (context *Context) LogQ() uint64 {
	return context.logQ
}

// T returns the plaintext modulus of the target context.
func (context *Context) T() uint64 {
	return context.t
}

// Sigma returns sigma, which is the variance used for the gaussian sampling of the target context.
func (context *Context) Sigma() float64 {
	return context.sigma
}

// ContextT returns the polynomial (ring) context of the plaintext modulus, of the target context.
func (context *Context) ContextT() *ring.Context {
	return context.contextT
}

// ContextQ returns the polynomial (ring) context of the ciphertext modulus, of the target context.
func (context *Context) ContextQ1() *ring.Context {
	return context.contextQ1
}

// ContextP returns the polynomial (ring) context of the secondary ciphertext modulus, of the target context.
func (context *Context) ContextQ2() *ring.Context {
	return context.contextQ2
}

// ContextKeys returns the polynomial (ring) context used for the key-generation.
func (context *Context) ContextKeys() *ring.Context {
	return context.contextQ1P
}

// ContextPKeys returns the ring Context of the KeySwitchPrimes.
func (context *Context) ContextP() *ring.Context {
	return context.contextP
}

// Alpha returns #Pi.
func (context *Context) Alpha() uint64 {
	return context.alpha
}

// Beta returns ceil(#Qi/#Pi).
func (context *Context) Beta() uint64 {
	return context.beta
}

// GaussianSampler returns the context's gaussian sampler instance
func (context *Context) GaussianSampler() *ring.KYSampler {
	return context.gaussianSampler
}
