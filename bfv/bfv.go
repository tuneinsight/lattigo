// Package bfv implements a RNS-accelerated Fan-Vercauteren version of Brakerski's scale invariant homomorphic encryption scheme. It provides modular arithmetic over the integers.
package bfv

import (
	"github.com/ldsec/lattigo/ring"
	"math/bits"
)

// Context is a struct which contains all the elements required to instantiate the BFV Scheme. This includes the parameters (N, plaintext modulus, ciphertext modulus,
// sampling, polynomial contexts and other parameters required for the homomorphic operations).
type Context struct {

	// Polynomial degree
	n uint64

	// Plaintext Modulus
	t uint64

	logQ uint64
	logP uint64

	// floor(Q/T) mod each Qi in montgomeryform
	deltaMont []uint64
	delta     []uint64

	// Ternary and Gaussian samplers
	sigma           float64
	gaussianSampler *ring.KYSampler
	ternarySampler  *ring.TernarySampler

	// Maximum bitLength of the modulies
	maxBit uint64

	// Polynomial contexts
	contextT  *ring.Context
	contextQ  *ring.Context
	contextP  *ring.Context
	contextQP *ring.Context

	// Galois elements used to permute the batched plaintext in the encrypted domain
	gen    uint64
	genInv uint64

	galElRotRow      uint64
	galElRotColLeft  []uint64
	galElRotColRight []uint64
}

// NewContext creates a new empty Context.
func NewContext() *Context {
	return new(Context)
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
// - ModuliQ : the ciphertext modulus composed primes congruent to 1 mod 2N.
//
// - ModuliP : the secondary ciphertext modulus used during the multiplication, composed of primes congruent to 1 mod 2N. Must be bigger than ModuliQ by a margin of ~20 bits.
//
// - sigma    : the variance of the gaussian sampling.
func NewContextWithParam(params *Parameters) (newContext *Context, err error) {
	newContext = new(Context)
	if err := newContext.SetParameters(params); err != nil {
		return nil, err
	}
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
// - ModuliQ : the ciphertext modulus composed primes congruent to 1 mod 2N.
//
// - ModuliP : the secondary ciphertext modulus used during the multiplication, composed of primes congruent to 1 mod 2N. Must be bigger than ModuliQ by a margin of ~20 bits.
//
// - sigma    : the variance of the gaussian sampling.
func (context *Context) SetParameters(params *Parameters) (err error) {

	contextT := ring.NewContext()
	contextQ := ring.NewContext()
	contextP := ring.NewContext()
	contextQP := ring.NewContext()

	N := params.N
	t := params.T
	ModuliQ := params.Qi
	ModuliP := params.Pi
	sigma := params.Sigma

	// Plaintext NTT Parameters
	// We do not check for an error since the plaintext NTT is optional
	// it will still compute the other relevant parameters
	contextT.SetParameters(N, []uint64{t})
	contextT.GenNTTParams()
	// ========================

	if err := contextQ.SetParameters(N, ModuliQ); err != nil {
		return err
	}

	if err := contextQ.GenNTTParams(); err != nil {
		return err
	}

	if err := contextP.SetParameters(N, ModuliP); err != nil {
		return err
	}

	if err := contextP.GenNTTParams(); err != nil {
		return err
	}

	if err := contextQP.Merge(contextQ, contextP); err != nil {
		return err
	}

	context.n = N

	context.t = t

	context.logQ = uint64(contextQ.ModulusBigint.Value.BitLen())
	context.logP = uint64(contextP.ModulusBigint.Value.BitLen())

	delta := ring.NewUint(1).Div(contextQ.ModulusBigint, ring.NewUint(t))
	tmpBig := ring.NewUint(1)
	context.deltaMont = make([]uint64, len(ModuliQ))
	context.delta = make([]uint64, len(ModuliQ))
	for i, Qi := range ModuliQ {
		context.delta[i] = tmpBig.Mod(delta, ring.NewUint(Qi)).Uint64()
		context.deltaMont[i] = ring.MForm(context.delta[i], Qi, contextQ.GetBredParams()[i])
	}

	context.sigma = sigma

	context.gaussianSampler = contextQ.NewKYSampler(sigma, int(6*sigma))
	context.ternarySampler = contextQ.NewTernarySampler()

	context.maxBit = 0

	for i := range ModuliQ {
		tmp := uint64(bits.Len64(ModuliQ[i]))
		if tmp > context.maxBit {
			context.maxBit = tmp
		}
	}

	context.contextT = contextT
	context.contextQ = contextQ
	context.contextP = contextP
	context.contextQP = contextQP

	context.gen = 5
	context.genInv = ring.ModExp(context.gen, (N<<1)-1, N<<1)

	mask := (N << 1) - 1

	context.galElRotColLeft = make([]uint64, N>>1)
	context.galElRotColRight = make([]uint64, N>>1)

	context.galElRotColRight[0] = 1
	context.galElRotColLeft[0] = 1

	for i := uint64(1); i < N>>1; i++ {
		context.galElRotColLeft[i] = (context.galElRotColLeft[i-1] * context.gen) & mask
		context.galElRotColRight[i] = (context.galElRotColRight[i-1] * context.genInv) & mask

	}

	context.galElRotRow = (N << 1) - 1

	return nil
}

// N returns N which is the degree of the ring, of the target Context.
func (context *Context) N() uint64 {
	return context.n
}

// LogQ returns logQ which is the total bitzise of the ciphertext modulus.
func (context *Context) LogQ() uint64 {
	return context.logQ
}

// LogP returns logQ which is the total bitzise of the secondary ciphertext modulus.
func (context *Context) LogP() uint64 {
	return context.logP
}

// T returns the plaintext modulus of the target Context.
func (context *Context) T() uint64 {
	return context.t
}

// Delta returns t/Q, modulo each Qi, where t is the plaintext modulus, Q is the product of all the Qi, of the target Context.
func (context *Context) Delta() []uint64 {
	return context.delta
}

// Sigma returns sigma, which is the variance used for the gaussian sampling of the target Context.
func (context *Context) Sigma() float64 {
	return context.sigma
}

// ContextT returns the polynomial (ring) context of the plaintext modulus, of the target Context.
func (context *Context) ContextT() *ring.Context {
	return context.contextT
}

// ContextQ returns the polynomial (ring) context of the ciphertext modulus, of the target Context.
func (context *Context) ContextQ() *ring.Context {
	return context.contextQ
}

// ContextP returns the polynomial (ring) context of the secondary ciphertext modulus, of the target Context.
func (context *Context) ContextP() *ring.Context {
	return context.contextP
}

// ContextQP returns the polynomial (ring) context of the extended ciphertext modulus, of the target Context.
func (context *Context) ContextQP() *ring.Context {
	return context.contextQP
}
