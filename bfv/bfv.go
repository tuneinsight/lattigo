// Package bfv implements a RNS-accelerated Fan-Vercauteren version of Brakerski's scale invariant homomorphic encryption scheme. It provides modular arithmetic over the integers.
package bfv

import (
	"github.com/ldsec/lattigo/ring"
	"math"
)

// BfvContext is a struct which contains all the elements required to instantiate the BFV Scheme. This includes the parameters (N, plaintext modulus, ciphertext modulus,
// sampling, polynomial contexts and other parameters required for the homomorphic operations).
type BfvContext struct {

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

	// Polynomial contexts
	contextT  *ring.Context
	contextQ  *ring.Context
	contextP  *ring.Context
	contextQP *ring.Context

	contextKeys       *ring.Context
	specialprimes     []uint64
	alpha             uint64
	beta              uint64
	rescaleParamsKeys []uint64 // (P^-1) mod each qi

	// Galois elements used to permute the batched plaintext in the encrypted domain
	gen    uint64
	genInv uint64

	galElRotRow      uint64
	galElRotColLeft  []uint64
	galElRotColRight []uint64
}

// NewBfvContext creates a new empty BfvContext.
func NewBfvContext() *BfvContext {
	return new(BfvContext)
}

// NewBfvContextWithParam creates a new BfvContext with the given parameters. Returns an error if one of the parameters would not ensure the
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
func NewBfvContextWithParam(params *Parameters) (newbfvcontext *BfvContext, err error) {
	newbfvcontext = new(BfvContext)
	if err := newbfvcontext.SetParameters(params); err != nil {
		return nil, err
	}
	return
}

// SetParameters populates a new BfvContext with the given parameters. Returns an error if one of the parameters would not ensure the
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
func (bfvContext *BfvContext) SetParameters(params *Parameters) (err error) {

	bfvContext.contextT = ring.NewContext()
	bfvContext.contextQ = ring.NewContext()
	bfvContext.contextP = ring.NewContext()
	bfvContext.contextQP = ring.NewContext()
	bfvContext.contextKeys = ring.NewContext()

	N := params.N
	t := params.T
	ModuliQ := params.Qi
	ModuliP := params.Pi
	sigma := params.Sigma

	bfvContext.n = N
	bfvContext.t = t

	// Plaintext NTT Parameters
	// We do not check for an error since the plaintext NTT is optional
	// it will still compute the other relevant parameters
	bfvContext.contextT.SetParameters(N, []uint64{t})
	bfvContext.contextT.GenNTTParams()
	// ========================

	if err := bfvContext.contextQ.SetParameters(N, ModuliQ); err != nil {
		return err
	}

	if err := bfvContext.contextQ.GenNTTParams(); err != nil {
		return err
	}

	if err := bfvContext.contextP.SetParameters(N, ModuliP); err != nil {
		return err
	}

	if err := bfvContext.contextP.GenNTTParams(); err != nil {
		return err
	}

	if err := bfvContext.contextQP.Merge(bfvContext.contextQ, bfvContext.contextP); err != nil {
		return err
	}

	if err = bfvContext.contextKeys.SetParameters(N, append(ModuliQ, params.KeySwitchPrimes...)); err != nil {
		return err
	}

	if err = bfvContext.contextKeys.GenNTTParams(); err != nil {
		return err
	}

	bfvContext.specialprimes = make([]uint64, len(params.KeySwitchPrimes))
	for i := range params.KeySwitchPrimes {
		bfvContext.specialprimes[i] = params.KeySwitchPrimes[i]
	}

	bfvContext.rescaleParamsKeys = make([]uint64, len(ModuliQ))

	PBig := ring.NewUint(1)
	for _, pj := range bfvContext.specialprimes {
		PBig.Mul(PBig, ring.NewUint(pj))
	}

	bfvContext.alpha = uint64(len(bfvContext.specialprimes))
	bfvContext.beta = uint64(math.Ceil(float64(len(ModuliQ)) / float64(bfvContext.alpha)))

	tmp := ring.NewUint(0)
	bredParams := bfvContext.contextQ.GetBredParams()
	for i, Qi := range ModuliQ {
		tmp.Mod(PBig, ring.NewUint(Qi))
		bfvContext.rescaleParamsKeys[i] = ring.MForm(ring.ModExp(ring.BRedAdd(tmp.Uint64(), Qi, bredParams[i]), Qi-2, Qi), Qi, bredParams[i])
	}

	bfvContext.logQ = uint64(bfvContext.contextKeys.ModulusBigint.Value.BitLen())
	bfvContext.logP = uint64(bfvContext.contextP.ModulusBigint.Value.BitLen())

	delta := ring.NewUint(1).Div(bfvContext.contextQ.ModulusBigint, ring.NewUint(t))
	tmpBig := ring.NewUint(1)
	bfvContext.deltaMont = make([]uint64, len(ModuliQ))
	bfvContext.delta = make([]uint64, len(ModuliQ))
	for i, Qi := range ModuliQ {
		bfvContext.delta[i] = tmpBig.Mod(delta, ring.NewUint(Qi)).Uint64()
		bfvContext.deltaMont[i] = ring.MForm(bfvContext.delta[i], Qi, bfvContext.contextQ.GetBredParams()[i])
	}

	bfvContext.sigma = sigma

	bfvContext.gaussianSampler = bfvContext.contextKeys.NewKYSampler(sigma, int(6*sigma))
	bfvContext.ternarySampler = bfvContext.contextKeys.NewTernarySampler()

	bfvContext.gen = 5
	bfvContext.genInv = ring.ModExp(bfvContext.gen, (N<<1)-1, N<<1)

	mask := (N << 1) - 1

	bfvContext.galElRotColLeft = make([]uint64, N>>1)
	bfvContext.galElRotColRight = make([]uint64, N>>1)

	bfvContext.galElRotColRight[0] = 1
	bfvContext.galElRotColLeft[0] = 1

	for i := uint64(1); i < N>>1; i++ {
		bfvContext.galElRotColLeft[i] = (bfvContext.galElRotColLeft[i-1] * bfvContext.gen) & mask
		bfvContext.galElRotColRight[i] = (bfvContext.galElRotColRight[i-1] * bfvContext.genInv) & mask

	}

	bfvContext.galElRotRow = (N << 1) - 1

	return nil
}

// N returns N which is the degree of the ring, of the target bfvcontext.
func (bfvContext *BfvContext) N() uint64 {
	return bfvContext.n
}

// LogQ returns logQ which is the total bitzise of the ciphertext modulus.
func (bfvContext *BfvContext) LogQ() uint64 {
	return bfvContext.logQ
}

// LogP returns logQ which is the total bitzise of the secondary ciphertext modulus.
func (bfvContext *BfvContext) LogP() uint64 {
	return bfvContext.logP
}

// PlaintextModulus returns the plaintext modulus of the target bfvcontext.
func (bfvContext *BfvContext) T() uint64 {
	return bfvContext.t
}

// Delta returns Q/t, modulo each Qi, where t is the plaintext modulus, Q is the product of all the Qi, of the target bfvcontext.
func (bfvContext *BfvContext) Delta() []uint64 {
	return bfvContext.delta
}

// Delta returns Q/t, modulo each Qi, where t is the plaintext modulus, Q is the product of all the Qi, of the target bfvcontext.
func (bfvContext *BfvContext) DeltaMont() []uint64 {
	return bfvContext.deltaMont
}

// Sigma returns sigma, which is the variance used for the gaussian sampling of the target bfvcontext.
func (bfvContext *BfvContext) Sigma() float64 {
	return bfvContext.sigma
}

// ContextT returns the polynomial (ring) context of the plaintext modulus, of the target bfvcontext.
func (bfvContext *BfvContext) ContextT() *ring.Context {
	return bfvContext.contextT
}

// ContextQ returns the polynomial (ring) context of the ciphertext modulus, of the target bfvcontext.
func (bfvContext *BfvContext) ContextQ() *ring.Context {
	return bfvContext.contextQ
}

// ContextP returns the polynomial (ring) context of the secondary ciphertext modulus, of the target bfvcontext.
func (bfvContext *BfvContext) ContextP() *ring.Context {
	return bfvContext.contextP
}

// ContextQP returns the polynomial (ring) context of the extended ciphertext modulus, of the target bfvcontext.
func (bfvContext *BfvContext) ContextQP() *ring.Context {
	return bfvContext.contextQP
}

// ContextKeys returns the polynomial (ring) context used for the key-generation.
func (bfvContext *BfvContext) ContextKeys() *ring.Context {
	return bfvContext.contextKeys
}

func (bfvcontext *BfvContext) KeySwitchPrimes() []uint64 {
	return bfvcontext.specialprimes
}

func (bfvcontext *BfvContext) Beta() uint64 {
	return bfvcontext.beta
}

// GaussianSampler returns the context's gaussian sampler instance
func (bfvContext *BfvContext) GaussianSampler() *ring.KYSampler {
	return bfvContext.gaussianSampler
}

// TernarySampler returns the context's ternary sampler instance
func (bfvcontext *BfvContext) TernarySampler() *ring.TernarySampler {
	return bfvcontext.ternarySampler
}
