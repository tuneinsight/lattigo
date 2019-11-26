// Package bfv implements a RNS-accelerated Fan-Vercauteren version of Brakerski's scale invariant homomorphic encryption scheme. It provides modular arithmetic over the integers.
package bfv

import (
	"github.com/ldsec/lattigo/ring"
	"math"
	"math/big"
)

const GaloisGen uint64 = 5

// Context is a struct which contains all the elements required to instantiate the BFV Scheme. This includes the parameters (N, plaintext modulus, ciphertext modulus,
// sampling, polynomial contexts and other parameters required for the homomorphic operations).
type Context struct {

	// Polynomial degree
	n uint64

	// Plaintext Modulus
	t uint64

	logQ uint64
	logP uint64

	// floor(Q/T) mod each Qi in Montgomery form
	deltaMont []uint64
	delta     []uint64

	// Ternary and Gaussian samplers
	sigma           float64
	gaussianSampler *ring.KYSampler

	// Polynomial contexts
	contextT *ring.Context
	contextQ *ring.Context
	contextP *ring.Context

	QHalf *big.Int
	PHalf *big.Int

	rescaleParamsMul []uint64

	contextKeys       *ring.Context
	contextPKeys      *ring.Context
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

// NewContext creates a new empty context.
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
func NewContextWithParam(params *Parameters) (newContext *Context) {
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
// - ModuliQ : the ciphertext modulus composed primes congruent to 1 mod 2N.
//
// - ModuliP : the secondary ciphertext modulus used during the multiplication, composed of primes congruent to 1 mod 2N. Must be bigger than ModuliQ by a margin of ~20 bits.
//
// - sigma    : the variance of the gaussian sampling.
func (context *Context) SetParameters(params *Parameters) {

	context.contextT = ring.NewContext()
	context.contextQ = ring.NewContext()
	context.contextP = ring.NewContext()
	context.contextKeys = ring.NewContext()
	context.contextPKeys = ring.NewContext()

	N := params.N
	t := params.T
	ModuliQ := params.Qi
	ModuliP := params.Pi
	sigma := params.Sigma

	context.n = N
	context.t = t

	// Plaintext NTT Parameters
	// We do not check for an error since the plaintext NTT is optional
	// it will still compute the other relevant parameters
	context.contextT.SetParameters(N, []uint64{t})
	if err := context.contextT.GenNTTParams(); err != nil {
		panic(err)
	}
	// ========================

	context.contextQ.SetParameters(N, ModuliQ)

	if err := context.contextQ.GenNTTParams(); err != nil {
		panic(err)
	}

	context.contextP.SetParameters(N, ModuliP)

	if err := context.contextP.GenNTTParams(); err != nil {
		panic(err)
	}

	context.contextKeys.SetParameters(N, append(ModuliQ, params.KeySwitchPrimes...))

	if err := context.contextKeys.GenNTTParams(); err != nil {
		panic(err)
	}

	context.contextPKeys.SetParameters(N, params.KeySwitchPrimes)

	if err := context.contextPKeys.GenNTTParams(); err != nil {
		panic(err)
	}

	context.specialprimes = make([]uint64, len(params.KeySwitchPrimes))
	for i := range params.KeySwitchPrimes {
		context.specialprimes[i] = params.KeySwitchPrimes[i]
	}

	context.rescaleParamsKeys = make([]uint64, len(ModuliQ))

	PBig := ring.NewUint(1)
	for _, pj := range context.specialprimes {
		PBig.Mul(PBig, ring.NewUint(pj))
	}

	context.alpha = uint64(len(context.specialprimes))
	context.beta = uint64(math.Ceil(float64(len(ModuliQ)) / float64(context.alpha)))

	tmp := new(big.Int)
	bredParams := context.contextQ.GetBredParams()
	for i, Qi := range ModuliQ {
		tmp.Mod(PBig, ring.NewUint(Qi))
		context.rescaleParamsKeys[i] = ring.MForm(ring.ModExp(ring.BRedAdd(tmp.Uint64(), Qi, bredParams[i]), Qi-2, Qi), Qi, bredParams[i])
	}

	context.rescaleParamsMul = make([]uint64, len(context.contextP.Modulus))

	bredParams = context.contextP.GetBredParams()
	for i, Pi := range context.contextP.Modulus {
		tmp.Mod(context.contextQ.ModulusBigint, ring.NewUint(Pi))
		context.rescaleParamsMul[i] = ring.MForm(ring.ModExp(ring.BRedAdd(tmp.Uint64(), Pi, bredParams[i]), Pi-2, Pi), Pi, bredParams[i])
	}

	context.QHalf = new(big.Int).Rsh(context.contextQ.ModulusBigint, 1)
	context.PHalf = new(big.Int).Rsh(context.contextP.ModulusBigint, 1)

	context.logQ = uint64(context.contextKeys.ModulusBigint.BitLen())
	context.logP = uint64(context.contextP.ModulusBigint.BitLen())

	delta := new(big.Int).Quo(context.contextQ.ModulusBigint, ring.NewUint(t))
	tmpBig := new(big.Int)
	context.deltaMont = make([]uint64, len(ModuliQ))
	context.delta = make([]uint64, len(ModuliQ))
	for i, Qi := range ModuliQ {
		context.delta[i] = tmpBig.Mod(delta, ring.NewUint(Qi)).Uint64()
		context.deltaMont[i] = ring.MForm(context.delta[i], Qi, context.contextQ.GetBredParams()[i])
	}

	context.sigma = sigma

	context.gaussianSampler = context.contextKeys.NewKYSampler(sigma, int(6*sigma))

	context.gen = GaloisGen
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
}

// N returns N which is the degree of the ring, of the target context.
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

// T returns the plaintext modulus of the target context.
func (context *Context) T() uint64 {
	return context.t
}

// Delta returns Q/t, modulo each Qi, where t is the plaintext modulus, Q is the product of all the Qi, of the target context.
func (context *Context) Delta() []uint64 {
	return context.delta
}

// DeltaMont returns Q/t, modulo each Qi, where t is the plaintext modulus, Q is the product of all the Qi, of the target context.
func (context *Context) DeltaMont() []uint64 {
	return context.deltaMont
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
func (context *Context) ContextQ() *ring.Context {
	return context.contextQ
}

// ContextP returns the polynomial (ring) context of the secondary ciphertext modulus, of the target context.
func (context *Context) ContextP() *ring.Context {
	return context.contextP
}

// ContextKeys returns the polynomial (ring) context used for the key-generation.
func (context *Context) ContextKeys() *ring.Context {
	return context.contextKeys
}

// ContextPKeys returns the ring Context of the KeySwitchPrimes.
func (context *Context) ContextPKeys() *ring.Context {
	return context.contextPKeys
}

// KeySwitchPrimes returns the extended P moduli used for the KeySwitching operation.
func (context *Context) KeySwitchPrimes() []uint64 {
	return context.specialprimes
}

// Alpha returns #Pi.
func (context *Context) Alpha() uint64 {
	return context.alpha
}

// Beta returns ceil(#Qi/#Pi).
func (context *Context) Beta() uint64 {
	return context.beta
}

// RescaleParamsKeys returns the rescaling parameters for the P moduli.
func (context *Context) RescaleParamsKeys() []uint64 {
	return context.rescaleParamsKeys
}

// GaussianSampler returns the context's gaussian sampler instance
func (context *Context) GaussianSampler() *ring.KYSampler {
	return context.gaussianSampler
}
