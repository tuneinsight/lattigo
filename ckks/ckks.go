//Package ckks implements a RNS-accelerated version of the Homomorphic Encryption for Arithmetic for Approximate Numbers
//(HEAAN, a.k.a. CKKS) scheme. It provides approximate arithmetic over the complex numbers.
package ckks

import (
	"github.com/ldsec/lattigo/ring"
	"math"
	"math/big"
)

const GaloisGen uint64 = 5

// Context is a struct which contains all the elements required to instantiate the CKKS Scheme. This includes the parameters (N, ciphertext modulus,
// sampling, polynomial contexts and other parameters required for the homomorphic operations).
type Context struct {

	// Context parameters
	logN     uint64
	logQ     uint64
	scale    float64
	n        uint64
	maxSlots uint64

	// Number of available levels
	levels uint64

	bigintChain []*big.Int

	// Contexts
	alpha     uint64
	beta      uint64
	contextQ  *ring.Context
	contextP  *ring.Context
	contextQP *ring.Context

	// Sampling variance
	sigma float64

	// Samplers
	gaussianSampler *ring.KYSampler

	// Rotation params
	galElConjugate   uint64
	galElRotColLeft  []uint64
	galElRotColRight []uint64
}

// NewContext creates a new Context with the given parameters. Returns an error if one of the parameters would not ensure the
// correctness of the scheme (however it doesn't check for security).
func NewContext(params *Parameters) (ckkscontext *Context) {
	var err error

	ckkscontext = new(Context)

	ckkscontext.logN = uint64(params.LogN)
	ckkscontext.n = 1 << uint64(params.LogN)
	ckkscontext.maxSlots = 1 << (uint64(params.LogN) - 1)
	ckkscontext.scale = params.Scale
	ckkscontext.sigma = params.Sigma

	ckkscontext.levels = uint64(len(params.Q))

	ckkscontext.alpha = uint64(len(params.P))
	ckkscontext.beta = uint64(math.Ceil(float64(ckkscontext.levels) / float64(ckkscontext.alpha)))

	Q, P := GenModuli(params)
	N := ckkscontext.n

	ckkscontext.bigintChain = genBigIntChain(Q)

	if ckkscontext.contextQ, err = ring.NewContextWithParams(N, Q); err != nil {
		panic(err)
	}

	if ckkscontext.contextP, err = ring.NewContextWithParams(N, P); err != nil {
		panic(err)
	}

	if ckkscontext.contextQP, err = ring.NewContextWithParams(N, append(Q, P...)); err != nil {
		panic(err)
	}

	ckkscontext.logQ = uint64(ckkscontext.contextQP.ModulusBigint.BitLen())

	ckkscontext.gaussianSampler = ckkscontext.contextQP.NewKYSampler(params.Sigma, int(6*params.Sigma))

	ckkscontext.galElRotColLeft = ring.GenGaloisParams(N, GaloisGen)
	ckkscontext.galElRotColRight = ring.GenGaloisParams(N, ring.ModExp(GaloisGen, 2*N-1, 2*N))
	ckkscontext.galElConjugate = 2*N - 1

	return ckkscontext

}

// N returns logN of the Context.
func (ckkscontext *Context) N() uint64 {
	return 1 << ckkscontext.logN
}

// LogN returns logN of the Context.
func (ckkscontext *Context) LogN() uint64 {
	return ckkscontext.logN
}

// LogQ returns the log_2(prod(modulie)) of the Context.
func (ckkscontext *Context) LogQ() uint64 {
	return ckkscontext.logQ
}

// Moduli returns the moduli of the Context.
func (ckkscontext *Context) Moduli() []uint64 {
	return ckkscontext.contextQ.Modulus
}

// BigintChain returns the moduli chain in big.Int.
func (ckkscontext *Context) BigintChain() []*big.Int {
	return ckkscontext.bigintChain
}

// KeySwitchPrimes returns the extra moduli used for the KeySwitching operation.
func (ckkscontext *Context) KeySwitchPrimes() []uint64 {
	return ckkscontext.contextP.Modulus
}

// Alpha returns #Pi.
func (ckkscontext *Context) Alpha() uint64 {
	return ckkscontext.alpha
}

// Beta returns ceil(#Qi/#Pi)
func (ckkscontext *Context) Beta() uint64 {
	return ckkscontext.beta
}

// Levels returns the number of levels of the Context.
func (ckkscontext *Context) Levels() uint64 {
	return ckkscontext.levels
}

// Scale returns the default scalt of the Context.
func (ckkscontext *Context) Scale() float64 {
	return ckkscontext.scale
}

// ContextQ returns the rint context of the Ciphertexts.
func (ckkscontext *Context) ContextQ() *ring.Context {
	return ckkscontext.contextQ
}

// ContextP returns the rint context of the KeySwitchPrimes.
func (ckkscontext *Context) ContextP() *ring.Context {
	return ckkscontext.contextP
}

// contextQP returns the ring context of the keys.
func (ckkscontext *Context) ContextQP() *ring.Context {
	return ckkscontext.contextQP
}

// Slots returns the number of slots that the scheme can encrypt at the same time.
func (ckkscontext *Context) Slots() uint64 {
	return ckkscontext.maxSlots
}

// Sigma returns the variance used by the target context to sample gaussian polynomials.
func (ckkscontext *Context) Sigma() float64 {
	return ckkscontext.sigma
}

// GaussianSampler returns the context's gaussian sampler instance
func (ckkscontext *Context) GaussianSampler() *ring.KYSampler {
	return ckkscontext.gaussianSampler
}
