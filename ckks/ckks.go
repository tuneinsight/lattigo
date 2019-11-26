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

	// Moduli chain
	moduli      []uint64
	scalechain  []float64
	bigintChain []*big.Int

	// Contexts
	specialprimes []uint64
	alpha         uint64
	beta          uint64
	contextQ      *ring.Context
	contextP      *ring.Context
	contextKeys   *ring.Context

	// Pre-computed values for the rescaling
	rescaleParamsKeys []uint64 // (P^-1) mod each qi

	// Sampling variance
	sigma float64

	// Samplers
	gaussianSampler *ring.KYSampler

	// Galoi generator for the rotations, encoding and decoding params
	gen    uint64
	genInv uint64

	// Rotation params
	galElRotRow      uint64
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

	ckkscontext.levels = uint64(len(params.Modulichain))

	ckkscontext.alpha = uint64(len(params.P))
	ckkscontext.beta = uint64(math.Ceil(float64(ckkscontext.levels) / float64(ckkscontext.alpha)))

	// ========== START < PRIMES GENERATION > START ===============
	ckkscontext.scalechain = make([]float64, len(params.Modulichain))

	// Extracts all the different primes bit size and maps their number
	primesbitlen := make(map[uint64]uint64)
	for i, qi := range params.Modulichain {

		primesbitlen[uint64(qi)]++

		if uint64(params.Modulichain[i]) > 60 {
			panic("provided moduli must be smaller than 61")
		}
	}

	for _, pj := range params.P {
		primesbitlen[uint64(pj)]++

		if uint64(pj) > 60 {
			panic("provided P must be smaller than 61")
		}
	}

	// For each bitsize, finds that many primes
	primes := make(map[uint64][]uint64)
	for key, value := range primesbitlen {
		primes[key] = GenerateCKKSPrimes(key, uint64(params.LogN), value)
	}

	// Assigns the primes to the ckks moduli chain
	ckkscontext.moduli = make([]uint64, len(params.Modulichain))
	for i, qi := range params.Modulichain {
		ckkscontext.moduli[i] = primes[uint64(params.Modulichain[i])][0]
		primes[uint64(qi)] = primes[uint64(qi)][1:]

		ckkscontext.scalechain[i] = float64(ckkscontext.moduli[i])
	}

	// Assigns the primes to the special primes list for the the keyscontext
	ckkscontext.specialprimes = make([]uint64, len(params.P))
	for i, pj := range params.P {
		ckkscontext.specialprimes[i] = primes[uint64(pj)][0]
		primes[uint64(pj)] = primes[uint64(pj)][1:]
	}

	ckkscontext.bigintChain = make([]*big.Int, ckkscontext.levels)

	ckkscontext.bigintChain[0] = ring.NewUint(ckkscontext.moduli[0])
	for i := uint64(1); i < ckkscontext.levels; i++ {
		ckkscontext.bigintChain[i] = ring.NewUint(ckkscontext.moduli[i])
		ckkscontext.bigintChain[i].Mul(ckkscontext.bigintChain[i], ckkscontext.bigintChain[i-1])
	}

	//Contexts
	ckkscontext.contextQ = ring.NewContext()
	ckkscontext.contextQ.SetParameters(1<<ckkscontext.logN, ckkscontext.moduli)

	if err = ckkscontext.contextQ.GenNTTParams(); err != nil {
		panic(err)
	}

	ckkscontext.contextP = ring.NewContext()
	ckkscontext.contextP.SetParameters(1<<ckkscontext.logN, ckkscontext.specialprimes)

	if err = ckkscontext.contextP.GenNTTParams(); err != nil {
		panic(err)
	}

	ckkscontext.contextKeys = ring.NewContext()
	ckkscontext.contextKeys.SetParameters(1<<ckkscontext.logN, append(ckkscontext.moduli, ckkscontext.specialprimes...))

	if err = ckkscontext.contextKeys.GenNTTParams(); err != nil {
		panic(err)
	}

	ckkscontext.logQ = uint64(ckkscontext.contextKeys.ModulusBigint.BitLen())

	var Qi uint64

	bredParams := ckkscontext.contextQ.GetBredParams()

	ckkscontext.rescaleParamsKeys = make([]uint64, ckkscontext.levels)

	PBig := ring.NewUint(1)
	for _, pj := range ckkscontext.specialprimes {
		PBig.Mul(PBig, ring.NewUint(pj))
	}

	tmp := ring.NewUint(0)

	for i := uint64(0); i < ckkscontext.levels; i++ {

		Qi = ckkscontext.moduli[i]

		tmp.Mod(PBig, ring.NewUint(Qi))

		ckkscontext.rescaleParamsKeys[i] = ring.MForm(ring.ModExp(ring.BRedAdd(tmp.Uint64(), Qi, bredParams[i]), Qi-2, Qi), Qi, bredParams[i])
	}

	ckkscontext.gaussianSampler = ckkscontext.contextKeys.NewKYSampler(params.Sigma, int(6*params.Sigma))

	var m, mask uint64

	m = ckkscontext.n << 1

	mask = m - 1

	ckkscontext.gen = 5 // Any integer equal to 1 mod 4 and comprime to 2N will do fine
	ckkscontext.genInv = ring.ModExp(ckkscontext.gen, mask, m)

	ckkscontext.galElRotColLeft = make([]uint64, ckkscontext.maxSlots)
	ckkscontext.galElRotColRight = make([]uint64, ckkscontext.maxSlots)

	ckkscontext.galElRotColRight[0] = 1
	ckkscontext.galElRotColLeft[0] = 1

	for i := uint64(1); i < ckkscontext.maxSlots; i++ {
		ckkscontext.galElRotColLeft[i] = (ckkscontext.galElRotColLeft[i-1] * ckkscontext.gen) & mask
		ckkscontext.galElRotColRight[i] = (ckkscontext.galElRotColRight[i-1] * ckkscontext.genInv) & mask
	}

	ckkscontext.galElRotRow = mask

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
	return ckkscontext.moduli
}

// BigintChain returns the moduli chain in big.Int.
func (ckkscontext *Context) BigintChain() []*big.Int {
	return ckkscontext.bigintChain
}

// KeySwitchPrimes returns the extra moduli used for the KeySwitching operation.
func (ckkscontext *Context) KeySwitchPrimes() []uint64 {
	return ckkscontext.specialprimes
}

// Alpha returns #Pi.
func (ckkscontext *Context) Alpha() uint64 {
	return ckkscontext.alpha
}

// Beta returns ceil(#Qi/#Pi)
func (ckkscontext *Context) Beta() uint64 {
	return ckkscontext.beta
}

// RescaleParamsKeys returns the rescaling parameters for the KeySwitching operation.
func (ckkscontext *Context) RescaleParamsKeys() []uint64 {
	return ckkscontext.rescaleParamsKeys
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

// ContextKeys returns the ring context of the keys.
func (ckkscontext *Context) ContextKeys() *ring.Context {
	return ckkscontext.contextKeys
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
