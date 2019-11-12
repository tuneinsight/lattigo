//Package ckks implements a RNS-accelerated version of the Homomorphic Encryption for Arithmetic for Approximate Numbers
//(HEAAN, a.k.a. CKKS) scheme. It provides approximate arithmetic over the complex numbers.
package ckks

import (
	"errors"
	"github.com/ldsec/lattigo/ring"
	"math"
)

const GaloisGen uint64 = 5

// CkksContext is a struct which contains all the elements required to instantiate the CKKS Scheme. This includes the parameters (N, ciphertext modulus,
// sampling, polynomial contexts and other parameters required for the homomorphic operations).
type CkksContext struct {

	// Context parameters
	logN     uint64
	logQ     uint64
	scale    float64
	n        uint64
	maxSlots uint64

	// Number of avaliable levels
	levels uint64

	// Moduli chain
	moduli      []uint64
	scalechain  []float64
	bigintChain []*ring.Int

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
	ternarySampler  *ring.TernarySampler

	// Galoi generator for the rotations, encoding and decoding params
	gen    uint64
	genInv uint64

	// Rotation params
	galElRotRow      uint64
	galElRotColLeft  []uint64
	galElRotColRight []uint64
}

// NewCkksContext creates a new CkksContext with the given parameters. Returns an error if one of the parameters would not ensure the
// correctness of the scheme (however it doesn't check for security).
func NewCkksContext(params *Parameters) (ckkscontext *CkksContext, err error) {
	ckkscontext = new(CkksContext)

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

		primesbitlen[uint64(qi)] += 1

		if uint64(params.Modulichain[i]) > 60 {
			return nil, errors.New("error : provided moduli must be smaller than 61")
		}
	}

	for _, pj := range params.P {
		primesbitlen[uint64(pj)] += 1

		if uint64(pj) > 60 {
			return nil, errors.New("error : provided P must be smaller than 61")
		}
	}

	// For each bitsize, finds that many primes
	primes := make(map[uint64][]uint64)
	for key, value := range primesbitlen {
		primes[key], _ = GenerateCKKSPrimes(key, uint64(params.LogN), value)
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

	ckkscontext.bigintChain = make([]*ring.Int, ckkscontext.levels)

	ckkscontext.bigintChain[0] = ring.NewUint(ckkscontext.moduli[0])
	for i := uint64(1); i < ckkscontext.levels; i++ {
		ckkscontext.bigintChain[i] = ring.NewUint(ckkscontext.moduli[i])
		ckkscontext.bigintChain[i].Mul(ckkscontext.bigintChain[i], ckkscontext.bigintChain[i-1])
	}

	//Contexts
	ckkscontext.contextQ = ring.NewContext()
	if err = ckkscontext.contextQ.SetParameters(1<<ckkscontext.logN, ckkscontext.moduli); err != nil {
		return nil, err
	}
	if err = ckkscontext.contextQ.GenNTTParams(); err != nil {
		return nil, err
	}

	ckkscontext.contextP = ring.NewContext()
	if err = ckkscontext.contextP.SetParameters(1<<ckkscontext.logN, ckkscontext.specialprimes); err != nil {
		return nil, err
	}

	if err = ckkscontext.contextP.GenNTTParams(); err != nil {
		return nil, err
	}

	ckkscontext.contextKeys = ring.NewContext()
	if err = ckkscontext.contextKeys.SetParameters(1<<ckkscontext.logN, append(ckkscontext.moduli, ckkscontext.specialprimes...)); err != nil {
		return nil, err
	}

	if err = ckkscontext.contextKeys.GenNTTParams(); err != nil {
		return nil, err
	}

	ckkscontext.logQ = uint64(ckkscontext.contextKeys.ModulusBigint.Value.BitLen())

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
	ckkscontext.ternarySampler = ckkscontext.contextKeys.NewTernarySampler()

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

	return ckkscontext, nil

}

// LogN returns logN of the ckksContext.
func (ckksContext *CkksContext) N() uint64 {
	return 1 << ckksContext.logN
}

// LogN returns logN of the ckksContext.
func (ckksContext *CkksContext) LogN() uint64 {
	return ckksContext.logN
}

// LogQ returns the log_2(prod(modulie)) of the ckksContext.
func (ckksContext *CkksContext) LogQ() uint64 {
	return ckksContext.logQ
}

// Moduli returns the moduli of the ckksContext.
func (ckksContext *CkksContext) Moduli() []uint64 {
	return ckksContext.moduli
}

func (ckksContext *CkksContext) BigintChain() []*ring.Int {
	return ckksContext.bigintChain
}

func (ckksContext *CkksContext) KeySwitchPrimes() []uint64 {
	return ckksContext.specialprimes
}

func (ckksContext *CkksContext) Alpha() uint64 {
	return ckksContext.alpha
}

func (ckksContext *CkksContext) Beta() uint64 {
	return ckksContext.beta
}

func (ckksContext *CkksContext) RescaleParamsKeys() []uint64 {
	return ckksContext.rescaleParamsKeys
}

// Levels returns the number of levels of the ckksContext.
func (ckksContext *CkksContext) Levels() uint64 {
	return ckksContext.levels
}

// Scale returns the default scalt of the ckksContext.
func (ckksContext *CkksContext) Scale() float64 {
	return ckksContext.scale
}

func (ckksContext *CkksContext) ContextQ() *ring.Context {
	return ckksContext.contextQ
}

func (ckksContext *CkksContext) ContextP() *ring.Context {
	return ckksContext.contextP
}

// ContextKeys returns the ring context under which the keys are created.
func (ckksContext *CkksContext) ContextKeys() *ring.Context {
	return ckksContext.contextKeys
}

// Slots returns the number of slots that the scheme can encrypt at the same time.
func (ckksContext *CkksContext) Slots() uint64 {
	return ckksContext.maxSlots
}

// Sigma returns the variance used by the target context to sample gaussian polynomials.
func (ckksContext *CkksContext) Sigma() float64 {
	return ckksContext.sigma
}

// GaussianSampler returns the context's gaussian sampler instance
func (ckksContext *CkksContext) GaussianSampler() *ring.KYSampler {
	return ckksContext.gaussianSampler
}

// TernarySampler returns the context's ternary sampler instance
func (ckksContext *CkksContext) TernarySampler() *ring.TernarySampler {
	return ckksContext.ternarySampler
}
