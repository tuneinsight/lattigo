//Package ckks implements a RNS-accelerated version of the Homomorphic Encryption for Arithmetic for Approximate Numbers
//(HEAAN, a.k.a. CKKS) scheme. It provides approximate arithmetic over the complex numbers.
package ckks

import (
	"errors"
	"github.com/ldsec/lattigo/ring"
	"math"
)

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

	// Modulie chain
	moduli     []uint64
	scalechain []float64

	// Contexts chain
	contextLevel []*ring.Context

	// Keys' contexts
	specialprimes []uint64
	alpha         uint64
	beta          uint64
	keyscontext   []*ring.Context
	contexts      []*ring.Context

	// Pre-computed values for the rescaling
	rescaleParams     [][]uint64
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

	specialprimecontext := ring.NewContext()

	if err = specialprimecontext.SetParameters(1<<ckkscontext.logN, ckkscontext.specialprimes); err != nil {
		return nil, err
	}

	if err = specialprimecontext.GenNTTParams(); err != nil {
		return nil, err
	}

	ckkscontext.contextLevel = make([]*ring.Context, ckkscontext.levels)
	ckkscontext.keyscontext = make([]*ring.Context, ckkscontext.levels)

	ckkscontext.contextLevel[0] = ring.NewContext()
	ckkscontext.keyscontext[0] = ring.NewContext()

	if err = ckkscontext.contextLevel[0].SetParameters(1<<ckkscontext.logN, ckkscontext.moduli[:1]); err != nil {
		return nil, err
	}

	if err = ckkscontext.contextLevel[0].GenNTTParams(); err != nil {
		return nil, err
	}

	// Each key context is equal to the same context + a special prime
	if err = ckkscontext.keyscontext[0].Merge(ckkscontext.contextLevel[0], specialprimecontext); err != nil {
		return nil, err
	}

	for i := uint64(1); i < ckkscontext.levels; i++ {

		ckkscontext.contextLevel[i] = ring.NewContext()

		if err = ckkscontext.contextLevel[i].SetParameters(1<<ckkscontext.logN, ckkscontext.moduli[i:i+1]); err != nil {
			return nil, err
		}

		if err = ckkscontext.contextLevel[i].GenNTTParams(); err != nil {
			return nil, err
		}

		// Instead of recomputing and storing redundant context, subsequent contexts are a merge of previous contexts.
		if err = ckkscontext.contextLevel[i].Merge(ckkscontext.contextLevel[i-1], ckkscontext.contextLevel[i]); err != nil {
			return nil, err
		}

		// Special context used for the generation of the keys
		// Each key context is equal to the same context + a special prime
		ckkscontext.keyscontext[i] = ckkscontext.contextLevel[i].CopyNew()
		if err = ckkscontext.keyscontext[i].Merge(ckkscontext.keyscontext[i], specialprimecontext); err != nil {
			return nil, err
		}
	}
	// ========== END < CONTEXTS CHAIN > END ===============

	ckkscontext.logQ = uint64(ckkscontext.keyscontext[ckkscontext.levels-1].ModulusBigint.Value.BitLen())

	var Qi uint64

	bredParams := ckkscontext.contextLevel[ckkscontext.levels-1].GetBredParams()

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

	var Ql uint64

	ckkscontext.rescaleParams = make([][]uint64, ckkscontext.levels-1)

	for j := uint64(ckkscontext.levels) - 1; j > 0; j-- {

		ckkscontext.rescaleParams[j-1] = make([]uint64, j)

		Ql = ckkscontext.moduli[j]

		bredParams := ckkscontext.contextLevel[j-1].GetBredParams()

		for i := uint64(0); i < j; i++ {

			Qi = ckkscontext.moduli[i]

			ckkscontext.rescaleParams[j-1][i] = ring.MForm(ring.ModExp(Ql, Qi-2, Qi), Qi, bredParams[i])
		}
	}

	ckkscontext.gaussianSampler = ckkscontext.keyscontext[ckkscontext.levels-1].NewKYSampler(params.Sigma, int(6*params.Sigma))
	ckkscontext.ternarySampler = ckkscontext.keyscontext[ckkscontext.levels-1].NewTernarySampler()

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

func (ckksContext *CkksContext) KeySwitchPrimes() []uint64 {
	return ckksContext.specialprimes
}

// Levels returns the number of levels of the ckksContext.
func (ckksContext *CkksContext) Levels() uint64 {
	return ckksContext.levels
}

// Scale returns the default scalt of the ckksContext.
func (ckksContext *CkksContext) Scale() float64 {
	return ckksContext.scale
}

func (ckksContext *CkksContext) Context(level uint64) *ring.Context {
	return ckksContext.contextLevel[level]
}

// ContextKeys returns the ring context under which the keys are created.
func (ckksContext *CkksContext) ContextKey(level uint64) *ring.Context {
	return ckksContext.keyscontext[level]
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