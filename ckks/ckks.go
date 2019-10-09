//Package ckks implements a RNS-accelerated version of the Homomorphic Encryption for Arithmetic for Approximate Numbers
//(HEAAN, a.k.a. CKKS) scheme. It provides approximate arithmetic over the complex numbers.
package ckks

import (
	"errors"
	"github.com/ldsec/lattigo/ring"
	"math/bits"
)

// Context is a struct which contains all the elements required to instantiate the CKKS Scheme. This includes the parameters (N, ciphertext modulus,
// sampling, polynomial contexts and other parameters required for the homomorphic operations).
type Context struct {

	// Context parameters
	logN         uint64
	logQ         uint64
	logPrecision uint64
	logScale     uint64
	n            uint64
	slots        uint64

	// Uperbound in bits of the moduli
	maxBit uint64

	// Number of avaliable levels
	levels uint64

	// moduli chain
	moduli []uint64

	// Contexts chain
	contextLevel []*ring.Context

	// Keys' context
	keyscontext *ring.Context

	// Pre-computed values for the rescaling
	scalechain   []uint64
	rescalParams [][]uint64

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

// NewContext creates a new Context with the given parameters. Returns an error if one of the parameters would not ensure the
// correctness of the scheme (however it doesn't check for security).
func NewContext(params *Parameters) (context *Context, err error) {

	context = new(Context)

	context.logN = uint64(params.LogN)
	context.n = 1 << uint64(params.LogN)
	context.slots = 1 << (uint64(params.LogN) - 1)
	context.logScale = uint64(params.Logscale)
	context.levels = uint64(len(params.Modulichain))

	context.maxBit = 60 // The first prime is always 60 bits to ensure some bits of precision for integers.

	// ========== START < PRIMES GENERATION > START ===============
	context.scalechain = make([]uint64, context.levels)

	primesbitlen := make(map[uint64]uint64)

	for i := range params.Modulichain {

		primesbitlen[uint64(params.Modulichain[i])]++

		context.scalechain[i] = uint64(params.Modulichain[i])

		if uint64(params.Modulichain[i]) > 60 {
			return nil, errors.New("error : provided moduli must be smaller than 60")
		}
	}

	primes := make(map[uint64][]uint64)

	for key, value := range primesbitlen {
		primes[key], _ = GenerateCKKSPrimes(key, uint64(params.LogN), value)
	}

	context.moduli = make([]uint64, context.levels)

	for i := range params.Modulichain {
		context.moduli[i] = primes[uint64(params.Modulichain[i])][0]
		primes[uint64(params.Modulichain[i])] = primes[uint64(params.Modulichain[i])][1:]

		if uint64(bits.Len64(context.moduli[i])) > context.maxBit {
			context.maxBit = uint64(bits.Len64(context.moduli[i]))
		}
	}

	// ========== END < PRIMES GENERATION > END ===============

	// ========== START < CONTEXTS CHAIN > START ===============
	context.contextLevel = make([]*ring.Context, context.levels)

	context.contextLevel[0] = ring.NewContext()

	if err = context.contextLevel[0].SetParameters(1<<uint64(params.LogN), context.moduli[:1]); err != nil {
		return nil, err
	}

	if err = context.contextLevel[0].GenNTTParams(); err != nil {
		return nil, err
	}

	for i := uint64(1); i < context.levels; i++ {

		context.contextLevel[i] = ring.NewContext()

		if err = context.contextLevel[i].SetParameters(1<<uint64(params.LogN), context.moduli[i:i+1]); err != nil {
			return nil, err
		}

		if err = context.contextLevel[i].GenNTTParams(); err != nil {
			return nil, err
		}

		// Instead of recomputing and storing redundant context, subsequent contexts are a merge of previous contexts.
		if err = context.contextLevel[i].Merge(context.contextLevel[i-1], context.contextLevel[i]); err != nil {
			return nil, err
		}
	}
	// ========== END < CONTEXTS CHAIN > END ===============

	// Context used for the generation of the keys
	context.keyscontext = context.contextLevel[context.levels-1]
	context.logQ = uint64(context.contextLevel[context.levels-1].ModulusBigint.Value.BitLen())

	// ========== START < RESCALE PRE-COMPUATION PARAMETERS > START ===============
	var Qi, Ql uint64

	context.rescalParams = make([][]uint64, context.levels-1)

	for j := uint64(context.levels) - 1; j > 0; j-- {

		context.rescalParams[j-1] = make([]uint64, j)

		Ql = context.moduli[j]

		bredParams := context.contextLevel[j-1].GetBredParams()

		for i := uint64(0); i < j; i++ {

			Qi = context.moduli[i]

			context.rescalParams[j-1][i] = ring.MForm(ring.ModExp(Ql, Qi-2, Qi), Qi, bredParams[i])
		}
	}
	// ========== END < RESCALE PRE-COMPUATION PARAMETERS > END ===============

	// default variance
	context.sigma = params.Sigma

	context.gaussianSampler = context.keyscontext.NewKYSampler(params.Sigma, int(6*params.Sigma))
	context.ternarySampler = context.keyscontext.NewTernarySampler()

	// ========== START < ROTATION ELEMENTS > START ===============
	var m, mask uint64

	m = context.n << 1

	mask = m - 1

	context.gen = 5 // Any integer equal to 1 mod 4 and comprime to 2N will do fine
	context.genInv = ring.ModExp(context.gen, mask, m)

	context.galElRotColLeft = make([]uint64, context.slots)
	context.galElRotColRight = make([]uint64, context.slots)

	context.galElRotColRight[0] = 1
	context.galElRotColLeft[0] = 1

	for i := uint64(1); i < context.slots; i++ {
		context.galElRotColLeft[i] = (context.galElRotColLeft[i-1] * context.gen) & mask
		context.galElRotColRight[i] = (context.galElRotColRight[i-1] * context.genInv) & mask
	}

	context.galElRotRow = mask
	// ============ END < ROTATION ELEMENTS > END ================

	return context, nil
}

// LogN returns logN of the context.
func (context *Context) LogN() uint64 {
	return context.logN
}

// LogQ returns the log_2(prod(modulie)) of the context.
func (context *Context) LogQ() uint64 {
	return context.logQ
}

// Moduli returns the moduli of the context.
func (context *Context) Moduli() []uint64 {
	return context.moduli
}

// Levels returns the number of levels of the context.
func (context *Context) Levels() uint64 {
	return context.levels
}

// Scale returns the default scalt of the context.
func (context *Context) Scale() uint64 {
	return context.logScale
}

// ContextKeys returns the ring context under which the keys are created.
func (context *Context) ContextKeys() *ring.Context {
	return context.keyscontext
}

// Slots returns the number of slots that the scheme can encrypt at the same time.
func (context *Context) Slots() uint64 {
	return (1 << (context.logN - 1))
}

// Sigma returns the variance used by the target context to sample gaussian polynomials.
func (context *Context) Sigma() float64 {
	return context.sigma
}
