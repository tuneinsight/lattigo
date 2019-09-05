//Package ckks implements a RNS-accelerated version of the Homomorphic Encryption for Arithmetic for Approximate Numbers
//(HEAAN, a.k.a. CKKS) scheme. It provides approximate arithmetic over the complex numbers.
package ckks

import (
	"errors"
	"github.com/lca1/lattigo/ring"
	"math"
	"math/bits"
)

// CkksContext is a struct which contains all the elements required to instantiate the CKKS Scheme. This includes the parameters (N, ciphertext modulus,
// sampling, polynomial contexts and other parameters required for the homomorphic operations).
type CkksContext struct {

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

	// Encoding and Decoding params
	indexMatrix []uint64
	gap         uint64
	roots       []complex128
	inv_roots   []complex128
}

// NewCkksContext creates a new CkksContext with the given parameters. Returns an error if one of the parameters would not ensure the
// correctness of the scheme (however it doesn't check for security).
func NewCkksContext(params *Parameters) (ckkscontext *CkksContext, err error) {

	ckkscontext = new(CkksContext)

	ckkscontext.logN = params.logN
	ckkscontext.n = 1 << params.logN
	ckkscontext.slots = 1 << (params.logN - 1)
	ckkscontext.logScale = params.logScale
	ckkscontext.levels = uint64(len(params.modulichain))

	ckkscontext.maxBit = 60 // The first prime is always 60 bits to ensure some bits of precision for integers.

	// ========== START < PRIMES GENERATION > START ===============
	ckkscontext.scalechain = make([]uint64, ckkscontext.levels)

	primesbitlen := make(map[uint64]uint64)

	for i := range params.modulichain {

		primesbitlen[params.modulichain[i]] += 1

		ckkscontext.scalechain[i] = params.modulichain[i]

		if params.modulichain[i] > 60 {
			return nil, errors.New("error : provided moduli must be smaller than 60")
		}
	}

	primes := make(map[uint64][]uint64)

	for key, value := range primesbitlen {
		primes[key], _ = GenerateCKKSPrimes(key, params.logN, value)
	}

	ckkscontext.moduli = make([]uint64, ckkscontext.levels)

	for i := range params.modulichain {
		ckkscontext.moduli[i] = primes[params.modulichain[i]][0]
		primes[params.modulichain[i]] = primes[params.modulichain[i]][1:]

		if uint64(bits.Len64(ckkscontext.moduli[i])) > ckkscontext.maxBit {
			ckkscontext.maxBit = uint64(bits.Len64(ckkscontext.moduli[i]))
		}
	}

	// ========== END < PRIMES GENERATION > END ===============

	// ========== START < CONTEXTS CHAIN > START ===============
	ckkscontext.contextLevel = make([]*ring.Context, ckkscontext.levels)

	ckkscontext.contextLevel[0] = ring.NewContext()

	if err = ckkscontext.contextLevel[0].SetParameters(1<<params.logN, ckkscontext.moduli[:1]); err != nil {
		return nil, err
	}

	if err = ckkscontext.contextLevel[0].ValidateParameters(); err != nil {
		return nil, err
	}

	for i := uint64(1); i < ckkscontext.levels; i++ {

		ckkscontext.contextLevel[i] = ring.NewContext()

		if err = ckkscontext.contextLevel[i].SetParameters(1<<params.logN, ckkscontext.moduli[i:i+1]); err != nil {
			return nil, err
		}

		if err = ckkscontext.contextLevel[i].ValidateParameters(); err != nil {
			return nil, err
		}

		// Instead of recomputing and storing redundant context, subsequent contexts are a merge of previous contexts.
		if err = ckkscontext.contextLevel[i].Merge(ckkscontext.contextLevel[i-1], ckkscontext.contextLevel[i]); err != nil {
			return nil, err
		}
	}
	// ========== END < CONTEXTS CHAIN > END ===============

	// Context used for the generation of the keys
	ckkscontext.keyscontext = ckkscontext.contextLevel[ckkscontext.levels-1]
	ckkscontext.logQ = uint64(ckkscontext.contextLevel[ckkscontext.levels-1].ModulusBigint.Value.BitLen())

	// ========== START < RESCALE PRE-COMPUATION PARAMETERS > START ===============
	var Qi, Ql uint64

	ckkscontext.rescalParams = make([][]uint64, ckkscontext.levels-1)

	for j := uint64(ckkscontext.levels) - 1; j > 0; j-- {

		ckkscontext.rescalParams[j-1] = make([]uint64, j)

		Ql = ckkscontext.moduli[j]

		bredParams := ckkscontext.contextLevel[j-1].GetBredParams()

		for i := uint64(0); i < j; i++ {

			Qi = ckkscontext.moduli[i]

			ckkscontext.rescalParams[j-1][i] = ring.MForm(modexp(Ql, Qi-2, Qi), Qi, bredParams[i])
		}
	}
	// ========== END < RESCALE PRE-COMPUATION PARAMETERS > END ===============

	// default variance
	ckkscontext.sigma = params.sigma

	ckkscontext.gaussianSampler = ckkscontext.keyscontext.NewKYSampler(params.sigma, int(6*params.sigma))
	ckkscontext.ternarySampler = ckkscontext.keyscontext.NewTernarySampler()

	// ========== START < ROTATION ELEMENTS > START ===============
	var m, mask uint64

	m = ckkscontext.n << 1

	mask = m - 1

	ckkscontext.gen = 5 // Any integer equal to 1 mod 4 and comprime to 2N will do fine
	ckkscontext.genInv = modexp(ckkscontext.gen, mask, m)

	ckkscontext.galElRotColLeft = make([]uint64, ckkscontext.slots)
	ckkscontext.galElRotColRight = make([]uint64, ckkscontext.slots)

	ckkscontext.galElRotColRight[0] = 1
	ckkscontext.galElRotColLeft[0] = 1

	for i := uint64(1); i < ckkscontext.slots; i++ {
		ckkscontext.galElRotColLeft[i] = (ckkscontext.galElRotColLeft[i-1] * ckkscontext.gen) & mask
		ckkscontext.galElRotColRight[i] = (ckkscontext.galElRotColRight[i-1] * ckkscontext.genInv) & mask
	}

	ckkscontext.galElRotRow = mask
	// ============ END < ROTATION ELEMENTS > END ================

	// ========== START < ENCODER PARAMETERS > START ================
	var pos, index1, index2 uint64

	ckkscontext.gap = 1 //gap-1 is the gap between each slot, here 1 means no gap.

	ckkscontext.indexMatrix = make([]uint64, ckkscontext.n)

	pos = 1

	for i := uint64(0); i < ckkscontext.slots; i++ {

		index1 = (pos - 1) >> 1
		index2 = (m - pos - 1) >> 1

		ckkscontext.indexMatrix[i] = bitReverse64(index1, ckkscontext.logN)
		ckkscontext.indexMatrix[i|ckkscontext.slots] = bitReverse64(index2, ckkscontext.logN)

		pos *= ckkscontext.gen
		pos &= mask
	}

	ckkscontext.roots = make([]complex128, ckkscontext.n)
	ckkscontext.inv_roots = make([]complex128, ckkscontext.n)

	angle := 6.283185307179586 / float64(m)
	psi := complex(math.Cos(angle), math.Sin(angle))
	psiInv := complex(1, 0) / psi

	ckkscontext.roots[0] = 1
	ckkscontext.inv_roots[0] = 1

	for j := uint64(1); j < ckkscontext.n; j++ {

		indexReversePrev := bitReverse64(j-1, ckkscontext.logN)
		indexReverseNext := bitReverse64(j, ckkscontext.logN)

		ckkscontext.roots[indexReverseNext] = ckkscontext.roots[indexReversePrev] * psi
		ckkscontext.inv_roots[indexReverseNext] = ckkscontext.inv_roots[indexReversePrev] * psiInv
	}
	// ========== END < ENCODER PARAMETERS > END ================

	return ckkscontext, nil
}

// LogN returns logN of the ckkscontext.
func (ckkscontext *CkksContext) LogN() uint64 {
	return ckkscontext.logN
}

// LogQ returns the log_2(prod(modulie)) of the ckkscontext.
func (ckkscontext *CkksContext) LogQ() uint64 {
	return ckkscontext.logQ
}

// Moduli returns the moduli of the ckkscontext.
func (ckkscontext *CkksContext) Moduli() []uint64 {
	return ckkscontext.moduli
}

// Levels returns the number of levels of the ckkscontext.
func (ckkscontext *CkksContext) Levels() uint64 {
	return ckkscontext.levels
}

// Scale returns the default scalt of the ckkscontext.
func (ckkscontext *CkksContext) Scale() uint64 {
	return ckkscontext.logScale
}

// ContextKeys returns the ring context under which the keys are created.
func (ckkscontext *CkksContext) ContextKeys() *ring.Context {
	return ckkscontext.keyscontext
}

// Slots returns the number of slots that the scheme can encrypt at the same time.
func (ckkscontext *CkksContext) Slots() uint64 {
	return (1 << (ckkscontext.logN - 1))
}

// Sigma returns the variance used by the target context to sample gaussian polynomials.
func (ckkscontext *CkksContext) Sigma() float64 {
	return ckkscontext.sigma
}
