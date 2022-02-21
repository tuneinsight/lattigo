package dckks

import (
	"math"
	"math/bits"

	"github.com/tuneinsight/lattigo/v3/ckks"
	"github.com/tuneinsight/lattigo/v3/ring"
	"github.com/tuneinsight/lattigo/v3/rlwe"
)

// GetMinimumLevelForBootstrapping takes the security parameter lambda, the ciphertext scale, the number of parties and the moduli chain
// and returns the minimum level at which the collective refresh can be called with a security of at least 128-bits.
// It returns 3 parameters :
// minLevel : the minimum level at which the collective refresh must be called to ensure correctness
// logBound : the bit length of the masks to be sampled to mask the plaintext and ensure 128-bits of statistical indistinguishability
// ok 		: a boolean flag, which is set to false if no such instance exist
func GetMinimumLevelForBootstrapping(lambda int, scale float64, nParties int, moduli []uint64) (minLevel, logBound int, ok bool) {
	logBound = lambda + int(math.Ceil(math.Log2(scale)))
	maxBound := logBound + bits.Len64(uint64(nParties))
	minLevel = -1
	logQ := 0
	for i := 0; logQ < maxBound; i++ {
		if i >= len(moduli) {
			return 0, 0, false
		}

		logQ += bits.Len64(moduli[i])
		minLevel++
	}
	if len(moduli) < minLevel {
		return 0, 0, false
	}

	return minLevel, logBound, true
}

// NewAdditiveShareBigint instantiates a new additive share struct composed of "n" big.Int elements
func NewAdditiveShareBigint(params ckks.Parameters, logSlots int) *rlwe.AdditiveShareBigint {
	dslots := 1 << logSlots
	if params.RingType() == ring.Standard {
		dslots *= 2
	}
	return rlwe.NewAdditiveShareBigint(params.Parameters, dslots)
}
