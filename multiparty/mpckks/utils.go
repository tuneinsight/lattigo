package mpckks

import (
	"math"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
)

// GetMinimumLevelForRefresh takes the security parameter lambda, the ciphertext scale, the number of parties and the moduli chain
// and returns the minimum level at which the collective refresh can be called with a security of at least 128-bits.
// It returns 3 parameters:
//
//   - minLevel : the minimum level at which the collective refresh must be called to ensure correctness
//   - logBound : the bit length of the masks to be sampled to mask the plaintext and ensure 128-bits of statistical indistinguishability
//   - ok 		: a boolean flag, which is set to false if no such instance exist
func GetMinimumLevelForRefresh(lambda int, scale rlwe.Scale, nParties int, moduli []uint64) (minLevel int, logBound uint, ok bool) {
	logBound = uint(lambda + int(math.Ceil(math.Log2(scale.Float64()))))
	maxBound := math.Ceil(float64(logBound) + math.Log2(float64(nParties)))

	minLevel = -1
	logQ := 0.0

	for i := 0; logQ < maxBound; i++ {

		if i >= len(moduli) {
			return 0, 0, false
		}

		logQ += math.Log2(float64(moduli[i]))
		minLevel++
	}

	if len(moduli) < minLevel {
		return 0, 0, false
	}

	return minLevel, logBound, true
}
