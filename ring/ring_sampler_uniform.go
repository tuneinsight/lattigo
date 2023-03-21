package ring

import (
	"encoding/binary"

	"github.com/tuneinsight/lattigo/v4/utils/sampling"
)

// UniformSampler wraps a util.PRNG and represents the state of a sampler of uniform polynomials.
type UniformSampler struct {
	baseSampler
	randomBufferN []byte
}

// NewUniformSampler creates a new instance of UniformSampler from a PRNG and ring definition.
func NewUniformSampler(prng sampling.PRNG, baseRing *Ring) *UniformSampler {
	uniformSampler := new(UniformSampler)
	uniformSampler.baseRing = baseRing
	uniformSampler.prng = prng
	uniformSampler.randomBufferN = make([]byte, baseRing.N())
	return uniformSampler
}

// AtLevel returns an instance of the target UniformSampler that operates at the target level.
// This instance is not thread safe and cannot be used concurrently to the base instance.
func (u *UniformSampler) AtLevel(level int) *UniformSampler {
	return &UniformSampler{
		baseSampler:   u.baseSampler.AtLevel(level),
		randomBufferN: u.randomBufferN,
	}
}

// Read generates a new polynomial with coefficients following a uniform distribution over [0, Qi-1].
func (u *UniformSampler) Read(pol *Poly) {

	var randomUint, mask, qi uint64
	var ptr int

	u.prng.Read(u.randomBufferN)

	N := u.baseRing.N()

	buffer := u.randomBufferN

	for j := 0; j < u.baseRing.level+1; j++ {

		qi = u.baseRing.SubRings[j].Modulus

		// Starts by computing the mask
		mask = u.baseRing.SubRings[j].Mask

		ptmp := pol.Coeffs[j]

		// Iterates for each modulus over each coefficient
		for i := 0; i < N; i++ {

			// Samples an integer between [0, qi-1]
			for {

				// Refills the buff if it runs empty
				if ptr == N {
					u.prng.Read(buffer)
					ptr = 0
				}

				// Reads bytes from the buff
				randomUint = binary.LittleEndian.Uint64(buffer[ptr:ptr+8]) & mask
				ptr += 8

				// If the integer is between [0, qi-1], breaks the loop
				if randomUint < qi {
					break
				}
			}

			ptmp[i] = randomUint
		}
	}
}

// ReadNew generates a new polynomial with coefficients following a uniform distribution over [0, Qi-1].
// Polynomial is created at the max level.
func (u *UniformSampler) ReadNew() (Pol *Poly) {
	Pol = u.baseRing.NewPoly()
	u.Read(Pol)
	return
}

func (u *UniformSampler) WithPRNG(prng sampling.PRNG) *UniformSampler {
	return &UniformSampler{baseSampler: baseSampler{prng: prng, baseRing: u.baseRing}, randomBufferN: u.randomBufferN}
}

// RandUniform samples a uniform randomInt variable in the range [0, mask] until randomInt is in the range [0, v-1].
// mask needs to be of the form 2^n -1.
func RandUniform(prng sampling.PRNG, v uint64, mask uint64) (randomInt uint64) {
	for {
		randomInt = randInt64(prng, mask)
		if randomInt < v {
			return randomInt
		}
	}
}

// randInt32 samples a uniform variable in the range [0, mask], where mask is of the form 2^n-1, with n in [0, 32].
func randInt32(prng sampling.PRNG, mask uint64) uint64 {

	// generate random 4 bytes
	randomBytes := make([]byte, 4)
	prng.Read(randomBytes)

	// convert 4 bytes to a uint32
	randomUint32 := uint64(binary.LittleEndian.Uint32(randomBytes))

	// return required bits
	return mask & randomUint32
}

// randInt64 samples a uniform variable in the range [0, mask], where mask is of the form 2^n-1, with n in [0, 64].
func randInt64(prng sampling.PRNG, mask uint64) uint64 {

	// generate random 8 bytes
	randomBytes := make([]byte, 8)
	prng.Read(randomBytes)

	// convert 8 bytes to a uint64
	randomUint64 := binary.LittleEndian.Uint64(randomBytes)

	// return required bits
	return mask & randomUint64
}
