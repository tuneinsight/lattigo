package ring

import (
	"encoding/binary"

	"github.com/tuneinsight/lattigo/v6/utils/sampling"
)

// UniformSampler wraps a util.PRNG and represents the state of a sampler of uniform polynomials.
type UniformSampler struct {
	*baseSampler
	*randomBuffer
}

// NewUniformSampler creates a new instance of UniformSampler from a PRNG and ring definition.
func NewUniformSampler(prng sampling.PRNG, baseRing *Ring) (u *UniformSampler) {
	u = new(UniformSampler)
	u.baseSampler = &baseSampler{}
	u.baseRing = baseRing
	u.prng = prng
	u.randomBuffer = newRandomBuffer()
	return
}

// AtLevel returns an instance of the target UniformSampler to sample at the given level.
// The returned sampler cannot be used concurrently to the original sampler.
func (u *UniformSampler) AtLevel(level int) Sampler {
	return &UniformSampler{
		baseSampler:  u.baseSampler.AtLevel(level),
		randomBuffer: u.randomBuffer,
	}
}

func (u *UniformSampler) Read(pol Poly) {
	u.read(pol, func(a, b, c uint64) uint64 {
		return b
	})
}

func (u *UniformSampler) ReadAndAdd(pol Poly) {
	u.read(pol, func(a, b, c uint64) uint64 {
		return CRed(a+b, c)
	})
}

func (u *UniformSampler) read(pol Poly, f func(a, b, c uint64) uint64) {

	level := u.baseRing.Level()

	var randomUint, mask, qi uint64

	prng := u.prng
	N := u.baseRing.N()
	byteArrayLength := len(u.randomBufferN)

	var ptr int
	if ptr = u.ptr; ptr == 0 || ptr == byteArrayLength {
		if _, err := prng.Read(u.randomBufferN); err != nil {
			// Sanity check, this error should not happen.
			panic(err)
		}
		ptr = 0 // for the case where ptr == byteArrayLength
	}

	buffer := u.randomBufferN

	for j := 0; j < level+1; j++ {

		qi = u.baseRing.SubRings[j].Modulus

		// Starts by computing the mask
		mask = u.baseRing.SubRings[j].Mask

		coeffs := pol.Coeffs[j]

		// Iterates for each modulus over each coefficient
		for i := 0; i < N; i++ {

			// Samples an integer between [0, qi-1]
			for {

				// Refills the buff if it runs empty
				if ptr == byteArrayLength {
					if _, err := u.prng.Read(buffer); err != nil {
						// Sanity check, this error should not happen.
						panic(err)
					}
					ptr = 0
				}

				// Reads bytes from the buff
				randomUint = binary.BigEndian.Uint64(buffer[ptr:ptr+8]) & mask
				ptr += 8

				// If the integer is between [0, qi-1], breaks the loop
				if randomUint < qi {
					break
				}
			}

			coeffs[i] = f(coeffs[i], randomUint, qi)
		}
	}

	u.ptr = ptr
}

// ReadNew generates a new polynomial with coefficients following a uniform distribution over [0, Qi-1].
// Polynomial is created at the max level.
func (u *UniformSampler) ReadNew() (pol Poly) {
	pol = u.baseRing.NewPoly()
	u.Read(pol)
	return
}

func (u *UniformSampler) WithPRNG(prng sampling.PRNG) *UniformSampler {
	return &UniformSampler{
		baseSampler: &baseSampler{
			prng:     prng,
			baseRing: u.baseRing,
		},
		randomBuffer: newRandomBuffer(),
	}
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
	if _, err := prng.Read(randomBytes); err != nil {
		// Sanity check, this error should not happen.
		panic(err)
	}

	// convert 4 bytes to a uint32
	randomUint32 := uint64(binary.BigEndian.Uint32(randomBytes))

	// return required bits
	return mask & randomUint32
}

// randInt64 samples a uniform variable in the range [0, mask], where mask is of the form 2^n-1, with n in [0, 64].
func randInt64(prng sampling.PRNG, mask uint64) uint64 {

	// generate random 8 bytes
	randomBytes := make([]byte, 8)
	if _, err := prng.Read(randomBytes); err != nil {
		// Sanity check, this error should not happen.
		panic(err)
	}

	// convert 8 bytes to a uint64
	randomUint64 := binary.BigEndian.Uint64(randomBytes)

	// return required bits
	return mask & randomUint64
}
