package ring

import (
	"encoding/binary"

	"github.com/tuneinsight/lattigo/v4/utils"
)

// UniformSampler wraps a util.PRNG and represents the state of a sampler of uniform polynomials.
type UniformSampler struct {
	baseSampler
	randomBufferN []byte
	ptr           int
}

// NewUniformSampler creates a new instance of UniformSampler from a PRNG and ring definition.
func NewUniformSampler(prng utils.PRNG, baseRing *Ring) (u *UniformSampler) {
	u = new(UniformSampler)
	u.baseRing = baseRing
	u.prng = prng
	u.randomBufferN = make([]byte, baseRing.N())
	return
}

// Read generates a new polynomial with coefficients following a uniform distribution over [0, Qi-1].
func (u *UniformSampler) Read(pol *Poly) {
	u.ReadLvl(pol.Level(), pol)
}

// ReadLvl generates a new polynomial with coefficients following a uniform distribution over [0, Qi-1].
func (u *UniformSampler) ReadLvl(level int, pol *Poly) {

	var randomUint, mask, qi uint64

	prng := u.prng
	N := u.baseRing.N()

	var ptr int
	if ptr = u.ptr; ptr == 0 || ptr == N {
		prng.Read(u.randomBufferN)
	}

	buffer := u.randomBufferN

	for j := 0; j < level+1; j++ {

		qi = u.baseRing.Tables[j].Modulus

		// Starts by computing the mask
		mask = u.baseRing.Tables[j].Mask

		coeffs := pol.Coeffs[j]

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
				randomUint = binary.BigEndian.Uint64(buffer[ptr:ptr+8]) & mask
				ptr += 8

				// If the integer is between [0, qi-1], breaks the loop
				if randomUint < qi {
					break
				}
			}

			coeffs[i] = randomUint
		}
	}

	u.ptr = ptr
}

func (u *UniformSampler) ReadAndAddLvl(level int, pol *Poly) {
	panic("UniformSampler.ReadAndAddLvl is not implemented")
}

// ReadNew generates a new polynomial with coefficients following a uniform distribution over [0, Qi-1].
// Polynomial is created at the max level.
func (u *UniformSampler) ReadNew() (pol *Poly) {
	pol = u.baseRing.NewPoly()
	u.Read(pol)
	return
}

// ReadLvlNew generates a new polynomial with coefficients following a uniform distribution over [0, Qi-1].
// Polynomial is created at the specified level.
func (u *UniformSampler) ReadLvlNew(level int) (pol *Poly) {
	pol = u.baseRing.NewPolyLvl(level)
	u.ReadLvl(level, pol)
	return
}

func (u *UniformSampler) WithPRNG(prng utils.PRNG) *UniformSampler {
	return &UniformSampler{baseSampler: baseSampler{prng: prng, baseRing: u.baseRing}, randomBufferN: u.randomBufferN}
}

// RandUniform samples a uniform randomInt variable in the range [0, mask] until randomInt is in the range [0, v-1].
// mask needs to be of the form 2^n -1.
func RandUniform(prng utils.PRNG, v uint64, mask uint64) (randomInt uint64) {
	for {
		randomInt = randInt64(prng, mask)
		if randomInt < v {
			return randomInt
		}
	}
}

// randInt32 samples a uniform variable in the range [0, mask], where mask is of the form 2^n-1, with n in [0, 32].
func randInt32(prng utils.PRNG, mask uint64) uint64 {

	// generate random 4 bytes
	randomBytes := make([]byte, 4)
	prng.Read(randomBytes)

	// convert 4 bytes to a uint32
	randomUint32 := uint64(binary.BigEndian.Uint32(randomBytes))

	// return required bits
	return mask & randomUint32
}

// randInt64 samples a uniform variable in the range [0, mask], where mask is of the form 2^n-1, with n in [0, 64].
func randInt64(prng utils.PRNG, mask uint64) uint64 {

	// generate random 8 bytes
	randomBytes := make([]byte, 8)
	prng.Read(randomBytes)

	// convert 8 bytes to a uint64
	randomUint64 := binary.BigEndian.Uint64(randomBytes)

	// return required bits
	return mask & randomUint64
}
