package ring

import (
	"crypto/rand"
	"encoding/binary"
)

// UniformPoly generates a new polynomial with coefficients following a uniform distribution over [0, Qi-1].
func (context *Context) UniformPoly(Pol *Poly) {
	context.UniformPolyLvl(uint64(len(Pol.Coeffs)-1), Pol)
}

// UniformPoly generates a new polynomial with coefficients following a uniform distribution over [0, Qi-1].
func (context *Context) UniformPolyLvl(level uint64, Pol *Poly) {

	var randomUint, mask, qi uint64
	var ptr uint64

	randomBytes := make([]byte, context.N)
	if _, err := rand.Read(randomBytes); err != nil {
		panic("crypto rand error")
	}

	for j := uint64(0); j < level; j++ {

		qi = context.Modulus[j]

		// Starts by computing the mask
		mask = context.mask[j]

		ptmp := Pol.Coeffs[j]

		// Iterates for each modulus over each coefficient
		for i := uint64(0); i < context.N; i++ {

			// Samples an integer between [0, qi-1]
			for {

				// Replenishes the pool if it runs empty
				if ptr == context.N {
					if _, err := rand.Read(randomBytes); err != nil {
						panic("crypto rand error")
					}
					ptr = 0
				}

				// Reads bytes from the pool
				randomUint = binary.BigEndian.Uint64(randomBytes[ptr:ptr+8]) & mask
				ptr += 8

				// If the integer is between [0, qi-1], breaks the loop
				if randomUint < qi {
					break
				}
			}

			ptmp[i] = randomUint
		}
	}

	return
}

// NewUniformPoly generates a new polynomial with coefficients following a uniform distribution over [0, Qi-1].
func (context *Context) NewUniformPoly() (Pol *Poly) {

	Pol = context.NewPoly()

	context.UniformPoly(Pol)

	return
}

// NewUniformPolyLvl generates a new polynomial with coefficients following a uniform distribution over [0, Qi-1].
func (context *Context) NewUniformPolyLvl(level uint64) (Pol *Poly) {

	Pol = context.NewPolyLvl(level)

	context.UniformPoly(Pol)

	return
}

// RandUniform samples a uniform randomInt variable in the range [0, mask] until randomInt is in the range [0, v-1].
// mask needs to be of the form 2^n -1.
func RandUniform(v uint64, mask uint64) (randomInt uint64) {
	for {
		randomInt = randInt64(mask)
		if randomInt < v {
			return randomInt
		}
	}
}

// randInt32 samples a uniform variable in the range [0, mask], where mask is of the form 2^n-1, with n in [0, 32].
func randInt32(mask uint64) uint64 {

	// generate random 4 bytes
	randomBytes := make([]byte, 4)
	_, err := rand.Read(randomBytes)
	if err != nil {
		panic("crypto rand error")
	}

	// convert 4 bytes to a uint32
	randomUint32 := uint64(binary.BigEndian.Uint32(randomBytes))

	// return required bits
	return mask & randomUint32
}

// randInt64 samples a uniform variable in the range [0, mask], where mask is of the form 2^n-1, with n in [0, 64].
func randInt64(mask uint64) uint64 {

	// generate random 8 bytes
	randomBytes := make([]byte, 8)
	_, err := rand.Read(randomBytes)
	if err != nil {
		panic("crypto rand error")
	}

	// convert 8 bytes to a uint64
	randomUint64 := binary.BigEndian.Uint64(randomBytes)

	// return required bits
	return mask & randomUint64
}
