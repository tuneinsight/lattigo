package ring

import (
	"encoding/binary"
	"github.com/ldsec/lattigo/utils"
	"math/bits"
)

// CRPGenerator is the structure storing the parameters for deterministicaly securely
// generating random polynomials using the structure PRNG.
type CRPGenerator struct {
	prng    *utils.PRNG
	context *Context
	masks   []uint64
}

// NewCRPGenerator creates a new CRPGenerator, that will deterministicaly and securely generate uniform polynomials
// in the input context using the hash function blake2b. The PRNG can be instantiated with a key on top
// of the public seed. If not key is used, set key=nil.
func NewCRPGenerator(key []byte, context *Context) (*CRPGenerator, error) {
	var err error
	crpgenerator := new(CRPGenerator)
	crpgenerator.prng, err = utils.NewPRNG(key)
	crpgenerator.context = context
	crpgenerator.masks = make([]uint64, len(context.Modulus))
	for i, qi := range context.Modulus {
		crpgenerator.masks[i] = (1 << uint64(bits.Len64(qi))) - 1
	}
	return crpgenerator, err
}

// GetClock returns the current clock of the CRPGenerator.
func (crpgenerator *CRPGenerator) GetClock() uint64 {
	return crpgenerator.prng.GetClock()
}

// Seed resets the CRPGenerator and instantiate it with a new seed. Does not change the key.
func (crpgenerator *CRPGenerator) Seed(seed []byte) {
	crpgenerator.prng.Seed(seed)
}

// GetSeed returns the seed of the CRPGenerator.
func (crpgenerator *CRPGenerator) GetSeed() []byte {
	return crpgenerator.prng.GetSeed()[:]
}

// SetClock sets the clock of the CRPGenerator to the given input by clocking it until the
// clock cycle reach the desired number. If the given input is smaller than the current clock,
// it will return an error.
func (crpgenerator *CRPGenerator) SetClock(n uint64) error {
	if err := crpgenerator.prng.SetClock(n); err != nil {
		return err
	}

	return nil
}

// Clock generates and returns a new uniform polynomial. Also increases the clock cycle by 1.
func (crpgenerator *CRPGenerator) Clock() *Poly {
	var coeff uint64
	var randomBytes []byte

	crp := crpgenerator.context.NewPoly()

	// Starts with 32 random bytes from the prng
	randomBytes = crpgenerator.prng.Clock()

	for i := uint64(0); i < crpgenerator.context.N; i++ {

		for j, qi := range crpgenerator.context.Modulus {

			for {
				// If not enough randomBytes, generate 32 more random bytes from the PRNG
				if len(randomBytes) < 8 {
					randomBytes = crpgenerator.prng.Clock()
				}

				// Converts 4 randomBytes into an uint64 of at most 60 bits (maximum size of the modulus)
				coeff = binary.BigEndian.Uint64(randomBytes[:8]) & crpgenerator.masks[j]

				// Drops the used bytes
				randomBytes = randomBytes[8:]

				// If coeff is smaller than qi, breaks
				if coeff < qi {
					break
				}
			}

			// Assigns the coeff to the polynomial
			crp.Coeffs[j][i] = coeff
		}
	}

	return crp
}
