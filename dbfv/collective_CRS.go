package dbfv

import (
	"encoding/binary"
	"errors"
	"github.com/lca1/lattigo/ring"
	"golang.org/x/crypto/blake2b"
	"hash"
	"math/bits"
)

// PRNG is a structure storing the parameters used to securely and deterministicaly generate shared
// random bytes among different parties using the hash function blake2b.
type PRNG struct {
	clock uint64
	seed  []byte
	hash  hash.Hash
}

// NewPRNG creates a new instance of PRNG.
// Accepts an optional key, else set key=nil.
func NewPRNG(key []byte) (*PRNG, error) {
	var err error
	prng := new(PRNG)
	prng.clock = 0
	prng.hash, err = blake2b.New512(key)
	return prng, err
}

// GetClock returns the value of the clock cycle of the PRNG.
func (prng *PRNG) GetClock() uint64 {
	return prng.clock
}

// Seed resets the current state of the PRNG (without changing the
// optional key) and seeds it with the given bytes.
// Seed will also reset the clock cycle to 0.
func (prng *PRNG) Seed(seed []byte) {
	prng.hash.Reset()
	prng.seed = seed[:]
	prng.hash.Write(seed)
	prng.clock = 0
}

func (prng *PRNG) GetSeed() []byte {
	return prng.seed[:]
}

// Clock returns the right 32 bytes of the digest value
// of the current PRNG state and reseeds the PRNG with the
// left 32 bytes of the digest value. Also increases the clock cycle by 1.
func (prng *PRNG) Clock() []byte {
	tmp := prng.hash.Sum(nil)
	prng.hash.Write(tmp[:32])
	prng.clock += 1
	return tmp[32:]
}

// SetClock sets the clock cycle of the PRNG to a given number by calling Clock until
// the clock cycle reaches the desired number. Returns an error if the target clock
// cycle is smaller than the current clock cycle.
func (prng *PRNG) SetClock(n uint64) error {
	if prng.clock > n {
		return errors.New("error : cannot set prng clock to a previous state")
	}
	var tmp []byte
	for prng.clock != n {
		tmp = prng.hash.Sum(nil)
		prng.hash.Write(tmp[:32])
		prng.clock += 1
	}
	return nil
}

// CRPGenerator is the structure storing the parameters for deterministicaly securely
// generating random polynomials using the structure PRNG.
type CRPGenerator struct {
	prng    *PRNG
	context *ring.Context
	masks   []uint64
}

// NewCRPGenerator creates a new CRPGenerator, that will deterministicaly and securely generate uniform polynomials
// in the input context using the hash function blake2b. The PRNG can be instantiated with a key on top
// of the public seed. If not key is used, set key=nil.
func NewCRPGenerator(key []byte, context *ring.Context) (*CRPGenerator, error) {
	var err error
	crpgenerator := new(CRPGenerator)
	crpgenerator.prng, err = NewPRNG(key)
	crpgenerator.context = context
	crpgenerator.masks = make([]uint64, len(context.Modulus))
	for i, qi := range context.Modulus {
		crpgenerator.masks[i] = (1 << uint64(bits.Len64(qi))) - 1
	}
	return crpgenerator, err
}

// GetClock returns the current clock of the CRPGenerator.
func (crpgenerator *CRPGenerator) GetClock() uint64 {
	return crpgenerator.prng.clock
}

// Seed resets the CRPGenerator and instantiate it with a new seed. Does not change the key.
func (crpgenerator *CRPGenerator) Seed(seed []byte) {
	crpgenerator.prng.Seed(seed)
}

// GetSeed returns the seed of the CRPGenerator.
func (crpgenerator *CRPGenerator) GetSeed() []byte {
	return crpgenerator.prng.seed[:]
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
func (crpgenerator *CRPGenerator) Clock() *ring.Poly {
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
