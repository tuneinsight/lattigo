package dbfv

import (
	"encoding/binary"
	"errors"
	"github.com/ldsec/lattigo/ring"
	"golang.org/x/crypto/blake2b"
	"hash"
	"math/bits"
)

type PRNG struct {
	clock uint64
	seed  []byte
	hash  hash.Hash
}

// NewPRNG creates a new instance of PRNG
// Accepts an optional key, if no key, call with key=nill
func NewPRNG(key []byte) (*PRNG, error) {
	var err error
	prng := new(PRNG)
	prng.clock = 0
	prng.hash, err = blake2b.New512(key)
	return prng, err
}

// GetClock returns the value of the clock cycle of the PRNG
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

// Clock returns the right part of the digest value
// of the current PRNG state and reseeds the PRNG with the
// left part of the digest value. Increases the clock cycle by 1.
func (prng *PRNG) Clock() []byte {
	tmp := prng.hash.Sum(nil)
	prng.hash.Write(tmp[:32])
	prng.clock += 1
	return tmp[32:]
}

// Sets the clock cycle of the PRNG to a given number. Returns an error if
// the target clock cycle is smaller than the current clock cycle.
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

// Common Reference Polynomial Generator
type CRPGenerator struct {
	prng    *PRNG
	context *ring.Context
	masks   []uint64
}

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

func (crpgenerator *CRPGenerator) GetClock() uint64 {
	return crpgenerator.prng.clock
}

func (crpgenerator *CRPGenerator) Seed(seed []byte) {
	crpgenerator.prng.Seed(seed)
}

func (crpgenerator *CRPGenerator) GetSeed() []byte {
	return crpgenerator.prng.seed[:]
}

func (crpgenerator *CRPGenerator) SetClock(n uint64) error {
	if err := crpgenerator.prng.SetClock(n); err != nil {
		return err
	}

	return nil
}

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

			// Assigns coeff to the ring crp coefficient
			crp.Coeffs[j][i] = coeff
		}
	}

	return crp
}
