package utils

import (
	"crypto/rand"
	"errors"
	"golang.org/x/crypto/blake2b"
)

type PRNG interface {
	Clock(sum []byte)
}

// OSrandPRNG is an empty struct used as a method receiver for crypto/rand Read
type OSrandPRNG struct{}

// Clock reads bytes from crypto/rand on sum.
func (prng *OSrandPRNG) Clock(sum []byte) {
	if _, err := rand.Read(sum); err != nil {
		panic("crypto rand error")
	}
}

// seededPRNG is a structure storing the parameters used to securely and deterministically generate shared
// sequences of random bytes among different parties using the hash function blake2b. Backward sequence
// security (given the digest i, compute the digest i-1) is ensured by default, however forward sequence
// security (given the digest i, compute the digest i+1) is only ensured if the seededPRNG is keyed.
type SeededPRNG struct {
	clock uint64
	seed  []byte
	xof   blake2b.XOF
}

// NewSeededPRNG creates a new instance of seededPRNG.
// Accepts an optional key, else set key=nil.
func NewSeededPRNG(key []byte) (*SeededPRNG, error) {
	var err error
	prng := new(SeededPRNG)
	prng.clock = 0
	prng.xof, err = blake2b.NewXOF(blake2b.OutputLengthUnknown, key)
	return prng, err
}

// GetClock returns the value of the clock cycle of the seededPRNG.
func (prng *SeededPRNG) GetClock() uint64 {
	return prng.clock
}

// Seed resets the current state of the seededPRNG (without changing the
// optional key) and seeds it with the given bytes.
// Seed will also reset the clock cycle to 0.
func (prng *SeededPRNG) Seed(seed []byte) {
	prng.xof.Reset()
	prng.seed = make([]byte, len(seed))
	copy(prng.seed, seed)
	prng.xof.Write(seed)
	prng.clock = 0
}

// GetSeed returns the current seed of the seededPRNG.
func (prng *SeededPRNG) GetSeed() []byte {
	return prng.seed[:]
}

// Clock reads bytes from the seededPRNG on sum.
func (prng *SeededPRNG) Clock(sum []byte) {
	if _, err := prng.xof.Read(sum); err != nil {
		panic(err)
	}
	prng.clock++
}

// SetClock sets the clock cycle of the seededPRNG to a given number by calling Clock until
// the clock cycle reaches the desired number. Returns an error if the target clock
// cycle is smaller than the current clock cycle.
func (prng *SeededPRNG) SetClock(sum []byte, n uint64) error {
	if prng.clock > n {
		return errors.New("error : cannot set seededPRNG clock to a previous state")
	}
	for prng.clock != n {
		if _, err := prng.xof.Read(sum); err != nil {
			panic(err)
		}
		prng.clock++
	}
	return nil
}
