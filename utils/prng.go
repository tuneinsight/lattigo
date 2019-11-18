package utils

import (
	"errors"
	"golang.org/x/crypto/blake2b"
	"hash"
)

// PRNG is a structure storing the parameters used to securely and deterministicaly generate shared
// sequences of random bytes among different parties using the hash function blake2b. Backward sequence
// security (given the digest i, compute the digest i-1) is ensured by default, however forward sequence
// security (given the digest i, compute the digest i+1) is only ensured if the PRNG is given a key.
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

// Clock returns the right 64 bytes digest value of the current
// PRNG state and reseeds the PRNG with those same 64 bytes.
// Also increases the clock cycle by 1.
func (prng *PRNG) Clock() []byte {
	tmp := prng.hash.Sum(nil)
	prng.hash.Write(tmp)
	prng.clock += 1
	return tmp
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
		prng.hash.Write(tmp)
		prng.clock += 1
	}
	return nil
}
