package utils

import (
	"errors"
	"golang.org/x/crypto/blake2b"
)

type PRNG interface {
	GetClock() uint64
	Seed(seed []byte)
	GetSeed() []byte
	Clock(sum []byte)
	SetClock(sum []byte, n uint64) error
}

// prng is a structure storing the parameters used to securely and deterministically generate shared
// sequences of random bytes among different parties using the hash function blake2b. Backward sequence
// security (given the digest i, compute the digest i-1) is ensured by default, however forward sequence
// security (given the digest i, compute the digest i+1) is only ensured if the prng is keyed.
type prng struct {
	clock uint64
	seed  []byte
	xof   blake2b.XOF
}

// NewPRNG creates a new instance of prng.
// Accepts an optional key, else set key=nil.
func NewPRNG(key []byte) (PRNG, error) {
	var err error
	prng := new(prng)
	prng.clock = 0
	prng.xof, err = blake2b.NewXOF(blake2b.OutputLengthUnknown, key)
	return prng, err
}

// GetClock returns the value of the clock cycle of the prng.
func (prng *prng) GetClock() uint64 {
	return prng.clock
}

// Seed resets the current state of the prng (without changing the
// optional key) and seeds it with the given bytes.
// Seed will also reset the clock cycle to 0.
func (prng *prng) Seed(seed []byte) {
	prng.xof.Reset()
	prng.seed = make([]byte, len(seed))
	copy(prng.seed, seed)
	prng.xof.Write(seed)
	prng.clock = 0
}

// GetSeed returns the current seed of the prng.
func (prng *prng) GetSeed() []byte {
	return prng.seed[:]
}

// Clock reads bytes from the prng on sum.
func (prng *prng) Clock(sum []byte) {
	if _, err := prng.xof.Read(sum); err != nil {
		panic(err)
	}
	prng.clock++
}

// SetClock sets the clock cycle of the prng to a given number by calling Clock until
// the clock cycle reaches the desired number. Returns an error if the target clock
// cycle is smaller than the current clock cycle.
func (prng *prng) SetClock(sum []byte, n uint64) error {
	if prng.clock > n {
		return errors.New("error : cannot set prng clock to a previous state")
	}
	for prng.clock != n {
		if _, err := prng.xof.Read(sum); err != nil {
			panic(err)
		}
		prng.clock++
	}
	return nil
}
