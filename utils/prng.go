package utils

import (
	"crypto/rand"
	"errors"

	"golang.org/x/crypto/blake2b"
)

// PRNG is an interface for secure (keyed) deterministic generation of random bytes
type PRNG interface {
	Clock(sum []byte)
	GetClock() uint64
	SetClock(sum []byte, n uint64) error
}

// KeyedPRNG is a structure storing the parameters used to securely and deterministically generate shared
// sequences of random bytes among different parties using the hash function blake2b. Backward sequence
// security (given the digest i, compute the digest i-1) is ensured by default, however forward sequence
// security (given the digest i, compute the digest i+1) is only ensured if the KeyedPRNG is keyed.
type KeyedPRNG struct {
	clock uint64
	xof   blake2b.XOF
}

// NewKeyedPRNG creates a new instance of KeyedPRNG.
// Accepts an optional key, else set key=nil which is treated as key=[]byte{}
// WARNING: A PRNG INITIALISED WITH key=nil IS INSECURE!
func NewKeyedPRNG(key []byte) (*KeyedPRNG, error) {
	var err error
	prng := new(KeyedPRNG)
	prng.clock = 0
	prng.xof, err = blake2b.NewXOF(blake2b.OutputLengthUnknown, key)
	return prng, err
}

// NewPRNG creates KeyedPRNG keyed from rand.Read for instances were no key should be provided by the user
func NewPRNG() (*KeyedPRNG, error) {
	var err error
	prng := new(KeyedPRNG)
	prng.clock = 0
	randomBytes := make([]byte, 64)
	if _, err := rand.Read(randomBytes); err != nil {
		panic("crypto rand error")
	}
	prng.xof, err = blake2b.NewXOF(blake2b.OutputLengthUnknown, randomBytes)
	return prng, err
}

// GetClock returns the value of the clock cycle of the KeyedPRNG.
func (prng *KeyedPRNG) GetClock() uint64 {
	return prng.clock
}

// Clock reads bytes from the KeyedPRNG on sum.
func (prng *KeyedPRNG) Clock(sum []byte) {
	if _, err := prng.xof.Read(sum); err != nil {
		panic(err)
	}
	prng.clock++
}

// SetClock sets the clock cycle of the KeyedPRNG to a given number by calling Clock until
// the clock cycle reaches the desired number. Returns an error if the target clock
// cycle is smaller than the current clock cycle.
func (prng *KeyedPRNG) SetClock(sum []byte, n uint64) error {
	if prng.clock > n {
		return errors.New("error: cannot set KeyedPRNG clock to a previous state")
	}
	for prng.clock != n {
		if _, err := prng.xof.Read(sum); err != nil {
			panic(err)
		}
		prng.clock++
	}
	return nil
}
