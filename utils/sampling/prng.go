package sampling

import (
	"crypto/rand"
	"io"

	"golang.org/x/crypto/blake2b"
)

// PRNG is an interface for secure (keyed) deterministic generation of random bytes
type PRNG interface {
	io.Reader
}

// KeyedPRNG is a structure storing the parameters used to securely and deterministically generate shared
// sequences of random bytes among different parties using the hash function blake2b. Backward sequence
// security (given the digest i, compute the digest i-1) is ensured by default, however forward sequence
// security (given the digest i, compute the digest i+1) is only ensured if the KeyedPRNG is keyed.
type KeyedPRNG struct {
	key []byte
	xof blake2b.XOF
}

// NewKeyedPRNG creates a new instance of KeyedPRNG.
// Accepts an optional key, else set key=nil which is treated as key=[]byte{}
// WARNING: A PRNG INITIALISED WITH key=nil IS INSECURE!
func NewKeyedPRNG(key []byte) (*KeyedPRNG, error) {
	var err error
	prng := new(KeyedPRNG)
	prng.xof, err = blake2b.NewXOF(blake2b.OutputLengthUnknown, key)
	return prng, err
}

// NewPRNG creates KeyedPRNG keyed from rand.Read for instances were no key should be provided by the user
func NewPRNG() (*KeyedPRNG, error) {
	var err error
	prng := new(KeyedPRNG)
	key := make([]byte, 64)
	if _, err := rand.Read(key); err != nil {
		panic("crypto rand error")
	}
	prng.key = key
	prng.xof, err = blake2b.NewXOF(blake2b.OutputLengthUnknown, key)
	return prng, err
}

// Read reads bytes from the KeyedPRNG on sum.
func (prng *KeyedPRNG) Read(sum []byte) (n int, err error) {
	if n, err = prng.xof.Read(sum); err != nil {
		panic(err)
	}
	return n, nil
}

// Reset resets the PRNG to its initial state.
func (prng *KeyedPRNG) Reset() {
	prng.xof.Reset()
}
