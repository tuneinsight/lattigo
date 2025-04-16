package sampling

import (
	"crypto/rand"
	"io"
	"sync"

	"golang.org/x/crypto/blake2b"
)

// PRNG is an interface for secure generation of random bytes
type PRNG interface {
	io.Reader
}

type ThreadSafePRNG struct {
}

// NewPRNG returns a new PRNG that is thread-safe
func NewPRNG() (*ThreadSafePRNG, error) {
	return &ThreadSafePRNG{}, nil
}

// Read reads bytes from the KeyedPRNG on sum.
func (prng *ThreadSafePRNG) Read(sum []byte) (n int, err error) {
	return rand.Read(sum)
}

// KeyedPRNG is a structure storing the parameters used to securely and *deterministically* generate shared
// sequences of random bytes among different parties using the hash function blake2b. Backward sequence
// security (given the digest i, compute the digest i-1) is ensured by default, however forward sequence
// security (given the digest i, compute the digest i+1) is only ensured if the KeyedPRNG is keyed.
// WARNING: KeyedPRNG should NOT be called by multiple threads. It does not make sense to do so as the resulting
// sequence will not be deterministic for a given key. For a PRNG securely seeded with a private keyuse [ThreadSafePRNG].
type KeyedPRNG struct {
	mutex sync.Mutex
	key   []byte
	xof   blake2b.XOF
}

// NewKeyedPRNG creates a new instance of KeyedPRNG.
// Accepts an optional key, else set key=nil which is treated as key=[]byte{}
// WARNING: A PRNG INITIALISED WITH key=nil IS INSECURE!
// WARNING: KeyedPRNG should NOT be called by multiple threads. If that occurs, the generated sequence will not be deterministic.
func NewKeyedPRNG(key []byte) (*KeyedPRNG, error) {
	var err error
	prng := new(KeyedPRNG)
	prng.xof, err = blake2b.NewXOF(blake2b.OutputLengthUnknown, key)
	return prng, err
}

// Key returns a copy of the key used to seed the PRNG.
// This value can be used with `NewKeyedPRNG` to instantiate
// a new PRNG that will produce the same stream of bytes.
func (prng *KeyedPRNG) Key() (key []byte) {
	key = make([]byte, len(prng.key))
	copy(key, prng.key)
	return
}

// Read reads bytes from the KeyedPRNG on sum.
// WARNING: Read() should NOT be called concurrently by multiple threads. If that occurs, the generated sequence will not be deterministic.
func (prng *KeyedPRNG) Read(sum []byte) (n int, err error) {
	prng.mutex.Lock()
	defer prng.mutex.Unlock()
	return prng.xof.Read(sum)
}

// Reset resets the PRNG to its initial state.
// WARNING: KeyedPRNG's methods should not be called concurrently. If that occurs, the generated sequence will not be deterministic.
func (prng *KeyedPRNG) Reset() {
	prng.mutex.Lock()
	defer prng.mutex.Unlock()
	prng.xof.Reset()
}
