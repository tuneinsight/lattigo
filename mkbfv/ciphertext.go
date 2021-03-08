package mkbfv

import "github.com/ldsec/lattigo/v2/ring"

// MKCiphertext is type for a bfv ciphertext in a multi key setting
// it has a list of bfv ciphertexts along with the list of participants
type MKCiphertext struct {
	ciphertexts []*ring.Poly
	peerIDs     []uint64
}
