package mkbfv

import (
	"github.com/ldsec/lattigo/v2/bfv"
)

// MKCiphertext is type for a bfv ciphertext in a multi key setting
// it has a list of bfv ciphertexts along with the list of participants
type MKCiphertext struct {
	bfvCiphertexts []bfv.Ciphertext
	peerIDs        []uint64
}
