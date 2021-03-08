package mkbfv

import (
	"github.com/ldsec/lattigo/v2/bfv"
	"github.com/ldsec/lattigo/v2/ring"
)

// PartDec computes a partial decription key for the ciphertext component of a given participant
// for participant i, ski and cti must be used
func PartDec(ct *ring.Poly, sk *MKSecretKey, out *ring.Poly) {
	// TODO: implement partDec section 4.3 Chen et al.
}

// MergeDec merges the partial decription parts and returns the plaintext
func MergeDec(ct *ring.Poly, partialKeys []*ring.Poly, out *bfv.Plaintext) {
	// TODO: implement merge as described in section 4.3 Chen et al.
}
