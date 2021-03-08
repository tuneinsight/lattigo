package mkbfv

import "github.com/ldsec/lattigo/v2/ring"

// Add adds the cyphertexts component wise and expend their
func Add(c1 *MKCiphertext, c2 *MKCiphertext, res *MKCiphertext, ringQP *ring.Ring) {

	length := len(c1.ciphertexts)

	for c := uint64(0); c < uint64(length); c++ {
		ringQP.Add(c1.ciphertexts[c], c2.ciphertexts[c], res.ciphertexts[c])
	}

	MergeSlices(c1.peerIDs, c2.peerIDs, res.peerIDs)
}

// Mult will compute the homomorphic multiplication and relinearize the resulting cyphertext
func Mult() {
	// TODO: implement multiplication
}
