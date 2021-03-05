package mkbfv

import (
	"github.com/ldsec/lattigo/v2/bfv"
	"github.com/ldsec/lattigo/v2/ring"
)

// Convert creates a switching key K_ij from the evaluation key of the i-th peer and the public key of the j-th peer.
func Convert(D *MKEvaluationKey, publicKey *MKPublicKey, params *bfv.Parameters) *MKSwitchingKey {
	res := new(MKSwitchingKey)
	d0 := D.key[0]
	d1 := D.key[1]
	beta := params.Beta() // size of the decomposition
	k0 := new(MKDecomposedPoly)
	k1 := new(MKDecomposedPoly)
	k2 := new(MKDecomposedPoly)
	ring := GetRing(params)

	for l := uint64(0); l < beta; l++ {
		gInv := GInverse(publicKey.key[1].poly[l])
		Dot(gInv, d0, k0.poly[l], ring, beta) // TODO : g^-1 ( b_j[l]) dot d_i[0]
		Dot(gInv, d1, k1.poly[l], ring, beta) // TODO : g^-1 ( b_j[l]) dot d_i[1]

	}
	copy(k2.poly, D.key[2].poly) // TODO : We could also use the ring_poly.Copy method, if this does not do a deep copy
	res.key = [3]*MKDecomposedPoly{k0, k1, k2}

	return res
}

// CreateSharedRelinearizationKey generates a shared relinearization key containing the switching key for all pair of participants.
func CreateSharedRelinearizationKey(params *bfv.Parameters) *MKRelinearizationKey {
	res := new(MKRelinearizationKey)

	return res
}

// GInverse is a method that returns g^(-1). TODO : See what team says regarding if the algo already exists somewhere.
func GInverse(ring *ring.Poly) *MKDecomposedPoly {
	return new(MKDecomposedPoly)
}
