package mkrlwe

import (
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/rlwe"
)

// MKCiphertext is type for a rlwe ciphertext in a multi key setting
// it contains a rlwe Element along with the list of participants in the right order
type MKCiphertext struct {
	Value   []*ring.Poly
	PeerIDs []uint64
}

// NewMKCiphertext returns a new MKciphertext corresponding to the given slice of peers with values initialized to 0
func NewMKCiphertext(peerIDs []uint64, r *ring.Ring, params *rlwe.Parameters) *MKCiphertext {

	poly := make([]*ring.Poly, len(peerIDs)+1)

	for i := 0; i < len(peerIDs); i++ {
		poly[i] = r.NewPoly()
	}

	res := new(MKCiphertext)
	res.Value = poly
	res.PeerIDs = peerIDs

	return res
}

// PadCiphers pad two ciphertext corresponding to a different set of parties
// to make their dimension match. ci = 0 if the participant i is not involved in the ciphertext
// the peerIDs list is also updated to match the new dimension
func PadCiphers(c1, c2 *MKCiphertext, params *rlwe.Parameters) (c1Out, c2Out *MKCiphertext) {

	ringQ := GetRingQ(params)
	allPeers := MergeSlices(c1.PeerIDs, c2.PeerIDs)
	k := len(allPeers)

	res1 := make([]*ring.Poly, k+1) // + 1 for c0
	res2 := make([]*ring.Poly, k+1)

	// put c0 in
	res1[0] = c1.Value[0].CopyNew()
	res2[0] = c2.Value[0].CopyNew()

	// copy ciphertext values if participant involved
	// else put a 0 polynomial
	for i, peer := range allPeers {
		index := i + 1

		index1 := Contains(c1.PeerIDs, peer)
		index2 := Contains(c2.PeerIDs, peer)

		if index1 >= 0 {
			res1[index] = c1.Value[index1+1].CopyNew()
		} else {
			res1[index] = nil
		}

		if index2 >= 0 {
			res2[index] = c2.Value[index2+1].CopyNew()
		} else {
			res2[index] = nil
		}
	}

	c1out := NewMKCiphertext(allPeers, ringQ, params)
	c2out := NewMKCiphertext(allPeers, ringQ, params)

	c1out.Value = res1
	c2out.Value = res2

	return c1out, c2out
}
