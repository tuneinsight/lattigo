package mkckks

import (
	"github.com/ldsec/lattigo/v2/ckks"
	"github.com/ldsec/lattigo/v2/ring"
)

// MKCiphertext is type for a ckks ciphertext in a multi key setting
// it contains a ckks ciphertexts along with the list of participants in the right order
type MKCiphertext struct {
	ciphertexts *ckks.Ciphertext
	peerIDs     []uint64
}

// PadCiphers pad two ciphertext corresponding to a different set of parties
// to make their dimension match. ci = 0 if the participant i is not involved in the ciphertext
// the peerIDs list is also updated to match the new dimension
func PadCiphers(c1, c2 *MKCiphertext, params *ckks.Parameters) (c1Out, c2Out *MKCiphertext) {

	ringQ := GetRingQ(params)
	allPeers := MergeSlices(c1.peerIDs, c2.peerIDs)
	k := len(allPeers)

	res1 := make([]*ring.Poly, k+1) // + 1 for c0
	res2 := make([]*ring.Poly, k+1)

	// put c0 in
	res1[0] = c1.ciphertexts.Element.Value()[0].CopyNew()
	res2[0] = c2.ciphertexts.Element.Value()[0].CopyNew()

	// copy ciphertext values if participant involved
	// else put a 0 polynomial
	for i, peer := range allPeers {
		index := i + 1

		index1 := Contains(c1.peerIDs, peer)
		index2 := Contains(c2.peerIDs, peer)

		if index1 >= 0 {
			res1[index] = c1.ciphertexts.Element.Value()[index1+1].CopyNew()
		} else {
			res1[index] = ringQ.NewPoly()
		}

		if index2 >= 0 {
			res2[index] = c2.ciphertexts.Element.Value()[index2+1].CopyNew()
		} else {
			res2[index] = ringQ.NewPoly()
		}
	}

	c1out := NewMKCiphertext(allPeers, ringQ, params, c1.ciphertexts.Level())
	c2out := NewMKCiphertext(allPeers, ringQ, params, c2.ciphertexts.Level())

	c1out.ciphertexts.SetValue(res1)
	c2out.ciphertexts.SetValue(res2)

	return c1out, c2out
}

// NewMKCiphertext returns a new MKciphertext corresponding to the given slice of peers
func NewMKCiphertext(peerIDs []uint64, r *ring.Ring, params *ckks.Parameters, level uint64) *MKCiphertext {

	res := new(MKCiphertext)
	res.ciphertexts = ckks.NewCiphertext(params, uint64(len(peerIDs)), level, params.Scale())
	res.peerIDs = peerIDs

	return res
}
