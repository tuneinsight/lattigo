package mkbfv

import (
	"github.com/ldsec/lattigo/v2/bfv"
	"github.com/ldsec/lattigo/v2/mkrlwe"
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/rlwe"
)

// MKCiphertext is type for a rlwe ciphertext in a multi key setting
// it contains a rlwe Element along with the list of participants in the right order
type MKCiphertext struct {
	Ciphertexts *bfv.Ciphertext
	PeerID      []uint64
}

// NewMKCiphertext returns a new MKciphertext corresponding to the given slice of peers
func NewMKCiphertext(peerIDs []uint64, r *ring.Ring, params *bfv.Parameters) *MKCiphertext {

	res := new(MKCiphertext)
	res.Ciphertexts = bfv.NewCiphertext(*params, uint64(len(peerIDs)))
	res.PeerID = peerIDs

	return res
}

// PadCiphers call the padding method from mkrlwe and put the result into new multi key ciphertexts
func PadCiphers(c1, c2 *MKCiphertext, params *bfv.Parameters) (c1Out, c2Out *MKCiphertext) {

	p1, p2, ids := mkrlwe.PadCiphers(&mkrlwe.MKCiphertext{Value: c1.Ciphertexts.Value, PeerIDs: c1.PeerID}, &mkrlwe.MKCiphertext{Value: c2.Ciphertexts.Value, PeerIDs: c2.PeerID}, &params.Parameters)

	res1 := new(bfv.Ciphertext)
	res1.Element = new(rlwe.Element)
	res1.IsNTT = false

	res2 := new(bfv.Ciphertext)
	res2.Element = new(rlwe.Element)
	res2.IsNTT = false

	res1.SetValue(p1)
	res2.SetValue(p2)

	return &MKCiphertext{res1, ids}, &MKCiphertext{res2, ids}
}
