package mkckks

import (
	"github.com/ldsec/lattigo/v2/ckks"
	"github.com/ldsec/lattigo/v2/mkrlwe"
	"github.com/ldsec/lattigo/v2/ring"
)

// MKCiphertext is type for a rlwe ciphertext in a multi key setting
// it contains a rlwe Element along with the list of participants in the right order
type MKCiphertext struct {
	Ciphertexts *ckks.Ciphertext
	PeerID      []uint64
}

// NewMKCiphertext returns a new MKciphertext corresponding to the given slice of peers
func NewMKCiphertext(peerIDs []uint64, r *ring.Ring, params *ckks.Parameters, level uint64, scale float64) *MKCiphertext {

	res := new(MKCiphertext)
	res.Ciphertexts = ckks.NewCiphertext(*params, uint64(len(peerIDs)), level, scale)
	res.PeerID = peerIDs

	return res
}

// PadCiphers call the padding method from mkrlwe and put the result into new multi key ciphertexts
func PadCiphers(c1, c2 *MKCiphertext, params *ckks.Parameters) (c1Out, c2Out *MKCiphertext) {

	p1, p2, ids := mkrlwe.PadCiphers(&mkrlwe.MKCiphertext{Value: c1.Ciphertexts.Value, PeerIDs: c1.PeerID}, &mkrlwe.MKCiphertext{Value: c2.Ciphertexts.Value, PeerIDs: c2.PeerID}, &params.Parameters)

	res1 := new(ckks.Ciphertext)
	res1.Element = new(ckks.Element)
	res1.SetScale(c1.Ciphertexts.Scale())
	res1.IsNTT = true

	res2 := new(ckks.Ciphertext)
	res2.Element = new(ckks.Element)
	res2.SetScale(c2.Ciphertexts.Scale())
	res2.IsNTT = true

	res1.SetValue(p1)
	res2.SetValue(p2)

	return &MKCiphertext{res1, ids}, &MKCiphertext{res2, ids}
}
