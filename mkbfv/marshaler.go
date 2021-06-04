package mkbfv

import (
	"encoding/binary"

	"github.com/ldsec/lattigo/v2/bfv"
	"github.com/ldsec/lattigo/v2/mkrlwe"
)

// MarshalBinary encodes a MKCiphertext in a byte slice.
func (ciphertext *MKCiphertext) MarshalBinary() (data []byte, err error) {

	d1, err := ciphertext.Ciphertexts.MarshalBinary()

	if err != nil {
		panic("Error while marshaling ciphertext")
	}

	d2 := mkrlwe.MarshallPeerID(ciphertext.PeerID)
	data = mkrlwe.MergeTwoByteSlices(d1, d2)

	return data, nil
}

// UnmarshalBinary decodes a previously marshaled MKCiphertext in the target MKCiphertext.
func (ciphertext *MKCiphertext) UnmarshalBinary(data []byte) (err error) {

	ciphertext.Ciphertexts = new(bfv.Ciphertext)
	separator := binary.LittleEndian.Uint64(data[:8])
	ciphertext.Ciphertexts.UnmarshalBinary(data[8 : separator+8])
	ciphertext.PeerID = mkrlwe.UnMarshallPeerID(data[separator+8:])

	return nil
}
