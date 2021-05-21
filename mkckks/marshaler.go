package mkckks

import (
	"github.com/ldsec/lattigo/v2/ckks"
	"github.com/ldsec/lattigo/v2/mkrlwe"
)

// MarshalBinary encodes a MKCiphertext in a byte slice.
func (ciphertext *MKCiphertext) MarshalBinary() (data [2][]byte) {

	var err error
	data[0], err = ciphertext.Ciphertexts.MarshalBinary()

	if err != nil {
		panic("Error while marshaling ciphertext")
	}

	data[1] = mkrlwe.MarshallPeerID(ciphertext.PeerID)

	return data
}

// UnmarshalBinary decodes a previously marshaled MKCiphertext in the target MKCiphertext.
func (ciphertext *MKCiphertext) UnmarshalBinary(data [2][]byte) (err error) {

	ciphertext.Ciphertexts = new(ckks.Ciphertext)
	ciphertext.Ciphertexts.UnmarshalBinary(data[0])
	ciphertext.PeerID = mkrlwe.UnMarshallPeerID(data[1])

	return nil
}
