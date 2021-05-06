package mkckks

import (
	"github.com/ldsec/lattigo/v2/ckks"
	"github.com/ldsec/lattigo/v2/mkrlwe"
)

// MKEncryptor is an interface wrapping the ckks.Encryptor with the ring used for encryption
type MKEncryptor interface {
	EncryptMK(plaintext *ckks.Plaintext) *mkrlwe.MKCiphertext
}

// mkEncryptor is a struct wrapping the ckks.Encryptor with the ring used for encryption
type mkEncryptor struct {
	ckksEncryptor ckks.Encryptor
	peerID        uint64
	params        *ckks.Parameters
}

// NewMKEncryptor creates a new ckks encryptor fromm the given MKPublicKey and the ckks parameters
func NewMKEncryptor(pk *mkrlwe.MKPublicKey, params *ckks.Parameters, id uint64) MKEncryptor {

	ckksPublicKey := new(ckks.PublicKey)
	ckksPublicKey.Value[0] = pk.Key[0].Poly[0] // b[0]
	ckksPublicKey.Value[1] = pk.Key[1].Poly[0] // a[0]

	return &mkEncryptor{ckks.NewEncryptorFromPk(params, ckksPublicKey), pk.PeerID, params}
}

// EncryptMK encrypt the plaintext and put id in the ciphertext's peerIds
func (encryptor *mkEncryptor) EncryptMK(plaintext *ckks.Plaintext) *mkrlwe.MKCiphertext {

	mkCiphertext := new(mkrlwe.MKCiphertext)

	if encryptor.params.PiCount() != 0 {
		mkCiphertext.Ciphertexts = encryptor.ckksEncryptor.EncryptNew(plaintext)

	} else {
		mkCiphertext.Ciphertexts = encryptor.ckksEncryptor.EncryptFastNew(plaintext)
	}

	mkCiphertext.PeerIDs = []uint64{encryptor.peerID}

	return mkCiphertext
}
