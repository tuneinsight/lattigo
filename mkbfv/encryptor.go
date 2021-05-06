package mkbfv

import (
	"github.com/ldsec/lattigo/v2/bfv"
	"github.com/ldsec/lattigo/v2/mkrlwe"
	"github.com/ldsec/lattigo/v2/rlwe"
)

// MKEncryptor is an interface wrapping the bfv.Encryptor with the ring used for encryption
type MKEncryptor interface {
	EncryptMK(plaintext *bfv.Plaintext) *MKCiphertext
}

// mkEncryptor is a struct wrapping the bfv.Encryptor with the ring used for encryption
type mkEncryptor struct {
	bfvEncryptor bfv.Encryptor
	peerID       uint64
	params       *bfv.Parameters
}

// NewMKEncryptor creates a new bfv encryptor fromm the given MKPublicKey and the bfv parameters
func NewMKEncryptor(pk *mkrlwe.MKPublicKey, params *bfv.Parameters, id uint64) MKEncryptor {

	bfvPublicKey := new(rlwe.PublicKey)
	bfvPublicKey.Value[0] = pk.Key[0].Poly[0] // b[0]
	bfvPublicKey.Value[1] = pk.Key[1].Poly[0] // a[0]

	return &mkEncryptor{bfv.NewEncryptorFromPk(*params, bfvPublicKey), pk.PeerID, params}
}

// EncryptMK encrypt the plaintext and put id in the ciphertext's peerIds
func (encryptor *mkEncryptor) EncryptMK(plaintext *bfv.Plaintext) *MKCiphertext {

	mkCiphertext := new(MKCiphertext)

	if encryptor.params.PCount() != 0 {
		mkCiphertext.Ciphertexts = encryptor.bfvEncryptor.EncryptNew(plaintext)
	} else {
		mkCiphertext.Ciphertexts = encryptor.bfvEncryptor.EncryptFastNew(plaintext)
	}

	mkCiphertext.PeerID = []uint64{encryptor.peerID}

	return mkCiphertext
}
