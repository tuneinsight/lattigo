package mkbfv

import (
	"github.com/ldsec/lattigo/v2/bfv"
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

// NewEncryptor creates a new bfv encryptor fromm the given MKPublicKey and the bfv parameters
func NewMKEncryptor(pk *MKPublicKey, params *bfv.Parameters, id uint64) MKEncryptor {

	bfvPublicKey := new(bfv.PublicKey)
	bfvPublicKey.Value[0] = pk.key[0].poly[0] // b[0]
	bfvPublicKey.Value[1] = pk.key[1].poly[0] // a[0]

	return &mkEncryptor{bfv.NewEncryptorFromPk(params, bfvPublicKey), pk.peerID, params}
}

// EncryptMK encrypt the plaintext and put id in the ciphertext's peerIds
func (encryptor *mkEncryptor) EncryptMK(plaintext *bfv.Plaintext) *MKCiphertext {

	mkCiphertext := new(MKCiphertext)

	if encryptor.params.PiCount() != 0 {
		mkCiphertext.ciphertexts = encryptor.bfvEncryptor.EncryptNew(plaintext)
	} else {
		mkCiphertext.ciphertexts = encryptor.bfvEncryptor.EncryptFastNew(plaintext)
	}

	mkCiphertext.peerIDs = []uint64{encryptor.peerID}

	return mkCiphertext
}
