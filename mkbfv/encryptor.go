package mkbfv

import (
	"github.com/ldsec/lattigo/v2/bfv"
	"github.com/ldsec/lattigo/v2/mkrlwe"
	"github.com/ldsec/lattigo/v2/rlwe"
)

// MKEncryptor is an interface wrapping the bfv.Encryptor with the ring used for encryption
type MKEncryptor interface {
	Encrypt(plaintext *bfv.Plaintext) *bfv.Ciphertext
}

// mkEncryptor is a struct wrapping the bfv.Encryptor with the ring used for encryption
type mkEncryptor struct {
	bfvEncryptor bfv.Encryptor
	params       *bfv.Parameters
}

// NewMKEncryptor creates a new bfv encryptor fromm the given MKPublicKey and the bfv parameters
func NewMKEncryptor(pk *mkrlwe.MKPublicKey, params *bfv.Parameters) MKEncryptor {

	bfvPublicKey := new(rlwe.PublicKey)
	bfvPublicKey.Value[0] = pk.Key[0].Poly[0] // b[0]
	bfvPublicKey.Value[1] = pk.Key[1].Poly[0] // a[0]

	return &mkEncryptor{bfv.NewEncryptorFromPk(*params, bfvPublicKey), params}
}

// EncryptMK encrypt the plaintext and put id in the ciphertext's peerIds
func (encryptor *mkEncryptor) Encrypt(plaintext *bfv.Plaintext) *bfv.Ciphertext {

	if encryptor.params.PCount() != 0 {
		return encryptor.bfvEncryptor.EncryptNew(plaintext)
	}

	return encryptor.bfvEncryptor.EncryptFastNew(plaintext)
}
