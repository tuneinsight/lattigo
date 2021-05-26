package mkckks

import (
	"github.com/ldsec/lattigo/v2/ckks"
	"github.com/ldsec/lattigo/v2/mkrlwe"
	"github.com/ldsec/lattigo/v2/rlwe"
)

// MKEncryptor is an interface wrapping the ckks.Encryptor with the ring used for encryption
type MKEncryptor interface {
	Encrypt(plaintext *ckks.Plaintext) *ckks.Ciphertext
}

// mkEncryptor is a struct wrapping the ckks.Encryptor with the ring used for encryption
type mkEncryptor struct {
	ckksEncryptor ckks.Encryptor
	params        *ckks.Parameters
}

// NewMKEncryptor creates a new ckks encryptor fromm the given MKPublicKey and the ckks parameters
func NewMKEncryptor(pk *mkrlwe.MKPublicKey, params *ckks.Parameters) MKEncryptor {

	ckksPublicKey := new(rlwe.PublicKey)
	ckksPublicKey.Value[0] = pk.Key[0].Poly[0] // b[0]
	ckksPublicKey.Value[1] = pk.Key[1].Poly[0] // a[0]

	return &mkEncryptor{ckks.NewEncryptorFromPk(*params, ckksPublicKey), params}
}

// EncryptMK encrypt the plaintext and put id in the ciphertext's peerIds
func (encryptor *mkEncryptor) Encrypt(plaintext *ckks.Plaintext) *ckks.Ciphertext {

	if encryptor.params.PCount() != 0 {
		return encryptor.ckksEncryptor.EncryptNew(plaintext)
	}

	return encryptor.ckksEncryptor.EncryptFastNew(plaintext)

}
