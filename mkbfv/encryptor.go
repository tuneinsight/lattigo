package mkbfv

import (
	"github.com/ldsec/lattigo/v2/bfv"
	"github.com/ldsec/lattigo/v2/ring"
)

// MKEncryptor is an interface wrapping the bfv.Encryptor with the ring used for encryption
type MKEncryptor interface {
	EncryptMK(plaintext *bfv.Plaintext, id uint64) *MKCiphertext
}

// mkEncryptor is a struct wrapping the bfv.Encryptor with the ring used for encryption
type mkEncryptor struct {
	bfvEncryptor bfv.Encryptor
	ringQ        *ring.Ring
}

// NewEncryptor creates a new bfv encryptor fromm the given MKPublicKey and the bfv parameters
func NewMKEncryptor(pk *MKPublicKey, params *bfv.Parameters) MKEncryptor {

	bfvPublicKey := new(bfv.PublicKey)
	bfvPublicKey.Value[0] = pk.key[0].poly[0]
	bfvPublicKey.Value[1] = pk.key[1].poly[0]

	return &mkEncryptor{bfv.NewEncryptorFromPk(params, bfvPublicKey), GetRingQ(params)}
}

// EncryptMK encrypt the plaintext and put id in the ciphertext's peerIds
func (encryptor *mkEncryptor) EncryptMK(plaintext *bfv.Plaintext, id uint64) *MKCiphertext {

	ciphertext := new(bfv.Element)
	val := make([]*ring.Poly, 2) // ciphertext has 2 components (c0, c1)
	val[0] = encryptor.ringQ.NewPoly()
	val[1] = encryptor.ringQ.NewPoly()
	ciphertext.SetValue(val)

	encryptor.bfvEncryptor.Encrypt(plaintext, ciphertext.Ciphertext())

	mkCiphertext := new(MKCiphertext)
	mkCiphertext.ciphertexts = ciphertext.Ciphertext()
	mkCiphertext.peerIDs = []uint64{id}

	return mkCiphertext
}
