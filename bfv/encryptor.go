package bfv

import (
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/rlwe"
)

// Encryptor in an interface for encryptors
type Encryptor struct {
	rlwe.Encryptor
}

func NewEncryptor(params Parameters, key interface{}) *Encryptor {
	switch key := key.(type) {
	case *rlwe.PublicKey:
		return &Encryptor{rlwe.NewEncryptorFromPk(params.Parameters, key)}
	case *rlwe.SecretKey:
		return &Encryptor{rlwe.NewEncryptorFromSk(params.Parameters, key)}
	default:
		panic("key must be either rlwe.PublicKey or rlwe.SecretKey")
	}
}

func (encryptor *Encryptor) EncryptNew(plaintext *Plaintext) *Ciphertext {
	return &Ciphertext{Element: encryptor.Encryptor.EncryptNew(&rlwe.Plaintext{Value: plaintext.Value})}
}

func (encryptor *Encryptor) Encrypt(plaintext *Plaintext, ciphertext *Ciphertext) {
	encryptor.Encryptor.Encrypt(&rlwe.Plaintext{Value: plaintext.Value}, &rlwe.Element{ciphertext.Value})
}

func (encryptor *Encryptor) EncryptFastNew(plaintext *Plaintext) *Ciphertext {
	return &Ciphertext{Element: encryptor.Encryptor.EncryptFastNew(&rlwe.Plaintext{Value: plaintext.Value})}
}

func (encryptor *Encryptor) EncryptFast(plaintext *Plaintext, ciphertext *Ciphertext) {
	encryptor.Encryptor.EncryptFast(&rlwe.Plaintext{Value: plaintext.Value}, &rlwe.Element{ciphertext.Value})
}

func (encryptor *Encryptor) EncryptFromCRP(plaintext *Plaintext, ciphertext *Ciphertext, crp *ring.Poly) {
	encryptor.Encryptor.EncryptFromCRP(&rlwe.Plaintext{Value: plaintext.Value}, &rlwe.Element{ciphertext.Value}, crp)
}

func (encryptor *Encryptor) EncryptFromCRPNew(plaintext *Plaintext, crp *ring.Poly) *Ciphertext {
	return &Ciphertext{Element: encryptor.Encryptor.EncryptFromCRPNew(&rlwe.Plaintext{Value: plaintext.Value}, crp)}
}
