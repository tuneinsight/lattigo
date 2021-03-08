package mkbfv

import "github.com/ldsec/lattigo/v2/bfv"

// NewEncryptor creates a new bfv encryptor fromm the given MKPublicKey and the bfv parameters
func NewEncryptor(pk *MKPublicKey, params *bfv.Parameters) bfv.Encryptor {

	bfvPublicKey := new(bfv.PublicKey)
	bfvPublicKey.Value[0] = pk.key[0].poly[0]
	bfvPublicKey.Value[1] = pk.key[1].poly[0]

	return bfv.NewEncryptorFromPk(params, bfvPublicKey)
}
