package mkbfv

import (
	"github.com/ldsec/lattigo/v2/bfv"
	"github.com/ldsec/lattigo/v2/ring"
)

// KeyGen generated a secret key, a public key and a relinearization key
// given BFV paramters and the peer id.
func KeyGen(params *bfv.Parameters, peerID uint64) *MKKeys {

	generator := bfv.NewKeyGenerator(params)

	keyBag := new(MKKeys)

	// generate private and public BFV keys

	keyBag.secretKey.key = generator.GenSecretKey()
	keyBag.secretKey.peerID = peerID

	keyBag.publicKey.key = generator.GenPublicKey(keyBag.secretKey.key)
	keyBag.publicKey.peerID = peerID

	// generate relinearization key TODO: create uniENC function

	return keyBag
}

// Symmetric encryption of a single ring element.
func uniEnc(plaintext ring.Poly, sk *MKSecretKey) {

}
