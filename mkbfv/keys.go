package mkbfv

import (
	"github.com/ldsec/lattigo/v2/bfv"
	"github.com/ldsec/lattigo/v2/ring"
)

// MKSecretKey is a type for BFV secret keys in a multi key context.
type MKSecretKey struct {
	key    *bfv.SecretKey
	peerID uint64
}

// MKPublicKey is a type for BFV public keys and ID in a multi key context. key[1] = a and key[0] = -s * a + e mod q
type MKPublicKey struct {
	key    [2]*MKDecomposedPoly
	peerID uint64
}

// MKDecomposedPoly is a type for vectors decomposed in a basis w (belong to Rq^d)(gadget decomposition)
type MKDecomposedPoly struct {
	poly []*ring.Poly
}

// MKEvaluationKey is a type for BFV evaluation keys in a multy key context.
type MKEvaluationKey struct {
	key    [3]*MKDecomposedPoly
	peerID uint64
}

// MKSwitchingKey is a type for BFV switching keys in a multy key context.
type MKSwitchingKey struct {
	key [3]*MKDecomposedPoly
	//peerID uint64 // Commented because in relinkey_gen.Convert we might not need a peerID, or might need multiple
}

// MKRelinearizationKey is a type for BFV relinearization keys in a multy key context.
type MKRelinearizationKey struct {
	key [][]*MKSwitchingKey
}

// MKKeys is a type that contains all keys necessary for the multy key protocol.
type MKKeys struct {
	secretKey MKSecretKey
	publicKey MKPublicKey
	evalKey   MKEvaluationKey
	relinKey  MKRelinearizationKey
}
