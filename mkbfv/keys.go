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

// MKPublicKey is a type for BFV public keys and ID in a multi key context.
type MKPublicKey struct {
	key    *bfv.PublicKey
	peerID uint64
}

// MKRelinearizationKey is a type for BFV relinearization keys in a multy key context.
type MKRelinearizationKey struct {
	key    [][3]*ring.Poly
	peerID uint64
}

// MKEvaluationKey is a type for BFV evaluation keys in a multy key context.
type MKEvaluationKey struct {
	//TODO: ask if same as in bfv or follow the paper's implementation (Hao Chen)
	peerID uint64
}

// MKRotationKey is a type for BFV evaluation keys in a multy key context.
type MKRotationKey struct {
	//TODO: ask if same as in bfv or follow the paper's implementation (Hao Chen)
	peerID uint64
}

// MKKeys is a type that contains all keys necessary for the multy key protocol.
type MKKeys struct {
	secretKey MKSecretKey
	publicKey MKPublicKey
	relinKey  MKRelinearizationKey
	evalKey   MKEvaluationKey
	rotKey    MKRotationKey
}
