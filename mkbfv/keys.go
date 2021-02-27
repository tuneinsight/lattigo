package mkbfv

import (
	"github.com/ldsec/lattigo/v2/bfv"
)

// MKSecretKey is a type for BFV secret keys in a multi key context.
type MKSecretKey struct {
	bfv.SecretKey
	peerID uint64
}

// MKPublicKey is a type for BFV public keys and ID in a multi key context.
type MKPublicKey struct {
	key    bfv.PublicKey
	peerID uint64
}

// MKRelinearizationKey is a type for BFV relinearization keys in a multy key context.
type MKRelinearizationKey struct {
	//TODO: ask if relin key same as bfv or follow the paper's implementation (Hao Chen)
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
