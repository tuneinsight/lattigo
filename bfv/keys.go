package bfv

import "github.com/ldsec/lattigo/v2/rlwe"

// SecretKey is a type for BFV secret keys.
type SecretKey struct{ rlwe.SecretKey }

// PublicKey is a type for BFV public keys.
type PublicKey struct{ rlwe.PublicKey }

// SwitchingKey is a type for BFV public switching keys.
type SwitchingKey struct{ rlwe.SwitchingKey }

// RelinearizationKey is a type for BFV public relinearization keys.
type RelinearizationKey struct{ rlwe.RelinearizationKey }

// RotationKeySet is a type for storing BFV public rotation keys.
type RotationKeySet struct{ rlwe.RotationKeySet }

// EvaluationKey is a type composing the relinearization and rotation keys into an evaluation
// key that can be used to initialize bfv.Evaluator types.
type EvaluationKey struct {
	Rlk  *RelinearizationKey
	Rtks *RotationKeySet
}

// NewSecretKey returns an allocated BFV secret key with zero values.
func NewSecretKey(params *Parameters) (sk *SecretKey) {
	return &SecretKey{*rlwe.NewSecretKey(params.N(), params.QPiCount())}
}

// NewPublicKey returns an allocated BFV public with zero values.
func NewPublicKey(params *Parameters) (pk *PublicKey) {
	return &PublicKey{*rlwe.NewPublicKey(params.N(), params.QPiCount())}
}

// NewSwitchingKey returns an allocated BFV public switching key with zero values.
func NewSwitchingKey(params *Parameters) *SwitchingKey {
	return &SwitchingKey{*rlwe.NewSwitchingKey(params.N(), params.QPiCount(), params.Beta())}
}

// NewRelinearizationKey returns an allocated BFV public relinearization key with zero value for each degree in [2 < maxRelinDegree].
func NewRelinearizationKey(params *Parameters, maxRelinDegree int) *RelinearizationKey {
	return &RelinearizationKey{*rlwe.NewRelinKey(maxRelinDegree, params.N(), params.QPiCount(), params.Beta())}
}

// NewRotationKeySet return an allocated set of BFV public relineariation keys with zero values for each galois element
// (i.e., for each supported rotation).
func NewRotationKeySet(params *Parameters, galoisElements []uint64) *RotationKeySet {
	return &RotationKeySet{*rlwe.NewRotationKeySet(galoisElements, params.N(), params.QPiCount(), params.Beta())}
}
