package ckks

import "github.com/ldsec/lattigo/v2/rlwe"

// SecretKey is a type for CKKS secret keys.
type SecretKey struct{ rlwe.SecretKey }

// PublicKey is a type for CKKS public keys.
type PublicKey struct{ rlwe.PublicKey }

// SwitchingKey is a type for CKKS public switching keys.
type SwitchingKey struct{ rlwe.SwitchingKey }

// RelinearizationKey is a type for CKKS public relinearization keys.
type RelinearizationKey struct{ rlwe.RelinearizationKey }

// RotationKeySet is a type for storing CKKS public rotation keys.
type RotationKeySet struct{ rlwe.RotationKeySet }

// EvaluationKey is a type composing the relinearization and rotation keys into an evaluation
// key that can be used to initialize bfv.Evaluator types.
type EvaluationKey struct {
	Rlk  *RelinearizationKey
	Rtks *RotationKeySet
}

// BootstrappingKey is a type for a CKKS bootstrapping key, wich regroups the necessary public relinearization
// and rotation keys (i.e., an EvaluationKey).
type BootstrappingKey EvaluationKey

// NewSecretKey returns an allocated CKKS secret key with zero values.
func NewSecretKey(params *Parameters) (sk *SecretKey) {
	return &SecretKey{*rlwe.NewSecretKey(params.N(), params.QPiCount())}
}

// NewPublicKey returns an allocated CKKS public with zero values.
func NewPublicKey(params *Parameters) (pk *PublicKey) {
	return &PublicKey{*rlwe.NewPublicKey(params.N(), params.QPiCount())}
}

// NewSwitchingKey returns an allocated CKKS public switching key with zero values.
func NewSwitchingKey(params *Parameters) *SwitchingKey {
	return &SwitchingKey{*rlwe.NewSwitchingKey(params.N(), params.QPiCount(), params.Beta())}
}

// NewRelinearizationKey returns an allocated CKKS public relinearization key with zero value.
func NewRelinearizationKey(params *Parameters) *RelinearizationKey {
	return &RelinearizationKey{*rlwe.NewRelinKey(2, params.N(), params.QPiCount(), params.Beta())}
}

// NewRotationKeySet return an allocated set of CKKS public relineariation keys with zero values for each galois element
// (i.e., for each supported rotation).
func NewRotationKeySet(params *Parameters, galoisElements []uint64) *RotationKeySet {
	return &RotationKeySet{*rlwe.NewRotationKeySet(galoisElements, params.N(), params.QPiCount(), params.Beta())}
}
