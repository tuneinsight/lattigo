package ckks

import "github.com/ldsec/lattigo/v2/rlwe"

// BootstrappingKey is a type for a CKKS bootstrapping key, wich regroups the necessary public relinearization
// and rotation keys (i.e., an EvaluationKey).
type BootstrappingKey rlwe.EvaluationKey

// NewSecretKey returns an allocated CKKS secret key with zero values.
func NewSecretKey(params *Parameters) (sk *rlwe.SecretKey) {
	return rlwe.NewSecretKey(params.RLWEParameters())
}

// NewPublicKey returns an allocated CKKS public with zero values.
func NewPublicKey(params *Parameters) (pk *rlwe.PublicKey) {
	return rlwe.NewPublicKey(params.RLWEParameters())
}

// NewSwitchingKey returns an allocated CKKS public switching key with zero values.
func NewSwitchingKey(params *Parameters) *rlwe.SwitchingKey {
	return rlwe.NewSwitchingKey(params.RLWEParameters())
}

// NewRelinearizationKey returns an allocated CKKS public relinearization key with zero value.
func NewRelinearizationKey(params *Parameters) *rlwe.RelinearizationKey {
	return rlwe.NewRelinKey(params.RLWEParameters(), 2)
}

// NewRotationKeySet return an allocated set of CKKS public rotation keys with zero values for each galois element
// (i.e., for each supported rotation).
func NewRotationKeySet(params *Parameters, galoisElements []uint64) *rlwe.RotationKeySet {
	return rlwe.NewRotationKeySet(params.RLWEParameters(), galoisElements)
}
