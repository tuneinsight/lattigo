package ckks

import "github.com/tuneinsight/lattigo/v3/rlwe"

// KeyGenerator is an interface for the generation of CKKS keys.
type KeyGenerator interface {
	rlwe.KeyGenerator
	GenSwitchingKeysForBridge(skCKKS, skCI *rlwe.SecretKey) (*SwkComplexToReal, *SwkRealToComplex)
}

// SwkComplexToReal is a SwitchingKey to switch from the standard domain to the conjugate invariant domain.
type SwkComplexToReal struct {
	rlwe.SwitchingKey
}

// SwkRealToComplex is a SwitchingKey to switch from the conjugate invariant domain to the standard domain.
type SwkRealToComplex struct {
	rlwe.SwitchingKey
}

type keyGenerator struct {
	rlwe.KeyGenerator

	params *Parameters
}

// GenSwitchingKeysForBridge generates the necessary switching keys to switch from the standard domain to the conjugate invariant domain
// and vice-versa.
func (keygen *keyGenerator) GenSwitchingKeysForBridge(skStd, skConjugateInvariant *rlwe.SecretKey) (*SwkComplexToReal, *SwkRealToComplex) {
	swkStdToCi, swkCitoStd := keygen.GenSwitchingKeysForRingSwap(skStd, skConjugateInvariant)
	return &SwkComplexToReal{*swkStdToCi}, &SwkRealToComplex{*swkCitoStd}
}

// NewKeyGenerator creates a rlwe.KeyGenerator instance from the CKKS parameters.
func NewKeyGenerator(params Parameters) KeyGenerator {
	return &keyGenerator{rlwe.NewKeyGenerator(params.Parameters), &params}
}

// NewSecretKey returns an allocated CKKS secret key with zero values.
func NewSecretKey(params Parameters) (sk *rlwe.SecretKey) {
	return rlwe.NewSecretKey(params.Parameters)
}

// NewPublicKey returns an allocated CKKS public with zero values.
func NewPublicKey(params Parameters) (pk *rlwe.PublicKey) {
	return rlwe.NewPublicKey(params.Parameters)
}

// NewSwitchingKey returns an allocated CKKS public switching key with zero values.
func NewSwitchingKey(params Parameters) *rlwe.SwitchingKey {
	return rlwe.NewSwitchingKey(params.Parameters, params.QCount()-1, params.PCount()-1)
}

// NewRelinearizationKey returns an allocated CKKS public relinearization key with zero value.
func NewRelinearizationKey(params Parameters) *rlwe.RelinearizationKey {
	return rlwe.NewRelinKey(params.Parameters, 2)
}

// NewRotationKeySet returns an allocated set of CKKS public rotation keys with zero values for each galois element
// (i.e., for each supported rotation).
func NewRotationKeySet(params Parameters, galoisElements []uint64) *rlwe.RotationKeySet {
	return rlwe.NewRotationKeySet(params.Parameters, galoisElements)
}
