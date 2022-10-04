package bgv

import "github.com/tuneinsight/lattigo/v4/rlwe"

// NewKeyGenerator creates a rlwe.KeyGenerator instance from the BGV parameters.
func NewKeyGenerator(params Parameters) rlwe.KeyGenerator {
	return rlwe.NewKeyGenerator(params.Parameters)
}

// NewSecretKey returns an allocated BGV secret key with zero values.
func NewSecretKey(params Parameters) (sk *rlwe.SecretKey) {
	return rlwe.NewSecretKey(params.Parameters)
}

// NewPublicKey returns an allocated BGV public key with zero values.
func NewPublicKey(params Parameters) (pk *rlwe.PublicKey) {
	return rlwe.NewPublicKey(params.Parameters)
}

// NewSwitchingKey returns an allocated BGV public switching key with zero values.
func NewSwitchingKey(params Parameters) *rlwe.SwitchingKey {
	return rlwe.NewSwitchingKey(params.Parameters, params.QCount()-1, params.PCount()-1)
}

// NewRelinearizationKey returns an allocated BGV public relinearization key with zero values for each degree in [2 < maxRelinDegree].
func NewRelinearizationKey(params Parameters, maxRelinDegree int) *rlwe.RelinearizationKey {
	return rlwe.NewRelinKey(params.Parameters, maxRelinDegree)
}

// NewRotationKeySet returns an allocated set of BGV public rotation keys with zero values for each Galois element
// (i.e., for each supported rotation).
func NewRotationKeySet(params Parameters, galoisElements []uint64) *rlwe.RotationKeySet {
	return rlwe.NewRotationKeySet(params.Parameters, galoisElements)
}
