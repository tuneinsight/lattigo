package bfv

import "github.com/ldsec/lattigo/v2/rlwe"

// NewKeyGenerator creates a rlwe.KeyGenerator instance from the BFV parameters.
func NewKeyGenerator(params Parameters) rlwe.KeyGenerator {
	return rlwe.NewKeyGenerator(params.Parameters)
}
