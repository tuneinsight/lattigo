package bfv

import (
	"github.com/ldsec/lattigo/v2/rlwe"
)

// KeyGenerator is an interface wrapping an rlwe.KeyGenerator.
type KeyGenerator interface {
	rlwe.KeyGenerator
}

type keyGenerator struct {
	rlwe.KeyGenerator
}

// NewKeyGenerator creates a new KeyGenerator, from which the secret and public keys, as well as the evaluation,
// rotation and switching keys can be generated.
func NewKeyGenerator(params Parameters) KeyGenerator {
	return &keyGenerator{KeyGenerator: rlwe.NewKeyGenerator(params.Parameters)}
}
