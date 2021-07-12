package ckks

import "github.com/ldsec/lattigo/v2/rlwe"

func NewKeyGenerator(params Parameters) rlwe.KeyGenerator {
	return rlwe.NewKeyGenerator(params.Parameters)
}
