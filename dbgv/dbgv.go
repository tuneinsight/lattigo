// Package dbgv implements a distributed (or threshold) version of the
// unified RNS-accelerated version of the Fan-Vercauteren version of
// Brakerski's scale invariant homomorphic encryption scheme (BFV)
// and Brakerski-Gentry-Vaikuntanathan (BGV) homomorphic encryption scheme.
// It provides modular arithmetic over the integers and enables secure
// multiparty computation solutions.
// See `drlwe/README.md` for additional information on multiparty schemes.
package dbgv

import (
	"github.com/tuneinsight/lattigo/v4/bgv"
	"github.com/tuneinsight/lattigo/v4/drlwe"
	"github.com/tuneinsight/lattigo/v4/ring"
)

// NewPublicKeyGenProtocol creates a new drlwe.PublicKeyGenProtocol instance from the BGV parameters.
// The returned protocol instance is generic and can be used in other multiparty schemes.
func NewPublicKeyGenProtocol(params bgv.Parameters) drlwe.PublicKeyGenProtocol {
	return drlwe.NewPublicKeyGenProtocol(params.Parameters)
}

// NewRelinearizationKeyGenProtocol creates a new drlwe.RKGProtocol instance from the BGV parameters.
// The returned protocol instance is generic and can be used in other multiparty schemes.
func NewRelinearizationKeyGenProtocol(params bgv.Parameters) drlwe.RelinearizationKeyGenProtocol {
	return drlwe.NewRelinearizationKeyGenProtocol(params.Parameters)
}

// NewGaloisKeyGenProtocol creates a new drlwe.GaloisKeyGenProtocol instance from the BGV parameters.
// The returned protocol instance is generic and can be used in other multiparty schemes.
func NewGaloisKeyGenProtocol(params bgv.Parameters) drlwe.GaloisKeyGenProtocol {
	return drlwe.NewGaloisKeyGenProtocol(params.Parameters)
}

// NewKeySwitchProtocol creates a new drlwe.KeySwitchProtocol instance from the BGV parameters.
// The returned protocol instance is generic and can be used in other multiparty schemes.
func NewKeySwitchProtocol(params bgv.Parameters, noiseFlooding ring.DistributionParameters) (drlwe.KeySwitchProtocol, error) {
	return drlwe.NewKeySwitchProtocol(params.Parameters, noiseFlooding)
}

// NewPublicKeySwitchProtocol creates a new drlwe.PublicKeySwitchProtocol instance from the BGV paramters.
// The returned protocol instance is generic and can be used in other multiparty schemes.
func NewPublicKeySwitchProtocol(params bgv.Parameters, noiseFlooding ring.DistributionParameters) (drlwe.PublicKeySwitchProtocol, error) {
	return drlwe.NewPublicKeySwitchProtocol(params.Parameters, noiseFlooding)
}
