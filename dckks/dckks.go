// Package dckks implements a distributed (or threshold) version of the CKKS scheme that
// enables secure multiparty computation solutions.
// See `drlwe/README.md` for additional information on multiparty schemes.
package dckks

import (
	"github.com/tuneinsight/lattigo/v4/ckks"
	"github.com/tuneinsight/lattigo/v4/drlwe"
	"github.com/tuneinsight/lattigo/v4/ring"
)

// NewPublicKeyGenProtocol creates a new drlwe.PublicKeyGenProtocol instance from the CKKS parameters.
// The returned protocol instance is generic and can be used in other multiparty schemes.
func NewPublicKeyGenProtocol(params ckks.Parameters) drlwe.PublicKeyGenProtocol {
	return drlwe.NewPublicKeyGenProtocol(params.Parameters)
}

// NewRelinearizationKeyGenProtocol creates a new drlwe.RelinearizationKeyGenProtocol instance from the CKKS parameters.
// The returned protocol instance is generic and can be used in other multiparty schemes.
func NewRelinearizationKeyGenProtocol(params ckks.Parameters) drlwe.RelinearizationKeyGenProtocol {
	return drlwe.NewRelinearizationKeyGenProtocol(params.Parameters)
}

// NewGaloisKeyGenProtocol creates a new drlwe.GaloisKeyGenProtocol instance from the CKKS parameters.
// The returned protocol instance is generic and can be used in other multiparty schemes.
func NewGaloisKeyGenProtocol(params ckks.Parameters) drlwe.GaloisKeyGenProtocol {
	return drlwe.NewGaloisKeyGenProtocol(params.Parameters)
}

// NewKeySwitchProtocol creates a new drlwe.KeySwitchProtocol instance from the CKKS parameters.
// The returned protocol instance is generic and can be used in other multiparty schemes.
func NewKeySwitchProtocol(params ckks.Parameters, noise ring.DistributionParameters) (drlwe.KeySwitchProtocol, error) {
	return drlwe.NewKeySwitchProtocol(params.Parameters, noise)
}

// NewPublicKeySwitchProtocol creates a new drlwe.PublicKeySwitchProtocol instance from the CKKS paramters.
// The returned protocol instance is generic and can be used in other multiparty schemes.
func NewPublicKeySwitchProtocol(params ckks.Parameters, noise ring.DistributionParameters) (drlwe.PublicKeySwitchProtocol, error) {
	return drlwe.NewPublicKeySwitchProtocol(params.Parameters, noise)
}
