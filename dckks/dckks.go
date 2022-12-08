// Package dckks implements a distributed (or threshold) version of the CKKS scheme that
// enables secure multiparty computation solutions.
// See `drlwe/README.md` for additional information on multiparty schemes.
package dckks

import (
	"github.com/tuneinsight/lattigo/v4/ckks"
	"github.com/tuneinsight/lattigo/v4/drlwe"
	"github.com/tuneinsight/lattigo/v4/ring"
)

// NewCKGProtocol creates a new drlwe.CKGProtocol instance from the CKKS parameters.
// The returned protocol instance is generic and can be used in other multiparty schemes.
func NewCKGProtocol(params ckks.Parameters) *drlwe.CKGProtocol {
	return drlwe.NewCKGProtocol(params.Parameters)
}

// NewRKGProtocol creates a new drlwe.RKGProtocol instance from the CKKS parameters.
// The returned protocol instance is generic and can be used in other multiparty schemes.
func NewRKGProtocol(params ckks.Parameters) *drlwe.RKGProtocol {
	return drlwe.NewRKGProtocol(params.Parameters)
}

// NewGKGProtocol creates a new drlwe.GKGProtocol instance from the CKKS parameters.
// The returned protocol instance is generic and can be used in other multiparty schemes.
func NewGKGProtocol(params ckks.Parameters) *drlwe.GKGProtocol {
	return drlwe.NewGKGProtocol(params.Parameters)
}

// NewCKSProtocol creates a new drlwe.CKSProtocol instance from the CKKS parameters.
// The returned protocol instance is generic and can be used in other multiparty schemes.
func NewCKSProtocol(params ckks.Parameters, noise ring.Distribution) *drlwe.CKSProtocol {
	return drlwe.NewCKSProtocol(params.Parameters, noise)
}

// NewPCKSProtocol creates a new drlwe.PCKSProtocol instance from the CKKS paramters.
// The returned protocol instance is generic and can be used in other multiparty schemes.
func NewPCKSProtocol(params ckks.Parameters, noise ring.Distribution) *drlwe.PCKSProtocol {
	return drlwe.NewPCKSProtocol(params.Parameters, noise)
}
