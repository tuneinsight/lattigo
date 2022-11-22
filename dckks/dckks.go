// Package dckks implements a distributed (or threshold) version of the CKKS scheme that
// enables secure multiparty computation solutions.
// See `drlwe/README.md` for additional information on multiparty schemes.
package dckks

import (
	"github.com/tuneinsight/lattigo/v4/ckks"
	"github.com/tuneinsight/lattigo/v4/drlwe"
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

// NewRTGProtocol creates a new drlwe.RTGProtocol instance from the CKKS parameters.
// The returned protocol instance is generic and can be used in other multiparty schemes.
func NewRTGProtocol(params ckks.Parameters) *drlwe.RTGProtocol {
	return drlwe.NewRTGProtocol(params.Parameters)
}

// NewCKSProtocol creates a new drlwe.CKSProtocol instance from the CKKS parameters.
// The returned protocol instance is generic and can be used in other multiparty schemes.
func NewCKSProtocol(params ckks.Parameters, sigmaSmudging float64) *drlwe.CKSProtocol {
	return drlwe.NewCKSProtocol(params.Parameters, sigmaSmudging)
}

// NewPCKSProtocol creates a new drlwe.PCKSProtocol instance from the CKKS paramters.
// The returned protocol instance is generic and can be used in other multiparty schemes.
func NewPCKSProtocol(params ckks.Parameters, sigmaSmudging float64) *drlwe.PCKSProtocol {
	return drlwe.NewPCKSProtocol(params.Parameters, sigmaSmudging)
}
