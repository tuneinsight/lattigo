// Package dbgv implements a distributed (or threshold) version of the BGV scheme that
// enables secure multiparty computation solutions.
// See `drlwe/README.md` for additional information on multiparty schemes.
package dbgv

import (
	"github.com/tuneinsight/lattigo/v4/bgv"
	"github.com/tuneinsight/lattigo/v4/drlwe"
)

// NewCKGProtocol creates a new drlwe.CKGProtocol instance from the BGV parameters.
// The returned protocol instance is generic and can be used in other multiparty schemes.
func NewCKGProtocol(params bgv.Parameters) *drlwe.CKGProtocol {
	return drlwe.NewCKGProtocol(params.Parameters)
}

// NewRKGProtocol creates a new drlwe.RKGProtocol instance from the BGV parameters.
// The returned protocol instance is generic and can be used in other multiparty schemes.
func NewRKGProtocol(params bgv.Parameters) *drlwe.RKGProtocol {
	return drlwe.NewRKGProtocol(params.Parameters)
}

// NewRTGProtocol creates a new drlwe.RTGProtocol instance from the BGV parameters.
// The returned protocol instance is generic and can be used in other multiparty schemes.
func NewRTGProtocol(params bgv.Parameters) *drlwe.RTGProtocol {
	return drlwe.NewRTGProtocol(params.Parameters)
}

// NewCKSProtocol creates a new drlwe.CKSProtocol instance from the BGV parameters.
// The returned protocol instance is generic and can be used in other multiparty schemes.
func NewCKSProtocol(params bgv.Parameters, sigmaSmudging float64) *drlwe.CKSProtocol {
	return drlwe.NewCKSProtocol(params.Parameters, sigmaSmudging)
}

// NewPCKSProtocol creates a new drlwe.PCKSProtocol instance from the BGV paramters.
// The returned protocol instance is generic and can be used in other multiparty schemes.
func NewPCKSProtocol(params bgv.Parameters, sigmaSmudging float64) *drlwe.PCKSProtocol {
	return drlwe.NewPCKSProtocol(params.Parameters, sigmaSmudging)
}
