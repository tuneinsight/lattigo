// Package dbfv implements a distributed (or threshold) version of the BFV scheme that
// enables secure multiparty computation solutions.
// See `drlwe/README.md` for additional information on multiparty schemes.
package dbfv

import (
	"github.com/tuneinsight/lattigo/v4/bfv"
	"github.com/tuneinsight/lattigo/v4/drlwe"
)

// NewCKGProtocol creates a new drlwe.CKGProtocol instance from the BFV parameters.
// The returned protocol instance is generic and can be used in other multiparty schemes.
func NewCKGProtocol(params bfv.Parameters) *drlwe.CKGProtocol {
	return drlwe.NewCKGProtocol(params.Parameters)
}

// NewRKGProtocol creates a new drlwe.RKGProtocol instance from the BFV parameters.
// The returned protocol instance is generic and can be used in other multiparty schemes.
func NewRKGProtocol(params bfv.Parameters) *drlwe.RKGProtocol {
	return drlwe.NewRKGProtocol(params.Parameters)
}

// NewGKGProtocol creates a new drlwe.GKGProtocol instance from the BFV parameters.
// The returned protocol instance is generic and can be used in other multiparty schemes.
func NewGKGProtocol(params bfv.Parameters) *drlwe.GKGProtocol {
	return drlwe.NewGKGProtocol(params.Parameters)
}

// NewCKSProtocol creates a new drlwe.CKSProtocol instance from the BFV parameters.
// The returned protocol instance is generic and can be used in other multiparty schemes.
func NewCKSProtocol(params bfv.Parameters, sigmaSmudging float64) *drlwe.CKSProtocol {
	return drlwe.NewCKSProtocol(params.Parameters, sigmaSmudging)
}

// NewPCKSProtocol creates a new drlwe.PCKSProtocol instance from the BFV paramters.
// The returned protocol instance is generic and can be used in other multiparty schemes.
func NewPCKSProtocol(params bfv.Parameters, sigmaSmudging float64) *drlwe.PCKSProtocol {
	return drlwe.NewPCKSProtocol(params.Parameters, sigmaSmudging)
}
