// Package dbfv implements a distributed (or threshold) version of the BFV scheme that
// enables secure multiparty computation solutions.
// See `drlwe/README.md` for additional information on multiparty schemes.
package dbfv

import (
	"github.com/tuneinsight/lattigo/v4/bfv"
	"github.com/tuneinsight/lattigo/v4/dbgv"
	"github.com/tuneinsight/lattigo/v4/drlwe"
	"github.com/tuneinsight/lattigo/v4/ring/distribution"
)

// NewCKGProtocol creates a new drlwe.CKGProtocol instance from the BFV parameters.
// The returned protocol instance is generic and can be used in other multiparty schemes.
func NewCKGProtocol(params bfv.Parameters) *drlwe.CKGProtocol {
	return drlwe.NewCKGProtocol(params.Parameters.Parameters)
}

// NewRKGProtocol creates a new drlwe.RKGProtocol instance from the BFV parameters.
// The returned protocol instance is generic and can be used in other multiparty schemes.
func NewRKGProtocol(params bfv.Parameters) *drlwe.RKGProtocol {
	return drlwe.NewRKGProtocol(params.Parameters.Parameters)
}

// NewGKGProtocol creates a new drlwe.GKGProtocol instance from the BFV parameters.
// The returned protocol instance is generic and can be used in other multiparty schemes.
func NewGKGProtocol(params bfv.Parameters) *drlwe.GKGProtocol {
	return drlwe.NewGKGProtocol(params.Parameters.Parameters)
}

// NewCKSProtocol creates a new drlwe.CKSProtocol instance from the BFV parameters.
// The returned protocol instance is generic and can be used in other multiparty schemes.
func NewCKSProtocol(params bfv.Parameters, noise distribution.Distribution) *drlwe.CKSProtocol {
	return drlwe.NewCKSProtocol(params.Parameters.Parameters, noise)
}

// NewPCKSProtocol creates a new drlwe.PCKSProtocol instance from the BFV paramters.
// The returned protocol instance is generic and can be used in other multiparty schemes.
func NewPCKSProtocol(params bfv.Parameters, noise distribution.Distribution) *drlwe.PCKSProtocol {
	return drlwe.NewPCKSProtocol(params.Parameters.Parameters, noise)
}

// NewRefreshProtocol creates a new instance of the RefreshProtocol.
func NewRefreshProtocol(params bfv.Parameters, noise distribution.Distribution) (rft *dbgv.RefreshProtocol) {
	return dbgv.NewRefreshProtocol(params.Parameters, noise)
}

// NewE2SProtocol creates a new instance of the E2SProtocol.
func NewE2SProtocol(params bfv.Parameters, noise distribution.Distribution) (e2s *dbgv.E2SProtocol) {
	return dbgv.NewE2SProtocol(params.Parameters, noise)
}

// NewS2EProtocol creates a new instance of the S2EProtocol.
func NewS2EProtocol(params bfv.Parameters, noise distribution.Distribution) (e2s *dbgv.S2EProtocol) {
	return dbgv.NewS2EProtocol(params.Parameters, noise)
}

// NewMaskedTransformProtocol creates a new instance of the MaskedTransformProtocol.
func NewMaskedTransformProtocol(paramsIn, paramsOut bfv.Parameters, noise distribution.Distribution) (rfp *dbgv.MaskedTransformProtocol, err error) {
	return dbgv.NewMaskedTransformProtocol(paramsIn.Parameters, paramsOut.Parameters, noise)
}
