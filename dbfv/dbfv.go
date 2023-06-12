// Package dbfv implements a distributed (or threshold) version of the BFV scheme that
// enables secure multiparty computation solutions.
// See `drlwe/README.md` for additional information on multiparty schemes.
package dbfv

import (
	"github.com/tuneinsight/lattigo/v4/bfv"
	"github.com/tuneinsight/lattigo/v4/dbgv"
	"github.com/tuneinsight/lattigo/v4/drlwe"
	"github.com/tuneinsight/lattigo/v4/ring"
)

// NewPublicKeyGenProtocol creates a new drlwe.PublicKeyGenProtocol instance from the BFV parameters.
// The returned protocol instance is generic and can be used in other multiparty schemes.
func NewPublicKeyGenProtocol(params bfv.Parameters) *drlwe.PublicKeyGenProtocol {
	return drlwe.NewPublicKeyGenProtocol(params.Parameters.Parameters)
}

// NewRelinKeyGenProtocol creates a new drlwe.RelinKeyGenProtocol instance from the BFV parameters.
// The returned protocol instance is generic and can be used in other multiparty schemes.
func NewRelinKeyGenProtocol(params bfv.Parameters) *drlwe.RelinKeyGenProtocol {
	return drlwe.NewRelinKeyGenProtocol(params.Parameters.Parameters)
}

// NewGaloisKeyGenProtocol creates a new drlwe.RelinKeyGenProtocol instance from the BFV parameters.
// The returned protocol instance is generic and can be used in other multiparty schemes.
func NewGaloisKeyGenProtocol(params bfv.Parameters) *drlwe.GaloisKeyGenProtocol {
	return drlwe.NewGaloisKeyGenProtocol(params.Parameters.Parameters)
}

// NewKeySwitchProtocol creates a new drlwe.KeySwitchProtocol instance from the BFV parameters.
// The returned protocol instance is generic and can be used in other multiparty schemes.
func NewKeySwitchProtocol(params bfv.Parameters, noiseFlooding ring.DistributionParameters) *drlwe.KeySwitchProtocol {
	return drlwe.NewKeySwitchProtocol(params.Parameters.Parameters, noiseFlooding)
}

// NewPublicKeySwitchProtocol creates a new drlwe.PublicKeySwitchProtocol instance from the BFV paramters.
// The returned protocol instance is generic and can be used in other multiparty schemes.
func NewPublicKeySwitchProtocol(params bfv.Parameters, noiseFlooding ring.DistributionParameters) *drlwe.PublicKeySwitchProtocol {
	return drlwe.NewPublicKeySwitchProtocol(params.Parameters.Parameters, noiseFlooding)
}

// NewRefreshProtocol creates a new instance of the RefreshProtocol.
func NewRefreshProtocol(params bfv.Parameters, noiseFlooding ring.DistributionParameters) (rft *dbgv.RefreshProtocol) {
	return dbgv.NewRefreshProtocol(params.Parameters, noiseFlooding)
}

// NewEncToShareProtocol creates a new instance of the EncToShareProtocol.
func NewEncToShareProtocol(params bfv.Parameters, noiseFlooding ring.DistributionParameters) (e2s *dbgv.EncToShareProtocol) {
	return dbgv.NewEncToShareProtocol(params.Parameters, noiseFlooding)
}

// NewShareToEncProtocol creates a new instance of the ShareToEncProtocol.
func NewShareToEncProtocol(params bfv.Parameters, noiseFlooding ring.DistributionParameters) (e2s *dbgv.ShareToEncProtocol) {
	return dbgv.NewShareToEncProtocol(params.Parameters, noiseFlooding)
}

// NewMaskedTransformProtocol creates a new instance of the MaskedTransformProtocol.
func NewMaskedTransformProtocol(paramsIn, paramsOut bfv.Parameters, noiseFlooding ring.DistributionParameters) (rfp *dbgv.MaskedTransformProtocol, err error) {
	return dbgv.NewMaskedTransformProtocol(paramsIn.Parameters, paramsOut.Parameters, noiseFlooding)
}
