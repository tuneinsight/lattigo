//Package dckks implements a distributed (or threshold) version of the CKKS scheme that enables secure multiparty computation solutions with secret-shared secret keys.
package dckks

import (
	"github.com/ldsec/lattigo/v2/ckks"
	"github.com/ldsec/lattigo/v2/drlwe"
	"github.com/ldsec/lattigo/v2/ring"
)

// CKGProtocol is the structure storing the parameters and state for a party in the collective key generation protocol.
type CKGProtocol struct {
	drlwe.CKGProtocol
}

// NewCKGProtocol creates a new CKGProtocol instance
func NewCKGProtocol(params *ckks.Parameters) *CKGProtocol {

	ckg := new(CKGProtocol)
	ckg.CKGProtocol = *drlwe.NewCKGProtocol(params.N(), params.Qi(), params.Pi(), params.Sigma())
	return ckg
}

// GenCKKSPublicKey return the current aggregation of the received shares as a ckks.PublicKey.
func (ckg *CKGProtocol) GenCKKSPublicKey(roundShare *drlwe.CKGShare, crs *ring.Poly, pubkey *ckks.PublicKey) {
	ckg.CKGProtocol.GenPublicKey(roundShare, crs, &pubkey.PublicKey)
}
