//Package dbfv implements a distributed (or threshold) version of the BFV scheme that enables secure multiparty computation solutions with secret-shared secret keys.
package dbfv

import (
	"github.com/ldsec/lattigo/v2/bfv"
	"github.com/ldsec/lattigo/v2/drlwe"
	"github.com/ldsec/lattigo/v2/ring"
)

// CKGProtocol is the structure storing the parameters and state for a party in the collective key generation protocol.
type CKGProtocol struct {
	drlwe.CKGProtocol
}

// NewCKGProtocol creates a new CKGProtocol instance
func NewCKGProtocol(params *bfv.Parameters) *CKGProtocol {
	ckg := new(CKGProtocol)
	ckg.CKGProtocol = *drlwe.NewCKGProtocol(params.N(), params.Qi(), params.Pi(), params.Sigma())
	return ckg
}

// AllocateShares allocates the CKG shares.
// func (ckg *CKGProtocol) AllocateShares() *drlwe.CKGShare {
// 	return ckg.CKGProtocol.AllocateShares()
// }

// // GenShare generates the party's public key share from its secret key as:
// //
// // crs*s_i + e_i
// //
// // for the receiver protocol. Has no effect is the share was already generated.
// func (ckg *CKGProtocol) GenShare(sk *ring.Poly, crs *ring.Poly, shareOut *drlwe.CKGShare) {
// 	ckg.CKGProtocol.GenShare(sk, crs, shareOut)
// }

// // // AggregateShares aggregates a new share to the aggregate key
// func (ckg *CKGProtocol) AggregateShares(share1, share2, shareOut *drlwe.CKGShare) {
// 	ckg.CKGProtocol.AggregateShares(share1, share2, shareOut)
// }

// GenPublicKey return the current aggregation of the received shares as a bfv.PublicKey.
func (ckg *CKGProtocol) GenPublicKey(roundShare *drlwe.CKGShare, crs *ring.Poly, pubkey *bfv.PublicKey) {
	pubkey.Set([2]*ring.Poly{roundShare.Poly, crs})
}
