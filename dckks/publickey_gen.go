//Package dckks implements a distributed (or threshold) version of the CKKS scheme that enables secure multiparty computation solutions with secret-shared secret keys.
package dckks

import (
	"github.com/ldsec/lattigo/ckks"
	"github.com/ldsec/lattigo/ring"
)

// CKGProtocol is the structure storing the parameters and state for a party in the collective key generation protocol.
type CKGProtocol struct {
	ringContext     *ring.Context
	gaussianSampler *ring.KYSampler
}

// CKGShare is a struct storing the CKG protocol's share.
type CKGShare *ring.Poly

// NewCKGProtocol creates a new CKGProtocol instance
func NewCKGProtocol(ckksCtx *ckks.Context) *CKGProtocol {
	ckg := new(CKGProtocol)
	ckg.ringContext = ckksCtx.ContextKeys()
	ckg.gaussianSampler = ckksCtx.GaussianSampler()
	return ckg
}

// AllocateShares allocates the share of the CKG protocol.
func (ckg *CKGProtocol) AllocateShares() CKGShare {
	return ckg.ringContext.NewPoly()
}

// GenShare generates the party's public key share from its secret key as:
//
// crs*s_i + e_i
//
// for the receiver protocol. Has no effect is the share was already generated.
func (ckg *CKGProtocol) GenShare(sk *ring.Poly, crs *ring.Poly, shareOut CKGShare) {
	ckg.gaussianSampler.SampleNTT(shareOut)
	ckg.ringContext.MulCoeffsMontgomeryAndSub(sk, crs, shareOut)
}

// AggregateShares aggregates a new share to the aggregate key
func (ckg *CKGProtocol) AggregateShares(share1, share2, shareOut CKGShare) {
	ckg.ringContext.Add(share1, share2, shareOut)
}

// GenPublicKey return the current aggregation of the received shares as a bfv.PublicKey.
func (ckg *CKGProtocol) GenPublicKey(roundShare CKGShare, crs *ring.Poly, pubkey *ckks.PublicKey) {
	pubkey.Set([2]*ring.Poly{roundShare, crs})
}
