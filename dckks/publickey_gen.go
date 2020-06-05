//Package dckks implements a distributed (or threshold) version of the CKKS scheme that enables secure multiparty computation solutions with secret-shared secret keys.
package dckks

import (
	"github.com/ldsec/lattigo/ckks"
	"github.com/ldsec/lattigo/ring"
	"github.com/ldsec/lattigo/utils"
)

// CKGProtocol is the structure storing the parameters and state for a party in the collective key generation protocol.
type CKGProtocol struct {
	dckksContext *dckksContext
}

// CKGShare is a struct storing the CKG protocol's share.
type CKGShare *ring.Poly

// NewCKGProtocol creates a new CKGProtocol instance
func NewCKGProtocol(params *ckks.Parameters) *CKGProtocol {

	if !params.IsValid() {
		panic("cannot NewCKGProtocol : params not valid (check if they where generated properly)")
	}

	ckg := new(CKGProtocol)
	ckg.dckksContext = newDckksContext(params)
	return ckg
}

// AllocateShares allocates the share of the CKG protocol.
func (ckg *CKGProtocol) AllocateShares() CKGShare {
	return ckg.dckksContext.contextQP.NewPoly()
}

// GenShare generates the party's public key share from its secret key as:
//
// crs*s_i + e_i
//
// for the receiver protocol. Has no effect is the share was already generated.
func (ckg *CKGProtocol) GenShare(sk *ring.Poly, crs *ring.Poly, shareOut CKGShare) {
	contextQP := ckg.dckksContext.contextQP
	prng, err := utils.NewPRNG()
	if err != nil {
		panic(err)
	}
	gaussianSampler := ring.NewGaussianSampler(prng, contextQP)

	gaussianSampler.SampleNTTLvl(uint64(len(contextQP.Modulus)-1), shareOut, ckg.dckksContext.params.Sigma, uint64(6*ckg.dckksContext.params.Sigma))
	contextQP.MulCoeffsMontgomeryAndSub(sk, crs, shareOut)
}

// AggregateShares aggregates a new share to the aggregate key
func (ckg *CKGProtocol) AggregateShares(share1, share2, shareOut CKGShare) {
	ckg.dckksContext.contextQP.Add(share1, share2, shareOut)
}

// GenPublicKey return the current aggregation of the received shares as a bfv.PublicKey.
func (ckg *CKGProtocol) GenPublicKey(roundShare CKGShare, crs *ring.Poly, pubkey *ckks.PublicKey) {
	pubkey.Set([2]*ring.Poly{roundShare, crs})
}
