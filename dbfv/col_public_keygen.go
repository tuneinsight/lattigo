//Package dbfv implements a distributed (or threshold) version of the BFV scheme that enables secure multiparty computation solutions with secret-shared secret keys.
package dbfv

import (
	"github.com/ldsec/lattigo/bfv"
	"github.com/ldsec/lattigo/ring"
)

// ckgProtocolState is the structure storing the parameters and state for a party in the collective key generation protocol.
type ckgProtocolState struct {
	context         *ring.Context
	gaussianSampler *ring.KYSampler
}

type ckgShare *ring.Poly

// NewCKGProtocol creates a new ckgProtocolState instance
func NewCKGProtocol(bfvCtx *bfv.BfvContext) *ckgProtocolState {
	ckg := new(ckgProtocolState)
	ckg.context = bfvCtx.ContextQ()
	ckg.gaussianSampler = bfvCtx.GaussianSampler()
	return ckg
}

func (ckg *ckgProtocolState) AllocateShares() ckgShare {
	return ckg.context.NewPoly()
}

// GenShare generates the party's public key share from its secret key as:
//
// crs*s_i + e_i
//
// for the receiver protocol. Has no effect is the share was already generated.
func (ckg *ckgProtocolState) GenShare(sk *ring.Poly,  crs *ring.Poly, shareOut ckgShare) {
	ckg.gaussianSampler.SampleNTT(shareOut)
	ckg.context.MulCoeffsMontgomeryAndSub(sk, crs, shareOut)
}

// AggregateShare aggregates a new share to the aggregate key
func (ckg *ckgProtocolState) AggregateShare(share1, share2, shareOut ckgShare) {
	ckg.context.Add(share1, share2, shareOut)
}

// GetAggregatedKey return the current aggregation of the received shares as a bfv.PublicKey.
func (ckg *ckgProtocolState) GetAggregatedKey(roundShare ckgShare, crs *ring.Poly, pubkey  *bfv.PublicKey) {
	pubkey.Set([2]*ring.Poly{roundShare, crs})
}
