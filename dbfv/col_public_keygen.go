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

	share *ring.Poly
	cpk   [2]*ring.Poly
}

// NewCKGProtocol creates a new ckgProtocolState instance
func NewCKGProtocol(bfvCtx *bfv.BfvContext, crs *ring.Poly) *ckgProtocolState {
	ckg := new(ckgProtocolState)
	ckg.context = bfvCtx.ContextQ()
	ckg.gaussianSampler = bfvCtx.GaussianSampler()
	ckg.cpk[0] = ckg.context.NewPoly()
	ckg.cpk[1] = crs.CopyNew()
	return ckg
}

// GenShare generates the party's public key share from its secret key as:
//
// crs*s_i + e_i
//
// for the receiver protocol. Has no effect is the share was already generated.
func (ckg *ckgProtocolState) GenShare(sk *ring.Poly) {
	if ckg.share == nil {
		// -(sk * crs) + e
		ckg.share = ckg.gaussianSampler.SampleNTTNew()
		ckg.context.MulCoeffsMontgomeryAndSub(sk, ckg.cpk[1], ckg.share)
	}
}

// AggregateShare aggregates a new share to the aggregate key
func (ckg *ckgProtocolState) AggregateShare(share *ring.Poly) {
	if share == nil {
		return
	}
	ckg.context.Add(ckg.cpk[0], share, ckg.cpk[0])
}

// AggregateShare aggregates several shares to the aggregate key
func (ckg *ckgProtocolState) AggregateShares(shares []*ring.Poly) {
	for _, share := range shares {
		ckg.AggregateShare(share)
	}
}

// GetShare returns the party's public share. If not generated using returns nil.
func (ckg *ckgProtocolState) GetShare() *ring.Poly {
	return ckg.share
}

// GetAggregatedShare return the current aggregation of the received shares.
func (ckg *ckgProtocolState) GetAggregatedShare() *ring.Poly {
	return ckg.cpk[0]
}

// GetAggregatedKey return the current aggregation of the received shares as a bfv.PublicKey.
func (ckg *ckgProtocolState) GetAggregatedKey() *bfv.PublicKey {
	collectivePk := new(bfv.PublicKey)
	collectivePk.Set(ckg.cpk)
	return collectivePk
}
