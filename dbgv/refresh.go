// Package dbgv implements a distributed (or threshold) version of the BGV scheme that enables secure multiparty computation solutions with secret-shared secret keys.
package dbgv

import (
	"github.com/tuneinsight/lattigo/v4/bgv"
	"github.com/tuneinsight/lattigo/v4/drlwe"
	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/rlwe"
)

// RefreshProtocol is a struct storing the relevant parameters for the Refresh protocol.
type RefreshProtocol struct {
	MaskedTransformProtocol
}

// ShallowCopy creates a shallow copy of RefreshProtocol in which all the read-only data-structures are
// shared with the receiver and the temporary buffers are reallocated. The receiver and the returned
// RefreshProtocol can be used concurrently.
func (rfp *RefreshProtocol) ShallowCopy() *RefreshProtocol {
	return &RefreshProtocol{*rfp.MaskedTransformProtocol.ShallowCopy()}
}

// RefreshShare is a struct storing a party's share in the Refresh protocol.
type RefreshShare struct {
	MaskedTransformShare
}

// NewRefreshProtocol creates a new Refresh protocol instance.
func NewRefreshProtocol(params bgv.Parameters, sigmaSmudging float64) (rfp *RefreshProtocol) {
	rfp = new(RefreshProtocol)
	mt, _ := NewMaskedTransformProtocol(params, params, sigmaSmudging)
	rfp.MaskedTransformProtocol = *mt
	return
}

// AllocateShare allocates the shares of the PermuteProtocol
func (rfp *RefreshProtocol) AllocateShare(inputLevel, outputLevel int) *RefreshShare {
	share := rfp.MaskedTransformProtocol.AllocateShare(inputLevel, outputLevel)
	return &RefreshShare{*share}
}

// GenShare generates a share for the Refresh protocol.
// ct1 is degree 1 element of a rlwe.Ciphertext, i.e. rlwe.Ciphertext.Value[1].
func (rfp *RefreshProtocol) GenShare(sk *rlwe.SecretKey, ct1 *ring.Poly, metadata rlwe.MetaData, crp drlwe.CKSCRP, shareOut *RefreshShare) {
	rfp.MaskedTransformProtocol.GenShare(sk, sk, ct1, metadata, crp, nil, &shareOut.MaskedTransformShare)
}

// AggregateShares aggregates two parties' shares in the Refresh protocol.
func (rfp *RefreshProtocol) AggregateShares(share1, share2, shareOut *RefreshShare) {
	rfp.MaskedTransformProtocol.AggregateShares(&share1.MaskedTransformShare, &share2.MaskedTransformShare, &shareOut.MaskedTransformShare)
}

// Finalize applies Decrypt, Recode and Recrypt on the input ciphertext.
func (rfp *RefreshProtocol) Finalize(ctIn *rlwe.Ciphertext, crp drlwe.CKSCRP, share *RefreshShare, ctOut *rlwe.Ciphertext) {
	rfp.MaskedTransformProtocol.Transform(ctIn, nil, crp, &share.MaskedTransformShare, ctOut)
}
