package mpbgv

import (
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/multiparty"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/schemes/bgv"
)

// RefreshProtocol is a struct storing the relevant parameters for the Refresh protocol.
type RefreshProtocol struct {
	MaskedTransformProtocol
}

// ShallowCopy creates a shallow copy of [RefreshProtocol] in which all the read-only data-structures are
// shared with the receiver and the temporary buffers are reallocated. The receiver and the returned
// [RefreshProtocol] can be used concurrently.
func (rfp *RefreshProtocol) ShallowCopy() RefreshProtocol {
	return RefreshProtocol{rfp.MaskedTransformProtocol.ShallowCopy()}
}

// NewRefreshProtocol creates a new Refresh protocol instance.
func NewRefreshProtocol(params bgv.Parameters, noiseFlooding ring.DistributionParameters) (rfp RefreshProtocol, err error) {
	rfp = RefreshProtocol{}
	mt, err := NewMaskedTransformProtocol(params, params, noiseFlooding)
	rfp.MaskedTransformProtocol = mt
	return rfp, err
}

// AllocateShare allocates the shares of the PermuteProtocol
func (rfp RefreshProtocol) AllocateShare(inputLevel, outputLevel int) multiparty.RefreshShare {
	return rfp.MaskedTransformProtocol.AllocateShare(inputLevel, outputLevel)
}

// GenShare generates a share for the Refresh protocol.
// ct1 is degree 1 element of a rlwe.Ciphertext, i.e. rlwe.Ciphertext.Value[1].
func (rfp RefreshProtocol) GenShare(sk *rlwe.SecretKey, ct *rlwe.Ciphertext, crp multiparty.KeySwitchCRP, shareOut *multiparty.RefreshShare) (err error) {
	return rfp.MaskedTransformProtocol.GenShare(sk, sk, ct, crp, nil, shareOut)
}

// AggregateShares aggregates two parties' shares in the Refresh protocol.
func (rfp RefreshProtocol) AggregateShares(share1, share2 multiparty.RefreshShare, shareOut *multiparty.RefreshShare) (err error) {
	return rfp.MaskedTransformProtocol.AggregateShares(share1, share2, shareOut)
}

// Finalize applies Decrypt, Recode and Recrypt on the input ciphertext.
func (rfp RefreshProtocol) Finalize(ctIn *rlwe.Ciphertext, crp multiparty.KeySwitchCRP, share multiparty.RefreshShare, opOut *rlwe.Ciphertext) (err error) {
	return rfp.MaskedTransformProtocol.Transform(ctIn, nil, crp, share, opOut)
}
