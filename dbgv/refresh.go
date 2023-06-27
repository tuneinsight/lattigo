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
func (rfp *RefreshProtocol) ShallowCopy() RefreshProtocol {
	return RefreshProtocol{rfp.MaskedTransformProtocol.ShallowCopy()}
}

// NewRefreshProtocol creates a new Refresh protocol instance.
func NewRefreshProtocol(params bgv.Parameters, noiseFlooding ring.DistributionParameters) (rfp RefreshProtocol) {
	rfp = RefreshProtocol{}
	mt, _ := NewMaskedTransformProtocol(params, params, noiseFlooding)
	rfp.MaskedTransformProtocol = mt
	return
}

// AllocateShare allocates the shares of the PermuteProtocol
func (rfp RefreshProtocol) AllocateShare(inputLevel, outputLevel int) drlwe.RefreshShare {
	return rfp.MaskedTransformProtocol.AllocateShare(inputLevel, outputLevel)
}

// GenShare generates a share for the Refresh protocol.
// ct1 is degree 1 element of a rlwe.Ciphertext, i.e. rlwe.Ciphertext.Value[1].
func (rfp RefreshProtocol) GenShare(sk *rlwe.SecretKey, ct *rlwe.Ciphertext, scale rlwe.Scale, crp drlwe.KeySwitchCRP, shareOut *drlwe.RefreshShare) {
	rfp.MaskedTransformProtocol.GenShare(sk, sk, ct, scale, crp, nil, shareOut)
}

// AggregateShares aggregates two parties' shares in the Refresh protocol.
func (rfp RefreshProtocol) AggregateShares(share1, share2 drlwe.RefreshShare, shareOut *drlwe.RefreshShare) {
	rfp.MaskedTransformProtocol.AggregateShares(share1, share2, shareOut)
}

// Finalize applies Decrypt, Recode and Recrypt on the input ciphertext.
func (rfp RefreshProtocol) Finalize(ctIn *rlwe.Ciphertext, crp drlwe.KeySwitchCRP, share drlwe.RefreshShare, ctOut *rlwe.Ciphertext) {
	rfp.MaskedTransformProtocol.Transform(ctIn, nil, crp, share, ctOut)
}
