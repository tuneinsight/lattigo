package mheint

import (
	"github.com/tuneinsight/lattigo/v5/he/heint"
	"github.com/tuneinsight/lattigo/v5/mhe"
	"github.com/tuneinsight/lattigo/v5/ring"

	"github.com/tuneinsight/lattigo/v5/core/rlwe"
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
func NewRefreshProtocol(params heint.Parameters, noiseFlooding ring.DistributionParameters) (rfp RefreshProtocol, err error) {
	rfp = RefreshProtocol{}
	mt, err := NewMaskedTransformProtocol(params, params, noiseFlooding)
	rfp.MaskedTransformProtocol = mt
	return rfp, err
}

// AllocateShare allocates the shares of the PermuteProtocol
func (rfp RefreshProtocol) AllocateShare(inputLevel, outputLevel int) mhe.RefreshShare {
	return rfp.MaskedTransformProtocol.AllocateShare(inputLevel, outputLevel)
}

// GenShare generates a share for the Refresh protocol.
// ct1 is degree 1 element of a rlwe.Ciphertext, i.e. rlwe.Ciphertext.Value[1].
func (rfp RefreshProtocol) GenShare(sk *rlwe.SecretKey, ct *rlwe.Ciphertext, crp mhe.KeySwitchCRP, shareOut *mhe.RefreshShare) (err error) {
	return rfp.MaskedTransformProtocol.GenShare(sk, sk, ct, crp, nil, shareOut)
}

// AggregateShares aggregates two parties' shares in the Refresh protocol.
func (rfp RefreshProtocol) AggregateShares(share1, share2 mhe.RefreshShare, shareOut *mhe.RefreshShare) (err error) {
	return rfp.MaskedTransformProtocol.AggregateShares(share1, share2, shareOut)
}

// Finalize applies Decrypt, Recode and Recrypt on the input ciphertext.
func (rfp RefreshProtocol) Finalize(ctIn *rlwe.Ciphertext, crp mhe.KeySwitchCRP, share mhe.RefreshShare, opOut *rlwe.Ciphertext) (err error) {
	return rfp.MaskedTransformProtocol.Transform(ctIn, nil, crp, share, opOut)
}
