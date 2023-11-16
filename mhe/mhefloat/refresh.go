package mhefloat

import (
	"github.com/tuneinsight/lattigo/v5/he/hefloat"
	"github.com/tuneinsight/lattigo/v5/mhe"
	"github.com/tuneinsight/lattigo/v5/ring"

	"github.com/tuneinsight/lattigo/v5/core/rlwe"
)

// RefreshProtocol is a struct storing the relevant parameters for the Refresh protocol.
type RefreshProtocol struct {
	MaskedLinearTransformationProtocol
}

// NewRefreshProtocol creates a new Refresh protocol instance.
// prec : the log2 of decimal precision of the internal encoder.
func NewRefreshProtocol(params hefloat.Parameters, prec uint, noise ring.DistributionParameters) (rfp RefreshProtocol, err error) {
	rfp = RefreshProtocol{}
	mt, err := NewMaskedLinearTransformationProtocol(params, params, prec, noise)
	rfp.MaskedLinearTransformationProtocol = mt
	return rfp, err
}

// ShallowCopy creates a shallow copy of RefreshProtocol in which all the read-only data-structures are
// shared with the receiver and the temporary buffers are reallocated. The receiver and the returned
// RefreshProtocol can be used concurrently.
func (rfp RefreshProtocol) ShallowCopy() RefreshProtocol {
	return RefreshProtocol{rfp.MaskedLinearTransformationProtocol.ShallowCopy()}
}

// AllocateShare allocates the shares of the PermuteProtocol
func (rfp RefreshProtocol) AllocateShare(inputLevel, outputLevel int) mhe.RefreshShare {
	return rfp.MaskedLinearTransformationProtocol.AllocateShare(inputLevel, outputLevel)
}

// GenShare generates a share for the Refresh protocol.
// This protocol requires additional inputs which are :
// logBound : the bit length of the masks
// ct1      : the degree 1 element the ciphertext to refresh, i.e. ct1 = ckk.Ciphetext.Value[1].
// scale    : the scale of the ciphertext entering the refresh.
// The method "GetMinimumLevelForBootstrapping" should be used to get the minimum level at which the refresh can be called while still ensure 128-bits of security, as well as the
// value for logBound.
func (rfp RefreshProtocol) GenShare(sk *rlwe.SecretKey, logBound uint, ct *rlwe.Ciphertext, crs mhe.KeySwitchCRP, shareOut *mhe.RefreshShare) (err error) {
	return rfp.MaskedLinearTransformationProtocol.GenShare(sk, sk, logBound, ct, crs, nil, shareOut)
}

// AggregateShares aggregates two parties' shares in the Refresh protocol.
func (rfp RefreshProtocol) AggregateShares(share1, share2, shareOut *mhe.RefreshShare) (err error) {
	return rfp.MaskedLinearTransformationProtocol.AggregateShares(share1, share2, shareOut)
}

// Finalize applies Decrypt, Recode and Recrypt on the input ciphertext.
// The ciphertext scale is reset to the default scale.
func (rfp RefreshProtocol) Finalize(ctIn *rlwe.Ciphertext, crs mhe.KeySwitchCRP, share mhe.RefreshShare, opOut *rlwe.Ciphertext) (err error) {
	return rfp.MaskedLinearTransformationProtocol.Transform(ctIn, nil, crs, share, opOut)
}
