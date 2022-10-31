package dckks

import (
	"github.com/tuneinsight/lattigo/v4/ckks"
	"github.com/tuneinsight/lattigo/v4/drlwe"
	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/rlwe"
)

// RefreshProtocol is a struct storing the relevant parameters for the Refresh protocol.
type RefreshProtocol struct {
	MaskedTransformProtocol
}

// RefreshShare is a struct storing a party's share in the Refresh protocol.
type RefreshShare struct {
	MaskedTransformShare
}

// NewRefreshProtocol creates a new Refresh protocol instance.
// prec : the log2 of decimal precision of the internal encoder.
func NewRefreshProtocol(params ckks.Parameters, prec uint, sigmaSmudging float64) (rfp *RefreshProtocol) {
	rfp = new(RefreshProtocol)
	mt, _ := NewMaskedTransformProtocol(params, params, prec, sigmaSmudging)
	rfp.MaskedTransformProtocol = *mt
	return
}

// ShallowCopy creates a shallow copy of RefreshProtocol in which all the read-only data-structures are
// shared with the receiver and the temporary buffers are reallocated. The receiver and the returned
// RefreshProtocol can be used concurrently.
func (rfp *RefreshProtocol) ShallowCopy() *RefreshProtocol {
	return &RefreshProtocol{*rfp.MaskedTransformProtocol.ShallowCopy()}
}

// AllocateShare allocates the shares of the PermuteProtocol
func (rfp *RefreshProtocol) AllocateShare(inputLevel, outputLevel int) *RefreshShare {
	share := rfp.MaskedTransformProtocol.AllocateShare(inputLevel, outputLevel)
	return &RefreshShare{*share}
}

// GenShare generates a share for the Refresh protocol.
// This protocol requires additional inputs which are :
// logBound : the bit length of the masks
// logSlots : the bit length of the number of slots
// ct1      : the degree 1 element the ciphertext to refresh, i.e. ct1 = ckk.Ciphetext.Value[1].
// scale    : the scale of the ciphertext entering the refresh.
// The method "GetMinimumLevelForBootstrapping" should be used to get the minimum level at which the refresh can be called while still ensure 128-bits of security, as well as the
// value for logBound.
func (rfp *RefreshProtocol) GenShare(sk *rlwe.SecretKey, logBound uint, logSlots int, ct1 *ring.Poly, scale rlwe.Scale, crs drlwe.CKSCRP, shareOut *RefreshShare) {
	rfp.MaskedTransformProtocol.GenShare(sk, sk, logBound, logSlots, ct1, scale, crs, nil, &shareOut.MaskedTransformShare)
}

// AggregateShares aggregates two parties' shares in the Refresh protocol.
func (rfp *RefreshProtocol) AggregateShares(share1, share2, shareOut *RefreshShare) {
	rfp.MaskedTransformProtocol.AggregateShares(&share1.MaskedTransformShare, &share2.MaskedTransformShare, &shareOut.MaskedTransformShare)
}

// Finalize applies Decrypt, Recode and Recrypt on the input ciphertext.
// The ciphertext scale is reset to the default scale.
func (rfp *RefreshProtocol) Finalize(ctIn *rlwe.Ciphertext, logSlots int, crs drlwe.CKSCRP, share *RefreshShare, ctOut *rlwe.Ciphertext) {
	rfp.MaskedTransformProtocol.Transform(ctIn, logSlots, nil, crs, &share.MaskedTransformShare, ctOut)
}
