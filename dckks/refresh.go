package dckks

import (
	"github.com/ldsec/lattigo/v2/ckks"
	"github.com/ldsec/lattigo/v2/drlwe"
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/rlwe"
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
func NewRefreshProtocol(params ckks.Parameters, precision int, sigmaSmudging float64) (rfp *RefreshProtocol) {
	rfp = new(RefreshProtocol)
	rfp.MaskedTransformProtocol = *NewMaskedTransformProtocol(params, precision, sigmaSmudging)
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

// GenShares generates a share for the Refresh protocol.
// This protocol requires additional inputs which are :
// logBound : the bit length of the masks
// logSlots : the bit length of the number of slots
//
// The method "GetMinimumLevelForBootstrapping" should be used to get the minimum level at which the refresh can be called while still ensure 128-bits of security, as well as the
// value for logBound.
func (rfp *RefreshProtocol) GenShares(sk *rlwe.SecretKey, logBound, logSlots int, c1 *ring.Poly, scale float64, crs drlwe.CKSCRP, shareOut *RefreshShare) {
	rfp.MaskedTransformProtocol.GenShares(sk, logBound, logSlots, c1, scale, crs, nil, &shareOut.MaskedTransformShare)
}

// Aggregate aggregates two parties' shares in the Refresh protocol.
func (rfp *RefreshProtocol) Aggregate(share1, share2, shareOut *RefreshShare) {
	rfp.MaskedTransformProtocol.Aggregate(&share1.MaskedTransformShare, &share2.MaskedTransformShare, &shareOut.MaskedTransformShare)
}

// Finalize applies Decrypt, Recode and Recrypt on the input ciphertext.
func (rfp *RefreshProtocol) Finalize(ctIn *ckks.Ciphertext, logSlots int, crs drlwe.CKSCRP, share *RefreshShare, ctOut *ckks.Ciphertext) {
	rfp.MaskedTransformProtocol.Transform(ctIn, logSlots, nil, crs, &share.MaskedTransformShare, ctOut)
}
