package dbfv

import (
	"github.com/ldsec/lattigo/v2/bfv"
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/rlwe"
)

// RefreshProtocol is a struct storing the relevant parameters for the Refresh protocol.
type RefreshProtocol struct {
	MaskedTransformProtocol
}

type RefreshShare struct {
	MaskedTransformShare
}

// NewRefreshProtocol creates a new Refresh protocol instance.
func NewRefreshProtocol(params bfv.Parameters, sigmaSmudging float64) (rfp *RefreshProtocol) {
	rfp = new(RefreshProtocol)
	rfp.MaskedTransformProtocol = *NewMaskedTransformProtocol(params, sigmaSmudging)
	return
}

// AllocateShares allocates the shares of the PermuteProtocol
func (rfp *RefreshProtocol) AllocateShares() RefreshShare {
	return RefreshShare{rfp.MaskedTransformProtocol.AllocateShares()}
}

// GenShares generates a share for the Refresh protocol.
func (rfp *RefreshProtocol) GenShares(sk *rlwe.SecretKey, ciphertext *bfv.Ciphertext, crs *ring.Poly, shareOut RefreshShare) {
	rfp.MaskedTransformProtocol.GenShares(sk, ciphertext, crs, nil, shareOut.MaskedTransformShare)
}

func (rfp *RefreshProtocol) Aggregate(share1, share2, shareOut RefreshShare) {
	rfp.MaskedTransformProtocol.Aggregate(share1.MaskedTransformShare, share2.MaskedTransformShare, shareOut.MaskedTransformShare)
}

// Finalize applies Decrypt, Recode and Recrypt on the input ciphertext.
func (rfp *RefreshProtocol) Finalize(ciphertext *bfv.Ciphertext, crs *ring.Poly, share RefreshShare, ciphertextOut *bfv.Ciphertext) {
	rfp.MaskedTransformProtocol.Finalize(ciphertext, nil, crs, share.MaskedTransformShare, ciphertextOut)
}
