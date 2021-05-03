package dbfv

import (
	"github.com/ldsec/lattigo/v2/bfv"
	"github.com/ldsec/lattigo/v2/drlwe"
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/rlwe"
)

// RefreshProtocol is a struct storing the relevant parameters for the Refresh protocol.
type RefreshProtocol struct {
	e2s E2SProtocol
	s2e S2EProtocol

	encoder bfv.Encoder
	ringQ   *ring.Ring

	tmpPt   bfv.Plaintext
	tmpMask *ring.Poly
}

// RefreshShare is a struct storing the decryption and recryption shares.
type RefreshShare struct {
	e2sShare drlwe.CKSShare
	s2eShare drlwe.CKSShare
}

// MarshalBinary encodes a RefreshShare on a slice of bytes.
func (share *RefreshShare) MarshalBinary() ([]byte, error) {
	e2sData, err := share.e2sShare.MarshalBinary()
	if err != nil {
		return nil, err
	}
	s2eData, err := share.s2eShare.MarshalBinary()
	if err != nil {
		return nil, err
	}
	return append(e2sData, s2eData...), nil
}

// UnmarshalBinary decodes a marshaled RefreshShare on the target RefreshShare.
func (share *RefreshShare) UnmarshalBinary(data []byte) error {
	shareLen := len(data) >> 1
	share.e2sShare.UnmarshalBinary(data[:shareLen])
	share.s2eShare.UnmarshalBinary(data[shareLen:])
	return nil
}

// NewRefreshProtocol creates a new Refresh protocol instance.
func NewRefreshProtocol(params bfv.Parameters, sigmaSmudging float64) (rfp *RefreshProtocol) {

	rfp = new(RefreshProtocol)
	rfp.e2s = *NewE2SProtocol(params, sigmaSmudging)
	rfp.s2e = *NewS2EProtocol(params, sigmaSmudging)

	rfp.encoder = rfp.e2s.encoder

	rfp.ringQ = rfp.e2s.ringQ
	rfp.tmpPt = *bfv.NewPlaintext(params)
	rfp.tmpMask = rfp.ringQ.NewPoly()

	return
}

// AllocateShares allocates the shares of the Refresh protocol.
func (rfp *RefreshProtocol) AllocateShares() RefreshShare {
	return RefreshShare{*rfp.e2s.AllocateShare(), *rfp.s2e.AllocateShare()}
}

// GenShares generates a share for the Refresh protocol.
func (rfp *RefreshProtocol) GenShares(sk *rlwe.SecretKey, ciphertext *bfv.Ciphertext, crs *ring.Poly, shareOut RefreshShare) {
	rfp.e2s.GenShare(sk, ciphertext, AdditiveShare{*rfp.tmpMask}, &shareOut.e2sShare)
	rfp.s2e.GenShare(sk, crs, AdditiveShare{*rfp.tmpMask}, &shareOut.s2eShare)
}

// Aggregate sums share1 and share2 on shareOut.
func (rfp *RefreshProtocol) Aggregate(share1, share2, shareOut RefreshShare) {
	rfp.ringQ.Add(share1.e2sShare.Value, share2.e2sShare.Value, shareOut.e2sShare.Value)
	rfp.ringQ.Add(share1.s2eShare.Value, share2.s2eShare.Value, shareOut.s2eShare.Value)
}

// Finalize applies Decrypt, Recode and Recrypt on the input ciphertext.
func (rfp *RefreshProtocol) Finalize(ciphertext *bfv.Ciphertext, crs *ring.Poly, share RefreshShare, ciphertextOut *bfv.Ciphertext) {
	rfp.e2s.Finalize(nil, &share.e2sShare, ciphertext, &AdditiveShare{*rfp.tmpMask}) // tmpMask RingT(m - sum M_i)
	rfp.encoder.EncodeUint(rfp.tmpMask.Coeffs[0], &rfp.tmpPt)                        // tmpMask RingQ( Delta*(m - sum M_i))
	rfp.ringQ.Add(rfp.tmpPt.Value[0], share.s2eShare.Value, ciphertextOut.Value[0])
	rfp.s2e.Finalize(&drlwe.CKSShare{Value: ciphertextOut.Value[0]}, crs, *ciphertextOut)
}
