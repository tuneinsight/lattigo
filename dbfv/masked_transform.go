package dbfv

import (
	"github.com/ldsec/lattigo/v2/bfv"
	"github.com/ldsec/lattigo/v2/drlwe"
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/rlwe"
)

// MaskedTransformProtocol is a struct storing the parameters for the MaskedTransformProtocol protocol.
type MaskedTransformProtocol struct {
	e2s E2SProtocol
	s2e S2EProtocol

	encoder bfv.Encoder
	ringQ   *ring.Ring
	ringT   *ring.Ring

	tmpPt       bfv.Plaintext
	tmpPtRt     bfv.PlaintextRingT
	tmpMask     *ring.Poly
	tmpMaskPerm *ring.Poly
}

// MaskedTransformShare is a struct storing the decryption and recryption shares.
type MaskedTransformShare struct {
	e2sShare drlwe.CKSShare
	s2eShare drlwe.CKSShare
}

// MarshalBinary encodes a RefreshShare on a slice of bytes.
func (share *MaskedTransformShare) MarshalBinary() ([]byte, error) {
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
func (share *MaskedTransformShare) UnmarshalBinary(data []byte) error {
	shareLen := len(data) >> 1
	share.e2sShare.UnmarshalBinary(data[:shareLen])
	share.s2eShare.UnmarshalBinary(data[shareLen:])
	return nil
}

// NewMaskedTransformProtocol creates a new instance of the PermuteProtocol.
func NewMaskedTransformProtocol(params bfv.Parameters, sigmaSmudging float64) (rfp *MaskedTransformProtocol) {

	rfp = new(MaskedTransformProtocol)
	rfp.e2s = *NewE2SProtocol(params, sigmaSmudging)
	rfp.s2e = *NewS2EProtocol(params, sigmaSmudging)

	rfp.encoder = rfp.e2s.encoder

	rfp.ringQ = rfp.e2s.ringQ
	rfp.ringT = rfp.e2s.ringT
	rfp.tmpPt = *bfv.NewPlaintext(params)
	rfp.tmpPtRt = *bfv.NewPlaintextRingT(params)
	rfp.tmpMask = rfp.ringT.NewPoly()
	rfp.tmpMaskPerm = rfp.ringT.NewPoly()
	return
}

// AllocateShares allocates the shares of the PermuteProtocol
func (rfp *MaskedTransformProtocol) AllocateShares() MaskedTransformShare {
	return MaskedTransformShare{*rfp.e2s.AllocateShare(), *rfp.s2e.AllocateShare()}
}

// GenShares generates the shares of the PermuteProtocol
func (rfp *MaskedTransformProtocol) GenShares(sk *rlwe.SecretKey, ciphertext *bfv.Ciphertext, crs *ring.Poly, transform func(*ring.Poly, *ring.Poly), shareOut MaskedTransformShare) {
	rfp.e2s.GenShare(sk, ciphertext, AdditiveShare{*rfp.tmpMask}, &shareOut.e2sShare)
	mask := rfp.tmpMask
	if transform != nil {
		transform(rfp.tmpMask, rfp.tmpMaskPerm)
		mask = rfp.tmpMaskPerm
	}
	rfp.s2e.GenShare(sk, crs, AdditiveShare{*mask}, &shareOut.s2eShare)
}

// Aggregate sums share1 and share2 on shareOut.
func (rfp *MaskedTransformProtocol) Aggregate(share1, share2, shareOut MaskedTransformShare) {
	rfp.ringQ.Add(share1.e2sShare.Value, share2.e2sShare.Value, shareOut.e2sShare.Value)
	rfp.ringQ.Add(share1.s2eShare.Value, share2.s2eShare.Value, shareOut.s2eShare.Value)
}

// Finalize applies Decrypt, Recode and Recrypt on the input ciphertext.
func (rfp *MaskedTransformProtocol) Finalize(ciphertext *bfv.Ciphertext, transform func(*ring.Poly, *ring.Poly), crs *ring.Poly, share MaskedTransformShare, ciphertextOut *bfv.Ciphertext) {
	rfp.e2s.Finalize(nil, &share.e2sShare, ciphertext, &AdditiveShare{*rfp.tmpMask}) // tmpMask RingT(m - sum M_i)
	mask := rfp.tmpMask
	if transform != nil {
		transform(rfp.tmpMask, rfp.tmpMaskPerm)
		mask = rfp.tmpMaskPerm
	}
	rfp.tmpPtRt.Value[0].Copy(mask)
	rfp.encoder.ScaleUp(&rfp.tmpPtRt, &rfp.tmpPt)
	//rfp.encoder.EncodeUint(mask.Coeffs[0], &rfp.tmpPt) // tmpMask RingQ( Delta*(m - sum M_i))
	rfp.ringQ.Add(rfp.tmpPt.Value[0], share.s2eShare.Value, ciphertextOut.Value[0])
	rfp.s2e.Finalize(&drlwe.CKSShare{Value: ciphertextOut.Value[0]}, crs, *ciphertextOut)
}
