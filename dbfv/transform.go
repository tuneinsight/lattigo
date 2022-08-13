package dbfv

import (
	"github.com/tuneinsight/lattigo/v3/bfv"
	"github.com/tuneinsight/lattigo/v3/drlwe"
	"github.com/tuneinsight/lattigo/v3/ring"
	"github.com/tuneinsight/lattigo/v3/rlwe"
	"github.com/tuneinsight/lattigo/v3/utils"
)

// MaskedTransformProtocol is a struct storing the parameters for the MaskedTransformProtocol protocol.
type MaskedTransformProtocol struct {
	e2s E2SProtocol
	s2e S2EProtocol

	tmpPt       bfv.Plaintext
	tmpMask     *ring.Poly
	tmpMaskPerm *ring.Poly
}

// ShallowCopy creates a shallow copy of MaskedTransformProtocol in which all the read-only data-structures are
// shared with the receiver and the temporary buffers are reallocated. The receiver and the returned
// MaskedTransformProtocol can be used concurrently.
func (rfp *MaskedTransformProtocol) ShallowCopy() *MaskedTransformProtocol {
	params := rfp.e2s.params

	return &MaskedTransformProtocol{
		e2s:         *rfp.e2s.ShallowCopy(),
		s2e:         *rfp.s2e.ShallowCopy(),
		tmpPt:       *bfv.NewPlaintext(params),
		tmpMask:     params.RingT().NewPoly(),
		tmpMaskPerm: params.RingT().NewPoly(),
	}
}

// MaskedTransformFunc represents a user-defined in-place function that can be applied to masked BFV plaintexts, as a part of the
// Masked Transform Protocol.
// The function is called with a vector of integers modulo bfv.Parameters.T() of size bfv.Parameters.N() as input, and must write
// its output on the same buffer.
type MaskedTransformFunc func(coeffs []uint64)

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

	rfp.tmpPt = *bfv.NewPlaintext(params)
	rfp.tmpMask = params.RingT().NewPoly()
	rfp.tmpMaskPerm = params.RingT().NewPoly()
	return
}

// SampleCRP samples a common random polynomial to be used in the Masked-Transform protocol from the provided
// common reference string.
func (rfp *MaskedTransformProtocol) SampleCRP(level int, crs utils.PRNG) drlwe.CKSCRP {
	return rfp.e2s.SampleCRP(level, crs)
}

// AllocateShare allocates the shares of the PermuteProtocol.
func (rfp *MaskedTransformProtocol) AllocateShare() *MaskedTransformShare {
	return &MaskedTransformShare{*rfp.e2s.AllocateShare(), *rfp.s2e.AllocateShare()}
}

// GenShare generates the shares of the PermuteProtocol.
// ct1 is the degree 1 element of a bfv.Ciphertext, i.e. bfv.Ciphertext.Value[1].
func (rfp *MaskedTransformProtocol) GenShare(sk *rlwe.SecretKey, c1 *ring.Poly, crs drlwe.CKSCRP, transform MaskedTransformFunc, shareOut *MaskedTransformShare) {
	rfp.e2s.GenShare(sk, c1, &rlwe.AdditiveShare{Value: *rfp.tmpMask}, &shareOut.e2sShare)
	mask := rfp.tmpMask
	if transform != nil {
		coeffs := rfp.e2s.encoder.DecodeUintNew(&bfv.PlaintextRingT{Plaintext: &rlwe.Plaintext{Value: mask}})
		transform(coeffs)
		rfp.e2s.encoder.EncodeRingT(coeffs, &bfv.PlaintextRingT{Plaintext: &rlwe.Plaintext{Value: rfp.tmpMaskPerm}})
		mask = rfp.tmpMaskPerm
	}
	rfp.s2e.GenShare(sk, crs, &rlwe.AdditiveShare{Value: *mask}, &shareOut.s2eShare)
}

// Aggregate sums share1 and share2 on shareOut.
func (rfp *MaskedTransformProtocol) Aggregate(share1, share2, shareOut *MaskedTransformShare) {
	rfp.e2s.params.RingQ().Add(share1.e2sShare.Value, share2.e2sShare.Value, shareOut.e2sShare.Value)
	rfp.e2s.params.RingQ().Add(share1.s2eShare.Value, share2.s2eShare.Value, shareOut.s2eShare.Value)
}

// Transform applies Decrypt, Recode and Recrypt on the input ciphertext.
func (rfp *MaskedTransformProtocol) Transform(ciphertext *bfv.Ciphertext, transform MaskedTransformFunc, crs drlwe.CKSCRP, share *MaskedTransformShare, ciphertextOut *bfv.Ciphertext) {
	rfp.e2s.GetShare(nil, &share.e2sShare, ciphertext, &rlwe.AdditiveShare{Value: *rfp.tmpMask}) // tmpMask RingT(m - sum M_i)
	mask := rfp.tmpMask
	if transform != nil {
		coeffs := rfp.e2s.encoder.DecodeUintNew(&bfv.PlaintextRingT{Plaintext: &rlwe.Plaintext{Value: mask}})
		transform(coeffs)
		rfp.e2s.encoder.EncodeRingT(coeffs, &bfv.PlaintextRingT{Plaintext: &rlwe.Plaintext{Value: rfp.tmpMaskPerm}})
		mask = rfp.tmpMaskPerm
	}
	rfp.e2s.encoder.ScaleUp(&bfv.PlaintextRingT{Plaintext: &rlwe.Plaintext{Value: mask}}, &rfp.tmpPt)
	rfp.e2s.params.RingQ().Add(rfp.tmpPt.Value, share.s2eShare.Value, ciphertextOut.Value[0])
	rfp.s2e.GetEncryption(&drlwe.CKSShare{Value: ciphertextOut.Value[0]}, crs, ciphertextOut)
}
