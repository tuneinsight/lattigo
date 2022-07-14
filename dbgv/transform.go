package dbgv

import (
	"encoding/binary"

	"github.com/tuneinsight/lattigo/v3/bgv"
	"github.com/tuneinsight/lattigo/v3/drlwe"
	"github.com/tuneinsight/lattigo/v3/ring"
	"github.com/tuneinsight/lattigo/v3/rlwe"
	"github.com/tuneinsight/lattigo/v3/utils"
)

// MaskedTransformProtocol is a struct storing the parameters for the MaskedTransformProtocol protocol.
type MaskedTransformProtocol struct {
	e2s E2SProtocol
	s2e S2EProtocol

	tmpPt       *ring.Poly
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
		tmpPt:       params.RingQ().NewPoly(),
		tmpMask:     params.RingT().NewPoly(),
		tmpMaskPerm: params.RingT().NewPoly(),
	}
}

// MaskedTransformFunc represents a user-defined in-place function that can be applied to masked bgv plaintexts, as a part of the
// Masked Transform Protocol.
// The function is called with a vector of integers modulo bgv.Parameters.T() of size bgv.Parameters.N() as input, and must write
// its output on the same buffer.
type MaskedTransformFunc func(coeffs []uint64)

// MaskedTransformShare is a struct storing the decryption and recryption shares.
type MaskedTransformShare struct {
	e2sShare drlwe.CKSShare
	s2eShare drlwe.CKSShare
}

// MarshalBinary encodes a RefreshShare on a slice of bytes.
func (share *MaskedTransformShare) MarshalBinary() (data []byte, err error) {
	var e2sData, s2eData []byte
	if e2sData, err = share.e2sShare.MarshalBinary(); err != nil {
		return nil, err
	}
	if s2eData, err = share.s2eShare.MarshalBinary(); err != nil {
		return nil, err
	}
	data = make([]byte, 8)
	binary.LittleEndian.PutUint64(data, uint64(len(e2sData)))
	data = append(data, e2sData...)
	data = append(data, s2eData...)
	return data, nil
}

// UnmarshalBinary decodes a marshaled RefreshShare on the target RefreshShare.
func (share *MaskedTransformShare) UnmarshalBinary(data []byte) error {

	e2sDataLen := binary.LittleEndian.Uint64(data[:8])

	if err := share.e2sShare.UnmarshalBinary(data[8 : e2sDataLen+8]); err != nil {
		return err
	}
	if err := share.s2eShare.UnmarshalBinary(data[8+e2sDataLen:]); err != nil {
		return err
	}
	return nil
}

// NewMaskedTransformProtocol creates a new instance of the PermuteProtocol.
func NewMaskedTransformProtocol(params bgv.Parameters, sigmaSmudging float64) (rfp *MaskedTransformProtocol) {

	rfp = new(MaskedTransformProtocol)
	rfp.e2s = *NewE2SProtocol(params, sigmaSmudging)
	rfp.s2e = *NewS2EProtocol(params, sigmaSmudging)

	rfp.tmpPt = params.RingQ().NewPoly()
	rfp.tmpMask = params.RingT().NewPoly()
	rfp.tmpMaskPerm = params.RingT().NewPoly()
	return
}

// SampleCRP samples a common random polynomial to be used in the Masked-Transform protocol from the provided
// common reference string.
func (rfp *MaskedTransformProtocol) SampleCRP(level int, crs utils.PRNG) drlwe.CKSCRP {
	return rfp.e2s.SampleCRP(level, crs)
}

// AllocateShare allocates the shares of the PermuteProtocol
func (rfp *MaskedTransformProtocol) AllocateShare(levelDecrypt, levelRecrypt int) *MaskedTransformShare {
	return &MaskedTransformShare{*rfp.e2s.AllocateShare(levelDecrypt), *rfp.s2e.AllocateShare(levelRecrypt)}
}

// GenShare generates the shares of the PermuteProtocol.
// ct1 is the degree 1 element of a bgv.Ciphertext, i.e. bgv.Ciphertext.Value[1].
func (rfp *MaskedTransformProtocol) GenShare(sk *rlwe.SecretKey, ct1 *ring.Poly, scale uint64, crs drlwe.CKSCRP, transform MaskedTransformFunc, shareOut *MaskedTransformShare) {

	if ct1.Level() < shareOut.e2sShare.Value.Level() {
		panic("ct[1] level must be at least equal to e2sShare level")
	}

	if (*ring.Poly)(&crs).Level() != shareOut.s2eShare.Value.Level() {
		panic("crs level must be equal to s2eShare")
	}

	rfp.e2s.GenShare(sk, ct1, &rlwe.AdditiveShare{Value: *rfp.tmpMask}, &shareOut.e2sShare)
	mask := rfp.tmpMask
	if transform != nil {
		coeffs := make([]uint64, len(mask.Coeffs[0]))
		rfp.e2s.encoder.DecodeRingT(mask, scale, coeffs)
		transform(coeffs)
		rfp.e2s.encoder.EncodeRingT(coeffs, scale, rfp.tmpMaskPerm)
		mask = rfp.tmpMaskPerm
	}
	rfp.s2e.GenShare(sk, crs, &rlwe.AdditiveShare{Value: *mask}, &shareOut.s2eShare)
}

// AggregateShare sums share1 and share2 on shareOut.
func (rfp *MaskedTransformProtocol) AggregateShare(share1, share2, shareOut *MaskedTransformShare) {

	if share1.e2sShare.Value.Level() != share2.e2sShare.Value.Level() || share1.e2sShare.Value.Level() != shareOut.e2sShare.Value.Level() {
		panic("all e2s shares must be at the same level")
	}

	if share1.s2eShare.Value.Level() != share2.s2eShare.Value.Level() || share1.s2eShare.Value.Level() != shareOut.s2eShare.Value.Level() {
		panic("all s2e shares must be at the same level")
	}

	rfp.e2s.params.RingQ().AddLvl(share1.e2sShare.Value.Level(), share1.e2sShare.Value, share2.e2sShare.Value, shareOut.e2sShare.Value)
	rfp.e2s.params.RingQ().AddLvl(share1.s2eShare.Value.Level(), share1.s2eShare.Value, share2.s2eShare.Value, shareOut.s2eShare.Value)
}

// Transform applies Decrypt, Recode and Recrypt on the input ciphertext.
func (rfp *MaskedTransformProtocol) Transform(ct *bgv.Ciphertext, transform MaskedTransformFunc, crs drlwe.CKSCRP, share *MaskedTransformShare, ciphertextOut *bgv.Ciphertext) {

	if ct.Level() < share.e2sShare.Value.Level() {
		panic("input ciphertext level must be at least equal to e2s level")
	}

	maxLevel := (*ring.Poly)(&crs).Level()

	if maxLevel != share.s2eShare.Value.Level() {
		panic("crs level and s2e level must be the same")
	}

	rfp.e2s.GetShare(nil, &share.e2sShare, ct, &rlwe.AdditiveShare{Value: *rfp.tmpMask}) // tmpMask RingT(m - sum M_i)
	mask := rfp.tmpMask
	if transform != nil {
		coeffs := make([]uint64, len(mask.Coeffs[0]))
		rfp.e2s.encoder.DecodeRingT(mask, ct.Scale, coeffs)
		transform(coeffs)
		rfp.e2s.encoder.EncodeRingT(coeffs, ciphertextOut.Scale, rfp.tmpMaskPerm)
		mask = rfp.tmpMaskPerm
	}

	ciphertextOut.Resize(ciphertextOut.Degree(), maxLevel)

	rfp.e2s.encoder.RingT2Q(maxLevel, mask, rfp.tmpPt)
	rfp.e2s.params.RingQ().NTTLvl(maxLevel, rfp.tmpPt, rfp.tmpPt)
	rfp.e2s.params.RingQ().AddLvl(maxLevel, rfp.tmpPt, share.s2eShare.Value, ciphertextOut.Value[0])
	rfp.s2e.GetEncryption(&drlwe.CKSShare{Value: ciphertextOut.Value[0]}, crs, ciphertextOut)
}
