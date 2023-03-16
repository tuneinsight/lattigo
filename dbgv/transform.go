package dbgv

import (
	"encoding/binary"
	"fmt"

	"github.com/tuneinsight/lattigo/v4/bgv"
	"github.com/tuneinsight/lattigo/v4/drlwe"
	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/rlwe"
	"github.com/tuneinsight/lattigo/v4/utils"
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

// MaskedTransformFunc is a struct containing a user-defined in-place function that can be applied to masked bgv plaintexts, as a part of the
// Masked Transform Protocol.
// The function is called with a vector of integers modulo bgv.Parameters.T() of size bgv.Parameters.N() as input, and must write
// its output on the same buffer.
// Transform can be the identity.
// Decode: if true, then the masked BFV plaintext will be decoded before applying Transform.
// Recode: if true, then the masked BFV plaintext will be recoded after applying Transform.
// i.e. : Decode (true/false) -> Transform -> Recode (true/false).
type MaskedTransformFunc struct {
	Decode bool
	Func   func(coeffs []uint64)
	Encode bool
}

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

// UnmarshalBinary decodes a marshalled RefreshShare on the target RefreshShare.
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
func NewMaskedTransformProtocol(paramsIn, paramsOut bgv.Parameters, sigmaSmudging float64) (rfp *MaskedTransformProtocol, err error) {

	if paramsIn.N() > paramsOut.N() {
		return nil, fmt.Errorf("newMaskedTransformProtocol: paramsIn.N() != paramsOut.N()")
	}

	rfp = new(MaskedTransformProtocol)
	rfp.e2s = *NewE2SProtocol(paramsIn, sigmaSmudging)
	rfp.s2e = *NewS2EProtocol(paramsOut, sigmaSmudging)

	rfp.tmpPt = paramsOut.RingQ().NewPoly()
	rfp.tmpMask = paramsIn.RingT().NewPoly()
	rfp.tmpMaskPerm = paramsIn.RingT().NewPoly()
	return
}

// SampleCRP samples a common random polynomial to be used in the Masked-Transform protocol from the provided
// common reference string.
func (rfp *MaskedTransformProtocol) SampleCRP(level int, crs utils.PRNG) drlwe.CKSCRP {
	return rfp.s2e.SampleCRP(level, crs)
}

// AllocateShare allocates the shares of the PermuteProtocol
func (rfp *MaskedTransformProtocol) AllocateShare(levelDecrypt, levelRecrypt int) *MaskedTransformShare {
	return &MaskedTransformShare{*rfp.e2s.AllocateShare(levelDecrypt), *rfp.s2e.AllocateShare(levelRecrypt)}
}

// GenShare generates the shares of the PermuteProtocol.
// ct1 is the degree 1 element of a bgv.Ciphertext, i.e. bgv.Ciphertext.Value[1].
func (rfp *MaskedTransformProtocol) GenShare(skIn, skOut *rlwe.SecretKey, ct *rlwe.Ciphertext, scale rlwe.Scale, crs drlwe.CKSCRP, transform *MaskedTransformFunc, shareOut *MaskedTransformShare) {

	if ct.Level() < shareOut.e2sShare.Value.Level() {
		panic("cannot GenShare: ct[1] level must be at least equal to e2sShare level")
	}

	if (*ring.Poly)(&crs).Level() != shareOut.s2eShare.Value.Level() {
		panic("cannot GenShare: crs level must be equal to s2eShare")
	}

	rfp.e2s.GenShare(skIn, ct, &rlwe.AdditiveShare{Value: *rfp.tmpMask}, &shareOut.e2sShare)
	mask := rfp.tmpMask
	if transform != nil {
		coeffs := make([]uint64, len(mask.Coeffs[0]))

		if transform.Decode {
			rfp.e2s.encoder.DecodeRingT(mask, scale, coeffs)
		} else {
			copy(coeffs, mask.Coeffs[0])
		}

		transform.Func(coeffs)

		if transform.Encode {
			rfp.s2e.encoder.EncodeRingT(coeffs, scale, rfp.tmpMaskPerm)
		} else {
			copy(rfp.tmpMaskPerm.Coeffs[0], coeffs)
		}

		mask = rfp.tmpMaskPerm
	}
	rfp.s2e.GenShare(skOut, crs, &rlwe.AdditiveShare{Value: *mask}, &shareOut.s2eShare)
}

// AggregateShares sums share1 and share2 on shareOut.
func (rfp *MaskedTransformProtocol) AggregateShares(share1, share2, shareOut *MaskedTransformShare) {

	if share1.e2sShare.Value.Level() != share2.e2sShare.Value.Level() || share1.e2sShare.Value.Level() != shareOut.e2sShare.Value.Level() {
		panic("cannot AggregateShares: all e2s shares must be at the same level")
	}

	if share1.s2eShare.Value.Level() != share2.s2eShare.Value.Level() || share1.s2eShare.Value.Level() != shareOut.s2eShare.Value.Level() {
		panic("cannot AggregateShares: all s2e shares must be at the same level")
	}

	rfp.e2s.params.RingQ().AtLevel(share1.e2sShare.Value.Level()).Add(share1.e2sShare.Value, share2.e2sShare.Value, shareOut.e2sShare.Value)
	rfp.s2e.params.RingQ().AtLevel(share1.s2eShare.Value.Level()).Add(share1.s2eShare.Value, share2.s2eShare.Value, shareOut.s2eShare.Value)
}

// Transform applies Decrypt, Recode and Recrypt on the input ciphertext.
func (rfp *MaskedTransformProtocol) Transform(ct *rlwe.Ciphertext, transform *MaskedTransformFunc, crs drlwe.CKSCRP, share *MaskedTransformShare, ciphertextOut *rlwe.Ciphertext) {

	if ct.Level() < share.e2sShare.Value.Level() {
		panic("cannot Transform: input ciphertext level must be at least equal to e2s level")
	}

	maxLevel := (*ring.Poly)(&crs).Level()

	if maxLevel != share.s2eShare.Value.Level() {
		panic("cannot Transform: crs level and s2e level must be the same")
	}

	rfp.e2s.GetShare(nil, &share.e2sShare, ct, &rlwe.AdditiveShare{Value: *rfp.tmpMask}) // tmpMask RingT(m - sum M_i)
	mask := rfp.tmpMask
	if transform != nil {
		coeffs := make([]uint64, len(mask.Coeffs[0]))

		if transform.Decode {
			rfp.e2s.encoder.DecodeRingT(mask, ciphertextOut.Scale, coeffs)
		} else {
			copy(coeffs, mask.Coeffs[0])
		}

		transform.Func(coeffs)

		if transform.Encode {
			rfp.s2e.encoder.EncodeRingT(coeffs, ciphertextOut.Scale, rfp.tmpMaskPerm)
		} else {
			copy(rfp.tmpMaskPerm.Coeffs[0], coeffs)
		}

		mask = rfp.tmpMaskPerm
	}

	ciphertextOut.Resize(ciphertextOut.Degree(), maxLevel)

	rfp.s2e.encoder.RingT2Q(maxLevel, mask, rfp.tmpPt)
	rfp.s2e.encoder.ScaleUp(maxLevel, rfp.tmpPt, rfp.tmpPt)
	rfp.s2e.params.RingQ().AtLevel(maxLevel).NTT(rfp.tmpPt, rfp.tmpPt)
	rfp.s2e.params.RingQ().AtLevel(maxLevel).Add(rfp.tmpPt, share.s2eShare.Value, ciphertextOut.Value[0])
	rfp.s2e.GetEncryption(&drlwe.CKSShare{Value: ciphertextOut.Value[0]}, crs, ciphertextOut)
}
