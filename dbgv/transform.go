package dbgv

import (
	"fmt"

	"github.com/tuneinsight/lattigo/v4/bgv"
	"github.com/tuneinsight/lattigo/v4/drlwe"
	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/ring/distribution"
	"github.com/tuneinsight/lattigo/v4/rlwe"
	"github.com/tuneinsight/lattigo/v4/utils/sampling"
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

// NewMaskedTransformProtocol creates a new instance of the PermuteProtocol.
func NewMaskedTransformProtocol(paramsIn, paramsOut bgv.Parameters, noise distribution.Distribution) (rfp *MaskedTransformProtocol, err error) {

	if paramsIn.N() > paramsOut.N() {
		return nil, fmt.Errorf("newMaskedTransformProtocol: paramsIn.N() != paramsOut.N()")
	}

	rfp = new(MaskedTransformProtocol)
	rfp.e2s = *NewE2SProtocol(paramsIn, noise)
	rfp.s2e = *NewS2EProtocol(paramsOut, noise)

	rfp.tmpPt = paramsOut.RingQ().NewPoly()
	rfp.tmpMask = paramsIn.RingT().NewPoly()
	rfp.tmpMaskPerm = paramsIn.RingT().NewPoly()
	return
}

// SampleCRP samples a common random polynomial to be used in the Masked-Transform protocol from the provided
// common reference string.
func (rfp *MaskedTransformProtocol) SampleCRP(level int, crs sampling.PRNG) drlwe.CKSCRP {
	return rfp.s2e.SampleCRP(level, crs)
}

// AllocateShare allocates the shares of the PermuteProtocol
func (rfp *MaskedTransformProtocol) AllocateShare(levelDecrypt, levelRecrypt int) *drlwe.RefreshShare {
	return &drlwe.RefreshShare{E2SShare: *rfp.e2s.AllocateShare(levelDecrypt), S2EShare: *rfp.s2e.AllocateShare(levelRecrypt)}
}

// GenShare generates the shares of the PermuteProtocol.
// ct1 is the degree 1 element of a bgv.Ciphertext, i.e. bgv.Ciphertext.Value[1].
func (rfp *MaskedTransformProtocol) GenShare(skIn, skOut *rlwe.SecretKey, ct *rlwe.Ciphertext, scale rlwe.Scale, crs drlwe.CKSCRP, transform *MaskedTransformFunc, shareOut *drlwe.RefreshShare) {

	if ct.Level() < shareOut.E2SShare.Value.Level() {
		panic("cannot GenShare: ct[1] level must be at least equal to E2SShare level")
	}

	if crs.Value.Level() != shareOut.S2EShare.Value.Level() {
		panic("cannot GenShare: crs level must be equal to S2EShare")
	}

	rfp.e2s.GenShare(skIn, ct, &drlwe.AdditiveShare{Value: *rfp.tmpMask}, &shareOut.E2SShare)
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
	rfp.s2e.GenShare(skOut, crs, &drlwe.AdditiveShare{Value: *mask}, &shareOut.S2EShare)
}

// AggregateShares sums share1 and share2 on shareOut.
func (rfp *MaskedTransformProtocol) AggregateShares(share1, share2, shareOut *drlwe.RefreshShare) {

	if share1.E2SShare.Value.Level() != share2.E2SShare.Value.Level() || share1.E2SShare.Value.Level() != shareOut.E2SShare.Value.Level() {
		panic("cannot AggregateShares: all e2s shares must be at the same level")
	}

	if share1.S2EShare.Value.Level() != share2.S2EShare.Value.Level() || share1.S2EShare.Value.Level() != shareOut.S2EShare.Value.Level() {
		panic("cannot AggregateShares: all s2e shares must be at the same level")
	}

	rfp.e2s.params.RingQ().AtLevel(share1.E2SShare.Value.Level()).Add(share1.E2SShare.Value, share2.E2SShare.Value, shareOut.E2SShare.Value)
	rfp.s2e.params.RingQ().AtLevel(share1.S2EShare.Value.Level()).Add(share1.S2EShare.Value, share2.S2EShare.Value, shareOut.S2EShare.Value)
}

// Transform applies Decrypt, Recode and Recrypt on the input ciphertext.
func (rfp *MaskedTransformProtocol) Transform(ct *rlwe.Ciphertext, transform *MaskedTransformFunc, crs drlwe.CKSCRP, share *drlwe.RefreshShare, ciphertextOut *rlwe.Ciphertext) {

	if ct.Level() < share.E2SShare.Value.Level() {
		panic("cannot Transform: input ciphertext level must be at least equal to e2s level")
	}

	maxLevel := crs.Value.Level()

	if maxLevel != share.S2EShare.Value.Level() {
		panic("cannot Transform: crs level and s2e level must be the same")
	}

	rfp.e2s.GetShare(nil, &share.E2SShare, ct, &drlwe.AdditiveShare{Value: *rfp.tmpMask}) // tmpMask RingT(m - sum M_i)
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

	rfp.s2e.encoder.RingT2Q(maxLevel, true, mask, rfp.tmpPt)
	rfp.s2e.params.RingQ().AtLevel(maxLevel).NTT(rfp.tmpPt, rfp.tmpPt)
	rfp.s2e.params.RingQ().AtLevel(maxLevel).Add(rfp.tmpPt, share.S2EShare.Value, ciphertextOut.Value[0])
	rfp.s2e.GetEncryption(&drlwe.CKSShare{Value: ciphertextOut.Value[0]}, crs, ciphertextOut)
}
