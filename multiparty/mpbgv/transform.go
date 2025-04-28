package mpbgv

import (
	"fmt"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/multiparty"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/schemes/bgv"
	"github.com/tuneinsight/lattigo/v6/utils/sampling"
)

// MaskedTransformProtocol is a struct storing the parameters for the [MaskedTransformProtocol] protocol.
type MaskedTransformProtocol struct {
	e2s EncToShareProtocol
	s2e ShareToEncProtocol
}

// MaskedTransformFunc is a struct containing a user-defined in-place function that can be applied to masked integer plaintexts, as a part of the
// Masked Transform Protocol.
// The function is called with a vector of integers modulo bgv.Parameters.PlaintextModulus() of size bgv.Parameters.N() as input, and must write
// its output on the same buffer.
// Transform can be the identity.
//
//   - Decode: if true, then the masked BFV plaintext will be decoded before applying Transform.
//   - Recode: if true, then the masked BFV plaintext will be recoded after applying Transform.
//
// Decode (true/false) -> Transform -> Recode (true/false).
type MaskedTransformFunc struct {
	Decode bool
	Func   func(coeffs []uint64)
	Encode bool
}

// NewMaskedTransformProtocol creates a new instance of the PermuteProtocol.
func NewMaskedTransformProtocol(paramsIn, paramsOut bgv.Parameters, noiseFlooding ring.DistributionParameters) (rfp MaskedTransformProtocol, err error) {

	if paramsIn.N() > paramsOut.N() {
		return MaskedTransformProtocol{}, fmt.Errorf("newMaskedTransformProtocol: paramsIn.N() != paramsOut.N()")
	}

	rfp = MaskedTransformProtocol{}
	if rfp.e2s, err = NewEncToShareProtocol(paramsIn, noiseFlooding); err != nil {
		return
	}

	if rfp.s2e, err = NewShareToEncProtocol(paramsOut, noiseFlooding); err != nil {
		return
	}

	return
}

// SampleCRP samples a common random polynomial to be used in the Masked-Transform protocol from the provided
// common reference string.
func (rfp *MaskedTransformProtocol) SampleCRP(level int, crs sampling.PRNG) multiparty.KeySwitchCRP {
	return rfp.s2e.SampleCRP(level, crs)
}

// AllocateShare allocates the shares of the PermuteProtocol
func (rfp MaskedTransformProtocol) AllocateShare(levelDecrypt, levelRecrypt int) multiparty.RefreshShare {
	return multiparty.RefreshShare{EncToShareShare: rfp.e2s.AllocateShare(levelDecrypt), ShareToEncShare: rfp.s2e.AllocateShare(levelRecrypt)}
}

// GenShare generates the shares of the PermuteProtocol.
// ct1 is the degree 1 element of a rlwe.Ciphertext, i.e. rlwe.Ciphertext.Value[1].
func (rfp MaskedTransformProtocol) GenShare(skIn, skOut *rlwe.SecretKey, ct *rlwe.Ciphertext, crs multiparty.KeySwitchCRP, transform *MaskedTransformFunc, shareOut *multiparty.RefreshShare) (err error) {

	if ct.Level() < shareOut.EncToShareShare.Value.Level() {
		return fmt.Errorf("cannot GenShare: ct[1] level must be at least equal to EncToShareShare level")
	}

	if crs.Value.Level() != shareOut.ShareToEncShare.Value.Level() {
		return fmt.Errorf("cannot GenShare: crs level must be equal to ShareToEncShare")
	}

	tmpMask := rfp.e2s.params.RingT().NewPoly()
	rfp.e2s.GenShare(skIn, ct, &multiparty.AdditiveShare{Value: tmpMask}, &shareOut.EncToShareShare)
	mask := tmpMask
	if transform != nil {
		coeffs := make([]uint64, len(mask.Coeffs[0]))

		if transform.Decode {
			if err := rfp.e2s.encoder.DecodeRingT(mask, ct.Scale, coeffs); err != nil {
				return fmt.Errorf("cannot GenShare: %w", err)
			}
		} else {
			copy(coeffs, mask.Coeffs[0])
		}

		transform.Func(coeffs)

		tmpMaskPerm := rfp.e2s.params.RingT().NewPoly()
		if transform.Encode {
			if err := rfp.s2e.encoder.EncodeRingT(coeffs, ct.Scale, tmpMaskPerm); err != nil {
				return fmt.Errorf("cannot GenShare: %w", err)
			}
		} else {
			copy(tmpMaskPerm.Coeffs[0], coeffs)
		}

		mask = tmpMaskPerm
	}

	// Stores the ciphertext metadata on the public share
	shareOut.MetaData = *ct.MetaData

	return rfp.s2e.GenShare(skOut, crs, multiparty.AdditiveShare{Value: mask}, &shareOut.ShareToEncShare)
}

// AggregateShares sums share1 and share2 on shareOut.
func (rfp MaskedTransformProtocol) AggregateShares(share1, share2 multiparty.RefreshShare, shareOut *multiparty.RefreshShare) (err error) {

	if share1.EncToShareShare.Value.Level() != share2.EncToShareShare.Value.Level() || share1.EncToShareShare.Value.Level() != shareOut.EncToShareShare.Value.Level() {
		return fmt.Errorf("cannot AggregateShares: all e2s shares must be at the same level")
	}

	if share1.ShareToEncShare.Value.Level() != share2.ShareToEncShare.Value.Level() || share1.ShareToEncShare.Value.Level() != shareOut.ShareToEncShare.Value.Level() {
		return fmt.Errorf("cannot AggregateShares: all s2e shares must be at the same level")
	}

	rfp.e2s.params.RingQ().AtLevel(share1.EncToShareShare.Value.Level()).Add(share1.EncToShareShare.Value, share2.EncToShareShare.Value, shareOut.EncToShareShare.Value)
	rfp.s2e.params.RingQ().AtLevel(share1.ShareToEncShare.Value.Level()).Add(share1.ShareToEncShare.Value, share2.ShareToEncShare.Value, shareOut.ShareToEncShare.Value)

	return
}

// Transform applies Decrypt, Recode and Recrypt on the input ciphertext.
func (rfp MaskedTransformProtocol) Transform(ct *rlwe.Ciphertext, transform *MaskedTransformFunc, crs multiparty.KeySwitchCRP, share multiparty.RefreshShare, ciphertextOut *rlwe.Ciphertext) (err error) {

	if !ct.MetaData.Equal(&share.MetaData) {
		return fmt.Errorf("cannot Transform: input ct.MetaData != share.MetaData")
	}

	if ct.Level() < share.EncToShareShare.Value.Level() {
		return fmt.Errorf("cannot Transform: input ciphertext level must be at least equal to e2s level")
	}

	maxLevel := crs.Value.Level()

	if maxLevel != share.ShareToEncShare.Value.Level() {
		return fmt.Errorf("cannot Transform: crs level and s2e level must be the same")
	}

	tmpMask := rfp.e2s.params.RingT().NewPoly()
	rfp.e2s.GetShare(nil, share.EncToShareShare, ct, &multiparty.AdditiveShare{Value: tmpMask}) // tmpMask RingT(m - sum M_i)
	mask := tmpMask
	if transform != nil {
		coeffs := make([]uint64, len(mask.Coeffs[0]))

		if transform.Decode {
			if err := rfp.e2s.encoder.DecodeRingT(mask, ciphertextOut.Scale, coeffs); err != nil {
				return fmt.Errorf("cannot Transform: %w", err)
			}
		} else {
			copy(coeffs, mask.Coeffs[0])
		}

		transform.Func(coeffs)

		tmpMaskPerm := rfp.e2s.params.RingT().NewPoly()
		if transform.Encode {
			if err := rfp.s2e.encoder.EncodeRingT(coeffs, ciphertextOut.Scale, tmpMaskPerm); err != nil {
				return fmt.Errorf("cannot Transform: %w", err)
			}
		} else {
			copy(tmpMaskPerm.Coeffs[0], coeffs)
		}

		mask = tmpMaskPerm
	}

	ciphertextOut.Resize(ciphertextOut.Degree(), maxLevel)

	tmpPt := rfp.s2e.params.RingQ().AtLevel(maxLevel).NewPoly()
	rfp.s2e.encoder.RingT2Q(maxLevel, true, mask, tmpPt)
	rfp.s2e.params.RingQ().AtLevel(maxLevel).NTT(tmpPt, tmpPt)
	rfp.s2e.params.RingQ().AtLevel(maxLevel).Add(tmpPt, share.ShareToEncShare.Value, ciphertextOut.Value[0])

	return rfp.s2e.GetEncryption(multiparty.KeySwitchShare{Value: ciphertextOut.Value[0]}, crs, ciphertextOut)
}
