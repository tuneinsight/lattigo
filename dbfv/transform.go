package dbfv

import (
	"fmt"

	"github.com/tuneinsight/lattigo/v4/bfv"
	"github.com/tuneinsight/lattigo/v4/drlwe"
	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/rlwe"
	"github.com/tuneinsight/lattigo/v4/utils/sampling"
)

// MaskedTransformProtocol is a struct storing the parameters for the MaskedTransformProtocol protocol.
type MaskedTransformProtocol struct {
	e2s E2SProtocol
	s2e S2EProtocol

	tmpPt       *rlwe.Plaintext
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
		tmpPt:       bfv.NewPlaintext(params, params.MaxLevel()),
		tmpMask:     params.RingT().NewPoly(),
		tmpMaskPerm: params.RingT().NewPoly(),
	}
}

// MaskedTransformFunc is struct containing user defined in-place function that can be applied to masked BFV plaintexts, as a part of the
// Masked Transform Protocol.
// Transform is a function called with a vector of integers modulo bfv.Parameters.T() of size bfv.Parameters.N() as input, and must write
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
func NewMaskedTransformProtocol(paramsIn, paramsOut bfv.Parameters, sigmaSmudging float64) (rfp *MaskedTransformProtocol, err error) {

	if paramsIn.N() > paramsOut.N() {
		return nil, fmt.Errorf("newMaskedTransformProtocol: paramsIn.N() != paramsOut.N()")
	}

	rfp = new(MaskedTransformProtocol)

	rfp.e2s = *NewE2SProtocol(paramsIn, sigmaSmudging)
	rfp.s2e = *NewS2EProtocol(paramsOut, sigmaSmudging)

	rfp.tmpPt = bfv.NewPlaintext(paramsOut, paramsOut.MaxLevel())
	rfp.tmpMask = paramsIn.RingT().NewPoly()
	rfp.tmpMaskPerm = paramsIn.RingT().NewPoly()
	return
}

// SampleCRP samples a common random polynomial to be used in the Masked-Transform protocol from the provided
// common reference string.
func (rfp *MaskedTransformProtocol) SampleCRP(level int, crs sampling.PRNG) drlwe.CKSCRP {
	return rfp.s2e.SampleCRP(level, crs)
}

// AllocateShare allocates the shares of the PermuteProtocol.
func (rfp *MaskedTransformProtocol) AllocateShare(levelIn, levelOut int) *drlwe.RefreshShare {
	return &drlwe.RefreshShare{E2SShare: *rfp.e2s.AllocateShare(levelIn), S2EShare: *rfp.s2e.AllocateShare(levelOut)}
}

// GenShare generates the shares of the PermuteProtocol.
// ct1 is the degree 1 element of a bfv.Ciphertext, i.e. bfv.Ciphertext.Value[1].
func (rfp *MaskedTransformProtocol) GenShare(skIn, skOut *rlwe.SecretKey, ct *rlwe.Ciphertext, crs drlwe.CKSCRP, transform *MaskedTransformFunc, shareOut *drlwe.RefreshShare) {

	rfp.e2s.GenShare(skIn, ct, &drlwe.AdditiveShare{Value: *rfp.tmpMask}, &shareOut.E2SShare)

	mask := rfp.tmpMask
	if transform != nil {
		coeffs := make([]uint64, rfp.e2s.params.N())
		ecd := rfp.e2s.encoder
		ptT := &bfv.PlaintextRingT{Plaintext: &rlwe.Plaintext{Value: mask}}

		if transform.Decode {
			ecd.Decode(ptT, coeffs)
		} else {
			copy(coeffs, ptT.Value.Coeffs[0])
		}

		transform.Func(coeffs)

		if transform.Encode {
			ecd.EncodeRingT(coeffs, &bfv.PlaintextRingT{Plaintext: &rlwe.Plaintext{Value: rfp.tmpMaskPerm}})
		} else {
			copy(rfp.tmpMaskPerm.Coeffs[0], coeffs)
		}

		mask = rfp.tmpMaskPerm
	}

	rfp.s2e.GenShare(skOut, crs, &drlwe.AdditiveShare{Value: *mask}, &shareOut.S2EShare)
}

// AggregateShares sums share1 and share2 on shareOut.
func (rfp *MaskedTransformProtocol) AggregateShares(share1, share2, shareOut *drlwe.RefreshShare) {
	rfp.e2s.params.RingQ().Add(share1.E2SShare.Value, share2.E2SShare.Value, shareOut.E2SShare.Value)
	rfp.s2e.params.RingQ().Add(share1.S2EShare.Value, share2.S2EShare.Value, shareOut.S2EShare.Value)
}

// Transform applies Decrypt, Recode and Recrypt on the input ciphertext.
func (rfp *MaskedTransformProtocol) Transform(ciphertext *rlwe.Ciphertext, transform *MaskedTransformFunc, crs drlwe.CKSCRP, share *drlwe.RefreshShare, ciphertextOut *rlwe.Ciphertext) {

	rfp.e2s.GetShare(nil, &share.E2SShare, ciphertext, &drlwe.AdditiveShare{Value: *rfp.tmpMask}) // tmpMask RingT(m - sum M_i)

	mask := rfp.tmpMask

	if transform != nil {
		coeffs := make([]uint64, rfp.e2s.params.N())
		ecd := rfp.e2s.encoder
		ptT := &bfv.PlaintextRingT{Plaintext: &rlwe.Plaintext{Value: mask}}

		if transform.Decode {
			ecd.Decode(ptT, coeffs)
		} else {
			copy(coeffs, ptT.Value.Coeffs[0])
		}

		transform.Func(coeffs)

		if transform.Encode {
			ecd.EncodeRingT(coeffs, &bfv.PlaintextRingT{Plaintext: &rlwe.Plaintext{Value: rfp.tmpMaskPerm}})
		} else {
			copy(rfp.tmpMaskPerm.Coeffs[0], coeffs)
		}

		mask = rfp.tmpMaskPerm
	}

	ciphertextOut.Resize(1, rfp.s2e.params.MaxLevel())
	rfp.s2e.encoder.ScaleUp(&bfv.PlaintextRingT{Plaintext: &rlwe.Plaintext{Value: mask}}, rfp.tmpPt)
	rfp.s2e.params.RingQ().Add(rfp.tmpPt.Value, share.S2EShare.Value, ciphertextOut.Value[0])
	rfp.s2e.GetEncryption(&drlwe.CKSShare{Value: ciphertextOut.Value[0]}, crs, ciphertextOut)
}
