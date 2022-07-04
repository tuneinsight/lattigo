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

// MaskedTransformFunc represents a struct containing user devined in-place function that can be applied to masked BFV plaintexts, as a part of the
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
func (share *MaskedTransformShare) UnmarshalBinary(data []byte) (err error) {
	shareLen := len(data) >> 1
	if err = share.e2sShare.UnmarshalBinary(data[:shareLen]); err != nil {
		return
	}
	if err = share.s2eShare.UnmarshalBinary(data[shareLen:]); err != nil {
		return
	}
	return
}

// NewMaskedTransformProtocol creates a new instance of the PermuteProtocol.
func NewMaskedTransformProtocol(paramsIn, paramsOut bfv.Parameters, sigmaSmudging float64) (rfp *MaskedTransformProtocol) {

	rfp = new(MaskedTransformProtocol)

	rfp.e2s = *NewE2SProtocol(paramsIn, sigmaSmudging)
	rfp.s2e = *NewS2EProtocol(paramsOut, sigmaSmudging)

	rfp.tmpPt = *bfv.NewPlaintext(paramsOut)
	rfp.tmpMask = paramsIn.RingT().NewPoly()
	rfp.tmpMaskPerm = paramsIn.RingT().NewPoly()
	return
}

// SampleCRP samples a common random polynomial to be used in the Masked-Transform protocol from the provided
// common reference string.
func (rfp *MaskedTransformProtocol) SampleCRP(level int, crs utils.PRNG) drlwe.CKSCRP {
	return rfp.s2e.SampleCRP(level, crs)
}

// AllocateShare allocates the shares of the PermuteProtocol.
func (rfp *MaskedTransformProtocol) AllocateShare() *MaskedTransformShare {
	return &MaskedTransformShare{*rfp.e2s.AllocateShare(), *rfp.s2e.AllocateShare()}
}

// GenShare generates the shares of the PermuteProtocol.
// ct1 is the degree 1 element of a bfv.Ciphertext, i.e. bfv.Ciphertext.Value[1].
func (rfp *MaskedTransformProtocol) GenShare(sk *rlwe.SecretKey, c1 *ring.Poly, crs drlwe.CKSCRP, transform *MaskedTransformFunc, shareOut *MaskedTransformShare) {
	rfp.e2s.GenShare(sk, c1, &rlwe.AdditiveShare{Value: *rfp.tmpMask}, &shareOut.e2sShare)
	mask := rfp.tmpMask
	if transform != nil {
		coeffs := make([]uint64, rfp.e2s.params.N())
		ecd := rfp.e2s.encoder
		ptT := &bfv.PlaintextRingT{Plaintext: &rlwe.Plaintext{Value: mask}}

		if transform.Decode {
			ecd.DecodeRingT(ptT, coeffs)
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

	rfp.s2e.GenShare(skOut, crs, &rlwe.AdditiveShare{Value: *mask}, &shareOut.s2eShare)
}

// AggregateShare sums share1 and share2 on shareOut.
func (rfp *MaskedTransformProtocol) AggregateShare(share1, share2, shareOut *MaskedTransformShare) {
	rfp.e2s.params.RingQ().Add(share1.e2sShare.Value, share2.e2sShare.Value, shareOut.e2sShare.Value)
	rfp.s2e.params.RingQ().Add(share1.s2eShare.Value, share2.s2eShare.Value, shareOut.s2eShare.Value)
}

// Transform applies Decrypt, Recode and Recrypt on the input ciphertext.
func (rfp *MaskedTransformProtocol) Transform(ciphertext *bfv.Ciphertext, transform *MaskedTransformFunc, crs drlwe.CKSCRP, share *MaskedTransformShare, ciphertextOut *bfv.Ciphertext) {
	rfp.e2s.GetShare(nil, &share.e2sShare, ciphertext, &rlwe.AdditiveShare{Value: *rfp.tmpMask}) // tmpMask RingT(m - sum M_i)

	mask := rfp.tmpMask

	if transform != nil {
		coeffs := make([]uint64, rfp.e2s.params.N())
		ecd := rfp.e2s.encoder
		ptT := &bfv.PlaintextRingT{Plaintext: &rlwe.Plaintext{Value: mask}}

		if transform.Decode {
			ecd.DecodeRingT(ptT, coeffs)
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
	rfp.s2e.encoder.ScaleUp(&bfv.PlaintextRingT{Plaintext: &rlwe.Plaintext{Value: mask}}, &rfp.tmpPt)
	rfp.s2e.params.RingQ().Add(rfp.tmpPt.Value, share.s2eShare.Value, ciphertextOut.Value[0])
	rfp.s2e.GetEncryption(&drlwe.CKSShare{Value: ciphertextOut.Value[0]}, crs, ciphertextOut)
}
