package dckks

import (
	"math/big"

	"github.com/ldsec/lattigo/v2/ckks"
	"github.com/ldsec/lattigo/v2/drlwe"
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/rlwe"
)

// MaskedTransformProtocol is a struct storing the parameters for the MaskedTransformProtocol protocol.
type MaskedTransformProtocol struct {
	e2s             E2SProtocol
	s2e             S2EProtocol
	ringQ           *ring.Ring
	encoder         ckks.EncoderBigComplex
	tmpMask         *ring.Poly
	tmpMaskPerm     *ring.Poly
	maskBigint      []*big.Int
	maskFloat       []*big.Float
	maskComplex     []*ring.Complex
	gaussianSampler *ring.GaussianSampler
}

type MaskedTransformFunc func(coeffsIn, coeffsOut ckks.Plaintext)

// MaskedTransformShare is a struct storing the decryption and recryption shares.
type MaskedTransformShare struct {
	e2sShare drlwe.CKSShare
	s2eShare drlwe.CKSShare
}

func NewMaskedTransformProtocol(params ckks.Parameters, sigmaSmudging float64) (rfp *MaskedTransformProtocol) {

	rfp = new(MaskedTransformProtocol)
	rfp.e2s = *NewE2SProtocol(params, sigmaSmudging)
	rfp.s2e = *NewS2EProtocol(params, sigmaSmudging)

	rfp.ringQ = rfp.e2s.ringQ

	rfp.tmpMask = rfp.ringQ.NewPoly()
	rfp.tmpMaskPerm = rfp.ringQ.NewPoly()

	return
}

// AllocateShares allocates the shares of the PermuteProtocol
func (rfp *MaskedTransformProtocol) AllocateShares(levelDecrypt, levelRecrypt int) MaskedTransformShare {
	return MaskedTransformShare{*rfp.e2s.AllocateShareAtLevel(levelDecrypt), *rfp.s2e.AllocateShareAtLevel(levelRecrypt)}
}

// GenShares generates the shares of the PermuteProtocol
func (rfp *MaskedTransformProtocol) GenShares(sk *rlwe.SecretKey, nbParties int, ct *ckks.Ciphertext, crs *ring.Poly, transform MaskedTransformFunc, shareOut MaskedTransformShare) {

	if ct.Level() != shareOut.e2sShare.Value.Level() {
		panic("ciphertext level must be equal to e2sShare")
	}

	if crs.Level() != shareOut.s2eShare.Value.Level() {
		panic("crs level must be equal to s2eShare")
	}

	// Generates the enc 2 share share
	rfp.e2s.GenShare(sk, nbParties, ct, rlwe.AdditiveShare{Value: *rfp.tmpMask}, &shareOut.e2sShare)
	mask := rfp.tmpMask
	if transform != nil {
		transform(ckks.Plaintext{Plaintext : &rlwe.Plaintext{Value: mask, IsNTT: true}, Scale:ct.Scale}, ckks.Plaintext{Plaintext: &rlwe.Plaintext{Value:  rfp.tmpMaskPerm, IsNTT: true}, Scale:ct.Scale})
		mask = rfp.tmpMaskPerm
	}

	rfp.s2e.GenShare(sk, crs, rlwe.AdditiveShare{Value: *mask}, &shareOut.s2eShare)
}
