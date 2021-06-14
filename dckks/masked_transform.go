package dckks

import (
	"math/big"

	"github.com/ldsec/lattigo/v2/utils"
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
	tmpMask         *ring.Poly
	tmpMaskPerm     *ring.Poly
	vBigint      []*big.Int
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

	rfp.vBigint = make([]*big.Int, rfp.ringQ.N)

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

// Aggregate sums share1 and share2 on shareOut.
func (rfp *MaskedTransformProtocol) Aggregate(share1, share2, shareOut MaskedTransformShare) {

	if share1.e2sShare.Value.Level() != share2.e2sShare.Value.Level() || share1.e2sShare.Value.Level() != shareOut.e2sShare.Value.Level(){
		panic("all e2s shares must be at the same level")
	}

	if share1.s2eShare.Value.Level() != share2.s2eShare.Value.Level() || share1.s2eShare.Value.Level() != shareOut.s2eShare.Value.Level(){
		panic("all s2e shares must be at the same level")
	}

	rfp.ringQ.AddLvl(share1.e2sShare.Value.Level(), share1.e2sShare.Value, share2.e2sShare.Value, shareOut.e2sShare.Value)
	rfp.ringQ.AddLvl(share1.s2eShare.Value.Level(), share1.s2eShare.Value, share2.s2eShare.Value, shareOut.s2eShare.Value)
}

// Transform applies Decrypt, Recode and Recrypt on the input ciphertext.
func (rfp *MaskedTransformProtocol) Transform(ct *ckks.Ciphertext, transform MaskedTransformFunc, crs *ring.Poly, share MaskedTransformShare, ciphertextOut *ckks.Ciphertext) {

	if ct.Level() != share.e2sShare.Value.Level(){
		panic("ciphertext level and e2s level must be the same")
	}

	if crs.Level() != share.s2eShare.Value.Level(){
		panic("crs level and s2e level must be the same")
	}

	// Adds the aggregated decryption shares on the ciphertext (decrypts to LSSS)
	rfp.e2s.GetShare(nil, &share.e2sShare, ct, &rlwe.AdditiveShare{Value: *rfp.tmpMask})

	// Applies the transform on the LSSS
	mask := rfp.tmpMask
	tmp := rfp.tmpMaskPerm
	if transform != nil {
		transform(ckks.Plaintext{Plaintext : &rlwe.Plaintext{Value: mask, IsNTT: true}, Scale:ct.Scale}, ckks.Plaintext{Plaintext: &rlwe.Plaintext{Value:  rfp.tmpMaskPerm, IsNTT: true}, Scale:ct.Scale})
		mask, tmp = rfp.tmpMaskPerm, rfp.tmpMask
	}

	minLevel := utils.MinInt(crs.Level(), ct.Level())
	maxLevel := utils.MaxInt(crs.Level(), ct.Level())

	// Extends the level of the transformed LSSS to the level of the recryption shares
	centerAndExtendBasisLargeNorm(minLevel, maxLevel, rfp.ringQ, mask, rfp.vBigint, tmp)
	
	// Adds the aggregated recryption shares on the LSSS
	rfp.ringQ.AddLvl(maxLevel, tmp, share.s2eShare.Value, ciphertextOut.Value[0])

	// Copies the result on the out ciphertext
	rfp.s2e.GetEncryption(&drlwe.CKSShare{Value: ciphertextOut.Value[0]}, crs, *ciphertextOut)
}