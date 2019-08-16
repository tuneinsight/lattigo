package dckks

import (
	"github.com/lca1/lattigo-private/ckks"
	"github.com/lca1/lattigo-private/ring"
)

type CKG struct {
	context         *ring.Context
	gaussianSampler *ring.KYSampler

	share *ring.Poly
	cpk   [2]*ring.Poly
}

func NewCKG(context *ring.Context, crs *ring.Poly) *CKG {
	ckg := new(CKG)
	ckg.context = context
	ckg.gaussianSampler = context.NewKYSampler(3.19, 19)

	ckg.cpk[0] = context.NewPoly()
	ckg.cpk[1] = crs.CopyNew()

	return ckg
}

func (ckg *CKG) GenShare(sk *ring.Poly) error {

	// -(sk * crs) + e
	ckg.share = ckg.gaussianSampler.SampleNTTNew()
	ckg.context.MulCoeffsMontgomeryAndSub(sk, ckg.cpk[1], ckg.share)

	return nil
}

func (ckg *CKG) GetShare() *ring.Poly {
	return ckg.share
}

func (ckg *CKG) AggregateShares(shares []*ring.Poly) error {

	for i := 0; i < len(shares); i++ {
		ckg.context.Add(ckg.cpk[0], shares[i], ckg.cpk[0])
	}

	return nil
}

func (ckg *CKG) Finalize() (*ckks.PublicKey, error) {
	collectivePk := new(ckks.PublicKey)
	collectivePk.Set(ckg.cpk)
	return collectivePk, nil
}
