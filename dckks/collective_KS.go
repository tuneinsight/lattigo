package dckks

import (
	"github.com/ldsec/lattigo/ring"
)

type CKS struct {
	context *ring.Context

	sigmaSmudging   float64
	gaussianSampler *ring.KYSampler

	deltaSk *ring.Poly

	polypool *ring.Poly
}

func NewCKS(skInput, skOutput *ring.Poly, context *ring.Context, sigmaSmudging float64) *CKS {

	cks := new(CKS)
	cks.context = context

	cks.sigmaSmudging = sigmaSmudging
	cks.gaussianSampler = context.NewKYSampler(sigmaSmudging, int(6*sigmaSmudging))

	cks.deltaSk = cks.context.NewPoly()
	context.Sub(skInput, skOutput, cks.deltaSk)
	cks.polypool = context.NewPoly()
	return cks
}

func (cks *CKS) KeySwitch(c1 *ring.Poly) *ring.Poly {

	h := c1.CopyNew()

	cks.context.MulCoeffsMontgomery(h, cks.deltaSk, h)
	cks.gaussianSampler.SampleNTT(cks.polypool)
	cks.context.Add(h, cks.polypool, h)
	cks.polypool.Zero()

	return h
}

func (cks *CKS) Aggregate(c0 *ring.Poly, h []*ring.Poly) {

	for i := range h {
		cks.context.AddNoMod(c0, h[i], c0)

		if i&7 == 1 {
			cks.context.Reduce(c0, c0)
		}
	}

	if len(h)&7 != 7 {
		cks.context.Reduce(c0, c0)
	}
}
