package dbfv

import (
	"github.com/lca1/lattigo-private/ring"
	//"fmt"
)

type CKS struct {
	context *ring.Context

	sigmaSmudging   float64
	gaussianSampler *ring.KYSampler

	deltaSk *ring.Poly
}

func NewCKS(skInput, skOutput *ring.Poly, context *ring.Context, sigmaSmudging float64) *CKS {

	cks := new(CKS)
	cks.context = context

	cks.sigmaSmudging = sigmaSmudging
	cks.gaussianSampler = context.NewKYSampler(sigmaSmudging, int(6*sigmaSmudging))

	cks.deltaSk = cks.context.NewPoly()
	context.Sub(skInput, skOutput, cks.deltaSk)
	return cks
}

func (cks *CKS) KeySwitch(c1 *ring.Poly) *ring.Poly {

	h := c1.CopyNew()

	cks.context.NTT(h, h)
	cks.context.MulCoeffsMontgomery(h, cks.deltaSk, h)
	cks.context.InvNTT(h, h)
	cks.context.Add(h, cks.gaussianSampler.SampleNew(), h)

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
