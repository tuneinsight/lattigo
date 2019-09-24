package dckks

import (
	"github.com/ldsec/lattigo/ring"
)

type PCKS struct {
	context *ring.Context

	sigmaSmudging         float64
	gaussianSamplerSmudge *ring.KYSampler
	gaussianSampler       *ring.KYSampler
	ternarySampler        *ring.TernarySampler

	skInput  *ring.Poly
	pkOutput [2]*ring.Poly

	polypool *ring.Poly
}

func NewPCKS(skInput *ring.Poly, pkOutput [2]*ring.Poly, context *ring.Context, sigmaSmudging float64) *PCKS {

	pcks := new(PCKS)
	pcks.context = context
	pcks.sigmaSmudging = sigmaSmudging

	pcks.gaussianSamplerSmudge = context.NewKYSampler(sigmaSmudging, int(6*sigmaSmudging))
	pcks.gaussianSampler = context.NewKYSampler(3.19, 19)
	pcks.ternarySampler = context.NewTernarySampler()

	pcks.skInput = skInput
	pcks.pkOutput = pkOutput

	pcks.polypool = context.NewPoly()

	return pcks
}

func (pcks *PCKS) KeySwitch(ct1 *ring.Poly) (h [2]*ring.Poly) {

	// h_0 = pk_0 (NTT)
	h[0] = pcks.pkOutput[0].CopyNew()
	// h_1 = pk_1 (NTT)
	h[1] = pcks.pkOutput[1].CopyNew()

	//u_i
	pcks.ternarySampler.SampleMontgomeryNTT(pcks.polypool)

	// h_0 = u_i * pk_0 (NTT)
	pcks.context.MulCoeffsMontgomery(h[0], pcks.polypool, h[0])
	// h_1 = u_i * pk_1 (NTT)
	pcks.context.MulCoeffsMontgomery(h[1], pcks.polypool, h[1])

	// h0 = u_i * pk_0 + s_i*c_1 (NTT)
	pcks.context.MulCoeffsMontgomeryAndAdd(ct1, pcks.skInput, h[0])

	// h_0 = InvNTT(s_i*c_1 + u_i * pk_0) + e0
	pcks.gaussianSamplerSmudge.SampleNTT(pcks.polypool)
	pcks.context.Add(h[0], pcks.polypool, h[0])

	// h_1 = InvNTT(u_i * pk_1) + e1
	pcks.gaussianSampler.SampleNTT(pcks.polypool)
	pcks.context.Add(h[1], pcks.polypool, h[1])

	return h
}

func (pcks *PCKS) Aggregate(ct []*ring.Poly, h [][2]*ring.Poly) {

	pcks.context.AddNoMod(ct[0], h[0][0], ct[0])
	ct[1] = h[0][1].CopyNew()

	for i := 1; i < len(h); i++ {
		pcks.context.AddNoMod(ct[0], h[i][0], ct[0])
		pcks.context.AddNoMod(ct[1], h[i][1], ct[1])

		if i&7 == 7 {
			pcks.context.Reduce(ct[0], ct[0])
			pcks.context.Reduce(ct[1], ct[1])
		}
	}

	if len(h)&7 != 7 {
		pcks.context.Reduce(ct[0], ct[0])
		pcks.context.Reduce(ct[1], ct[1])
	}

}
