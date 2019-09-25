package dbfv

import (
	"github.com/ldsec/lattigo/ring"
)

// PCKS is the structure storing the parameters for the collective public key-switching.
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

// NewPCKS creates a new PCKS object and will be used to re-encrypt a ciphertext ctx encrypted under a secret-shared key mong j parties under a new
// collective public-key.
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

// GenShareRound3 is the first part of the unique round of the PCKS protocol. Each party computes the following :
//
// [s_i * ctx[0] + u_i * pk[0] + e_0i, u_i * pk[1] + e_1i]
//
// and broadcasts the result to the other j-1 parties.
func (pcks *PCKS) KeySwitch(ct1 *ring.Poly) (h [2]*ring.Poly) {

	// h_0 = pk_0 (NTT)
	h[0] = pcks.pkOutput[0].CopyNew()
	// h_1 = pk_1 (NTT)
	h[1] = pcks.pkOutput[1].CopyNew()

	//u_i
	pcks.ternarySampler.SampleMontgomeryNTT(0.5, pcks.polypool)

	// h_0 = u_i * pk_0 (NTT)
	pcks.context.MulCoeffsMontgomery(h[0], pcks.polypool, h[0])
	// h_1 = u_i * pk_1 (NTT)
	pcks.context.MulCoeffsMontgomery(h[1], pcks.polypool, h[1])

	// h0 = u_i * pk_0 + s_i*c_1 (NTT)
	pcks.context.NTT(ct1, pcks.polypool)
	pcks.context.MulCoeffsMontgomeryAndAdd(pcks.polypool, pcks.skInput, h[0])

	pcks.context.InvNTT(h[0], h[0])
	pcks.context.InvNTT(h[1], h[1])

	// h_0 = InvNTT(s_i*c_1 + u_i * pk_0) + e0
	pcks.gaussianSamplerSmudge.Sample(pcks.polypool)
	pcks.context.Add(h[0], pcks.polypool, h[0])

	// h_1 = InvNTT(u_i * pk_1) + e1
	pcks.gaussianSampler.Sample(pcks.polypool)
	pcks.context.Add(h[1], pcks.polypool, h[1])

	pcks.polypool.Zero()

	return h
}

// GenShareRoundTwo is the second part of the first and unique round of the PCKS protocol. Each party uppon receiving the j-1 elements from the
// other parties computes :
//
// [ctx[0] + sum(s_i * ctx[0] + u_i * pk[0] + e_0i), sum(u_i * pk[1] + e_1i)]
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
