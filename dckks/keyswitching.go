package dckks

import (
	"github.com/ldsec/lattigo/ckks"
	"github.com/ldsec/lattigo/ring"
)

// CKS is a structure storing the parameters for the collective key-switching protocol.
type CKS struct {
	ckksContext *ckks.CkksContext

	sigmaSmudging         float64
	gaussianSamplerSmudge *ring.KYSampler

	deltaSk *ring.Poly

	polypool *ring.Poly

	baseconverter *ring.FastBasisExtender
}

// NewCKS creates a new CKS that will be used to operate a collective key-switching on a ciphertext encrypted under a collective public-key, whose
// secret-shares are distributed among j parties, re-encrypting the ciphertext under an other public-key, whose secret-shares are also known to the
// parties.
func NewCKS(skInput, skOutput *ring.Poly, ckksContext *ckks.CkksContext, sigmaSmudging float64) *CKS {

	cks := new(CKS)

	cks.ckksContext = ckksContext

	cks.sigmaSmudging = sigmaSmudging
	cks.gaussianSamplerSmudge = ckksContext.ContextKeys().NewKYSampler(sigmaSmudging, int(6*sigmaSmudging))

	cks.deltaSk = ckksContext.ContextKeys().NewPoly()
	ckksContext.ContextKeys().Sub(skInput, skOutput, cks.deltaSk)

	cks.polypool = ckksContext.ContextKeys().NewPoly()

	cks.baseconverter = ring.NewFastBasisExtender(ckksContext.ContextQ().Modulus, ckksContext.KeySwitchPrimes())

	return cks
}

// KeySwitch is the first and unique round of the CKS protocol. Each party holding a ciphertext ctx encrypted under a collective publick-key musth
// compute the following :
//
// [(P*(skInput_i - skOutput_i) * ctx[0] + e_i)/P]
//
// Each party then broadcast the result of this computation to the other j-1 parties.
func (cks *CKS) KeySwitch(c1 *ring.Poly) (h *ring.Poly) {

	contextQ := cks.ckksContext.ContextQ()
	contextP := cks.ckksContext.ContextP()
	contextKeys := cks.ckksContext.ContextKeys()

	level := uint64(len(c1.Coeffs) - 1)

	h = contextKeys.NewPoly()
	contextQ.CopyLvl(level, c1, h)

	contextQ.MulCoeffsMontgomeryLvl(level, h, cks.deltaSk, h)

	for _, pj := range cks.ckksContext.KeySwitchPrimes() {
		contextQ.MulScalarLvl(level, h, pj, h)
	}

	cks.gaussianSamplerSmudge.SampleNTT(cks.polypool)
	contextQ.AddLvl(level, h, cks.polypool, h)

	hP := contextP.NewPoly()

	for x, i := 0, uint64(len(contextQ.Modulus)); i < uint64(len(contextKeys.Modulus)); x, i = x+1, i+1 {
		for j := uint64(0); j < contextKeys.N; j++ {
			hP.Coeffs[x][j] += cks.polypool.Coeffs[i][j]
		}
	}

	cks.baseconverter.ModDownSplitedNTT(contextQ, contextP, cks.ckksContext.RescaleParamsKeys(), level, h, hP, h, cks.polypool)

	h.Coeffs = h.Coeffs[:level+1]

	cks.polypool.Zero()

	return h
}

// Aggregate is the second part of the unique round of the CKS protocol. Uppon receiving the j-1 elements each party computes :
//
// [ctx[0] + sum((P*(skInput_i - skOutput_i)*ctx[0] + e_i)/P), ctx[1]]
func (cks *CKS) Aggregate(c0 *ring.Poly, h []*ring.Poly) {

	contextQ := cks.ckksContext.ContextQ()
	level := uint64(len(c0.Coeffs) - 1)

	for i := range h {
		contextQ.AddNoModLvl(level, c0, h[i], c0)

		if i&7 == 1 {
			contextQ.ReduceLvl(level, c0, c0)
		}
	}

	if len(h)&7 != 7 {
		contextQ.ReduceLvl(level, c0, c0)
	}
}
