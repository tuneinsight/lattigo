package dckks

import (
	"github.com/ldsec/lattigo/ckks"
	"github.com/ldsec/lattigo/ring"
)

// CKS is a structure storing the parameters for the collective key-switching protocol.
type CKS struct {
	contextCiphertexts *ring.Context
	contextKeys        *ring.Context
	keyswitchprimes    []uint64

	sigmaSmudging   float64
	gaussianSampler *ring.KYSampler

	deltaSk *ring.Poly

	polypool *ring.Poly

	baseconverter *ckks.FastBasisExtender

	rescaleParamsKeys []uint64
}

// NewCKS creates a new CKS that will be used to operate a collective key-switching on a ciphertext encrypted under a collective public-key, whose
// secret-shares are distributed among j parties, re-encrypting the ciphertext under an other public-key, whose secret-shares are also known to the
// parties.
func NewCKS(skInput, skOutput *ring.Poly, contextCiphertexts, contextKeys *ring.Context, keyswitchprimes []uint64, sigmaSmudging float64) *CKS {

	cks := new(CKS)

	cks.contextCiphertexts = contextCiphertexts
	cks.contextKeys = contextKeys

	cks.keyswitchprimes = make([]uint64, len(keyswitchprimes))
	for i := range keyswitchprimes {
		cks.keyswitchprimes[i] = keyswitchprimes[i]
	}

	cks.sigmaSmudging = sigmaSmudging
	cks.gaussianSampler = contextKeys.NewKYSampler(sigmaSmudging, int(6*sigmaSmudging))

	cks.deltaSk = cks.contextKeys.NewPoly()
	contextKeys.Sub(skInput, skOutput, cks.deltaSk)

	for _, pj := range cks.keyswitchprimes {
		cks.contextKeys.MulScalar(cks.deltaSk, pj, cks.deltaSk)
	}

	cks.polypool = contextKeys.NewPoly()

	cks.baseconverter = ckks.NewFastBasisExtender(contextCiphertexts.Modulus, cks.keyswitchprimes)

	cks.rescaleParamsKeys = make([]uint64, len(contextCiphertexts.Modulus))

	PBig := ring.NewUint(1)
	for _, pj := range cks.keyswitchprimes {
		PBig.Mul(PBig, ring.NewUint(pj))
	}

	tmp := ring.NewUint(0)
	bredParams := contextCiphertexts.GetBredParams()
	for i, Qi := range contextCiphertexts.Modulus {
		tmp.Mod(PBig, ring.NewUint(Qi))
		cks.rescaleParamsKeys[i] = ring.MForm(ring.ModExp(ring.BRedAdd(tmp.Uint64(), Qi, bredParams[i]), Qi-2, Qi), Qi, bredParams[i])
	}

	return cks
}

// KeySwitch is the first and unique round of the CKS protocol. Each party holding a ciphertext ctx encrypted under a collective publick-key musth
// compute the following :
//
// [(P*(skInput_i - skOutput_i) * ctx[0] + e_i)/P]
//
// Each party then broadcast the result of this computation to the other j-1 parties.
func (cks *CKS) KeySwitch(c1 *ring.Poly) (h *ring.Poly) {

	level := uint64(len(c1.Coeffs) - 1)

	h = cks.contextKeys.NewPoly()
	cks.contextCiphertexts.Copy(c1, h)

	cks.contextCiphertexts.InvNTT(h, h)
	cks.baseconverter.ModUp(level, h, h)
	cks.contextKeys.NTT(h, h)

	cks.contextKeys.MulCoeffsMontgomery(h, cks.deltaSk, h)

	cks.gaussianSampler.SampleNTT(cks.polypool)
	cks.contextKeys.Add(h, cks.polypool, h)

	cks.baseconverter.ModDown(cks.contextKeys, cks.rescaleParamsKeys, level, h, h, cks.polypool)

	h.Coeffs = h.Coeffs[:level+1]

	cks.polypool.Zero()

	return h
}

// Aggregate is the second part of the unique round of the CKS protocol. Uppon receiving the j-1 elements each party computes :
//
// [ctx[0] + sum((P*(skInput_i - skOutput_i)*ctx[0] + e_i)/P), ctx[1]]
func (cks *CKS) Aggregate(c0 *ring.Poly, h []*ring.Poly) {

	for i := range h {
		cks.contextCiphertexts.AddNoMod(c0, h[i], c0)

		if i&7 == 1 {
			cks.contextCiphertexts.Reduce(c0, c0)
		}
	}

	if len(h)&7 != 7 {
		cks.contextCiphertexts.Reduce(c0, c0)
	}
}

func modDownKeyscks(cks *CKS, level uint64, p0, p1 *ring.Poly, rescalepool []uint64) {

	var Qi uint64

	context := cks.contextKeys

	ring.InvNTT(p0.Coeffs[level], p0.Coeffs[level], context.N, context.GetNttPsiInv()[level], context.GetNttNInv()[level], context.Modulus[level], context.GetMredParams()[level])

	for i := uint64(0); i < level; i++ {

		ring.NTT(p0.Coeffs[level], rescalepool, context.N, context.GetNttPsi()[i], context.Modulus[i], context.GetMredParams()[i], context.GetBredParams()[i])

		Qi = cks.contextCiphertexts.Modulus[i]

		for j := uint64(0); j < cks.contextCiphertexts.N; j++ {
			p1.Coeffs[i][j] = p0.Coeffs[i][j] + Qi - ring.BRedAdd(rescalepool[j], Qi, context.GetBredParams()[i])  //  x[i] - x[-1]
			p1.Coeffs[i][j] = ring.MRed(p1.Coeffs[i][j], cks.rescaleParamsKeys[i], Qi, context.GetMredParams()[i]) // (x[i] - x[-1]) * InvQl
		}
	}
}
