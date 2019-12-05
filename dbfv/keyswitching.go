package dbfv

import (
	"github.com/ldsec/lattigo/bfv"
	"github.com/ldsec/lattigo/ring"
)

// CKSProtocol is a structure storing the parameters for the collective key-switching protocol.
type CKSProtocol struct {
	context *dbfvContext

	sigmaSmudging         float64
	gaussianSamplerSmudge *ring.KYSampler

	tmpNtt   *ring.Poly
	tmpDelta *ring.Poly
	hP       *ring.Poly

	baseconverter *ring.FastBasisExtender
}

// CKSShare is a type for the CKS protocol shares.
type CKSShare struct {
	*ring.Poly
}

// UnmarshalBinary decodes a previouls marshaled share on the target share.
func (share *CKSShare) UnmarshalBinary(data []byte) error {
	share.Poly = new(ring.Poly)
	err := share.Poly.UnmarshalBinary(data)
	return err

}

// NewCKSProtocol creates a new CKSProtocol that will be used to operate a collective key-switching on a ciphertext encrypted under a collective public-key, whose
// secret-shares are distributed among j parties, re-encrypting the ciphertext under another public-key, whose secret-shares are also known to the
// parties.
func NewCKSProtocol(params *bfv.Parameters, sigmaSmudging float64) *CKSProtocol {
	context := newDbfvContext(params)

	cks := new(CKSProtocol)

	cks.context = context

	cks.gaussianSamplerSmudge = context.contextQP.NewKYSampler(sigmaSmudging, int(6*sigmaSmudging))

	cks.tmpNtt = cks.context.contextQP.NewPoly()
	cks.tmpDelta = cks.context.contextQ.NewPoly()
	cks.hP = cks.context.contextP.NewPoly()

	cks.baseconverter = ring.NewFastBasisExtender(cks.context.contextQ, cks.context.contextP)

	return cks
}

// AllocateShare allocates the shares of the CKSProtocol
func (cks *CKSProtocol) AllocateShare() CKSShare {

	return CKSShare{cks.context.contextQ.NewPoly()}

}

// GenShare is the first and unique round of the CKSProtocol protocol. Each party holding a ciphertext ctx encrypted under a collective publick-key musth
// compute the following :
//
// [(skInput_i - skOutput_i) * ctx[0] + e_i]
//
// Each party then broadcast the result of this computation to the other j-1 parties.
func (cks *CKSProtocol) GenShare(skInput, skOutput *ring.Poly, ct *bfv.Ciphertext, shareOut CKSShare) {

	cks.context.contextQ.Sub(skInput, skOutput, cks.tmpDelta)

	cks.genShareDelta(cks.tmpDelta, ct, shareOut)
}

func (cks *CKSProtocol) genShareDelta(skDelta *ring.Poly, ct *bfv.Ciphertext, shareOut CKSShare) {

	level := uint64(len(ct.Value()[1].Coeffs) - 1)

	contextQ := cks.context.contextQ
	contextP := cks.context.contextP

	contextQ.NTT(ct.Value()[1], cks.tmpNtt)
	contextQ.MulCoeffsMontgomery(cks.tmpNtt, skDelta, shareOut.Poly)
	contextQ.MulScalarBigint(shareOut.Poly, contextP.ModulusBigint, shareOut.Poly)

	contextQ.InvNTT(shareOut.Poly, shareOut.Poly)

	cks.gaussianSamplerSmudge.Sample(cks.tmpNtt)
	contextQ.Add(shareOut.Poly, cks.tmpNtt, shareOut.Poly)

	for x, i := 0, uint64(len(contextQ.Modulus)); i < uint64(len(cks.context.contextQP.Modulus)); x, i = x+1, i+1 {
		tmphP := cks.hP.Coeffs[x]
		tmpNTT := cks.tmpNtt.Coeffs[i]
		for j := uint64(0); j < contextQ.N; j++ {
			tmphP[j] += tmpNTT[j]
		}
	}

	cks.baseconverter.ModDownSplitedPQ(level, shareOut.Poly, cks.hP, shareOut.Poly)

	cks.tmpNtt.Zero()
	cks.hP.Zero()
}

// AggregateShares is the second part of the unique round of the CKSProtocol protocol. Upon receiving the j-1 elements each party computes :
//
// [ctx[0] + sum((skInput_i - skOutput_i) * ctx[0] + e_i), ctx[1]]
func (cks *CKSProtocol) AggregateShares(share1, share2, shareOut CKSShare) {
	cks.context.contextQ.Add(share1.Poly, share2.Poly, shareOut.Poly)
}

// KeySwitch performs the actual keyswitching operation on a ciphertext ct and put the result in ctOut
func (cks *CKSProtocol) KeySwitch(combined CKSShare, ct *bfv.Ciphertext, ctOut *bfv.Ciphertext) {
	cks.context.contextQ.Add(ct.Value()[0], combined.Poly, ctOut.Value()[0])
	cks.context.contextQ.Copy(ct.Value()[1], ctOut.Value()[1])
}
