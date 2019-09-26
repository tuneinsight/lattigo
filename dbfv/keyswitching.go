package dbfv

import (
	"github.com/ldsec/lattigo/bfv"
	"github.com/ldsec/lattigo/ring"
)

// cksProtocol is a structure storing the parameters for the collective key-switching protocol.
type cksProtocol struct {
	ringContext     *ring.Context
	gaussianSampler *ring.KYSampler

	tmpNtt   *ring.Poly
	tmpDelta *ring.Poly
}

type cksShare *ring.Poly

// NewCKSProtocol creates a new cksProtocol that will be used to operate a collective key-switching on a ciphertext encrypted under a collective public-key, whose
// secret-shares are distributed among j parties, re-encrypting the ciphertext under an other public-key, whose secret-shares are also known to the
// parties.
func NewCKSProtocol(bfvContext *bfv.BfvContext, sigmaSmudging float64) *cksProtocol {

	cks := new(cksProtocol)
	cks.ringContext = bfvContext.ContextQ()
	cks.gaussianSampler = bfvContext.ContextQ().NewKYSampler(sigmaSmudging, int(6*sigmaSmudging))

	cks.tmpNtt = cks.ringContext.NewPoly()
	cks.tmpDelta = cks.ringContext.NewPoly()
	return cks
}

func (cks *cksProtocol) AllocateShare() cksShare {
	return cks.ringContext.NewPoly()
}

// GenShare is the first and unique round of the cksProtocol protocol. Each party holding a ciphertext ctx encrypted under a collective publick-key musth
// compute the following :
//
// [(skInput_i - skOutput_i) * ctx[0] + e_i]
//
// Each party then broadcast the result of this computation to the other j-1 parties.
func (cks *cksProtocol) GenShare(skInput, skOutput  *ring.Poly, ct *bfv.Ciphertext, shareOut cksShare) {
	cks.ringContext.Sub(skInput, skOutput, cks.tmpDelta)
	cks.GenShareDelta(cks.tmpDelta, ct, shareOut)
}



func (cks *cksProtocol) GenShareDelta(skDelta  *ring.Poly, ct *bfv.Ciphertext, shareOut cksShare) {

	cks.gaussianSampler.Sample(shareOut)

	cks.ringContext.NTT(ct.Value()[1], cks.tmpNtt)
	cks.ringContext.MulCoeffsMontgomery(cks.tmpNtt, skDelta, cks.tmpNtt)
	cks.ringContext.InvNTT(cks.tmpNtt, cks.tmpNtt)

	cks.ringContext.Add(cks.tmpNtt, shareOut, shareOut)
}


// AggregateShares is the second part of the unique round of the cksProtocol protocol. Uppon receiving the j-1 elements each party computes :
//
// [ctx[0] + sum((skInput_i - skOutput_i) * ctx[0] + e_i), ctx[1]]
func (cks *cksProtocol) AggregateShares(share1, share2, shareOut cksShare) {
	cks.ringContext.Add(share1, share2, shareOut)
}


// KeySwitch performs the actual keyswitching operation on a ciphertext ct and put the result in ctOut
func (cks *cksProtocol) KeySwitch(combined cksShare, ct *bfv.Ciphertext, ctOut *bfv.Ciphertext) {
	cks.ringContext.Add(ct.Value()[0], combined, ctOut.Value()[0])
	if ct != ctOut {
		_ = ct.Value()[1].Copy(ctOut.Value()[1])
	}
}