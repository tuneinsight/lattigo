package dbfv

import (
	"github.com/ldsec/lattigo/bfv"
	"github.com/ldsec/lattigo/ring"
)

// CKSProtocol is a structure storing the parameters for the collective key-switching protocol.
type CKSProtocol struct {

	contextCiphertexts *ring.Context
	contextKeys        *ring.Context

	sigmaSmudging         float64
	gaussianSampler *ring.KYSampler

	tmpNtt   *ring.Poly
	tmpDelta *ring.Poly

	baseconverter *bfv.FastBasisExtender
	rescaleParamsKeys []uint64
	keyswitchprimes   []uint64
	rescalepool       []uint64
}

type CKSShare *ring.Poly

// NewCKSProtocol creates a new CKSProtocol that will be used to operate a collective key-switching on a ciphertext encrypted under a collective public-key, whose
// secret-shares are distributed among j parties, re-encrypting the ciphertext under an other public-key, whose secret-shares are also known to the
// parties.
func NewCKSProtocol(bfvContext *bfv.BfvContext, sigmaSmudging float64) *CKSProtocol {

	cks := new(CKSProtocol)

	cks.contextCiphertexts = bfvContext.ContextQ()
	cks.contextKeys = bfvContext.ContextKeys()

	cks.gaussianSampler = bfvContext.ContextKeys().NewKYSampler(sigmaSmudging, int(6*sigmaSmudging))

	cks.tmpNtt = cks.contextKeys.NewPoly()
	cks.tmpDelta = cks.contextKeys.NewPoly()

	cks.keyswitchprimes = make([]uint64, len(bfvContext.KeySwitchPrimes()))
	for i, pi := range bfvContext.KeySwitchPrimes() {
		cks.keyswitchprimes[i] = pi
	}

	cks.baseconverter = bfv.NewFastBasisExtender(cks.contextCiphertexts.Modulus, cks.keyswitchprimes)

	cks.rescaleParamsKeys = make([]uint64, len(cks.contextCiphertexts.Modulus))

	PBig := ring.NewUint(1)
	for _, pj := range cks.keyswitchprimes {
		PBig.Mul(PBig, ring.NewUint(pj))
	}

	tmp := ring.NewUint(0)
	bredParams := cks.contextCiphertexts.GetBredParams()
	for i, Qi := range cks.contextCiphertexts.Modulus {
		tmp.Mod(PBig, ring.NewUint(Qi))
		cks.rescaleParamsKeys[i] = ring.MForm(ring.ModExp(ring.BRedAdd(tmp.Uint64(), Qi, bredParams[i]), Qi-2, Qi), Qi, bredParams[i])
	}

	return cks
}

func (cks *CKSProtocol) AllocateShare() CKSShare {
	return cks.contextKeys.NewPoly()
}

// GenShare is the first and unique round of the CKSProtocol protocol. Each party holding a ciphertext ctx encrypted under a collective publick-key musth
// compute the following :
//
// [(skInput_i - skOutput_i) * ctx[0] + e_i]
//
// Each party then broadcast the result of this computation to the other j-1 parties.
func (cks *CKSProtocol) GenShare(skInput, skOutput *ring.Poly, ct *bfv.Ciphertext, shareOut CKSShare) {

	cks.contextKeys.Sub(skInput, skOutput, cks.tmpDelta)

	cks.GenShareDelta(cks.tmpDelta, ct, shareOut)
}

func (cks *CKSProtocol) GenShareDelta(skDelta *ring.Poly, ct *bfv.Ciphertext, shareOut CKSShare) {

	level := uint64(len(ct.Value()[1].Coeffs) - 1)

	cks.contextCiphertexts.Copy(ct.Value()[1], cks.tmpNtt)
	cks.baseconverter.ModUp(level, cks.tmpNtt, cks.tmpNtt)

	cks.contextKeys.NTT(cks.tmpNtt, cks.tmpNtt)
	cks.contextKeys.MulCoeffsMontgomery(cks.tmpNtt, skDelta, shareOut)

	for _, pj := range cks.keyswitchprimes {
		cks.contextKeys.MulScalar(shareOut, pj, shareOut)
	}

	cks.contextCiphertexts.InvNTT(shareOut, shareOut)

	cks.gaussianSampler.Sample(cks.tmpNtt)
	cks.contextKeys.Add(shareOut, cks.tmpNtt, shareOut)

	cks.baseconverter.ModDown(cks.contextKeys, cks.rescaleParamsKeys, level, shareOut, shareOut, cks.tmpNtt)



	shareOut.Coeffs = shareOut.Coeffs[:level+1]
}

// AggregateShares is the second part of the unique round of the CKSProtocol protocol. Uppon receiving the j-1 elements each party computes :
//
// [ctx[0] + sum((skInput_i - skOutput_i) * ctx[0] + e_i), ctx[1]]
func (cks *CKSProtocol) AggregateShares(share1, share2, shareOut CKSShare) {
	cks.contextCiphertexts.Add(share1, share2, shareOut)
}

// KeySwitch performs the actual keyswitching operation on a ciphertext ct and put the result in ctOut
func (cks *CKSProtocol) KeySwitch(combined CKSShare, ct *bfv.Ciphertext, ctOut *bfv.Ciphertext) {
	cks.contextCiphertexts.Add(ct.Value()[0], combined, ctOut.Value()[0])
	cks.contextCiphertexts.Copy(ct.Value()[1], ctOut.Value()[1])
}
