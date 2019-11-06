package dbfv

import (
	"github.com/ldsec/lattigo/bfv"
	"github.com/ldsec/lattigo/ring"
	//"fmt"
)

type RefreshProtocol struct {
	bfvContext *bfv.BfvContext
	tmp        *ring.Poly
}

type RefreshShareDecrypt *ring.Poly
type RefreshShareRecrypt *ring.Poly
type RefreshShare struct {
	RefreshShareDecrypt RefreshShareDecrypt
	RefreshShareRecrypt RefreshShareRecrypt
}

func NewRefreshProtocol(bfvContext *bfv.BfvContext) (refreshProtocol *RefreshProtocol) {
	refreshProtocol = new(RefreshProtocol)
	refreshProtocol.bfvContext = bfvContext
	refreshProtocol.tmp = bfvContext.ContextQ().NewPoly()
	return
}

func (rfp *RefreshProtocol) AllocateShares() RefreshShare {
	return RefreshShare{rfp.bfvContext.ContextQ().NewPoly(),
		rfp.bfvContext.ContextQ().NewPoly()}
}

func (rfp *RefreshProtocol) GenShares(sk *ring.Poly, ciphertext *bfv.Ciphertext, crs *ring.Poly, share RefreshShare) {

	contextQ := rfp.bfvContext.ContextQ()
	contextT := rfp.bfvContext.ContextT()
	sampler := rfp.bfvContext.ContextQ().NewKYSampler(3.19, 19) // TODO : add smudging noise

	// h0 = s*ct[1]
	contextQ.NTT(ciphertext.Value()[1], rfp.tmp)
	contextQ.MulCoeffsMontgomeryAndAdd(sk, rfp.tmp, share.RefreshShareDecrypt)

	// h1 = -s*a
	contextQ.NTT(crs, rfp.tmp)
	contextQ.MulCoeffsMontgomeryAndSub(sk, rfp.tmp, share.RefreshShareRecrypt)

	contextQ.InvNTT(share.RefreshShareDecrypt, share.RefreshShareDecrypt)
	contextQ.InvNTT(share.RefreshShareRecrypt, share.RefreshShareRecrypt)

	// h0 = s*ct[1] + e
	sampler.Sample(rfp.tmp)
	contextQ.Add(share.RefreshShareDecrypt, rfp.tmp, share.RefreshShareDecrypt)

	// h1 = s*a + e'
	sampler.Sample(rfp.tmp)
	contextQ.Add(share.RefreshShareRecrypt, rfp.tmp, share.RefreshShareRecrypt)

	// mask = (uniform plaintext in [0, T-1]) * floor(Q/T)
	coeffs := contextT.NewUniformPoly()
	lift(coeffs, rfp.tmp, rfp.bfvContext)

	// h0 = s*ct[1] + mask
	contextQ.Add(share.RefreshShareDecrypt, rfp.tmp, share.RefreshShareDecrypt)

	// h0 = -s*a - mask
	contextQ.Sub(share.RefreshShareRecrypt, rfp.tmp, share.RefreshShareRecrypt)
}

func (rfp *RefreshProtocol) Aggregate(share1, share2, shareOut RefreshShare) {
	rfp.bfvContext.ContextQ().Add(share1.RefreshShareDecrypt, share2.RefreshShareDecrypt, shareOut.RefreshShareDecrypt)
	rfp.bfvContext.ContextQ().Add(share1.RefreshShareRecrypt, share2.RefreshShareRecrypt, shareOut.RefreshShareRecrypt)
}

func (rfp *RefreshProtocol) Decrypt(ciphertext *bfv.Ciphertext, shareDecrypt RefreshShareDecrypt, sharePlaintext *ring.Poly) {
	rfp.bfvContext.ContextQ().Add(ciphertext.Value()[0], shareDecrypt, sharePlaintext)
}

func (rfp *RefreshProtocol) Recode(sharePlaintext *ring.Poly, sharePlaintextOut *ring.Poly) {
	scaler := ring.NewSimpleScaler(rfp.bfvContext.T(), rfp.bfvContext.ContextQ())

	scaler.Scale(sharePlaintext, sharePlaintextOut)
	lift(sharePlaintextOut, sharePlaintextOut, rfp.bfvContext)
}

func (rfp *RefreshProtocol) Recrypt(sharePlaintext *ring.Poly, crs *ring.Poly, shareRecrypt RefreshShareRecrypt, ciphertextOut *bfv.Ciphertext) {

	rfp.bfvContext.ContextQ().Add(sharePlaintext, shareRecrypt, ciphertextOut.Value()[0])

	ciphertextOut.Value()[1].Copy(crs)
}

func (rfp *RefreshProtocol) Finalize(ciphertext *bfv.Ciphertext, crs *ring.Poly, share RefreshShare, ciphertextOut *bfv.Ciphertext) {
	rfp.Decrypt(ciphertext, share.RefreshShareDecrypt, rfp.tmp)
	rfp.Recode(rfp.tmp, rfp.tmp)
	rfp.Recrypt(rfp.tmp, crs, share.RefreshShareRecrypt, ciphertextOut)
}

func lift(p0, p1 *ring.Poly, bfvcontext *bfv.BfvContext) {
	for j := uint64(0); j < bfvcontext.N(); j++ {
		for i := len(bfvcontext.ContextQ().Modulus) - 1; i >= 0; i-- {
			p1.Coeffs[i][j] = ring.MRed(p0.Coeffs[0][j], bfvcontext.DeltaMont()[i], bfvcontext.ContextQ().Modulus[i], bfvcontext.ContextQ().GetMredParams()[i])
		}
	}
}
