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

func NewRefreshProtocol(bfvContext *bfv.BfvContext) (refreshProtocol *RefreshProtocol) {
	refreshProtocol = new(RefreshProtocol)
	refreshProtocol.bfvContext = bfvContext
	refreshProtocol.tmp = bfvContext.ContextQ().NewPoly()
	return
}

func (refreshProtocol *RefreshProtocol) AllocateShares() (RefreshShareDecrypt, RefreshShareRecrypt) {
	return refreshProtocol.bfvContext.ContextQ().NewPoly(), refreshProtocol.bfvContext.ContextQ().NewPoly()
}

func (refreshProtocol *RefreshProtocol) GenShares(sk *ring.Poly, ciphertext *bfv.Ciphertext, crs *ring.Poly, shareDecrypt RefreshShareDecrypt, shareRecrypt RefreshShareRecrypt) {

	contextQ := refreshProtocol.bfvContext.ContextQ()
	contextT := refreshProtocol.bfvContext.ContextT()
	sampler := refreshProtocol.bfvContext.ContextQ().NewKYSampler(3.19, 19) // TODO : add smudging noise

	// h0 = s*ct[1]
	contextQ.NTT(ciphertext.Value()[1], refreshProtocol.tmp)
	contextQ.MulCoeffsMontgomeryAndAdd(sk, refreshProtocol.tmp, shareDecrypt)

	// h1 = -s*a
	contextQ.NTT(crs, refreshProtocol.tmp)
	contextQ.MulCoeffsMontgomeryAndSub(sk, refreshProtocol.tmp, shareRecrypt)

	contextQ.InvNTT(shareDecrypt, shareDecrypt)
	contextQ.InvNTT(shareRecrypt, shareRecrypt)

	// h0 = s*ct[1] + e
	sampler.Sample(refreshProtocol.tmp)
	contextQ.Add(shareDecrypt, refreshProtocol.tmp, shareDecrypt)

	// h1 = s*a + e'
	sampler.Sample(refreshProtocol.tmp)
	contextQ.Add(shareRecrypt, refreshProtocol.tmp, shareRecrypt)

	// mask = (uniform plaintext in [0, T-1]) * floor(Q/T)
	coeffs := contextT.NewUniformPoly()
	lift(coeffs, refreshProtocol.tmp, refreshProtocol.bfvContext)

	// h0 = s*ct[1] + mask
	contextQ.Add(shareDecrypt, refreshProtocol.tmp, shareDecrypt)

	// h0 = -s*a - mask
	contextQ.Sub(shareRecrypt, refreshProtocol.tmp, shareRecrypt)

	return
}

func (refreshProtocol *RefreshProtocol) Aggregate(share1, share2, shareOut *ring.Poly) {
	refreshProtocol.bfvContext.ContextQ().Add(share1, share2, shareOut)
}

func (refreshProtocol *RefreshProtocol) Decrypt(ciphertext *bfv.Ciphertext, shareDecrypt RefreshShareDecrypt) {
	refreshProtocol.bfvContext.ContextQ().Add(ciphertext.Value()[0], shareDecrypt, ciphertext.Value()[0])
}

func (refreshProtocol *RefreshProtocol) Recode(ciphertext *bfv.Ciphertext) {
	scaler := ring.NewSimpleScaler(refreshProtocol.bfvContext.T(), refreshProtocol.bfvContext.ContextQ())

	scaler.Scale(ciphertext.Value()[0], ciphertext.Value()[0])
	lift(ciphertext.Value()[0], ciphertext.Value()[0], refreshProtocol.bfvContext)
}

func (refreshProtocol *RefreshProtocol) Recrypt(ciphertext *bfv.Ciphertext, crs *ring.Poly, shareRecrypt RefreshShareRecrypt) {

	refreshProtocol.bfvContext.ContextQ().Add(ciphertext.Value()[0], shareRecrypt, ciphertext.Value()[0])

	ciphertext.Value()[1] = crs.CopyNew()
}

func lift(p0, p1 *ring.Poly, bfvcontext *bfv.BfvContext) {
	for j := uint64(0); j < bfvcontext.N(); j++ {
		for i := len(bfvcontext.ContextQ().Modulus) - 1; i >= 0; i-- {
			p1.Coeffs[i][j] = ring.MRed(p0.Coeffs[0][j], bfvcontext.DeltaMont()[i], bfvcontext.ContextQ().Modulus[i], bfvcontext.ContextQ().GetMredParams()[i])
		}
	}
}
