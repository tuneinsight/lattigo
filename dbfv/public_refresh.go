package dbfv

import (
	"github.com/ldsec/lattigo/bfv"
	"github.com/ldsec/lattigo/ring"
	//"fmt"
)

type RefreshShares struct {
	h0 *ring.Poly
	h1 *ring.Poly
}

func GenRefreshShares(sk *bfv.SecretKey, ciphertext *bfv.Ciphertext, bfvcontext *bfv.BfvContext, crs *ring.Poly, encoder *bfv.BatchEncoder) (refreshShares *RefreshShares) {

	refreshShares = new(RefreshShares)

	polypool := bfvcontext.ContextQ().NewPoly()

	refreshShares.h0 = bfvcontext.ContextQ().NewPoly()
	refreshShares.h1 = bfvcontext.ContextQ().NewPoly()

	// h0 = s*ct[1]
	bfvcontext.ContextQ().NTT(ciphertext.Value()[1], polypool)
	bfvcontext.ContextQ().MulCoeffsMontgomeryAndAdd(sk.Get(), polypool, refreshShares.h0)

	// h1 = -s*a
	bfvcontext.ContextQ().NTT(crs, polypool)
	bfvcontext.ContextQ().MulCoeffsMontgomeryAndSub(sk.Get(), polypool, refreshShares.h1)

	bfvcontext.ContextQ().InvNTT(refreshShares.h0, refreshShares.h0)
	bfvcontext.ContextQ().InvNTT(refreshShares.h1, refreshShares.h1)

	// TODO : add smudging noise
	sampler := bfvcontext.ContextQ().NewKYSampler(3.19, 19)

	// h0 = s*ct[1] + e
	sampler.Sample(polypool)
	bfvcontext.ContextQ().Add(refreshShares.h0, polypool, refreshShares.h0)

	// h1 = s*a + e'
	sampler.Sample(polypool)
	bfvcontext.ContextQ().Add(refreshShares.h1, polypool, refreshShares.h1)

	// mask = (uniform plaintext in [0, T-1]) * floor(Q/T)
	coeffs := bfvcontext.ContextT().NewUniformPoly()
	lift(coeffs, polypool, bfvcontext)

	// h0 = s*ct[1] + mask
	bfvcontext.ContextQ().Add(refreshShares.h0, polypool, refreshShares.h0)

	// h0 = -s*a - mask
	bfvcontext.ContextQ().Sub(refreshShares.h1, polypool, refreshShares.h1)

	return
}

func Refresh(ciphertext *bfv.Ciphertext, sk *ring.Poly, refreshShares []*RefreshShares, bfvcontext *bfv.BfvContext, crs *ring.Poly, encoder *bfv.BatchEncoder) {

	scaler := ring.NewSimpleScaler(bfvcontext.T(), bfvcontext.ContextQ())

	// ct[0] += sum(h0_i)
	for i := range refreshShares {
		bfvcontext.ContextQ().Add(ciphertext.Value()[0], refreshShares[i].h0, ciphertext.Value()[0])
	}

	// (round(ct[0] * T / Q) % T) * floor(Q/T)
	scaler.Scale(ciphertext.Value()[0], ciphertext.Value()[0])
	lift(ciphertext.Value()[0], ciphertext.Value()[0], bfvcontext)

	// ct[0] += sum(h1_i)
	for i := range refreshShares {
		bfvcontext.ContextQ().Add(ciphertext.Value()[0], refreshShares[i].h1, ciphertext.Value()[0])
	}

	// ct[1] = a
	ciphertext.Value()[1] = crs.CopyNew()

}

func lift(p0, p1 *ring.Poly, bfvcontext *bfv.BfvContext) {
	for j := uint64(0); j < bfvcontext.N(); j++ {
		for i := len(bfvcontext.ContextQ().Modulus) - 1; i >= 0; i-- {
			p1.Coeffs[i][j] = ring.MRed(p0.Coeffs[0][j], bfvcontext.DeltaMont()[i], bfvcontext.ContextQ().Modulus[i], bfvcontext.ContextQ().GetMredParams()[i])
		}
	}
}
