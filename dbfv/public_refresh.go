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

	sampler := bfvcontext.ContextQ().NewKYSampler(3.19, 19)

	polypool := bfvcontext.ContextQ().NewPoly()

	coeffs := bfvcontext.ContextT().NewUniformPoly()
	mask := bfvcontext.NewPlaintext()
	encoder.EncodeUint(coeffs.Coeffs[0], mask)

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

	// h0 = s*ct[1] - mask + e
	sampler.Sample(polypool)
	bfvcontext.ContextQ().Add(refreshShares.h0, polypool, refreshShares.h0)

	// h1 = s*a - mask + e'
	sampler.Sample(polypool)
	bfvcontext.ContextQ().Add(refreshShares.h1, polypool, refreshShares.h1)

	// h0 = s*ct[1] + mask
	bfvcontext.ContextQ().Add(refreshShares.h0, mask.Value()[0], refreshShares.h0)

	// h0 = -s*a - mask
	bfvcontext.ContextQ().Sub(refreshShares.h1, mask.Value()[0], refreshShares.h1)

	return
}

func Refresh(ciphertext *bfv.Ciphertext, sk *ring.Poly, refreshShares []*RefreshShares, bfvcontext *bfv.BfvContext, crs *ring.Poly, encoder *bfv.BatchEncoder) {

	simplescaler, _ := ring.NewSimpleScaler(bfvcontext.T(), bfvcontext.ContextQ())

	// ct[0] += sum(h0_i)
	for i := range refreshShares {
		bfvcontext.ContextQ().Add(ciphertext.Value()[0], refreshShares[i].h0, ciphertext.Value()[0])
	}

	// (floor(ct[0] * T / Q) % T) * Q/T
	reencode(ciphertext.Value()[0], ciphertext.Value()[0], simplescaler, bfvcontext)

	// ct[0] += sum(h1_i)
	for i := range refreshShares {
		bfvcontext.ContextQ().Add(ciphertext.Value()[0], refreshShares[i].h1, ciphertext.Value()[0])
	}

	// ct[1] = a
	ciphertext.Value()[1] = crs.CopyNew()

}

func reencode(p0, p1 *ring.Poly, scaler *ring.SimpleScaler, bfvcontext *bfv.BfvContext) {

	scaler.Scale(p0, p1)

	for j := uint64(0); j < bfvcontext.N(); j++ {
		for i := len(bfvcontext.ContextQ().Modulus) - 1; i >= 0; i-- {
			p1.Coeffs[i][j] = ring.MRed(p1.Coeffs[0][j], bfvcontext.DeltaMont()[i], bfvcontext.ContextQ().Modulus[i], bfvcontext.ContextQ().GetMredParams()[i])
		}
	}

}
