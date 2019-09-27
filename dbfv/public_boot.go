package dbfv

import (
	"github.com/ldsec/lattigo/bfv"
	"github.com/ldsec/lattigo/ring"
	//"fmt"
)

type BootShares struct {
	h0 *ring.Poly
	h1 *ring.Poly
}

func GenBootShares(sk *bfv.SecretKey, ciphertext *bfv.Ciphertext, bfvcontext *bfv.BfvContext, crs *ring.Poly, encoder *bfv.BatchEncoder) (bootshares *BootShares) {

	bootshares = new(BootShares)

	sampler := bfvcontext.ContextQ().NewKYSampler(3.19, 19)

	polypool := bfvcontext.ContextQ().NewPoly()

	coeffs := bfvcontext.ContextT().NewUniformPoly()
	mask := bfvcontext.NewPlaintext()
	encoder.EncodeUint(coeffs.Coeffs[0], mask)

	bootshares.h0 = bfvcontext.ContextQ().NewPoly()
	bootshares.h1 = bfvcontext.ContextQ().NewPoly()

	// h0 = s*ct[1]
	bfvcontext.ContextQ().NTT(ciphertext.Value()[1], polypool)
	bfvcontext.ContextQ().MulCoeffsMontgomeryAndAdd(sk.Get(), polypool, bootshares.h0)

	// h1 = -s*a
	bfvcontext.ContextQ().NTT(crs, polypool)
	bfvcontext.ContextQ().MulCoeffsMontgomeryAndSub(sk.Get(), polypool, bootshares.h1)

	bfvcontext.ContextQ().InvNTT(bootshares.h0, bootshares.h0)
	bfvcontext.ContextQ().InvNTT(bootshares.h1, bootshares.h1)

	// h0 = s*ct[1] - mask + e
	sampler.Sample(polypool)
	bfvcontext.ContextQ().Add(bootshares.h0, polypool, bootshares.h0)

	// h1 = s*a - mask + e'
	sampler.Sample(polypool)
	bfvcontext.ContextQ().Add(bootshares.h1, polypool, bootshares.h1)

	// h0 = s*ct[1] + mask
	bfvcontext.ContextQ().Add(bootshares.h0, mask.Value()[0], bootshares.h0)

	// h0 = -s*a - mask
	bfvcontext.ContextQ().Sub(bootshares.h1, mask.Value()[0], bootshares.h1)

	return
}

func Bootstrapp(ciphertext *bfv.Ciphertext, sk *ring.Poly, bootshares []*BootShares, bfvcontext *bfv.BfvContext, crs *ring.Poly, encoder *bfv.BatchEncoder) {

	simplescaler, _ := ring.NewSimpleScaler(bfvcontext.T(), bfvcontext.ContextQ())

	// ct[0] += sum(h0_i)
	for i := range bootshares {
		bfvcontext.ContextQ().Add(ciphertext.Value()[0], bootshares[i].h0, ciphertext.Value()[0])
	}

	// (floor(ct[0] * T / Q) % T) * Q/T
	reencode(ciphertext.Value()[0], ciphertext.Value()[0], simplescaler, bfvcontext)

	// ct[0] += sum(h1_i)
	for i := range bootshares {
		bfvcontext.ContextQ().Add(ciphertext.Value()[0], bootshares[i].h1, ciphertext.Value()[0])
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
