package dbfv

import (
	"github.com/ldsec/lattigo/bfv"
	"github.com/ldsec/lattigo/ring"
)

type BootShares struct {
	h0 *ring.Poly
	h1 *ring.Poly
}

func GenBootShares(sk *bfv.SecretKey, ciphertext *bfv.Ciphertext, bfvcontext *bfv.BfvContext, crs *ring.Poly, encoder *bfv.BatchEncoder) (bootshares *BootShares) {

	bootshares = new(BootShares)

	sampler := bfvcontext.ContextQ().NewKYSampler(3.19, 19)

	coeffs := bfvcontext.ContextT().NewUniformPoly()
	mask := bfvcontext.NewPlaintext()
	encoder.EncodeUint(coeffs.Coeffs[0], mask)

	// h0, h1 = e, e'
	bootshares.h0 = sampler.SampleNTTNew()
	bootshares.h1 = sampler.SampleNTTNew()

	// h0 = e + mask
	// h1 = e' - mask
	bfvcontext.ContextQ().Add(bootshares.h0, mask.Value()[0], bootshares.h0)
	bfvcontext.ContextQ().Sub(bootshares.h1, mask.Value()[0], bootshares.h1)

	// h0 = s*ct[1] - mask + e
	bfvcontext.ContextQ().NTT(ciphertext.Value()[1], ciphertext.Value()[1])
	bfvcontext.ContextQ().MulCoeffsMontgomeryAndAdd(sk.Get(), ciphertext.Value()[1], bootshares.h0)
	bfvcontext.ContextQ().InvNTT(ciphertext.Value()[1], ciphertext.Value()[1])

	// h1 = s*a + mask + e'
	bfvcontext.ContextQ().NTT(crs, crs)
	bfvcontext.ContextQ().MulCoeffsMontgomeryAndSub(sk.Get(), crs, bootshares.h1)
	bfvcontext.ContextQ().InvNTT(crs, crs)

	bfvcontext.ContextQ().InvNTT(bootshares.h0, bootshares.h0)
	bfvcontext.ContextQ().InvNTT(bootshares.h1, bootshares.h1)

	return
}

func Bootstrapp(ciphertext *bfv.Ciphertext, sk *ring.Poly, bootshares []*BootShares, bfvcontext *bfv.BfvContext, crs *ring.Poly, encoder *bfv.BatchEncoder) {

	// ct[0] += sum(h0_i)
	for i := range bootshares {
		bfvcontext.ContextQ().Add(ciphertext.Value()[0], bootshares[i].h0, ciphertext.Value()[0])
	}

	simplescaler, _ := ring.NewSimpleScaler(bfvcontext.T(), bfvcontext.ContextQ())

	// floor(ct[0] * T / Q) % T
	simplescaler.Scale(ciphertext.Value()[0], ciphertext.Value()[0])

	// ct[0] *= Q/T
	for j := uint64(0); j < bfvcontext.N(); j++ {
		for i := len(bfvcontext.ContextQ().Modulus) - 1; i >= 0; i-- {
			ciphertext.Value()[0].Coeffs[i][j] = ring.MRed(ciphertext.Value()[0].Coeffs[0][j], bfvcontext.DeltaMont()[i], bfvcontext.ContextQ().Modulus[i], bfvcontext.ContextQ().GetMredParams()[i])
		}
	}

	// ct[0] += sum(h1_i)
	for i := range bootshares {
		bfvcontext.ContextQ().Add(ciphertext.Value()[0], bootshares[i].h1, ciphertext.Value()[0])
	}

	ciphertext.Value()[1] = crs.CopyNew()

}
