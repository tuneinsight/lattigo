package dckks

import (
	"github.com/ldsec/lattigo/ckks"
	"github.com/ldsec/lattigo/ring"
)

type RefreshShares struct {
	h0 *ring.Poly
	h1 *ring.Poly
}

func GenRefreshShares(sk *ckks.SecretKey, levelStart, nParties uint64, ckkscontext *ckks.CkksContext, c1, crs *ring.Poly) (refreshShares *RefreshShares) {

	refreshShares = new(RefreshShares)

	contextStart := ckkscontext.Context(levelStart)
	contextEnd := ckkscontext.Context(ckkscontext.Levels() - 1)

	bound := new(ring.Int)
	bound.SetBigInt(ckkscontext.Context(levelStart).ModulusBigint)
	bound.Div(bound, ring.NewUint(2*nParties))

	coeffs_bigint := make([]*ring.Int, 1<<ckkscontext.LogN())

	for i := range coeffs_bigint {
		coeffs_bigint[i] = ring.RandInt(bound)
		coeffs_bigint[i].Center(bound)
	}

	refreshShares.h0 = contextStart.NewPoly()
	refreshShares.h1 = contextEnd.NewPoly()

	// h0 = mask (at level min)
	contextStart.SetCoefficientsBigint(coeffs_bigint, refreshShares.h0)
	contextEnd.SetCoefficientsBigint(coeffs_bigint, refreshShares.h1)

	// h1 = mask (at level max)

	contextStart.NTT(refreshShares.h0, refreshShares.h0)
	contextEnd.NTT(refreshShares.h1, refreshShares.h1)

	// h0 = sk*c1 + mask
	contextStart.MulCoeffsMontgomeryAndAdd(sk.Get(), c1, refreshShares.h0)

	// h1 = sk*a + mask
	contextEnd.MulCoeffsMontgomeryAndAdd(sk.Get(), crs, refreshShares.h1)

	sampler := contextEnd.NewKYSampler(3.19, 19)

	// h0 = sk*c1 + mask + e0
	contextStart.Add(refreshShares.h0, sampler.SampleNTTNew(), refreshShares.h0)

	// h1 = sk*a + mask + e1
	contextEnd.Add(refreshShares.h1, sampler.SampleNTTNew(), refreshShares.h1)

	// h1 = -sk*c1 - mask - e0
	contextEnd.Neg(refreshShares.h1, refreshShares.h1)

	return

}

func Refresh(ciphertext *ckks.Ciphertext, refreshShares []*RefreshShares, ckkscontext *ckks.CkksContext, crs *ring.Poly, encoder *ckks.Encoder, decryptor *ckks.Decryptor) {

	// ct[0] += sum(h0_i)
	for i := range refreshShares {
		ckkscontext.Context(ciphertext.Level()).Add(ciphertext.Value()[0], refreshShares[i].h0, ciphertext.Value()[0])
	}

	ckkscontext.Context(ciphertext.Level()).InvNTT(ciphertext.Value()[0], ciphertext.Value()[0])

	coeffs_bigint := make([]*ring.Int, 1<<ckkscontext.LogN())

	ckkscontext.Context(ciphertext.Level()).PolyToBigint(ciphertext.Value()[0], coeffs_bigint)

	QStart := ckkscontext.Context(ciphertext.Level()).ModulusBigint
	QHalf := QStart.Copy()
	QHalf.Rsh(QHalf, 1)

	for ciphertext.Level() != ckkscontext.Levels()-1 {
		ciphertext.Value()[0].Coeffs = append(ciphertext.Value()[0].Coeffs, make([][]uint64, 1)...)
		ciphertext.Value()[0].Coeffs[ciphertext.Level()] = make([]uint64, 1<<ckkscontext.LogN())
	}

	var sign int
	for i := uint64(0); i < 1<<ckkscontext.LogN(); i++ {
		sign = coeffs_bigint[i].Compare(QHalf)
		if sign == 1 || sign == 0 {
			coeffs_bigint[i].Sub(coeffs_bigint[i], QStart)
		}
	}

	ckkscontext.Context(ciphertext.Level()).SetCoefficientsBigint(coeffs_bigint, ciphertext.Value()[0])

	ciphertext.SetCurrentModulus(ckkscontext.Context(ciphertext.Level()).ModulusBigint)

	ckkscontext.Context(ciphertext.Level()).NTT(ciphertext.Value()[0], ciphertext.Value()[0])

	// ct[0] += sum(h0_i)
	for i := range refreshShares {
		ckkscontext.Context(ckkscontext.Levels()-1).Add(ciphertext.Value()[0], refreshShares[i].h1, ciphertext.Value()[0])
	}

	ciphertext.Value()[1] = crs.CopyNew()
}
