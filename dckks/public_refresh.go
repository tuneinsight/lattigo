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

	context := ckkscontext.ContextQ()

	bound := new(ring.Int)
	bound.SetBigInt(ckkscontext.BigintChain()[levelStart])
	bound.Div(bound, ring.NewUint(2*nParties))

	coeffs_bigint := make([]*ring.Int, 1<<ckkscontext.LogN())

	for i := range coeffs_bigint {
		coeffs_bigint[i] = ring.RandInt(bound)
		coeffs_bigint[i].Center(bound)
	}

	refreshShares.h0 = context.NewPolyLvl(levelStart)
	refreshShares.h1 = context.NewPoly()

	// h0 = mask (at level min)
	context.SetCoefficientsBigintLvl(levelStart, coeffs_bigint, refreshShares.h0)
	context.SetCoefficientsBigint(coeffs_bigint, refreshShares.h1)

	// h1 = mask (at level max)

	context.NTTLvl(levelStart, refreshShares.h0, refreshShares.h0)
	context.NTT(refreshShares.h1, refreshShares.h1)

	// h0 = sk*c1 + mask
	context.MulCoeffsMontgomeryAndAddLvl(levelStart, sk.Get(), c1, refreshShares.h0)

	// h1 = sk*a + mask
	context.MulCoeffsMontgomeryAndAdd(sk.Get(), crs, refreshShares.h1)

	sampler := context.NewKYSampler(3.19, 19)

	// h0 = sk*c1 + mask + e0
	context.AddLvl(levelStart, refreshShares.h0, sampler.SampleNTTNew(), refreshShares.h0)

	// h1 = sk*a + mask + e1
	context.Add(refreshShares.h1, sampler.SampleNTTNew(), refreshShares.h1)

	// h1 = -sk*c1 - mask - e0
	context.Neg(refreshShares.h1, refreshShares.h1)

	return

}

func Refresh(ciphertext *ckks.Ciphertext, refreshShares []*RefreshShares, ckkscontext *ckks.CkksContext, crs *ring.Poly) {

	// ct[0] += sum(h0_i)
	for i := range refreshShares {
		ckkscontext.ContextQ().AddLvl(ciphertext.Level(), ciphertext.Value()[0], refreshShares[i].h0, ciphertext.Value()[0])
	}

	ckkscontext.ContextQ().InvNTTLvl(ciphertext.Level(), ciphertext.Value()[0], ciphertext.Value()[0])

	coeffs_bigint := make([]*ring.Int, 1<<ckkscontext.LogN())

	ckkscontext.ContextQ().PolyToBigint(ciphertext.Value()[0], coeffs_bigint)

	QStart := ckkscontext.BigintChain()[ciphertext.Level()]
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

	ckkscontext.ContextQ().SetCoefficientsBigintLvl(ciphertext.Level(), coeffs_bigint, ciphertext.Value()[0])

	ciphertext.SetCurrentModulus(ckkscontext.BigintChain()[ciphertext.Level()])

	ckkscontext.ContextQ().NTTLvl(ciphertext.Level(), ciphertext.Value()[0], ciphertext.Value()[0])

	// ct[0] += sum(h0_i)
	for i := range refreshShares {
		ckkscontext.ContextQ().AddLvl(ciphertext.Level(), ciphertext.Value()[0], refreshShares[i].h1, ciphertext.Value()[0])
	}

	ciphertext.Value()[1] = crs.CopyNew()
}
