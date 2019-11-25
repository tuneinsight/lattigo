package dckks

import (
	"github.com/ldsec/lattigo/ckks"
	"github.com/ldsec/lattigo/ring"
	"math/big"
)

type RefreshProtocol struct {
	ckksContext *ckks.Context
	tmp         *ring.Poly
	maskBigint  []*big.Int
}

type RefreshShareDecrypt *ring.Poly
type RefreshShareRecrypt *ring.Poly

func NewRefreshProtocol(ckksContext *ckks.Context) (refreshProtocol *RefreshProtocol) {
	refreshProtocol = new(RefreshProtocol)
	refreshProtocol.ckksContext = ckksContext
	refreshProtocol.tmp = ckksContext.ContextQ().NewPoly()
	refreshProtocol.maskBigint = make([]*big.Int, 1<<ckksContext.LogN())
	return
}

func (refreshProtocol *RefreshProtocol) AllocateShares(levelStart uint64) (RefreshShareDecrypt, RefreshShareRecrypt) {
	return refreshProtocol.ckksContext.ContextQ().NewPolyLvl(levelStart), refreshProtocol.ckksContext.ContextQ().NewPoly()
}

func (refreshProtocol *RefreshProtocol) GenShares(sk *ring.Poly, levelStart, nParties uint64, ciphertext *ckks.Ciphertext, crs *ring.Poly, shareDecrypt RefreshShareDecrypt, shareRecrypt RefreshShareRecrypt) {

	context := refreshProtocol.ckksContext.ContextQ()
	sampler := context.NewKYSampler(3.19, 19)

	bound := new(big.Int).Set(refreshProtocol.ckksContext.BigintChain()[levelStart])
	bound.Quo(bound, ring.NewUint(2*nParties))
	boundHalf := new(big.Int).Rsh(bound, 1)

	var sign int
	for i := range refreshProtocol.maskBigint {
		refreshProtocol.maskBigint[i] = ring.RandInt(bound)
		sign = refreshProtocol.maskBigint[i].Cmp(boundHalf)
		if sign == 1 || sign == 0 {
			refreshProtocol.maskBigint[i].Sub(refreshProtocol.maskBigint[i], bound)
		}
	}

	// h0 = mask (at level min)
	context.SetCoefficientsBigintLvl(levelStart, refreshProtocol.maskBigint, shareDecrypt)
	// h1 = mask (at level max)
	context.SetCoefficientsBigint(refreshProtocol.maskBigint, shareRecrypt)

	for i := range refreshProtocol.maskBigint {
		refreshProtocol.maskBigint[i] = new(big.Int)
	}

	context.NTTLvl(levelStart, shareDecrypt, shareDecrypt)
	context.NTT(shareRecrypt, shareRecrypt)

	// h0 = sk*c1 + mask
	context.MulCoeffsMontgomeryAndAddLvl(levelStart, sk, ciphertext.Value()[1], shareDecrypt)

	// h1 = sk*a + mask
	context.MulCoeffsMontgomeryAndAdd(sk, crs, shareRecrypt)

	// h0 = sk*c1 + mask + e0
	sampler.SampleNTT(refreshProtocol.tmp)
	context.AddLvl(levelStart, shareDecrypt, refreshProtocol.tmp, shareDecrypt)

	// h1 = sk*a + mask + e1
	sampler.SampleNTT(refreshProtocol.tmp)
	context.Add(shareRecrypt, refreshProtocol.tmp, shareRecrypt)

	// h1 = -sk*c1 - mask - e0
	context.Neg(shareRecrypt, shareRecrypt)

	refreshProtocol.tmp.Zero()
}

func (refreshProtocol *RefreshProtocol) Aggregate(share1, share2, shareOut *ring.Poly) {
	refreshProtocol.ckksContext.ContextQ().AddLvl(uint64(len(share1.Coeffs)-1), share1, share2, shareOut)
}

func (refreshProtocol *RefreshProtocol) Decrypt(ciphertext *ckks.Ciphertext, shareDecrypt RefreshShareDecrypt) {
	refreshProtocol.ckksContext.ContextQ().AddLvl(ciphertext.Level(), ciphertext.Value()[0], shareDecrypt, ciphertext.Value()[0])
}

func (refreshProtocol *RefreshProtocol) Recode(ciphertext *ckks.Ciphertext) {
	ckksContext := refreshProtocol.ckksContext
	contextQ := refreshProtocol.ckksContext.ContextQ()

	contextQ.InvNTTLvl(ciphertext.Level(), ciphertext.Value()[0], ciphertext.Value()[0])

	contextQ.PolyToBigint(ciphertext.Value()[0], refreshProtocol.maskBigint)

	QStart := ckksContext.BigintChain()[ciphertext.Level()]
	QHalf := new(big.Int).Rsh(QStart, 1)

	for ciphertext.Level() != ckksContext.Levels()-1 {
		ciphertext.Value()[0].Coeffs = append(ciphertext.Value()[0].Coeffs, make([][]uint64, 1)...)
		ciphertext.Value()[0].Coeffs[ciphertext.Level()] = make([]uint64, 1<<ckksContext.LogN())
	}

	var sign int
	for i := uint64(0); i < 1<<ckksContext.LogN(); i++ {
		sign = refreshProtocol.maskBigint[i].Cmp(QHalf)
		if sign == 1 || sign == 0 {
			refreshProtocol.maskBigint[i].Sub(refreshProtocol.maskBigint[i], QStart)
		}
	}

	contextQ.SetCoefficientsBigintLvl(ciphertext.Level(), refreshProtocol.maskBigint, ciphertext.Value()[0])

	ciphertext.SetCurrentModulus(ckksContext.BigintChain()[ciphertext.Level()])

	contextQ.NTTLvl(ciphertext.Level(), ciphertext.Value()[0], ciphertext.Value()[0])
}

func (refreshProtocol *RefreshProtocol) Recrypt(ciphertext *ckks.Ciphertext, crs *ring.Poly, shareRecrypt RefreshShareRecrypt) {

	refreshProtocol.ckksContext.ContextQ().Add(ciphertext.Value()[0], shareRecrypt, ciphertext.Value()[0])

	ciphertext.Value()[1] = crs.CopyNew()
}
