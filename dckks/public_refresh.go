package dckks

import (
	"github.com/ldsec/lattigo/ckks"
	"github.com/ldsec/lattigo/ring"
	"math/big"
)

// RefreshProtocol is a struct storing the parameters for the Refresh protocol.
type RefreshProtocol struct {
	ckksContext *ckks.Context
	tmp         *ring.Poly
	maskBigint  []*big.Int
}

// RefreshShareDecrypt is a struct storing the masked decryption share.
type RefreshShareDecrypt *ring.Poly

// RefreshShareRecrypt is a struct storing the masked recryption share.
type RefreshShareRecrypt *ring.Poly

// NewRefreshProtocol creates a new instance of the Refresh protocol.
func NewRefreshProtocol(ckksContext *ckks.Context) (refreshProtocol *RefreshProtocol) {
	refreshProtocol = new(RefreshProtocol)
	refreshProtocol.ckksContext = ckksContext
	refreshProtocol.tmp = ckksContext.ContextQ().NewPoly()
	refreshProtocol.maskBigint = make([]*big.Int, 1<<ckksContext.LogN())
	return
}

// AllocateShares allocates the shares of the Refresh protocol.
func (refreshProtocol *RefreshProtocol) AllocateShares(levelStart uint64) (RefreshShareDecrypt, RefreshShareRecrypt) {
	return refreshProtocol.ckksContext.ContextQ().NewPolyLvl(levelStart), refreshProtocol.ckksContext.ContextQ().NewPoly()
}

// GenShares generates the decryption and recryption shares of the Refresh protocol.
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

// Aggregate adds share1 with share2 on shareOut.
func (refreshProtocol *RefreshProtocol) Aggregate(share1, share2, shareOut *ring.Poly) {
	refreshProtocol.ckksContext.ContextQ().AddLvl(uint64(len(share1.Coeffs)-1), share1, share2, shareOut)
}

// Decrypt operates a masked decryption on the ciphertext with the given decryption share.
func (refreshProtocol *RefreshProtocol) Decrypt(ciphertext *ckks.Ciphertext, shareDecrypt RefreshShareDecrypt) {
	refreshProtocol.ckksContext.ContextQ().AddLvl(ciphertext.Level(), ciphertext.Value()[0], shareDecrypt, ciphertext.Value()[0])
}

// Recode takes a masked decrypted ciphertext at modulus Q_0 and returns the same masked decrypted ciphertext at modulus Q_L, with Q_0 << Q_L.
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

	contextQ.NTTLvl(ciphertext.Level(), ciphertext.Value()[0], ciphertext.Value()[0])
}

// Recrypt operates a masked recryption on the masked decrypted ciphertext.
func (refreshProtocol *RefreshProtocol) Recrypt(ciphertext *ckks.Ciphertext, crs *ring.Poly, shareRecrypt RefreshShareRecrypt) {

	refreshProtocol.ckksContext.ContextQ().Add(ciphertext.Value()[0], shareRecrypt, ciphertext.Value()[0])

	ciphertext.Value()[1] = crs.CopyNew()
}
