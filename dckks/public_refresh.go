package dckks

import (
	"github.com/ldsec/lattigo/ckks"
	"github.com/ldsec/lattigo/ring"
	"math/big"
)

// RefreshProtocol is a struct storing the parameters for the Refresh protocol.
type RefreshProtocol struct {
	dckksContext *dckksContext
	tmp          *ring.Poly
	maskBigint   []*big.Int
}

// RefreshShareDecrypt is a struct storing the masked decryption share.
type RefreshShareDecrypt *ring.Poly

// RefreshShareRecrypt is a struct storing the masked recryption share.
type RefreshShareRecrypt *ring.Poly

// NewRefreshProtocol creates a new instance of the Refresh protocol.
func NewRefreshProtocol(params *ckks.Parameters) (refreshProtocol *RefreshProtocol) {

	if !params.IsValid() {
		panic("cannot NewRefreshProtocol : params not valid (check if they where generated properly)")
	}

	refreshProtocol = new(RefreshProtocol)
	dckksContext := newDckksContext(params)
	refreshProtocol.dckksContext = dckksContext
	refreshProtocol.tmp = dckksContext.contextQ.NewPoly()
	refreshProtocol.maskBigint = make([]*big.Int, dckksContext.n)
	return
}

// AllocateShares allocates the shares of the Refresh protocol.
func (refreshProtocol *RefreshProtocol) AllocateShares(levelStart uint64) (RefreshShareDecrypt, RefreshShareRecrypt) {
	return refreshProtocol.dckksContext.contextQ.NewPolyLvl(levelStart), refreshProtocol.dckksContext.contextQ.NewPoly()
}

// GenShares generates the decryption and recryption shares of the Refresh protocol.
func (refreshProtocol *RefreshProtocol) GenShares(sk *ring.Poly, levelStart, nParties uint64, ciphertext *ckks.Ciphertext, crs *ring.Poly, shareDecrypt RefreshShareDecrypt, shareRecrypt RefreshShareRecrypt) {

	context := refreshProtocol.dckksContext.contextQ
	sampler := context.NewSampler(3.19, 19)

	bound := ring.NewUint(context.Modulus[0])
	for i := uint64(1); i < levelStart+1; i++ {
		bound.Mul(bound, ring.NewUint(context.Modulus[i]))
	}

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
	refreshProtocol.dckksContext.contextQ.AddLvl(uint64(len(share1.Coeffs)-1), share1, share2, shareOut)
}

// Decrypt operates a masked decryption on the ciphertext with the given decryption share.
func (refreshProtocol *RefreshProtocol) Decrypt(ciphertext *ckks.Ciphertext, shareDecrypt RefreshShareDecrypt) {
	refreshProtocol.dckksContext.contextQ.AddLvl(ciphertext.Level(), ciphertext.Value()[0], shareDecrypt, ciphertext.Value()[0])
}

// Recode takes a masked decrypted ciphertext at modulus Q_0 and returns the same masked decrypted ciphertext at modulus Q_L, with Q_0 << Q_L.
func (refreshProtocol *RefreshProtocol) Recode(ciphertext *ckks.Ciphertext) {
	dckksContext := refreshProtocol.dckksContext
	context := refreshProtocol.dckksContext.contextQ

	context.InvNTTLvl(ciphertext.Level(), ciphertext.Value()[0], ciphertext.Value()[0])

	context.PolyToBigint(ciphertext.Value()[0], refreshProtocol.maskBigint)

	QStart := ring.NewUint(context.Modulus[0])
	for i := uint64(1); i < ciphertext.Level()+1; i++ {
		QStart.Mul(QStart, ring.NewUint(context.Modulus[i]))
	}

	QHalf := new(big.Int).Rsh(QStart, 1)

	for ciphertext.Level() != uint64(len(dckksContext.params.Qi)-1) {
		ciphertext.Value()[0].Coeffs = append(ciphertext.Value()[0].Coeffs, make([][]uint64, 1)...)
		ciphertext.Value()[0].Coeffs[ciphertext.Level()] = make([]uint64, dckksContext.n)
	}

	var sign int
	for i := uint64(0); i < dckksContext.n; i++ {
		sign = refreshProtocol.maskBigint[i].Cmp(QHalf)
		if sign == 1 || sign == 0 {
			refreshProtocol.maskBigint[i].Sub(refreshProtocol.maskBigint[i], QStart)
		}
	}

	context.SetCoefficientsBigintLvl(ciphertext.Level(), refreshProtocol.maskBigint, ciphertext.Value()[0])

	context.NTTLvl(ciphertext.Level(), ciphertext.Value()[0], ciphertext.Value()[0])
}

// Recrypt operates a masked recryption on the masked decrypted ciphertext.
func (refreshProtocol *RefreshProtocol) Recrypt(ciphertext *ckks.Ciphertext, crs *ring.Poly, shareRecrypt RefreshShareRecrypt) {

	refreshProtocol.dckksContext.contextQ.Add(ciphertext.Value()[0], shareRecrypt, ciphertext.Value()[0])

	ciphertext.Value()[1] = crs.CopyNew()
}
