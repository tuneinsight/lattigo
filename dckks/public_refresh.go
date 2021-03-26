package dckks

import (
	"math/big"

	"github.com/ldsec/lattigo/v2/ckks"
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/utils"
)

// RefreshProtocol is a struct storing the parameters for the Refresh protocol.
type RefreshProtocol struct {
	dckksContext    *dckksContext
	tmp             *ring.Poly
	maskBigint      []*big.Int
	gaussianSampler *ring.GaussianSampler
}

// RefreshShareDecrypt is a struct storing the masked decryption share.
type RefreshShareDecrypt *ring.Poly

// RefreshShareRecrypt is a struct storing the masked recryption share.
type RefreshShareRecrypt *ring.Poly

// NewRefreshProtocol creates a new instance of the Refresh protocol.
func NewRefreshProtocol(params *ckks.Parameters) (refreshProtocol *RefreshProtocol) {

	refreshProtocol = new(RefreshProtocol)
	dckksContext := newDckksContext(params)
	refreshProtocol.dckksContext = dckksContext
	refreshProtocol.tmp = dckksContext.ringQ.NewPoly()
	refreshProtocol.maskBigint = make([]*big.Int, dckksContext.n)
	prng, err := utils.NewPRNG()
	if err != nil {
		panic(err)
	}
	refreshProtocol.gaussianSampler = ring.NewGaussianSampler(prng)

	return
}

// AllocateShares allocates the shares of the Refresh protocol.
func (refreshProtocol *RefreshProtocol) AllocateShares(levelStart int) (RefreshShareDecrypt, RefreshShareRecrypt) {
	return refreshProtocol.dckksContext.ringQ.NewPolyLvl(levelStart), refreshProtocol.dckksContext.ringQ.NewPoly()
}

// GenShares generates the decryption and recryption shares of the Refresh protocol.
func (refreshProtocol *RefreshProtocol) GenShares(sk *ring.Poly, levelStart, nParties int, ciphertext *ckks.Ciphertext, targetScale float64, crs *ring.Poly, shareDecrypt RefreshShareDecrypt, shareRecrypt RefreshShareRecrypt) {

	ringQ := refreshProtocol.dckksContext.ringQ
	sigma := refreshProtocol.dckksContext.params.Sigma()

	bound := ring.NewUint(ringQ.Modulus[0])
	for i := 1; i < levelStart+1; i++ {
		bound.Mul(bound, ring.NewUint(ringQ.Modulus[i]))
	}

	bound.Quo(bound, ring.NewUint(uint64(2*nParties)))
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
	ringQ.SetCoefficientsBigintLvl(levelStart, refreshProtocol.maskBigint, shareDecrypt)

	inputScaleFlo := ring.NewFloat(ciphertext.Scale(), 256)
	outputScaleFlo := ring.NewFloat(targetScale, 256)

	inputScaleInt := new(big.Int)
	outputScaleInt := new(big.Int)

	inputScaleFlo.Int(inputScaleInt)
	outputScaleFlo.Int(outputScaleInt)

	// Scales the mask by the ratio between the two scales
	for i := range refreshProtocol.maskBigint {
		refreshProtocol.maskBigint[i].Mul(refreshProtocol.maskBigint[i], outputScaleInt)
		refreshProtocol.maskBigint[i].Quo(refreshProtocol.maskBigint[i], inputScaleInt)
	}

	// h1 = mask (at level max)
	ringQ.SetCoefficientsBigint(refreshProtocol.maskBigint, shareRecrypt)

	for i := range refreshProtocol.maskBigint {
		refreshProtocol.maskBigint[i].SetUint64(0)
	}

	ringQ.NTTLvl(levelStart, shareDecrypt, shareDecrypt)
	ringQ.NTT(shareRecrypt, shareRecrypt)

	// h0 = sk*c1 + mask
	ringQ.MulCoeffsMontgomeryAndAddLvl(levelStart, sk, ciphertext.Value()[1], shareDecrypt)

	// h1 = sk*a + mask
	ringQ.MulCoeffsMontgomeryAndAdd(sk, crs, shareRecrypt)

	// h0 = sk*c1 + mask + e0
	refreshProtocol.gaussianSampler.ReadLvl(levelStart, refreshProtocol.tmp, ringQ, sigma, int(6*sigma))
	ringQ.NTTLvl(levelStart, refreshProtocol.tmp, refreshProtocol.tmp)
	ringQ.AddLvl(levelStart, shareDecrypt, refreshProtocol.tmp, shareDecrypt)

	// h1 = sk*a + mask + e1
	refreshProtocol.gaussianSampler.Read(refreshProtocol.tmp, ringQ, sigma, int(6*sigma))
	ringQ.NTT(refreshProtocol.tmp, refreshProtocol.tmp)
	ringQ.Add(shareRecrypt, refreshProtocol.tmp, shareRecrypt)

	// h1 = -sk*c1 - mask - e0
	ringQ.Neg(shareRecrypt, shareRecrypt)

	refreshProtocol.tmp.Zero()
}

// Aggregate adds share1 with share2 on shareOut.
func (refreshProtocol *RefreshProtocol) Aggregate(share1, share2, shareOut *ring.Poly) {
	refreshProtocol.dckksContext.ringQ.AddLvl(len(share1.Coeffs)-1, share1, share2, shareOut)
}

// Decrypt operates a masked decryption on the ciphertext with the given decryption share.
func (refreshProtocol *RefreshProtocol) Decrypt(ciphertext *ckks.Ciphertext, shareDecrypt RefreshShareDecrypt) {
	refreshProtocol.dckksContext.ringQ.AddLvl(ciphertext.Level(), ciphertext.Value()[0], shareDecrypt, ciphertext.Value()[0])
}

// Recode takes a masked decrypted ciphertext at modulus Q_0 and returns the same masked decrypted ciphertext at modulus Q_L, with Q_0 << Q_L.
func (refreshProtocol *RefreshProtocol) Recode(ciphertext *ckks.Ciphertext, targetScale float64) {
	dckksContext := refreshProtocol.dckksContext
	ringQ := refreshProtocol.dckksContext.ringQ

	inputScaleFlo := ring.NewFloat(ciphertext.Scale(), 256)
	outputScaleFlo := ring.NewFloat(targetScale, 256)

	inputScaleInt := new(big.Int)
	outputScaleInt := new(big.Int)

	inputScaleFlo.Int(inputScaleInt)
	outputScaleFlo.Int(outputScaleInt)

	ringQ.InvNTTLvl(ciphertext.Level(), ciphertext.Value()[0], ciphertext.Value()[0])

	ringQ.PolyToBigint(ciphertext.Value()[0], refreshProtocol.maskBigint)

	QStart := ring.NewUint(ringQ.Modulus[0])
	for i := 1; i < ciphertext.Level()+1; i++ {
		QStart.Mul(QStart, ring.NewUint(ringQ.Modulus[i]))
	}

	QHalf := new(big.Int).Rsh(QStart, 1)

	for ciphertext.Level() != dckksContext.params.MaxLevel() {
		ciphertext.Value()[0].Coeffs = append(ciphertext.Value()[0].Coeffs, make([][]uint64, 1)...)
		ciphertext.Value()[0].Coeffs[ciphertext.Level()] = make([]uint64, dckksContext.n)
	}

	var sign int
	for i := 0; i < dckksContext.n; i++ {
		sign = refreshProtocol.maskBigint[i].Cmp(QHalf)
		if sign == 1 || sign == 0 {
			refreshProtocol.maskBigint[i].Sub(refreshProtocol.maskBigint[i], QStart)
		}

		refreshProtocol.maskBigint[i].Mul(refreshProtocol.maskBigint[i], outputScaleInt)
		refreshProtocol.maskBigint[i].Quo(refreshProtocol.maskBigint[i], inputScaleInt)
	}

	ringQ.SetCoefficientsBigintLvl(ciphertext.Level(), refreshProtocol.maskBigint, ciphertext.Value()[0])

	ringQ.NTTLvl(ciphertext.Level(), ciphertext.Value()[0], ciphertext.Value()[0])

	ciphertext.SetScale(targetScale)
}

// Recrypt operates a masked recryption on the masked decrypted ciphertext.
func (refreshProtocol *RefreshProtocol) Recrypt(ciphertext *ckks.Ciphertext, crs *ring.Poly, shareRecrypt RefreshShareRecrypt) {

	refreshProtocol.dckksContext.ringQ.Add(ciphertext.Value()[0], shareRecrypt, ciphertext.Value()[0])
	crs.Coeffs = crs.Coeffs[:ciphertext.Level()+1]
	ciphertext.Value()[1] = crs.CopyNew()
}
