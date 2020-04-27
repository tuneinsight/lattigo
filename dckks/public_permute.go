package dckks

import (
	"github.com/ldsec/lattigo/ckks"
	"github.com/ldsec/lattigo/ring"
	"math/big"
)

// RefreshProtocol is a struct storing the parameters for the Refresh protocol.
type PermuteProtocol struct {
	dckksContext *dckksContext
	encoder      ckks.EncoderBigComplex
	tmp          *ring.Poly
	maskBigint   []*big.Int
	maskFloat    []*big.Float
	maskComplex  []*ring.Complex
}

// NewRefreshProtocol creates a new instance of the Refresh protocol.
func NewPermuteProtocol(params *ckks.Parameters) (pp *PermuteProtocol) {

	if !params.IsValid() {
		panic("cannot NewRefreshProtocol : params not valid (check if they where generated properly)")
	}

	prec := uint64(256)

	pp = new(PermuteProtocol)
	pp.encoder = ckks.NewEncoderBigComplex(params, prec)
	dckksContext := newDckksContext(params)
	pp.dckksContext = dckksContext
	pp.tmp = dckksContext.contextQ.NewPoly()
	pp.maskBigint = make([]*big.Int, dckksContext.n)
	pp.maskFloat = make([]*big.Float, dckksContext.n)
	pp.maskComplex = make([]*ring.Complex, dckksContext.n>>1)

	for i := uint64(0); i < dckksContext.n>>1; i++ {
		pp.maskFloat[i] = new(big.Float)
		pp.maskFloat[i].SetPrec(uint(prec))

		pp.maskFloat[i+(dckksContext.n>>1)] = new(big.Float)
		pp.maskFloat[i+(dckksContext.n>>1)].SetPrec(uint(prec))

		pp.maskComplex[i] = new(ring.Complex)
	}

	return
}

// AllocateShares allocates the shares of the Refresh protocol.
func (pp *PermuteProtocol) AllocateShares(levelStart uint64) (RefreshShareDecrypt, RefreshShareRecrypt) {
	return pp.dckksContext.contextQ.NewPolyLvl(levelStart), pp.dckksContext.contextQ.NewPoly()
}

func (pp *PermuteProtocol) permuteWithIndex(permutation []uint64, values []*ring.Complex) {
	tmp := make([]*ring.Complex, len(values))

	for i := range values {
		tmp[i] = values[permutation[i]].Copy()
	}

	for i := range values {
		values[i] = tmp[i]
	}
}

// GenShares generates the decryption and recryption shares of the Refresh protocol.
func (pp *PermuteProtocol) GenShares(sk *ring.Poly, levelStart, nParties uint64, ciphertext *ckks.Ciphertext, crs *ring.Poly, slots uint64, permutation []uint64, shareDecrypt RefreshShareDecrypt, shareRecrypt RefreshShareRecrypt) {

	context := pp.dckksContext.contextQ

	bound := ring.NewUint(context.Modulus[0])
	for i := uint64(1); i < levelStart+1; i++ {
		bound.Mul(bound, ring.NewUint(context.Modulus[i]))
	}

	bound.Quo(bound, ring.NewUint(2*nParties))
	boundHalf := new(big.Int).Rsh(bound, 1)

	maxSlots := pp.dckksContext.n >> 1
	gap := maxSlots / slots

	// Samples the whole N coefficients for h0
	var sign int
	for i := uint64(0); i < 2*maxSlots; i++ {

		pp.maskBigint[i] = ring.RandInt(bound)

		sign = pp.maskBigint[i].Cmp(boundHalf)
		if sign == 1 || sign == 0 {
			pp.maskBigint[i].Sub(pp.maskBigint[i], bound)
		}
	}

	// h0 = mask (at level min)
	context.SetCoefficientsBigintLvl(levelStart, pp.maskBigint, shareDecrypt)
	context.NTTLvl(levelStart, shareDecrypt, shareDecrypt)
	// h0 = sk*c1 + mask
	context.MulCoeffsMontgomeryAndAddLvl(levelStart, sk, ciphertext.Value()[1], shareDecrypt)
	// h0 = sk*c1 + mask + e0
	context.SampleGaussianNTTLvl(uint64(len(context.Modulus)-1), pp.tmp, 3.19, 19)
	context.AddLvl(levelStart, shareDecrypt, pp.tmp, shareDecrypt)

	// Permutes only the (sparse) plaintext coefficients of h1
	for i, jdx, idx := uint64(0), maxSlots, uint64(0); i < slots; i, jdx, idx = i+1, jdx+gap, idx+gap {
		pp.maskFloat[idx].SetInt(pp.maskBigint[idx])
		pp.maskFloat[jdx].SetInt(pp.maskBigint[jdx])
		pp.maskComplex[idx][0] = pp.maskFloat[idx]
		pp.maskComplex[idx][1] = pp.maskFloat[jdx]
	}

	// h1 = pi(mask) (at level max)
	pp.encoder.FFT(pp.maskComplex, slots)
	pp.permuteWithIndex(permutation, pp.maskComplex)
	pp.encoder.InvFFT(pp.maskComplex, slots)

	for i, jdx, idx := uint64(0), maxSlots, uint64(0); i < slots; i, jdx, idx = i+1, jdx+gap, idx+gap {
		pp.maskComplex[i].Real().Int(pp.maskBigint[idx])
		pp.maskComplex[i].Imag().Int(pp.maskBigint[jdx])
	}

	context.SetCoefficientsBigint(pp.maskBigint, shareRecrypt)

	context.NTT(shareRecrypt, shareRecrypt)

	// h1 = sk*a + mask
	context.MulCoeffsMontgomeryAndAdd(sk, crs, shareRecrypt)

	// h1 = sk*a + mask + e1
	context.SampleGaussianNTTLvl(uint64(len(context.Modulus)-1), pp.tmp, 3.19, 19)
	context.Add(shareRecrypt, pp.tmp, shareRecrypt)

	// h1 = -sk*c1 - mask - e0
	context.Neg(shareRecrypt, shareRecrypt)

	pp.tmp.Zero()
}

// Aggregate adds share1 with share2 on shareOut.
func (pp *PermuteProtocol) Aggregate(share1, share2, shareOut *ring.Poly) {
	pp.dckksContext.contextQ.AddLvl(uint64(len(share1.Coeffs)-1), share1, share2, shareOut)
}

// Decrypt operates a masked decryption on the ciphertext with the given decryption share.
func (pp *PermuteProtocol) Decrypt(ciphertext *ckks.Ciphertext, shareDecrypt RefreshShareDecrypt) {
	pp.dckksContext.contextQ.AddLvl(ciphertext.Level(), ciphertext.Value()[0], shareDecrypt, ciphertext.Value()[0])
}

// Recode takes a masked decrypted ciphertext at modulus Q_0 and returns the same masked decrypted ciphertext at modulus Q_L, with Q_0 << Q_L.
func (pp *PermuteProtocol) Permute(ciphertext *ckks.Ciphertext, permutation []uint64, slots uint64) {
	dckksContext := pp.dckksContext
	context := pp.dckksContext.contextQ

	context.InvNTTLvl(ciphertext.Level(), ciphertext.Value()[0], ciphertext.Value()[0])

	context.PolyToBigint(ciphertext.Value()[0], pp.maskBigint)

	QStart := ring.NewUint(context.Modulus[0])
	for i := uint64(1); i < ciphertext.Level()+1; i++ {
		QStart.Mul(QStart, ring.NewUint(context.Modulus[i]))
	}
	QHalf := new(big.Int).Rsh(QStart, 1)

	maxSlots := pp.dckksContext.n >> 1
	gap := maxSlots / slots

	var sign int
	for i, idx := uint64(0), uint64(0); i < slots; i, idx = i+1, idx+gap {

		// Centers the value around the current modulus
		sign = pp.maskBigint[idx].Cmp(QHalf)
		if sign == 1 || sign == 0 {
			pp.maskBigint[idx].Sub(pp.maskBigint[idx], QStart)
		}

		// Centers the value around the current modulus
		sign = pp.maskBigint[idx+maxSlots].Cmp(QHalf)
		if sign == 1 || sign == 0 {
			pp.maskBigint[idx+maxSlots].Sub(pp.maskBigint[idx+maxSlots], QStart)
		}

		pp.maskComplex[i].Real().SetInt(pp.maskBigint[idx])
		pp.maskComplex[i].Imag().SetInt(pp.maskBigint[idx+maxSlots])
	}

	pp.encoder.FFT(pp.maskComplex, slots)

	pp.permuteWithIndex(permutation, pp.maskComplex)

	pp.encoder.InvFFT(pp.maskComplex, slots)

	for i, jdx, idx := uint64(0), maxSlots, uint64(0); i < slots; i, jdx, idx = i+1, jdx+gap, idx+gap {
		pp.maskComplex[i].Real().Int(pp.maskBigint[idx])
		pp.maskComplex[i].Imag().Int(pp.maskBigint[jdx])
	}

	for ciphertext.Level() != uint64(len(dckksContext.params.Qi)-1) {
		ciphertext.Value()[0].Coeffs = append(ciphertext.Value()[0].Coeffs, make([][]uint64, 1)...)
		ciphertext.Value()[0].Coeffs[ciphertext.Level()] = make([]uint64, dckksContext.n)
	}

	context.SetCoefficientsBigintLvl(ciphertext.Level(), pp.maskBigint, ciphertext.Value()[0])

	context.NTTLvl(ciphertext.Level(), ciphertext.Value()[0], ciphertext.Value()[0])
}

// Recrypt operates a masked recryption on the masked decrypted ciphertext.
func (pp *PermuteProtocol) Recrypt(ciphertext *ckks.Ciphertext, crs *ring.Poly, shareRecrypt RefreshShareRecrypt) {

	pp.dckksContext.contextQ.Add(ciphertext.Value()[0], shareRecrypt, ciphertext.Value()[0])

	ciphertext.Value()[1] = crs.CopyNew()
}
